# Reference: Julie Dudášová et al. (2026).
# "Annual immunogenicity measurements in randomized clinical trials help reveal heterogeneity in long-term vaccine efficacy: a pooled logistic regression modeling study"
# This program is released under the GNU GPL3 license. Copyright © 2026 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.

# Prepare environment ----------------------------------------------------------
rm(list=ls())

# Package loading
library(tidyverse)
library(ggplot2)
library(car)
library(vaxpmx) #  0.0.6, available at CRAN, run install.packages("vaxpmx")
library(nnet)
library(ggcorrplot) 
library(FactoMineR)
library(factoextra)
library(corrplot)
library(mstate)
library(ggridges)
library(gridExtra)
library(timereg)
library(kableExtra)
library(sandwich)

# Code sourcing
source('functions-manuscript-durability.R')

# Path to data specification
pathToData <- NA

pathToIMM <- paste0(pathToData,"/data/IMM.csv")
pathToEFF <- paste0(pathToData,"/data/EFF.csv")
pathToDEMO <- paste0(pathToData,"/data/DEMO.csv")
IMM <- read.csv(pathToIMM,header=TRUE)
EFF <- read.csv(pathToEFF,header=TRUE)
DEMO <- read.csv(pathToDEMO,header=TRUE)

# Clean data for analysis
Data <- ReshapeDataDurabilitySPS(IMM, EFF)

Data$result <- as.numeric(Data$result)

Data <- Data %>% 
  mutate(log2foldrise = log2(foldrise)) %>%
  mutate(log2foldrise_sq = log2(foldrise)*log2(foldrise)) %>%
  mutate(log2result = log2(result)) %>%
  mutate(log2result_sq = log2(result)*log2(result)) %>%
  mutate(log2foldrise_sqrt = sqrt(log2(foldrise))) %>%
  mutate(vaccine = TREATMNT=="Zoster Vaccine") %>%
  filter(time_1 > 0) %>% # shifting baseline to week 6 visit (day of foldrise measurement)
  filter(!is.na(log2foldrise)) %>%
  filter(!is.na(log2result)) %>%
  mutate(rel_y = rel_dy/365.25)
Data$vaccine <- as.numeric(Data$vaccine)

# Adding an additional record for diseased to be able to correctly assign weights for person-time at risk in logistic regression
duplicated_rows <- filter(Data, cens == 1) %>%
  mutate(cens = 2) %>% # cens value for added records for cases is 2
  mutate(time_1 = time_2)

Data <- rbind(Data, duplicated_rows) %>%
  arrange(patno, rel_dy) 

Data <- Data %>%
  mutate(cens=as.numeric(cens>1)) %>%
  # calculating person-days at risk for each titer measurement
  mutate(person_days = time_2-time_1) %>%
  # calculating person-years at risk for each titer measurement
  mutate(person_years = ifelse(cens==1, 1, person_days/365.25))

Data.control <- Data %>% filter(TREATMNT == "Placebo")
Data.vaccinated <- Data %>% filter(TREATMNT == "Zoster Vaccine")
Data.younger <- Data %>% filter(AGEN<70)
Data.older <- Data %>% filter(AGEN>=70)
Data.younger.vaccinated <- Data.younger %>% filter(TREATMNT == "Zoster Vaccine")
Data.younger.control <- Data.younger %>% filter(TREATMNT == "Placebo")
Data.older.vaccinated <- Data.older %>% filter(TREATMNT == "Zoster Vaccine")
Data.older.control <- Data.older %>% filter(TREATMNT == "Placebo")

sum(Data.younger$cens) # number of HZ cases in younger group
sum(Data.older$cens) # number of HZ cases in older group

Data.short <- Data %>%
  select(patno,rel_dy,log2result, AGEN, vaccine, f_mitt_04, t_hz_04, TREATMNT,
         time_1, time_2, cens, person_days, person_years, log2result_sq)

# Number of subjects with information needed for this analysis:
length(unique(Data$patno))
nrow(Data)

length(unique(Data.control$patno))
length(unique(Data.vaccinated$patno))

length(unique(Data.younger$patno))
length(unique(Data.older$patno))

# Number of subjects experiencing herpes zoster
Data.diseased <- Data %>% filter(f_mitt_04==1)
length(unique(Data.diseased$patno))

# Median time-to-HZ, range
summary(Data.diseased$t_hz_04)

# Age range
summary(Data$AGEN)
summary(Data.younger$AGEN)
summary(Data.older$AGEN)

sum(Data$person_years==0) # check for subjects with 0 person-years at risk

# Exploratory Data Analysis ----------------------------------------------------
## 1. Time profiles of ELISA titers by age group ----
p.pop.younger <- ggplot(data = Data.younger, aes(x = rel_y, y = log2result, color = TREATMNT, group = TREATMNT)) + 
  geom_point() +
  ylim(2.5,15)+
  scale_color_manual(values=c("darkslategray3", "gray20"))+
  geom_smooth(aes(fill = TREATMNT), method = "loess", span = 0.75, formula = y~x, se = TRUE) +  
  scale_fill_manual(values = c("darkslategray3", "gray20")) + 
  ggtitle("Younger")+
  xlab("Time (years)") +
  ylab("Log2 ELISA") + 
  theme_bw() 

p.pop.older <- ggplot(data = Data.older, aes(x = rel_y, y = log2result, color = TREATMNT, group = TREATMNT)) + 
  geom_point() +
  ylim(2.5,15)+
  scale_color_manual(values=c("darkslategray3", "gray20"))+
  geom_smooth(aes(fill = TREATMNT), method = "loess", span = 0.75, formula = y~x, se = TRUE) +  
  scale_fill_manual(values = c("darkslategray3", "gray20")) + 
  ggtitle("Older")+
  xlab("Time (years)") +
  ylab("Log2 ELISA") + 
  theme_bw() 

p.cases <- ggplot(data = Data[Data$f_mitt_04==1,], aes(x = rel_y, y = log2result, color = as.factor(cens), group = patno)) + 
  geom_point(aes(fill = as.factor(patno))) +
  geom_vline(aes(xintercept = t_hz_04/365.25)) + 
  scale_color_manual(values=c("darkgray", "firebrick"))+
  facet_wrap(~patno)+
  xlab("Time (years)") +
  ylab("Log2 ELISA") +
  theme_minimal() + 
  theme(
    panel.background = element_blank(), 
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "none"
  )

## 2. Effect of age group on immunogenicity ----

#### Log2 ELISA titer 
set.seed(1)
younger.boxplot.abs <- ggplot(data = Data.younger, aes(x=TREATMNT, y=log2result, fill=factor(cens))) +  
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitterdodge(), pch=21, aes(alpha = factor(cens)))+
  scale_fill_manual(values = c("1" = "firebrick",
                               "0" = "darkgray")) +
  scale_alpha_manual(values = c("1" = 1,
                                "0" = 0.2)) +
  ylim(2.5,15)+
  ggtitle("Younger") +
  xlab("Study arm") +
  ylab("Log2 ELISA") +
  theme_bw()

older.boxplot.abs <- ggplot(data = Data.older, aes(x=TREATMNT, y=log2result, fill=factor(cens))) +  
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitterdodge(), pch=21, aes(alpha = factor(cens)))+
  scale_fill_manual(values = c("1" = "firebrick",
                               "0" = "darkgray")) +
  scale_alpha_manual(values = c("1" = 1,
                                "0" = 0.2)) +
  ylim(2.5,15)+
  ggtitle("Older") +
  xlab("Study arm") +
  ylab("Log2 ELISA") +
  theme_bw()

make_plot <- function(df, title_text) {
  # Non-cases: subject-level means
  noncases_means <- df %>%
    filter(cens == 0) %>%
    group_by(patno, TREATMNT) %>%
    summarise(y = mean(log2result, na.rm = TRUE), .groups = "drop") %>%
    mutate(cens = 0L)
  
  # Cases: raw observations
  cases_raw <- df %>%
    filter(cens == 1) %>%
    transmute(patno, TREATMNT, y = log2result, cens = 1L)
  
  # Combine into a single plot-ready frame
  plot_df <- bind_rows(noncases_means, cases_raw) %>%
    mutate(cens = factor(cens, levels = c(0,1)))
  
  # N per arm & disease group (unique subjects)
  n_df <- plot_df %>%
    group_by(TREATMNT, cens) %>%
    summarise(n = n_distinct(patno), .groups = "drop")
  
  # Two separate label rows so they don't overlap
  n_df_non <- n_df %>% filter(cens == 0) %>% mutate(y = -0.5)
  n_df_case <- n_df %>% filter(cens == 1) %>% mutate(y = -1.0)
  
  ggplot(plot_df, aes(x = TREATMNT, y = y, fill = factor(cens))) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(aes(alpha = cens),
               position = position_jitterdodge(),
               pch = 21) +
    
    scale_fill_manual(values = c("0" = "darkgray", "1" = "firebrick")) +
    scale_alpha_manual(values = c("0" = 0.2, "1" = 1), guide = "none") +
    coord_cartesian(ylim = c(2.5,15)) +
    labs(title = title_text, x = "Study arm", y = "Log2 ELISA") +
    theme_bw()
}

older.boxplot.mean  <- make_plot(Data.older, "Older")
younger.boxplot.mean  <- make_plot(Data.younger, "Younger")

# Log2 ELISA titer at week 6 in vaccinated younger and older groups
Data.vaccinated.week6 <- Data %>% 
  filter(visit == "week 6") %>%
  filter(TREATMNT == "Zoster Vaccine") %>%
  mutate(age_group = as.factor(ifelse(AGEN<70, "Younger", "Older")))
Data.vaccinated.week6$age_group <- factor(Data.vaccinated.week6$age_group, 
                                          levels = c("Younger", "Older"))
overall.boxplot.abs <- ggplot(data = Data.vaccinated.week6, aes(x = age_group, y = log2result, fill = factor(TREATMNT))) + 
  geom_boxplot(outlier.shape = NA, aes(alpha = 0.5)) + 
  geom_point(position = position_jitterdodge(), pch = 21) +
  scale_fill_manual(values = c("Placebo" = "darkslategray3", 
                               "Zoster Vaccine" = "gray20")) +
  ylim(2.5, 15) +
  ggtitle("Overall") +
  xlab("Study arm") +
  ylab("Log2 ELISA") +
  theme_bw()

# Correlate of Risk and Correlate of Protection Analyses -----------------------

## A. Weighted pooled logistic regression (conventional, non-robust standard errors) ----

#### Log2 ELISA titer 
logisticFit.a.abs <- glm(cens ~ log2result, data = Data, family = binomial, weights = person_years) 
summary(logisticFit.a.abs) # model a, best-fitting shape

logisticFit.b.abs <- glm(cens ~ log2result + log2result_sq, data = Data, family = binomial, weights = person_years) 
summary(logisticFit.b.abs) # model b

logisticFit.c.abs <- glm(cens ~ log2result + AGEN, data = Data, family = binomial, weights = person_years) 
summary(logisticFit.c.abs) # model c

logisticFit.d.abs <- glm(cens ~ log2result * AGEN, data = Data, family = binomial, weights = person_years)  
summary(logisticFit.d.abs) # model d

logisticFit.e.CoP.abs <- glm(cens ~ log2result + vaccine, data = Data, family = binomial, weights = person_years) 
summary(logisticFit.e.CoP.abs) # model e, vaccine significant, best-fitting CoP model

logisticFit.f.CoP.abs <- glm(cens ~log2result * vaccine, data = Data, family = binomial, weights = person_years) 
summary(logisticFit.f.CoP.abs) # model f 

logisticFit.g.time.abs <- glm(cens ~ log2result * time_1 + vaccine, data = Data, family = binomial, weights = person_years) 
summary(logisticFit.g.time.abs) # model g, time-dependence of the coefficient for log2result

logisticFit.h.time.abs <- glm(cens ~ log2result + time_1 * vaccine, data = Data, family = binomial, weights = person_years) 
summary(logisticFit.h.time.abs) # model h, lowest AIC

model_m <- logisticFit.h.time.abs
model_j <- logisticFit.e.CoP.abs

par(mfrow=c(1,1))
plot(density(rstandard(model_m, type='deviance'))) 
lines(density(rstandard(model_j, type='deviance')), col='red')

par(mfrow=c(1,2))
scatter.smooth(sqrt(predict(model_m, type='response')), rstandard(model_m), col='gray')
scatter.smooth(sqrt(predict(model_j, type='response')), rstandard(model_j), col='gray')

# Plot deviance residuals
plot(resid(model_m, type = "deviance"), main = "Deviance Residuals")
plot(resid(model_j, type = "deviance"), main = "Deviance Residuals")

# Check for influential observations
influencePlot(model_m)
influencePlot(model_j)

par(mfrow=c(1,2))
plot(cooks.distance(model_m), type='h')
plot(cooks.distance(model_j), type='h')

cooksd_m <- cooks.distance(model_m)
cooksd_j <- cooks.distance(model_j)

length(cooksd_m[cooksd_m > mean(cooksd_m) * 2])
length(cooksd_j[cooksd_j > mean(cooksd_j) * 2])

which(influence.measures(model_m)$is.inf[,'cook.d'] )
which(influence.measures(model_j)$is.inf[,'cook.d'] )

# Diagnostic plots (similar to lm)
plot(model_m)
plot(model_j)

## B. Weighted pooled logistic regression with robust standard errors ----
robust_summary (model = logisticFit.a.abs,
                cluster = "patno") # model a

robust_summary(logisticFit.b.abs,
               cluster = "patno") # model b

robust_summary(logisticFit.c.abs,
               cluster = "patno") # model c

robust_summary(logisticFit.d.abs,
               cluster = "patno") # model d

robust_summary(logisticFit.e.CoP.abs,
               cluster = "patno") # model e

robust_summary(logisticFit.f.CoP.abs,
               cluster = "patno") # model f

robust_summary(logisticFit.g.time.abs,
               cluster = "patno") # model g

robust_summary(logisticFit.h.time.abs,
               cluster = "patno") # model h

## C. Cox proportional hazards regression ----

rows_to_modify <- which(Data$time_1 == Data$time_2 & Data$cens == 1)
Data.cox <- Data

# Shift the case status to the previous row
for (i in rows_to_modify) {
  if (i > 1) {  # Ensure we don't go out of bounds
    Data.cox$cens[i - 1] <- 1  # Assign case status to the previous row
  }
}

# Remove rows where time_1 equals time_2
Data.cox <- Data.cox[Data.cox$time_1 != Data.cox$time_2, ]

# Check the cleaned data
print(Data.cox)

coxFit.a.abs <- coxph(Surv(time_1, time_2, cens) ~ log2result, data = Data.cox) 
summary(coxFit.a.abs) # model a'
AIC(coxFit.a.abs)

coxFit.b.abs <- coxph(Surv(time_1, time_2, cens) ~ log2result + log2result_sq, data = Data.cox) 
summary(coxFit.b.abs) # model b'
AIC(coxFit.b.abs)

coxFit.c.abs <- coxph(Surv(time_1, time_2, cens) ~ log2result + AGEN, data = Data.cox)
summary(coxFit.c.abs) # model c'
AIC(coxFit.c.abs)

coxFit.d.abs <- coxph(Surv(time_1, time_2, cens) ~ log2result * AGEN, data = Data.cox)
summary(coxFit.d.abs) # model d'
AIC(coxFit.d.abs)

coxFit.e.CoP.abs <- coxph(Surv(time_1, time_2, cens) ~ log2result + vaccine, data = Data.cox)
summary(coxFit.e.CoP.abs) # model e', vaccine significant, best-fitting CoP model
AIC(coxFit.e.CoP.abs)

coxFit.f.CoP.abs <- coxph(Surv(time_1, time_2, cens) ~ log2result * vaccine, data = Data.cox)
summary(coxFit.f.CoP.abs) # model f'
AIC(coxFit.f.CoP.abs)

coxFit.g.time.abs <- coxph(Surv(time_1, time_2, cens) ~ log2result * time_1 + vaccine, data = Data.cox)
summary(coxFit.g.time.abs) # model g', time-dependence of the coefficient for log2result
AIC(coxFit.g.time.abs)

coxFit.h.time.abs <- coxph(Surv(time_1, time_2, cens) ~ log2result + time_1 * vaccine, data = Data.cox)
summary(coxFit.h.time.abs) # model h', lowest AIC, best-fitting model overall
AIC(coxFit.h.time.abs)

# Vaccine Efficacy Prediction --------------------------------------------------

#### Log2 ELISA titer 
set.seed(1)

# a: titers at all timepoints are used to estimate VE
logistic.abs <- glm(cens ~ log2result + vaccine, data = Data, family = binomial, weights = person_years) 
summary(logistic.abs)

VE.logistic.abs <- ve(logistic.abs, Data %>% filter(cens==0), nboot = 2000)
VE.logistic.younger.abs <- ve(logistic.abs, Data.younger %>% filter(cens==0), nboot = 2000)
VE.logistic.older.abs <- ve(logistic.abs, Data.older %>% filter(cens==0), nboot = 2000)

caseCountVEZostavax(Data)

# b: titers at peak (first timepoint) are used to estimate VE for the time period 0-1
Data.y0 <- Data %>% filter(visit == "week 6")
set.seed(1)

Data.control.y0 <- Data.y0 %>% filter(vaccine == 0)
Data.vaccinated.y0 <- Data.y0 %>% filter(vaccine == 1)
Data.younger.y0 <- Data.y0 %>% filter(AGEN<70)
Data.older.y0  <- Data.y0 %>% filter(AGEN>=70)
Data.younger.vaccinated.y0  <- Data.younger.y0 %>% filter(TREATMNT == "Zoster Vaccine")
Data.younger.control.y0  <- Data.younger.y0 %>% filter(TREATMNT == "Placebo")
Data.older.vaccinated.y0  <- Data.older.y0 %>% filter(TREATMNT == "Zoster Vaccine")
Data.older.control.y0  <- Data.older.y0 %>% filter(TREATMNT == "Placebo")

VE.logistic.y1.abs <- ve(logistic.abs, Data.y0 %>% filter(cens==0), nboot = 2000)
VE.logistic.younger.y1.abs <- ve(logistic.abs, Data.younger.y0 %>% filter(cens==0), nboot = 2000)
VE.logistic.older.y1.abs <- ve(logistic.abs, Data.older.y0 %>% filter(cens==0), nboot = 2000)

caseCountVEZostavax(Data.y0)

# c: titers at year 1  are used to estimate VE for the time period 1-2
Data.y1 <- Data %>% filter(visit=="year 1")
set.seed(1)

Data.control.y1 <- Data.y1 %>% filter(vaccine == 0)
Data.vaccinated.y1 <- Data.y1 %>% filter(vaccine == 1)
Data.younger.y1 <- Data.y1 %>% filter(AGEN<70)
Data.older.y1  <- Data.y1 %>% filter(AGEN>=70)
Data.younger.vaccinated.y1  <- Data.younger.y1 %>% filter(TREATMNT == "Zoster Vaccine")
Data.younger.control.y1  <- Data.younger.y1 %>% filter(TREATMNT == "Placebo")
Data.older.vaccinated.y1  <- Data.older.y1 %>% filter(TREATMNT == "Zoster Vaccine")
Data.older.control.y1  <- Data.older.y1 %>% filter(TREATMNT == "Placebo")

VE.logistic.y2.abs <- ve(logistic.abs, Data.y1 %>% filter(cens==0), nboot = 2000)
VE.logistic.younger.y2.abs <- ve(logistic.abs, Data.younger.y1 %>% filter(cens==0), nboot = 2000)
VE.logistic.older.y2.abs <- ve(logistic.abs, Data.older.y1 %>% filter(cens==0), nboot = 2000)

caseCountVEZostavax(Data.y1)

# d: titers at year 2  are used to estimate VE for the time period 2-3
Data.y2 <- Data %>% filter(visit=="year 2")
set.seed(1)

Data.control.y2 <- Data.y2 %>% filter(vaccine == 0)
Data.vaccinated.y2 <- Data.y2 %>% filter(vaccine == 1)
Data.younger.y2 <- Data.y2 %>% filter(AGEN<70)
Data.older.y2  <- Data.y2 %>% filter(AGEN>=70)
Data.younger.vaccinated.y2  <- Data.younger.y2 %>% filter(TREATMNT == "Zoster Vaccine")
Data.younger.control.y2  <- Data.younger.y2 %>% filter(TREATMNT == "Placebo")
Data.older.vaccinated.y2  <- Data.older.y2 %>% filter(TREATMNT == "Zoster Vaccine")
Data.older.control.y2  <- Data.older.y2 %>% filter(TREATMNT == "Placebo")

VE.logistic.y3.abs <- ve(logistic.abs, Data.y2 %>% filter(cens==0), nboot = 2000)
VE.logistic.younger.y3.abs <- ve(logistic.abs, Data.younger.y2 %>% filter(cens==0), nboot = 2000)
VE.logistic.older.y3.abs <- ve(logistic.abs, Data.older.y2 %>% filter(cens==0), nboot = 2000)

caseCountVEZostavax(Data.y2)

# e: titers at year 3  are used to estimate VE for the time period 3-4
Data.y3 <- Data %>% filter(visit=="year 3")
set.seed(1)

Data.control.y3 <- Data.y3 %>% filter(vaccine == 0)
Data.vaccinated.y3 <- Data.y3 %>% filter(vaccine == 1)
Data.younger.y3 <- Data.y3 %>% filter(AGEN<70)
Data.older.y3  <- Data.y3 %>% filter(AGEN>=70)
Data.younger.vaccinated.y3  <- Data.younger.y3 %>% filter(TREATMNT == "Zoster Vaccine")
Data.younger.control.y3  <- Data.younger.y3 %>% filter(TREATMNT == "Placebo")
Data.older.vaccinated.y3  <- Data.older.y3 %>% filter(TREATMNT == "Zoster Vaccine")
Data.older.control.y3  <- Data.older.y3 %>% filter(TREATMNT == "Placebo")

VE.logistic.y4.abs <- ve(logistic.abs, Data.y3 %>% filter(cens==0), nboot = 2000)
VE.logistic.younger.y4.abs <- ve(logistic.abs, Data.younger.y3 %>% filter(cens==0), nboot = 2000)
VE.logistic.older.y4.abs <- ve(logistic.abs, Data.older.y3 %>% filter(cens==0), nboot = 2000)

caseCountVEZostavax(Data.y3)

df.younger.abs <- data.frame("Age" = c("Younger","","","",""),
                             "Predictor" = c("Log2 ELISA","","","",""),
                             "Time" = c("Overall","Year 0 to 1", "Year 1 to 2", "Year 2 to 3", "Year 3 to 4"),
                             "VE" = c(round(VE.logistic.younger.abs$VE), 
                                      round(VE.logistic.younger.y1.abs$VE), 
                                      round(VE.logistic.younger.y2.abs$VE), 
                                      round(VE.logistic.younger.y3.abs$VE), 
                                      round(VE.logistic.younger.y4.abs$VE)),
                             "CI LB" = c(round(VE.logistic.younger.abs$CI$LB), 
                                         round(VE.logistic.younger.y1.abs$CI$LB), 
                                         round(VE.logistic.younger.y2.abs$CI$LB), 
                                         round(VE.logistic.younger.y3.abs$CI$LB), 
                                         round(VE.logistic.younger.y4.abs$CI$LB)),
                             "CI UB" = c(round(VE.logistic.younger.abs$CI$UB), 
                                         round(VE.logistic.younger.y1.abs$CI$UB), 
                                         round(VE.logistic.younger.y2.abs$CI$UB), 
                                         round(VE.logistic.younger.y3.abs$CI$UB), 
                                         round(VE.logistic.younger.y4.abs$CI$UB))
)

df.older.abs <- data.frame("Age" = c("Older","","","",""),
                           "Predictor" = c("Log2 ELISA","","","",""),
                           "Time" = c("Overall","Year 0 to 1", "Year 1 to 2", "Year 2 to 3", "Year 3 to 4"),
                           "VE" = c(round(VE.logistic.older.abs$VE), 
                                    round(VE.logistic.older.y1.abs$VE), 
                                    round(VE.logistic.older.y2.abs$VE), 
                                    round(VE.logistic.older.y3.abs$VE), 
                                    round(VE.logistic.older.y4.abs$VE)),
                           "CI LB" = c(round(VE.logistic.older.abs$CI$LB), 
                                       round(VE.logistic.older.y1.abs$CI$LB), 
                                       round(VE.logistic.older.y2.abs$CI$LB), 
                                       round(VE.logistic.older.y3.abs$CI$LB), 
                                       round(VE.logistic.older.y4.abs$CI$LB)),
                           "CI UB" = c(round(VE.logistic.older.abs$CI$UB), 
                                       round(VE.logistic.older.y1.abs$CI$UB), 
                                       round(VE.logistic.older.y2.abs$CI$UB), 
                                       round(VE.logistic.older.y3.abs$CI$UB), 
                                       round(VE.logistic.older.y4.abs$CI$UB))
)

# Time-varying Efficacy visualization ------------------------------------------

## A: final model (used to predict VE), overlay with observed proportions
set.seed(1)

logistic.abs <- glm(cens ~ log2result + vaccine, data = Data, family = binomial, weights = person_years) 
summary(logistic.abs) 

log2result <- seq(min(Data$log2result),max(Data$log2result),0.001)
vaccine <- rbinom(length(log2result),1,1/2)

df <- data.frame(log2result = log2result,
                 vaccine = vaccine)

r <- predict(logistic.abs, df, type = "link", se.fit = TRUE)

inverse_logit = function(x){
  exp(x)/(1+exp(x))
}

df$r <- inverse_logit(r$fit)
df$r.ub <- inverse_logit(r$fit + (r$se.fit*1.96))
df$r.lb <- inverse_logit(r$fit - (r$se.fit*1.96))
df$r.vaccinated <- df$r
df$r.vaccinated[df$vaccine==0] <- NA
df$r.lb.vaccinated <- df$r.lb
df$r.lb.vaccinated[df$vaccine==0] <- NA
df$r.ub.vaccinated <- df$r.ub
df$r.ub.vaccinated[df$vaccine==0] <- NA
df$r.control <- df$r
df$r.control[df$vaccine==1] <- NA
df$r.lb.control <- df$r.lb
df$r.lb.control[df$vaccine==1] <- NA
df$r.ub.control <- df$r.ub
df$r.ub.control[df$vaccine==1] <- NA

risk.zoster <- ggplot(df, aes(x = log2result)) +
  labs(x = 'log2 ELISA time-varying', y = 'Probability of disease') +
  theme_bw() +
  geom_line(data=df[!is.na(df$r.vaccinated),], aes(y = r.vaccinated), color = "gray20", linewidth=1) +
  geom_line(data=df[!is.na(df$r.control),], aes(y = r.control), color = "darkslategray3", linewidth=1) +
  geom_ribbon(data = df[!is.na(df$r.vaccinated),], aes(ymin = r.lb.vaccinated, ymax = r.ub.vaccinated), alpha = 0.15, fill = "gray20") +
  geom_ribbon(data = df[!is.na(df$r.control),], aes(ymin = r.lb.control, ymax = r.ub.control), alpha = 0.4, fill = "darkslategray3") +
  coord_cartesian(ylim = c(0, 0.25)) +
  xlim(min(data$log2result), 12.5) +
  ggtitle("PoD zoster")

OverlayObservedProportionsZostavax(plot_obj = risk.zoster,
                                   data = Data,
                                   predictor = "log2result",
                                   N = 12) 

## B: vaccine efficacy predictions over time

df.older.abs$point_color <- ifelse(
  df.older.abs$VE > 0 & df.older.abs$CI.LB > 0, "darkseagreen",
  ifelse(df.older.abs$VE > 0 & df.older.abs$CI.LB <= 0, "gray30",
         ifelse(df.older.abs$VE < 0 & df.older.abs$CI.LB < 0, "darksalmon", "gray30"))
)

pred.older <- ggplot(df.older.abs[2:5,], aes(x = Time, y = VE)) +
  labs(x = 'time', y = 'Predicted VE') +
  theme_bw() +
  geom_point(aes(group = 1, color = point_color), size=3) +
  geom_errorbar(aes(ymin = CI.LB, ymax = CI.UB), width = 0.1, color = "gray30") +
  coord_cartesian(ylim = c(-100, 100)) +
  ggtitle("predicted VE older") +
  scale_color_identity()  

df.younger.abs$point_color <- ifelse(
  df.younger.abs$VE > 0 & df.younger.abs$CI.LB > 0, "darkseagreen",
  ifelse(df.younger.abs$VE > 0 & df.younger.abs$CI.LB <= 0, "gray30",
         ifelse(df.younger.abs$VE < 0 & df.younger.abs$CI.LB < 0, "darksalmon", "gray30"))
)

pred.younger <- ggplot(df.younger.abs[2:5,], aes(x = Time, y = VE)) +
  labs(x = 'time', y = 'Predicted VE') +
  theme_bw() +
  geom_point(aes(group = 1, color = point_color), size=3) +
  geom_errorbar(aes(ymin = CI.LB, ymax = CI.UB), width = 0.1, color = "gray30") +
  coord_cartesian(ylim = c(-100, 100)) +
  ggtitle("predicted VE younger") +
  scale_color_identity()  

pred.older 
pred.younger 

# Risk Curve Visualization -----------------------------------------------------

#### Log2 ELISA titer 
set.seed(1)
logistic.abs <- glm(cens ~ log2result + vaccine, data = Data, family = binomial, weights = person_years)
summary(logistic.abs) 

log2result <- seq(2.5,15,0.001)
vaccine <- rbinom(length(log2result),1,1/2)

df <- data.frame(log2result = log2result,
                 vaccine = vaccine)

r <- predict(logistic.abs, df, type = "link", se.fit = TRUE)
inverse_logit = function(x){
  exp(x)/(1+exp(x))
}

df$r <- inverse_logit(r$fit)
df$r.ub <- inverse_logit(r$fit + (r$se.fit*1.96))
df$r.lb <- inverse_logit(r$fit - (r$se.fit*1.96))
df$r.vaccinated <- df$r
df$r.vaccinated[df$vaccine==0] <- NA
df$r.lb.vaccinated <- df$r.lb
df$r.lb.vaccinated[df$vaccine==0] <- NA
df$r.ub.vaccinated <- df$r.ub
df$r.ub.vaccinated[df$vaccine==0] <- NA
df$r.control <- df$r
df$r.control[df$vaccine==1] <- NA
df$r.lb.control <- df$r.lb
df$r.lb.control[df$vaccine==1] <- NA
df$r.ub.control <- df$r.ub
df$r.ub.control[df$vaccine==1] <- NA

### Younger
risk.younger.abs <- ggplot(df, aes(x = log2result)) +
  labs(x = 'log2 ELISA time-varying', y = 'Probability of disease') +
  theme_bw() +
  geom_line(data=df[!is.na(df$r.vaccinated),], aes(y = r.vaccinated), color = "gray20", linewidth=1) +
  geom_line(data = df[!is.na(df$r.vaccinated),], aes(y = r.lb.vaccinated), color = "gray20", linetype = "dotted", linewidth=1) +
  geom_line(data = df[!is.na(df$r.vaccinated),], aes(y = r.ub.vaccinated), color = "gray20", linetype = "dotted", linewidth=1) +
  geom_line(data = df[!is.na(df$r.control),], aes(y = r.lb.control), color = "darkslategray3", linetype = "dotted", linewidth=1) +
  geom_line(data = df[!is.na(df$r.control),], aes(y = r.ub.control), color = "darkslategray3", linetype = "dotted", linewidth=1) +
  geom_line(data=df[!is.na(df$r.control),], aes(y = r.control), color = "darkslategray3", linewidth=1) +
  geom_ribbon(data = df[!is.na(df$r.vaccinated),], aes(ymin = r.lb.vaccinated, ymax = r.ub.vaccinated), alpha = 0.1, fill = "gray20") +
  geom_ribbon(data = df[!is.na(df$r.control),], aes(ymin = r.lb.control, ymax = r.ub.control), alpha = 0.1, fill = "darkslategray3") +
  ggtitle("PoD younger")

immunogenicity.younger.y0.abs <- ggplot(Data.younger.y0, aes(x = log2result, y = TREATMNT)) +
  geom_density_ridges(
    jittered_points = TRUE, quantile_lines = TRUE, scale = 0.4, alpha = 0.7,
    vline_size = 2, vline_color = "darkgray",
    point_size = 0.5, point_alpha = 0.3,
    position = position_raincloud(adjust_vlines = TRUE)
  )+
  theme_ridges(grid = FALSE) + 
  theme(legend.position = "none",
        panel.spacing = unit(0.1, "lines"),
        strip.text.x = element_text(size = 8)) +
  xlim(2.5,15) +
  theme_bw() +
  theme(panel.background = element_blank()) + 
  xlab('Day 28') +
  ylab("")+
  theme(panel.background = element_blank())

immunogenicity.younger.y1.abs <- ggplot(Data.younger.y1, aes(x = log2result, y = TREATMNT)) +
  geom_density_ridges(
    jittered_points = TRUE, quantile_lines = TRUE, scale = 0.4, alpha = 0.7,
    vline_size = 2, vline_color = "darkgray",
    point_size = 0.5, point_alpha = 0.3,
    position = position_raincloud(adjust_vlines = TRUE)
  )+
  theme_ridges(grid = FALSE) + 
  theme(legend.position = "none",
        panel.spacing = unit(0.1, "lines"),
        strip.text.x = element_text(size = 8)) +
  xlim(2.5,15) +
  theme_bw() +
  theme(panel.background = element_blank()) + 
  xlab('Year 1') +
  ylab("")+
  theme(panel.background = element_blank()) 

immunogenicity.younger.y2.abs <- ggplot(Data.younger.y2, aes(x = log2result, y = TREATMNT)) +
  geom_density_ridges(
    jittered_points = TRUE, quantile_lines = TRUE, scale = 0.4, alpha = 0.7,
    vline_size = 2, vline_color = "darkgray",
    point_size = 0.5, point_alpha = 0.3,
    position = position_raincloud(adjust_vlines = TRUE)
  )+
  theme_ridges(grid = FALSE) + 
  theme(legend.position = "none",
        panel.spacing = unit(0.1, "lines"),
        strip.text.x = element_text(size = 8)) +
  xlim(2.5,15) +
  theme_bw() +
  theme(panel.background = element_blank()) + 
  xlab('Year 2') +
  ylab("")+
  theme(panel.background = element_blank()) 

immunogenicity.younger.y3.abs <- ggplot(Data.younger.y3, aes(x = log2result, y = TREATMNT)) +
  geom_density_ridges(
    jittered_points = TRUE, quantile_lines = TRUE, scale = 0.4, alpha = 0.7,
    vline_size = 2, vline_color = "darkgray",
    point_size = 0.5, point_alpha = 0.3,
    position = position_raincloud(adjust_vlines = TRUE)
  )+
  theme_ridges(grid = FALSE) + 
  theme(legend.position = "none",
        panel.spacing = unit(0.1, "lines"),
        strip.text.x = element_text(size = 8)) +
  xlim(2.5,15) +
  theme_bw() +
  theme(panel.background = element_blank()) + 
  xlab('Year 3') +
  ylab("")+
  theme(panel.background = element_blank()) 

### Older
risk.older.abs <- ggplot(df, aes(x = log2result)) +
  labs(x = 'log2 ELISA time-varying', y = 'Probability of disease') +
  theme_bw() +
  geom_line(data=df[!is.na(df$r.vaccinated),], aes(y = r.vaccinated), color = "gray20", linewidth=1) +
  geom_line(data = df[!is.na(df$r.vaccinated),], aes(y = r.lb.vaccinated), color = "gray20", linetype = "dotted", linewidth=1) +
  geom_line(data = df[!is.na(df$r.vaccinated),], aes(y = r.ub.vaccinated), color = "gray20", linetype = "dotted", linewidth=1) +
  geom_line(data = df[!is.na(df$r.control),], aes(y = r.lb.control), color = "darkslategray3", linetype = "dotted", linewidth=1) +
  geom_line(data = df[!is.na(df$r.control),], aes(y = r.ub.control), color = "darkslategray3", linetype = "dotted", linewidth=1) +
  geom_line(data=df[!is.na(df$r.control),], aes(y = r.control), color = "darkslategray3", linewidth=1) +
  geom_ribbon(data = df[!is.na(df$r.vaccinated),], aes(ymin = r.lb.vaccinated, ymax = r.ub.vaccinated), alpha = 0.1, fill = "gray20") +
  geom_ribbon(data = df[!is.na(df$r.control),], aes(ymin = r.lb.control, ymax = r.ub.control), alpha = 0.1, fill = "darkslategray3") +
  ggtitle("PoD older") 

immunogenicity.older.y0.abs <- ggplot(Data.older.y0, aes(x = log2result, y = TREATMNT)) +
  geom_density_ridges(
    jittered_points = TRUE, quantile_lines = TRUE, scale = 0.4, alpha = 0.7,
    vline_size = 2, vline_color = "darkgray",
    point_size = 0.5, point_alpha = 0.3,
    position = position_raincloud(adjust_vlines = TRUE)
  )+
  theme_ridges(grid = FALSE) + 
  theme(legend.position = "none",
        panel.spacing = unit(0.1, "lines"),
        strip.text.x = element_text(size = 8)) +
  xlim(2.5,15) +
  theme_bw() +
  theme(panel.background = element_blank()) + 
  xlab('Day 28') +
  ylab("")+
  theme(panel.background = element_blank()) 

immunogenicity.older.y1.abs <- ggplot(Data.older.y1, aes(x = log2result, y = TREATMNT)) +
  geom_density_ridges(
    jittered_points = TRUE, quantile_lines = TRUE, scale = 0.4, alpha = 0.7,
    vline_size = 2, vline_color = "darkgray",
    point_size = 0.5, point_alpha = 0.3,
    position = position_raincloud(adjust_vlines = TRUE)
  )+
  theme_ridges(grid = FALSE) + 
  theme(legend.position = "none",
        panel.spacing = unit(0.1, "lines"),
        strip.text.x = element_text(size = 8)) +
  xlim(2.5,15) +
  theme_bw() +
  theme(panel.background = element_blank()) + 
  xlab('Year 1') +
  ylab("")+
  theme(panel.background = element_blank()) 

immunogenicity.older.y2.abs <- ggplot(Data.older.y2, aes(x = log2result, y = TREATMNT)) +
  geom_density_ridges(
    jittered_points = TRUE, quantile_lines = TRUE, scale = 0.4, alpha = 0.7,
    vline_size = 2, vline_color = "darkgray",
    point_size = 0.5, point_alpha = 0.3,
    position = position_raincloud(adjust_vlines = TRUE)
  )+
  theme_ridges(grid = FALSE) + 
  theme(legend.position = "none",
        panel.spacing = unit(0.1, "lines"),
        strip.text.x = element_text(size = 8)) +
  xlim(2.5,15) +
  theme_bw() +
  theme(panel.background = element_blank()) + 
  xlab('Year 2') +
  ylab("")+
  theme(panel.background = element_blank()) 

immunogenicity.older.y3.abs <- ggplot(Data.older.y3, aes(x = log2result, y = TREATMNT)) +
  geom_density_ridges(
    jittered_points = TRUE, quantile_lines = TRUE, scale = 0.4, alpha = 0.7,
    vline_size = 2, vline_color = "darkgray",
    point_size = 0.5, point_alpha = 0.3,
    position = position_raincloud(adjust_vlines = TRUE)
  )+
  theme_ridges(grid = FALSE) + 
  theme(legend.position = "none",
        panel.spacing = unit(0.1, "lines"),
        strip.text.x = element_text(size = 8)) +
  xlim(2.5,15) +
  theme_bw() +
  theme(panel.background = element_blank()) + 
  xlab('Year 3') +
  ylab("")+
  theme(panel.background = element_blank()) 

grid.arrange(risk.younger.abs, risk.older.abs,
             immunogenicity.younger.y0.abs,immunogenicity.older.y0.abs,
             immunogenicity.younger.y1.abs,immunogenicity.older.y1.abs,
             immunogenicity.younger.y2.abs,immunogenicity.older.y2.abs,
             immunogenicity.younger.y3.abs,immunogenicity.older.y3.abs,
             ncol=2, nrow=5) 

# Supplementary Figure 3 -------------------------------------------------------

## Time-varying PoD model 
set.seed(1)

logistic.abs.time <- glm(cens ~ log2result + time_1 * vaccine, data = Data, family = binomial, weights = person_years) 
summary(logistic.abs) 

# Set sequence of log2result and time
log2_seq <- seq(min(Data$log2result), max(Data$log2result), by = 0.001)
time_vals <- c(28, 365.25, 730.5, 1095.75) # Corresponding to week 4, year 1, year 2, year 3

# Create prediction grid
df <- expand.grid(
  log2result = log2_seq,
  time_1 = time_vals,
  vaccine = c(0, 1)
)

# Predict
r <- predict(logistic.abs.time, newdata = df, type = "link", se.fit = TRUE)

inverse_logit = function(x){
  exp(x)/(1+exp(x))
}

# Add predictions to data frame
df$r <- inverse_logit(r$fit)
df$r.lb <- inverse_logit(r$fit - 1.96 * r$se.fit)
df$r.ub <- inverse_logit(r$fit + 1.96 * r$se.fit)

# Separate vaccinated and control for plotting
df$r.vaccinated <- ifelse(df$vaccine == 1, df$r, NA)
df$r.lb.vaccinated <- ifelse(df$vaccine == 1, df$r.lb, NA)
df$r.ub.vaccinated <- ifelse(df$vaccine == 1, df$r.ub, NA)

df$r.control <- ifelse(df$vaccine == 0, df$r, NA)
df$r.lb.control <- ifelse(df$vaccine == 0, df$r.lb, NA)
df$r.ub.control <- ifelse(df$vaccine == 0, df$r.ub, NA)

# Plot
risk.zoster.time <- ggplot(df, aes(x = log2result)) +
  labs(x = 'log2 ELISA time-varying', y = 'Probability of disease') +
  theme_bw() +
  facet_wrap(~ time_1, labeller = label_both) +
  
  # Vaccinated
  geom_line(data = df[!is.na(df$r.vaccinated),],
            aes(y = r.vaccinated), color = "gray20", linewidth = 1) +
  geom_ribbon(data = df[!is.na(df$r.vaccinated),],
              aes(ymin = r.lb.vaccinated, ymax = r.ub.vaccinated),
              alpha = 0.15, fill = "gray20") +
  
  # Control
  geom_line(data = df[!is.na(df$r.control),],
            aes(y = r.control), color = "darkslategray3", linewidth = 1) +
  geom_ribbon(data = df[!is.na(df$r.control),],
              aes(ymin = r.lb.control, ymax = r.ub.control),
              alpha = 0.4, fill = "darkslategray3") +
  
  coord_cartesian(ylim = c(0, 0.25)) +
  xlim(min(Data$log2result), 12.5) +
  ggtitle("PoD zoster by year and vaccine status")