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

# Data loading -----------------------------------------------------------------
subjectData <- read.csv("subjectData.csv")
symptomaticCases <- read.csv("symptomaticCases.csv")
titerData <- read.csv("titerData.csv")

# Number of events for each serotype
df <- symptomaticCases %>% filter(delayPD3>0)
NdiseaseDENV1 <- nrow(filter(df,serotype==1))
NdiseaseDENV2 <- nrow(filter(df,serotype==2))
NdiseaseDENV3 <- nrow(filter(df,serotype==3))
NdiseaseDENV4 <- nrow(filter(df,serotype==4))
NdiseaseDENVunknown <- nrow(filter(df,serotype==999))
df1 <- data.frame("Serotype" = c("DENV1", "DENV2", "DENV3", "DENV4"),
                  "Ndiseased" = c(NdiseaseDENV1, NdiseaseDENV2, NdiseaseDENV3,NdiseaseDENV4)) 

# Clean data for analysis
Data.all <- ReshapeDataDurability(subjectData, symptomaticCases, titerData)
Data <- Data.all$Data.all

# Adding an additional record for diseased with serotype 2 to be able to correctly assign weights for person-time at risk in logistic regression
duplicated_rows <- filter(Data, serotype_disease == 2 & cens == 1) %>%
  mutate(cens = 2) %>% # cens value for serotype 2, added records for cases is 2
  mutate(time_1 = time_2) 

# Preserving cens information for the other serotypes
Data <- Data %>%
  mutate(cens = ifelse(cens==1 & serotype_disease != 2, 3, cens)) # cens value for the other serotypes is 3

Data <- rbind(Data, duplicated_rows) %>%
  arrange(newID, time_NT) 

Data$time_disease[Data$cens==1] <- NA
Data$serotype_disease[Data$cens==1] <- NA
Data$disease[Data$cens==1] <- 0
Data$out[Data$cens==1] <- 0
Data$disease_any[Data$cens==1] <- 0
Data <- Data %>%
  mutate(cens=as.numeric(cens>1)) %>%
  # calculating person-days at risk for each titer measurement
  mutate(person_days_D2 = time_2-time_1) %>%
  # calculating person-years at risk for each titer measurement
  mutate(person_years_D2 = ifelse(cens==1 & serotype_disease==2, 1, person_days_D2/365.25))

# Serotype-specific disease status
Data$disease_D1 <- as.numeric(Data$serotype_disease==1 & Data$cens==1) 
Data$disease_D1[is.na(Data$disease_D1)] <- 0
Data$disease_D2 <- as.numeric(Data$serotype_disease==2 & Data$cens==1)
Data$disease_D2[is.na(Data$disease_D2)] <- 0
Data$disease_D3 <- as.numeric(Data$serotype_disease==3 & Data$cens==1)
Data$disease_D3[is.na(Data$disease_D3)] <- 0
Data$disease_D4 <- as.numeric(Data$serotype_disease==4 & Data$cens==1)
Data$disease_D4[is.na(Data$disease_D4)] <- 0
Data$serotype_disease <- Data$serotype_disease*Data$cens

# Average titers (not scaled and censored due to multiple timepoints)
Data$avg_NT <- (Data$D1_NT+Data$D2_NT+Data$D3_NT+Data$D4_NT)/4
Data$avg_NT_sq <- Data$avg_NT*Data$avg_NT

# Year of nAb titer measurement
Data <- Data %>% 
  mutate(year = 0)
Data$year[Data$time_NT>365 & Data$time_NT<=730] <- 1
Data$year[Data$time_NT>730 & Data$time_NT<=1095] <- 2
Data$year[Data$time_NT>1095 & Data$time_NT<=1460] <- 3
Data$year[Data$time_NT>1460 & Data$time_NT<=1825] <- 4
Data$year[Data$time_NT>1825 & Data$time_NT <= 2190] <- 5
Data$year[Data$time_NT>2190] <- 6
Data$year <- as.factor(Data$year)

# Serotype-specific nAb titers squared
Data <- Data %>%
  mutate(D1_NT_sq = D1_NT*D1_NT) %>%
  mutate(D2_NT_sq = D2_NT*D2_NT) %>%
  mutate(D3_NT_sq = D3_NT*D3_NT) %>%
  mutate(D4_NT_sq = D4_NT*D4_NT) 

# Additional minor data conversions
Data <- Data %>% 
  mutate(time_NT_y = time_NT/365.25) 

Data$serotype_disease[is.na(Data$serotype_disease)] <- 0
Data$serotype_disease[Data$serotype_disease==999] <- 0

# Sub-dataframes
Data.control <- Data %>% filter(vaccine == 0)
Data.vaccinated <- Data %>% filter(vaccine == 1)
Data.seropositive <- Data %>% filter(serostatus > 0)
Data.seronegative <- Data %>% filter(serostatus == 0)
Data.seropositive.vaccinated <- Data.vaccinated %>% filter(serostatus > 0)
Data.seropositive.control <- Data.control %>% filter(serostatus > 0)
Data.seronegative.vaccinated <- Data.vaccinated %>% filter(serostatus == 0)
Data.seronegative.control <- Data.control %>% filter(serostatus == 0)

sum(Data$person_years_D2==0) # check for subjects with 0 person-years at risk
sum(Data$person_years_D2==0)/length(Data$person_years_D2)
View(Data[Data$year==6,])

Data.short <- Data %>% select(newID, time_NT, D1_NT, D2_NT, D3_NT, D4_NT, serostatus, vaccine, disease, disease_D2, time_disease, time_event, treatment, time_1, time_2, person_days_D2, person_years_D2, D1_NT_sq, D2_NT_sq, D3_NT_sq, D4_NT_sq)
Data.manuscript <- Data %>% 
  select (newID, time_NT, time_1, time_2, D2_NT, disease_D2, serostatus, time_event, person_years_D2, vaccine) %>%
  group_by(newID) %>%
  mutate(disease_D2_newID = ifelse(sum(disease_D2) > 0, 1, 0)) %>%  # information about disease is copied to all rows per subject
  mutate(case_D2 = disease_D2) %>%
  select(-disease_D2) %>%
  filter(newID < 2)

new.line <- data.frame(newID = 0,
                       time_NT = 5,
                       time_1 = 5,
                       time_2 = 30,
                       D2_NT = 3.05,
                       serostatus = 0,
                       time_event = 305,
                       person_years_D2 = 0.068,
                       vaccine = 1,
                       disease_D2_newID = 1,
                       case_D2 = 0)
Data.manuscript <- rbind(new.line, Data.manuscript) 
Data.manuscript$color <- c(0,1,2,0,0,0,0,0,0,0)

ggplot(data = Data.manuscript, aes(x = time_1, y = D2_NT, color = as.factor(color), group = newID)) + 
  geom_line(color = "gray") +
  geom_point(aes(fill = as.factor(color), 
                 shape = as.factor(case_D2)), 
             size = 4) + 
  geom_vline(aes(xintercept = 305)) + 
  scale_color_manual(values = c("darkgray","firebrick", "coral3")) +
  scale_shape_manual(values = c(16, 1)) + 
  xlab("Time (days)") +
  ylab("Log-titer") +
  theme_bw() +
  theme(
    legend.position = "none"
  )

# Descriptive statistics of the data --------------------------------------------
# number of DENV-Any cases (first occurrences)
sum(Data$cens[Data$disease_any==1])
nrow(distinct(Data[Data$disease_any==1,], newID))

# number of DENV1 cases (first occurrences)
sum(Data$cens[Data$disease_D1==1])
nrow(distinct(Data[Data$disease_D1==1,], newID))
NdiseaseDENV1.comprisks <- sum(Data$disease_D1)

# number of DENV2 cases (first occurrences)
sum(Data$cens[Data$disease_D2==1])
nrow(distinct(Data[Data$disease_D2==1,], newID))
NdiseaseDENV2.comprisks <- sum(Data$disease_D2)

# number of DENV3 cases (first occurrences)
sum(Data$cens[Data$disease_D3==1])
nrow(distinct(Data[Data$disease_D3==1,], newID))
NdiseaseDENV3.comprisks <- sum(Data$disease_D3)

# number of DENV4 cases (first occurrences)
sum(Data$cens[Data$disease_D4==1])
nrow(distinct(Data[Data$disease_D4==1,], newID))
NdiseaseDENV4.comprisks <- sum(Data$disease_D4)

nrow(distinct(Data, newID)) # N subjects
nrow(distinct(Data.vaccinated, newID)) # N vaccinated
nrow(distinct(Data.control, newID)) # N control
nrow(distinct(Data.seropositive, newID)) # N seropositive
nrow(distinct(Data.seronegative, newID)) # N seronegative
summary(Data$time_event[Data$disease_D2==1]) # median time to DENV2 VCD

df1.comprisks <- data.frame("Serotype" = c("DENV1", "DENV2", "DENV3", "DENV4"),
                            "Ndiseased" = c(NdiseaseDENV1.comprisks, NdiseaseDENV2.comprisks, NdiseaseDENV3.comprisks,NdiseaseDENV4.comprisks)) 

# Exploratory Data Analysis ----------------------------------------------------
## 1. Time profiles of titers ----

# DENV2 titer

p.DENV2.pop.seropos <- ggplot(data = Data.seropositive %>% filter(disease_D2==0), aes(x = time_NT_y, y = D2_NT, color = treatment, group = treatment)) + 
  geom_point() +
  ylim(-1,12)+
  xlim(0,6.1) +
  scale_color_manual(values=c("darkslategray3", "gray20"))+
  geom_smooth(aes(fill = treatment), method = "loess", span = 0.75, formula = y~x, se = TRUE) +  
  scale_fill_manual(values = c("darkslategray3", "gray20")) + 
  ggtitle("Seropositive")+
  xlab("Time (years)") +
  ylab("Log2 PRNT50") + 
  theme_bw() 

p.DENV2.pop.seroneg <- ggplot(data = Data.seronegative %>% filter(disease_D2==0), aes(x = time_NT_y, y = D2_NT, color = treatment, group = treatment)) + 
  geom_point() +
  ylim(-1,12)+
  xlim(0,6.1) +
  scale_color_manual(values=c("darkslategray3", "gray20"))+
  geom_smooth(aes(fill = treatment), method = "loess", span = 0.75, formula = y~x, se = TRUE) +  
  scale_fill_manual(values = c("darkslategray3", "gray20")) + 
  ggtitle("Seronegative")+
  xlab("Time (years)") +
  ylab("Log2 PRNT50") + 
  theme_bw()

p.DENV2.cases <- ggplot(data = Data[Data$out==2,], aes(x = time_NT_y, y = D2_NT, color = as.factor(cens), group = newID)) + 
  geom_point(aes(fill = as.factor(newID))) +
  geom_vline(aes(xintercept = time_event/365.25)) + 
  scale_color_manual(values=c("darkgray", "firebrick"))+
  facet_wrap(~newID)+
  xlab("Time (years)") +
  ylab("Log2 PRNT50") +
  theme_minimal() + 
  theme(
    panel.background = element_blank(), 
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "none"
  )

## 2. Effect of serostatus on immunogenicity ----

set.seed(1)

DENV2.seronegative.boxplot <- ggplot(data = Data.seronegative, aes(x=treatment, y=D2_NT, fill=factor(disease_D2))) + 
  geom_boxplot(outlier.shape=NA) +
  geom_point(position=position_jitterdodge(), pch=21, aes(alpha = factor(disease_D2)))+
  scale_fill_manual(values = c("1" = "firebrick",
                               "0" = "darkgray")) +
  scale_alpha_manual(values = c("1" = 1,
                                "0" = 0.2)) +
  ylim(-0.1,11.5) +
  ggtitle("Seronegative") +
  xlab("Study arm") +
  ylab("Log2 PRNT50") +
  theme_bw()

DENV2.seropositive.boxplot <- ggplot(data = Data.seropositive, aes(x=treatment, y=D2_NT, fill=factor(disease_D2))) +  
  geom_boxplot(outlier.shape=NA) +
  geom_point(position=position_jitterdodge(), pch=21, aes(alpha = factor(disease_D2)))+
  scale_fill_manual(values = c("1" = "firebrick",
                               "0" = "darkgray")) +
  scale_alpha_manual(values = c("1" = 1,
                                "0" = 0.2)) +
  ylim(-0.1,11.5) +
  ggtitle("Seropositive") +
  xlab("Study arm") +
  ylab("Log2 PRNT50") +
  theme_bw()


make_plot <- function(df, title_text) {
  # Non-cases: subject-level means
  noncases_means <- df %>%
    filter(disease_D2 == 0) %>%
    group_by(newID, treatment) %>%
    summarise(y = mean(D2_NT, na.rm = TRUE), .groups = "drop") %>%
    mutate(disease_D2 = 0L)
  
  # Cases: raw observations
  cases_raw <- df %>%
    filter(disease_D2 == 1) %>%
    transmute(newID, treatment, y = D2_NT, disease_D2 = 1L)
  
  # Combine into a single plot-ready frame
  plot_df <- bind_rows(noncases_means, cases_raw) %>%
    mutate(disease_D2 = factor(disease_D2, levels = c(0,1)))
  
  # N per arm & disease group (unique subjects)
  n_df <- plot_df %>%
    group_by(treatment, disease_D2) %>%
    summarise(n = n_distinct(newID), .groups = "drop")
  
  # Two separate label rows so they don't overlap
  n_df_non <- n_df %>% filter(disease_D2 == 0) %>% mutate(y = -0.5)
  n_df_case <- n_df %>% filter(disease_D2 == 1) %>% mutate(y = -1.0)
  
  ggplot(plot_df, aes(x = treatment, y = y, fill = factor(disease_D2))) +
    geom_boxplot(outlier.shape = NA) +
    geom_point(aes(alpha = disease_D2),
               position = position_jitterdodge(),
               pch = 21) +

    scale_fill_manual(values = c("0" = "darkgray", "1" = "firebrick")) +
    scale_alpha_manual(values = c("0" = 0.2, "1" = 1), guide = "none") +
    coord_cartesian(ylim = c(-0.1, 11.5)) +
    labs(title = title_text, x = "Study arm", y = "Log2 PRNT50") +
    theme_bw()
}

DENV2.seronegative.boxplot.meanD2  <- make_plot(Data.seronegative, "Seronegative")
DENV2.seropositive.boxplot.meanD2  <- make_plot(Data.seropositive, "Seropositive")

Data.vaccinated.day28 <- Data %>% 
  filter(time_NT >0 & time_NT < 45) %>%
  filter(treatment == "vaccine") 
Data.vaccinated.day28$sero <- factor(Data.vaccinated.day28$sero, 
                                     levels = c("seropositive", "seronegative"))

overall.boxplot.abs <- ggplot(data = Data.vaccinated.day28, aes(x = sero, y = D2_NT, fill = factor(treatment))) + 
  geom_boxplot(outlier.shape = NA, aes(alpha = 0.5)) + 
  geom_point(position = position_jitterdodge(), pch = 21) +
  scale_fill_manual(values = c("placebo" = "darkslategray3", 
                               "vaccine" = "gray20")) +
  ylim(-0.1,11.5) +
  ggtitle("Overall") +
  xlab("Group") +
  ylab("Log2 PRNT50") +
  theme_bw()

# Correlate of Risk and Correlate of Protection Analyses ---------------------------------------------------

## A. Weighted pooled logistic regression (conventional, non-robust standard errors) ----

logisticFit.DENV2.A <- glm(disease_D2 ~ D2_NT, data = Data, family = binomial, weights = person_years_D2) 
summary(logisticFit.DENV2.A) # model A

logisticFit.DENV2.B <- glm(disease_D2 ~ D2_NT + D2_NT_sq, data = Data, family = binomial, weights = person_years_D2) 
summary(logisticFit.DENV2.B) # best-fitting shape, winning model (lowest AIC), model B

logisticFit.DENV2.C <- glm(disease_D2 ~ D2_NT + D2_NT_sq + serostatus, data = Data, family = binomial, weights = person_years_D2) 
summary(logisticFit.DENV2.C) # model C

logisticFit.DENV2.D <- glm(disease_D2 ~ D2_NT + D2_NT_sq * serostatus, data = Data, family = binomial, weights = person_years_D2) 
summary(logisticFit.DENV2.D) # model D

logisticFit.DENV2.E <- glm(disease_D2 ~ D2_NT * serostatus + D2_NT_sq, data = Data, family = binomial, weights = person_years_D2) 
summary(logisticFit.DENV2.E) # interaction only marginally significant, model E

logisticFit.DENV2.e2 <- glm(disease_D2 ~ D2_NT * serostatus + D2_NT_sq * serostatus, data = Data, family = binomial, weights = person_years_D2) 
summary(logisticFit.DENV2.e2) 

# Conditional independence of vaccination status with the best-fitting covariate models from CoR analysis
logisticFit.DENV2.F.CoP <- glm(disease_D2 ~ D2_NT + D2_NT_sq + vaccine, data = Data, family = binomial, weights = person_years_D2)
summary(logisticFit.DENV2.F.CoP) # vaccine insignificant, model F

logisticFit.DENV2.G.CoP <- glm(disease_D2 ~ D2_NT + D2_NT_sq * vaccine, data = Data, family = binomial, weights = person_years_D2)
summary(logisticFit.DENV2.G.CoP) # vaccine insignificant, model G

logisticFit.DENV2.H.CoP <- glm(disease_D2 ~ D2_NT * vaccine + D2_NT_sq, data = Data, family = binomial, weights = person_years_D2)
summary(logisticFit.DENV2.H.CoP) # vaccine insignificant, model H

# Check time-dependence of estimated coefficients with best-fitting model
logisticFit.DENV2.I.time <- glm(disease_D2 ~ D2_NT + D2_NT_sq * time_1, data = Data, family = binomial, weights = person_years_D2)
summary(logisticFit.DENV2.I.time) # time-dependence of the coefficient for D2_NT_sq, model I

logisticFit.DENV2.J.time <- glm(disease_D2 ~ D2_NT * time_1 + D2_NT_sq, data = Data, family = binomial, weights = person_years_D2)
summary(logisticFit.DENV2.J.time) # time-dependence of the coefficient for D2_NT_sq, model J

## B. Weighted pooled logistic regression with robust standard errors ----

robust_summary (model = logisticFit.DENV2.A,
                cluster = "newID") # model A

robust_summary(logisticFit.DENV2.B,
               cluster = "newID") # best-fitting shape, winning model (lowest AIC), model B

robust_summary(logisticFit.DENV2.C,
               cluster = "newID") # model C

robust_summary(logisticFit.DENV2.D,
               cluster = "newID") # model D

robust_summary(logisticFit.DENV2.E,
               cluster = "newID") # model E

robust_summary(logisticFit.DENV2.F.CoP,
               cluster = "newID") # model F

robust_summary(logisticFit.DENV2.G.CoP,
               cluster = "newID") # model G

robust_summary(logisticFit.DENV2.H.CoP,
               cluster = "newID") # model H

robust_summary(logisticFit.DENV2.I.time,
               cluster = "newID") # model I

robust_summary(logisticFit.DENV2.J.time,
               cluster = "newID") # model J

## C. Cox proportional hazards regression ----

rows_to_modify <- which(Data$time_1 == Data$time_2 & Data$disease_D2 == 1)
Data.cox <- Data

# Shift the case status to the previous row
for (i in rows_to_modify) {
  if (i > 1) {  # Ensure we don't go out of bounds
    Data.cox$disease_D2[i - 1] <- 1  # Assign case status to the previous row
  }
}

# Remove rows where time_1 equals time_2
Data.cox <- Data.cox[Data.cox$time_1 != Data.cox$time_2, ]

print(Data.cox)

coxFit.DENV2.A <- coxph(Surv(time_1, time_2, disease_D2) ~ D2_NT, data = Data.cox) 
summary(coxFit.DENV2.A) # model A'
AIC(coxFit.DENV2.A)

coxFit.DENV2.B <- coxph(Surv(time_1, time_2, disease_D2) ~ D2_NT + D2_NT_sq, data = Data.cox) 
summary(coxFit.DENV2.B) # model B'
AIC(coxFit.DENV2.B)

coxFit.DENV2.C <- coxph(Surv(time_1, time_2, disease_D2) ~ D2_NT + D2_NT_sq + serostatus, data = Data.cox) 
summary(coxFit.DENV2.C) # model C'
AIC(coxFit.DENV2.C)

coxFit.DENV2.D <- coxph(Surv(time_1, time_2, disease_D2) ~ D2_NT + D2_NT_sq * serostatus, data = Data.cox) 
summary(coxFit.DENV2.D) # model D'
AIC(coxFit.DENV2.D)

coxFit.DENV2.E <- coxph(Surv(time_1, time_2, disease_D2) ~ D2_NT * serostatus + D2_NT_sq, data = Data.cox) 
summary(coxFit.DENV2.E) # model E'
AIC(coxFit.DENV2.E)

coxFit.DENV2.F <- coxph(Surv(time_1, time_2, disease_D2) ~ D2_NT + D2_NT_sq + vaccine, data = Data.cox) 
summary(coxFit.DENV2.F) # model F'
AIC(coxFit.DENV2.F)

coxFit.DENV2.G <- coxph(Surv(time_1, time_2, disease_D2) ~ D2_NT + D2_NT_sq * vaccine, data = Data.cox) 
summary(coxFit.DENV2.G) # model G'
AIC(coxFit.DENV2.G)

coxFit.DENV2.H <- coxph(Surv(time_1, time_2, disease_D2) ~ D2_NT * vaccine + D2_NT_sq, data = Data.cox) 
summary(coxFit.DENV2.H) # model H'
AIC(coxFit.DENV2.H)

coxFit.DENV2.I <- coxph(Surv(time_1, time_2, disease_D2) ~ D2_NT + D2_NT_sq *time_1, data = Data.cox) 
summary(coxFit.DENV2.I) # model I'
AIC(coxFit.DENV2.I)

coxFit.DENV2.J <- coxph(Surv(time_1, time_2, disease_D2) ~ D2_NT * time_1 + D2_NT_sq, data = Data.cox) 
summary(coxFit.DENV2.J) # model J'
AIC(coxFit.DENV2.J)

# Vaccine Efficacy Prediction --------------------------------------------------
set.seed(1)

####### A: titers at all timepoints are used to estimate overall VE
logistic.DENV2 <- glm(disease_D2 ~ D2_NT + D2_NT_sq, data = Data, family = binomial, weights = person_years_D2) 
summary(logistic.DENV2) 

VE.logistic.DENV2 <- ve(logistic.DENV2, Data %>% filter(disease_D2==0), nboot = 2000) # we only use records of DENV2 non-cases to avoid duplicates of titer values (for DENV2 cases, we had added an additional observation to account correctly for person-years in PoD curve fitting)
VE.logistic.DENV2.seronegative <- ve(logistic.DENV2, Data.seronegative %>% filter(disease_D2==0), nboot = 2000)
VE.logistic.DENV2.seropositive <- ve(logistic.DENV2, Data.seropositive %>% filter(disease_D2==0), nboot = 2000)

caseCountVE(Data)

####### B: additional timepoints

## titers at peak (first timepoint) are used to estimate VE for the time period 0-1
Data.y0 <- Data %>% filter(year == 0)
set.seed(1)

Data.control.y0 <- Data.y0 %>% filter(vaccine == 0)
Data.vaccinated.y0 <- Data.y0 %>% filter(vaccine == 1)
Data.seropositive.y0 <- Data.y0 %>% filter(serostatus > 0)
Data.seronegative.y0 <- Data.y0 %>% filter(serostatus == 0)
Data.seropositive.vaccinated.y0 <- Data.vaccinated.y0 %>% filter(serostatus > 0)
Data.seropositive.control.y0 <- Data.control.y0 %>% filter(serostatus > 0)
Data.seronegative.vaccinated.y0 <- Data.vaccinated.y0 %>% filter(serostatus == 0)
Data.seronegative.control.y0 <- Data.control.y0 %>% filter(serostatus == 0)

VE.logistic.DENV2.y1 <- ve(logistic.DENV2, Data.y0 %>% filter(disease_D2==0), nboot = 2000)
VE.logistic.DENV2.seronegative.y1 <- ve(logistic.DENV2, Data.seronegative.y0 %>% filter(disease_D2==0), nboot = 2000)
VE.logistic.DENV2.seropositive.y1 <- ve(logistic.DENV2, Data.seropositive.y0 %>% filter(disease_D2==0), nboot = 2000)

caseCountVE(Data.y0)

## titers at year 1  are used to estimate VE for the time period 1-2
Data.y1 <- Data %>% filter(year == 1)
set.seed(1)

Data.control.y1 <- Data.y1 %>% filter(vaccine == 0)
Data.vaccinated.y1 <- Data.y1 %>% filter(vaccine == 1)
Data.seropositive.y1 <- Data.y1 %>% filter(serostatus > 0)
Data.seronegative.y1 <- Data.y1 %>% filter(serostatus == 0)
Data.seropositive.vaccinated.y1 <- Data.vaccinated.y1 %>% filter(serostatus > 0)
Data.seropositive.control.y1 <- Data.control.y1 %>% filter(serostatus > 0)
Data.seronegative.vaccinated.y1 <- Data.vaccinated.y1 %>% filter(serostatus == 0)
Data.seronegative.control.y1 <- Data.control.y1 %>% filter(serostatus == 0)

VE.logistic.DENV2.y2 <- ve(logistic.DENV2, Data.y1 %>% filter(disease_D2==0), nboot = 2000)
VE.logistic.DENV2.seronegative.y2 <- ve(logistic.DENV2, Data.seronegative.y1 %>% filter(disease_D2==0), nboot = 2000)
VE.logistic.DENV2.seropositive.y2 <- ve(logistic.DENV2, Data.seropositive.y1 %>% filter(disease_D2==0), nboot = 2000)

caseCountVE(Data.y1)

## titers at year 2  are used to estimate VE for the time period 2-3
Data.y2 <- Data %>% filter(year == 2)
set.seed(1)

Data.control.y2 <- Data.y2 %>% filter(vaccine == 0)
Data.vaccinated.y2 <- Data.y2 %>% filter(vaccine == 1)
Data.seropositive.y2 <- Data.y2 %>% filter(serostatus > 0)
Data.seronegative.y2 <- Data.y2 %>% filter(serostatus == 0)
Data.seropositive.vaccinated.y2 <- Data.vaccinated.y2 %>% filter(serostatus > 0)
Data.seropositive.control.y2 <- Data.control.y2 %>% filter(serostatus > 0)
Data.seronegative.vaccinated.y2 <- Data.vaccinated.y2 %>% filter(serostatus == 0)
Data.seronegative.control.y2 <- Data.control.y2 %>% filter(serostatus == 0)

VE.logistic.DENV2.y3 <- ve(logistic.DENV2, Data.y2 %>% filter(disease_D2==0), nboot = 2000)
VE.logistic.DENV2.seronegative.y3 <- ve(logistic.DENV2, Data.seronegative.y2 %>% filter(disease_D2==0), nboot = 2000)
VE.logistic.DENV2.seropositive.y3 <- ve(logistic.DENV2, Data.seropositive.y2 %>% filter(disease_D2==0), nboot = 2000)

caseCountVE(Data.y2)

## titers at year 3 (first timepoint) are used to estimate VE for the time period 3-4
Data.y3 <- Data %>% filter(year == 3)
set.seed(1)

Data.control.y3 <- Data.y3 %>% filter(vaccine == 0)
Data.vaccinated.y3 <- Data.y3 %>% filter(vaccine == 1)
Data.seropositive.y3 <- Data.y3 %>% filter(serostatus > 0)
Data.seronegative.y3 <- Data.y3 %>% filter(serostatus == 0)
Data.seropositive.vaccinated.y3 <- Data.vaccinated.y3 %>% filter(serostatus > 0)
Data.seropositive.control.y3 <- Data.control.y3 %>% filter(serostatus > 0)
Data.seronegative.vaccinated.y3 <- Data.vaccinated.y3 %>% filter(serostatus == 0)
Data.seronegative.control.y3 <- Data.control.y3 %>% filter(serostatus == 0)

VE.logistic.DENV2.y4 <- ve(logistic.DENV2, Data.y3 %>% filter(disease_D2==0), nboot = 2000)
VE.logistic.DENV2.seronegative.y4 <- ve(logistic.DENV2, Data.seronegative.y3 %>% filter(disease_D2==0), nboot = 2000)
VE.logistic.DENV2.seropositive.y4 <- ve(logistic.DENV2, Data.seropositive.y3 %>% filter(disease_D2==0), nboot = 2000)

caseCountVE(Data.y3)

## titers at year 4  are used to estimate VE for the time period 4-5
Data.y4 <- Data %>% filter(year == 4)
set.seed(1)

Data.control.y4 <- Data.y4 %>% filter(vaccine == 0)
Data.vaccinated.y4 <- Data.y4 %>% filter(vaccine == 1)
Data.seropositive.y4 <- Data.y4 %>% filter(serostatus > 0)
Data.seronegative.y4 <- Data.y4 %>% filter(serostatus == 0)
Data.seropositive.vaccinated.y4 <- Data.vaccinated.y4 %>% filter(serostatus > 0)
Data.seropositive.control.y4 <- Data.control.y4 %>% filter(serostatus > 0)
Data.seronegative.vaccinated.y4 <- Data.vaccinated.y4 %>% filter(serostatus == 0)
Data.seronegative.control.y4 <- Data.control.y4 %>% filter(serostatus == 0)

VE.logistic.DENV2.y5 <- ve(logistic.DENV2, Data.y4 %>% filter(disease_D2==0), nboot = 2000)
VE.logistic.DENV2.seronegative.y5 <- ve(logistic.DENV2, Data.seronegative.y4 %>% filter(disease_D2==0), nboot = 2000)
VE.logistic.DENV2.seropositive.y5 <- ve(logistic.DENV2, Data.seropositive.y4 %>% filter(disease_D2==0), nboot = 2000)

caseCountVE(Data.y4)

## titers at year 5  are used to estimate VE for the time period 5-6
Data.y5 <- Data %>% filter(year == 5)
set.seed(1)

Data.control.y5 <- Data.y5 %>% filter(vaccine == 0)
Data.vaccinated.y5 <- Data.y5 %>% filter(vaccine == 1)
Data.seropositive.y5 <- Data.y5 %>% filter(serostatus > 0)
Data.seronegative.y5 <- Data.y5 %>% filter(serostatus == 0)
Data.seropositive.vaccinated.y5 <- Data.vaccinated.y5 %>% filter(serostatus > 0)
Data.seropositive.control.y5 <- Data.control.y5 %>% filter(serostatus > 0)
Data.seronegative.vaccinated.y5 <- Data.vaccinated.y5 %>% filter(serostatus == 0)
Data.seronegative.control.y5 <- Data.control.y5 %>% filter(serostatus == 0)

VE.logistic.DENV2.y6 <- ve(logistic.DENV2, Data.y5 %>% filter(disease_D2==0), nboot = 2000)
VE.logistic.DENV2.seronegative.y6 <- ve(logistic.DENV2, Data.seronegative.y5 %>% filter(disease_D2==0), nboot = 2000)
VE.logistic.DENV2.seropositive.y6 <- ve(logistic.DENV2, Data.seropositive.y5 %>% filter(disease_D2==0), nboot = 2000)

caseCountVE(Data.y5)

## titers at year 6  are used to estimate VE for the time period 6-7
Data.y6 <- Data %>% filter(year == 6)
set.seed(1)

Data.control.y6 <- Data.y6 %>% filter(vaccine == 0)
Data.vaccinated.y6 <- Data.y6 %>% filter(vaccine == 1)
Data.seropositive.y6 <- Data.y6 %>% filter(serostatus > 0)
Data.seronegative.y6 <- Data.y6 %>% filter(serostatus == 0)
Data.seropositive.vaccinated.y6 <- Data.vaccinated.y6 %>% filter(serostatus > 0)
Data.seropositive.control.y6 <- Data.control.y6 %>% filter(serostatus > 0)
Data.seronegative.vaccinated.y6 <- Data.vaccinated.y6 %>% filter(serostatus == 0)
Data.seronegative.control.y6 <- Data.control.y6 %>% filter(serostatus == 0)

VE.logistic.DENV2.y7 <- ve(logistic.DENV2, Data.y6 %>% filter(disease_D2==0), nboot = 2000)
VE.logistic.DENV2.seronegative.y7 <- ve(logistic.DENV2, Data.seronegative.y6 %>% filter(disease_D2==0), nboot = 2000)
VE.logistic.DENV2.seropositive.y7 <- ve(logistic.DENV2, Data.seropositive.y6 %>% filter(disease_D2==0), nboot = 2000)

caseCountVE(Data.y6)

df.seropositive.DENV2 <- data.frame("Serotype"= c("DENV2", "", "","","","","",""),
                                    "Serostatus" = c("Seropositive", "", "","","","","",""),
                                    "Time" =  c("Overall","Year 0 to 1","Year 1 to 2","Year 2 to 3","Year 3 to 4","Year 4 to 5","Year 5 to 6", "Year 6 to 7"),
                                    "VE" = c(round(VE.logistic.DENV2.seropositive$VE),
                                             round(VE.logistic.DENV2.seropositive.y1$VE),
                                             round(VE.logistic.DENV2.seropositive.y2$VE),
                                             round(VE.logistic.DENV2.seropositive.y3$VE),
                                             round(VE.logistic.DENV2.seropositive.y4$VE),
                                             round(VE.logistic.DENV2.seropositive.y5$VE),
                                             round(VE.logistic.DENV2.seropositive.y6$VE),
                                             round(VE.logistic.DENV2.seropositive.y7$VE)),
                                    "CI LB" = c(round(VE.logistic.DENV2.seropositive$CI$LB),
                                                round(VE.logistic.DENV2.seropositive.y1$CI$LB),
                                                round(VE.logistic.DENV2.seropositive.y2$CI$LB),
                                                round(VE.logistic.DENV2.seropositive.y3$CI$LB),
                                                round(VE.logistic.DENV2.seropositive.y4$CI$LB),
                                                round(VE.logistic.DENV2.seropositive.y5$CI$LB),
                                                round(VE.logistic.DENV2.seropositive.y6$CI$LB),
                                                round(VE.logistic.DENV2.seropositive.y7$CI$LB)),
                                    "CI UB" = c(
                                      round(VE.logistic.DENV2.seropositive$CI$UB),
                                      round(VE.logistic.DENV2.seropositive.y1$CI$UB),
                                      round(VE.logistic.DENV2.seropositive.y2$CI$UB),
                                      round(VE.logistic.DENV2.seropositive.y3$CI$UB),
                                      round(VE.logistic.DENV2.seropositive.y4$CI$UB),
                                      round(VE.logistic.DENV2.seropositive.y5$CI$UB),
                                      round(VE.logistic.DENV2.seropositive.y6$CI$UB),
                                      round(VE.logistic.DENV2.seropositive.y7$CI$UB))
)

df.seronegative.DENV2 <- data.frame("Serotype"= c("DENV2", "", "","","","","",""),
                                    "Serostatus" = c("Seronegative", "", "","","","","",""),
                                    "Time" =  c("Overall","Year 0 to 1","Year 1 to 2","Year 2 to 3","Year 3 to 4","Year 4 to 5","Year 5 to 6", "Year 6 to 7"),
                                    "VE" = c(round(VE.logistic.DENV2.seronegative$VE),
                                             round(VE.logistic.DENV2.seronegative.y1$VE),
                                             round(VE.logistic.DENV2.seronegative.y2$VE),
                                             round(VE.logistic.DENV2.seronegative.y3$VE),
                                             round(VE.logistic.DENV2.seronegative.y4$VE),
                                             round(VE.logistic.DENV2.seronegative.y5$VE),
                                             round(VE.logistic.DENV2.seronegative.y6$VE),
                                             round(VE.logistic.DENV2.seronegative.y7$VE)),
                                    "CI LB" = c(round(VE.logistic.DENV2.seronegative$CI$LB),
                                                round(VE.logistic.DENV2.seronegative.y1$CI$LB),
                                                round(VE.logistic.DENV2.seronegative.y2$CI$LB),
                                                round(VE.logistic.DENV2.seronegative.y3$CI$LB),
                                                round(VE.logistic.DENV2.seronegative.y4$CI$LB),
                                                round(VE.logistic.DENV2.seronegative.y5$CI$LB),
                                                round(VE.logistic.DENV2.seronegative.y6$CI$LB),
                                                round(VE.logistic.DENV2.seronegative.y7$CI$LB)),
                                    "CI UB" = c(round(VE.logistic.DENV2.seronegative$CI$UB),
                                                round(VE.logistic.DENV2.seronegative.y1$CI$UB),
                                                round(VE.logistic.DENV2.seronegative.y2$CI$UB),
                                                round(VE.logistic.DENV2.seronegative.y3$CI$UB),
                                                round(VE.logistic.DENV2.seronegative.y4$CI$UB),
                                                round(VE.logistic.DENV2.seronegative.y5$CI$UB),
                                                round(VE.logistic.DENV2.seronegative.y6$CI$UB),
                                                round(VE.logistic.DENV2.seronegative.y7$CI$UB))
)

####### C: additional time points to compare against literature

## titers at year 4 and 5 are used to estimate VE for the time period 4-6
Data.y4.y5 <- Data %>% filter(year == 4 | year == 5)
set.seed(1)

Data.control.y4.y5 <- Data.y4.y5 %>% filter(vaccine == 0)
Data.vaccinated.y4.y5 <- Data.y4.y5 %>% filter(vaccine == 1)
Data.seropositive.y4.y5 <- Data.y4.y5 %>% filter(serostatus > 0)
Data.seronegative.y4.y5 <- Data.y4.y5 %>% filter(serostatus == 0)
Data.seropositive.vaccinated.y4.y5 <- Data.vaccinated.y4.y5 %>% filter(serostatus > 0)
Data.seropositive.control.y4.y5 <- Data.control.y4.y5 %>% filter(serostatus > 0)
Data.seronegative.vaccinated.y4.y5 <- Data.vaccinated.y4.y5 %>% filter(serostatus == 0)
Data.seronegative.control.y4.y5 <- Data.control.y4.y5 %>% filter(serostatus == 0)

VE.logistic.DENV2.y5.y6 <- ve(logistic.DENV2, Data.y4.y5 %>% filter(disease_D2==0), nboot = 2000)
VE.logistic.DENV2.seronegative.y5.y6 <- ve(logistic.DENV2, Data.seronegative.y4.y5 %>% filter(disease_D2==0), nboot = 2000)
VE.logistic.DENV2.seropositive.y5.y6 <- ve(logistic.DENV2, Data.seropositive.y4.y5 %>% filter(disease_D2==0), nboot = 2000)

# Time-varying Efficacy visualization ------------------------------------------
## Figure 3
## A: final model (used to predict VE), overlay with observed proportions
set.seed(1)

logistic.DENV2 <- glm(disease_D2 ~ D2_NT + D2_NT_sq, data = Data, family = binomial, weights = person_years_D2) 

summary(logistic.DENV2) 

d2_nt <- seq(0,7.5,0.001)
d2_nt_sq <- d2_nt*d2_nt
serostatus <- rbinom(length(d2_nt),1,1/2)

df <- data.frame(D2_NT = d2_nt,
                 D2_NT_sq = d2_nt_sq,
                 serostatus = serostatus)

r <- predict(logistic.DENV2, df, type = "link", se.fit = TRUE)
inverse_logit = function(x){
  exp(x)/(1+exp(x))
}

df$r <- inverse_logit(r$fit)
df$r.ub <- inverse_logit(r$fit + (r$se.fit*1.96))
df$r.lb <- inverse_logit(r$fit - (r$se.fit*1.96))

risk.DENV2 <- ggplot(df, aes(x = D2_NT)) +
  labs(x = 'log2 DENV2 TOE SN titer', y = 'Probability of disease') +
  theme_bw() +
  geom_line(aes(y = r), color = "black", linewidth=1) +
  geom_ribbon(aes(ymin = r.lb, ymax = r.ub),
              alpha = 0.1,
              fill = "black")+
  xlim(-0.1,7.5) +
  ylim(0,0.3)+
  ggtitle("DENV2 best fitting model") 

OverlayObservedProportions(plot_obj = risk.DENV2, data = Data, predictor = "D2_NT", N = 17) 

## B: vaccine efficacy predictions over time

df.seronegative.DENV2$point_color <- ifelse(
  df.seronegative.DENV2$VE > 0 & df.seronegative.DENV2$CI.LB > 0, "darkseagreen", 
  ifelse(df.seronegative.DENV2$VE > 0 & df.seronegative.DENV2$CI.LB <= 0, "gray30", 
         ifelse(df.seronegative.DENV2$VE < 0 & df.seronegative.DENV2$CI.LB < 0, "darksalmon", "gray30"))
)

pred.seronegative <- ggplot(df.seronegative.DENV2[2:8,], aes(x = Time, y = VE)) +
  labs(x = 'time', y = 'Predicted VE') +
  theme_bw() +
  geom_point(aes(group = 1, color = point_color), size=3) +
  geom_errorbar(aes(ymin = CI.LB, ymax = CI.UB), width = 0.1, color = "gray30") +
  coord_cartesian(ylim = c(-100, 100)) + 
  ggtitle("predicted VE seronegative") +
  scale_color_identity()  

df.seropositive.DENV2$point_color <- ifelse(
  df.seropositive.DENV2$VE > 0 & df.seropositive.DENV2$CI.LB > 0, "darkseagreen", 
  ifelse(df.seropositive.DENV2$VE > 0 & df.seropositive.DENV2$CI.LB <= 0, "gray30", 
         ifelse(df.seropositive.DENV2$VE < 0 & df.seropositive.DENV2$CI.LB < 0, "darksalmon", "gray30"))
)

pred.seropositive <- ggplot(df.seropositive.DENV2[2:8,], aes(x = Time, y = VE)) +
  labs(x = 'time', y = 'Predicted VE') +
  theme_bw() +
  geom_point(aes(group = 1, color = point_color), size=3) +
  geom_errorbar(aes(ymin = CI.LB, ymax = CI.UB), width = 0.1, color = "gray30") +
  coord_cartesian(ylim = c(-100, 100)) + 
  ggtitle("predicted VE seropositive") +
  scale_color_identity()  

pred.seronegative # save as png 600x600
pred.seropositive # save as png 600x600 

# Risk Curve Visualization -----------------------------------------------------

set.seed(1)

## Logistic, DENV2
logistic.DENV2 <- glm(disease_D2 ~ D2_NT + D2_NT_sq, data = Data, family = binomial, weights = person_years_D2) 

summary(logistic.DENV2) 

d2_nt <- seq(0,7.5,0.001)
d2_nt_sq <- d2_nt*d2_nt
serostatus <- rbinom(length(d2_nt),1,1/2)

df <- data.frame(D2_NT = d2_nt,
                 D2_NT_sq = d2_nt_sq,
                 serostatus = serostatus)

r <- predict(logistic.DENV2, df, type = "link", se.fit = TRUE)
inverse_logit = function(x){
  exp(x)/(1+exp(x))
}

df$r <- inverse_logit(r$fit)
df$r.ub <- inverse_logit(r$fit + (r$se.fit*1.96))
df$r.lb <- inverse_logit(r$fit - (r$se.fit*1.96))
df$r.seropositive <- df$r
df$r.seropositive[df$serostatus==0] <- NA
df$r.lb.seropositive <- df$r.lb
df$r.lb.seropositive[df$serostatus==0] <- NA
df$r.ub.seropositive <- df$r.ub
df$r.ub.seropositive[df$serostatus==0] <- NA
df$r.seronegative <- df$r
df$r.seronegative[df$serostatus==1] <- NA
df$r.lb.seronegative <- df$r.lb
df$r.lb.seronegative[df$serostatus==1] <- NA
df$r.ub.seronegative <- df$r.ub
df$r.ub.seronegative[df$serostatus==1] <- NA

### Seropositive
risk.DENV2.seropositive <- ggplot(df, aes(x = D2_NT)) +
  labs(x = 'log2 DENV2 TOE SN titer', y = 'Probability of disease') +
  theme_bw() +
  geom_line(data=df[!is.na(df$r.seropositive),], aes(y = r.seropositive), color = "black", linewidth=1) +
  geom_line(data=df[!is.na(df$r.seropositive),], aes(y = r.lb.seropositive), color = "black", linetype = "dotted", linewidth=1) +
  geom_line(data=df[!is.na(df$r.seropositive),], aes(y = r.ub.seropositive), color = "black", linetype = "dotted", linewidth=1) +
  geom_line(data=df[!is.na(df$r.seropositive),], aes(y = r.seropositive), color = "black", linewidth=1) +
  geom_ribbon(data = df[!is.na(df$r.seropositive),],
              aes(ymin = r.lb.seropositive, ymax = r.ub.seropositive),
              alpha = 0.1,
              fill = "black")+
  xlim(-0.1,7.5) +
  ggtitle("DENV2 Seropositive") 
risk.DENV2.seropositive

immunogenicity.DENV2.seropositive.y0 <- ggplot(Data.seropositive.y0, aes(x = D2_NT, y = treatment)) +
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
  xlim(0,7.5) +
  theme_bw() +
  theme(panel.background = element_blank()) + 
  xlab('Day 28') +
  ylab("")+
  theme(panel.background = element_blank()) 

immunogenicity.DENV2.seropositive.y1 <- ggplot(Data.seropositive.y1, aes(x = D2_NT, y = treatment)) +
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
  xlim(0,7.5) +
  theme_bw() +
  theme(panel.background = element_blank()) + 
  xlab('Year 1') +
  ylab("")+
  theme(panel.background = element_blank()) 

immunogenicity.DENV2.seropositive.y2 <- ggplot(Data.seropositive.y2, aes(x = D2_NT, y = treatment)) +
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
  xlim(0,7.5) +
  theme_bw() +
  theme(panel.background = element_blank()) + 
  xlab('Year 2') +
  ylab("")+
  theme(panel.background = element_blank()) 

immunogenicity.DENV2.seropositive.y3 <- ggplot(Data.seropositive.y3, aes(x = D2_NT, y = treatment)) +
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
  xlim(0,7.5) +
  theme_bw() +
  theme(panel.background = element_blank()) + 
  xlab('Year 3') +
  ylab("")+
  theme(panel.background = element_blank()) 

immunogenicity.DENV2.seropositive.y4 <- ggplot(Data.seropositive.y4, aes(x = D2_NT, y = treatment)) +
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
  xlim(0,7.5) +
  theme_bw() +
  theme(panel.background = element_blank()) + 
  xlab('Year 4') +
  ylab("")+
  theme(panel.background = element_blank()) 

immunogenicity.DENV2.seropositive.y5 <- ggplot(Data.seropositive.y5, aes(x = D2_NT, y = treatment)) +
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
  xlim(0,7.5) +
  theme_bw() +
  theme(panel.background = element_blank()) + 
  xlab('Year 5') +
  ylab("")+
  theme(panel.background = element_blank()) 

immunogenicity.DENV2.seropositive.y6 <- ggplot(Data.seropositive.y6, aes(x = D2_NT, y = treatment)) +
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
  xlim(0,7.5) +
  theme_bw() +
  theme(panel.background = element_blank()) + 
  xlab('Year 6') +
  ylab("")+
  theme(panel.background = element_blank()) 

### Seronegative
risk.DENV2.seronegative <- ggplot(df, aes(x = D2_NT)) +
  labs(x = 'log2 DENV2 TOE SN titer', y = 'Probability of disease') +
  theme_bw() +
  geom_line(data=df[!is.na(df$r.seronegative),], aes(y = r.seronegative), color = "black", linewidth=1) +
  geom_line(data=df[!is.na(df$r.seronegative),], aes(y = r.lb.seronegative), color = "black", linetype = "dotted", linewidth=1) +
  geom_line(data=df[!is.na(df$r.seropositive),], aes(y = r.ub.seropositive), color = "black", linetype = "dotted", linewidth=1) +
  geom_ribbon(data = df[!is.na(df$r.seronegative),],
              aes(ymin = r.lb.seronegative, ymax = r.ub.seronegative),
              alpha = 0.1,
              fill = "black")+
  xlim(-0.1,7.5) +
  ggtitle("DENV2 Seronegative") 
risk.DENV2.seronegative 

immunogenicity.DENV2.seronegative.y0 <- ggplot(Data.seronegative.y0, aes(x = D2_NT, y = treatment)) +
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
  xlim(0,7.5) +
  theme_bw() +
  theme(panel.background = element_blank()) + 
  xlab('Day 28') +
  ylab("")+
  theme(panel.background = element_blank()) 

immunogenicity.DENV2.seronegative.y1 <- ggplot(Data.seronegative.y1, aes(x = D2_NT, y = treatment)) +
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
  xlim(0,7.5) +
  theme_bw() +
  theme(panel.background = element_blank()) + 
  xlab('Year 1') +
  ylab("")+
  theme(panel.background = element_blank()) 

immunogenicity.DENV2.seronegative.y2 <- ggplot(Data.seronegative.y2, aes(x = D2_NT, y = treatment)) +
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
  xlim(0,7.5) +
  theme_bw() +
  theme(panel.background = element_blank()) + 
  xlab('Year 2') +
  ylab("")+
  theme(panel.background = element_blank()) 

immunogenicity.DENV2.seronegative.y3 <- ggplot(Data.seronegative.y3, aes(x = D2_NT, y = treatment)) +
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
  xlim(0,7.5) +
  theme_bw() +
  theme(panel.background = element_blank()) + 
  xlab('Year 3') +
  ylab("")+
  theme(panel.background = element_blank()) 

immunogenicity.DENV2.seronegative.y4 <- ggplot(Data.seronegative.y4, aes(x = D2_NT, y = treatment)) +
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
  xlim(0,7.5) +
  theme_bw() +
  theme(panel.background = element_blank()) + 
  xlab('Year 4') +
  ylab("")+
  theme(panel.background = element_blank()) 

immunogenicity.DENV2.seronegative.y5 <- ggplot(Data.seronegative.y5, aes(x = D2_NT, y = treatment)) +
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
  xlim(0,7.5) +
  theme_bw() +
  theme(panel.background = element_blank()) + 
  xlab('Year 5') +
  ylab("")+
  theme(panel.background = element_blank()) 

immunogenicity.DENV2.seronegative.y6 <- ggplot(Data.seronegative.y6, aes(x = D2_NT, y = treatment)) +
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
  xlim(0,7.5) +
  theme_bw() +
  theme(panel.background = element_blank()) + 
  xlab('Year 6') +
  ylab("")+
  theme(panel.background = element_blank()) 

grid.arrange(risk.DENV2.seronegative, risk.DENV2.seropositive,
             immunogenicity.DENV2.seronegative.y0,immunogenicity.DENV2.seropositive.y0,
             immunogenicity.DENV2.seronegative.y1,immunogenicity.DENV2.seropositive.y1,
             immunogenicity.DENV2.seronegative.y2,immunogenicity.DENV2.seropositive.y2,
             immunogenicity.DENV2.seronegative.y3,immunogenicity.DENV2.seropositive.y3,
             immunogenicity.DENV2.seronegative.y4, immunogenicity.DENV2.seropositive.y4,
             immunogenicity.DENV2.seronegative.y5,immunogenicity.DENV2.seropositive.y5,
             immunogenicity.DENV2.seronegative.y6,immunogenicity.DENV2.seropositive.y6,
             ncol=2, nrow=8) 