# Reference: Julie Dudášová, Zdeněk Valenta, Jeffrey R. Sachs (2024).
# "Improving precision of vaccine efficacy evaluation using immune correlate data in time-to-event models"
# This program is released under the GNU GPL3 license. Copyright © 2024 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.

# Prepare environment ----------------------------------------------------------
rm(list=ls())

# Package loading
library(tidyverse)
library(survival)
library(ggplot2)
library(car)
library(vaxpmx) # version 0.0.4
library(nnet)
library(ggcorrplot)
library(FactoMineR)
library(factoextra)
library(corrplot)
library(mstate)
library(cmprsk)
library(ggridges)
library(gridExtra)
library(car)

# Code sourcing
source('functions-reshapedata.R')

# Path to data specification
pathToData <- NA

pathToIMM <- paste0(pathToData,"/data/IMM.csv")
pathToEFF <- paste0(pathToData,"/data/EFF.csv")
pathToDEMO <- paste0(pathToData,"/data/DEMO.csv")
IMM <- read.csv(pathToIMM,header=TRUE)
EFF <- read.csv(pathToEFF,header=TRUE)
DEMO <- read.csv(pathToDEMO,header=TRUE)

# Clean data for analysis
Data <- ReshapeDataSPS(IMM, EFF)

# Exclude patients with duplicated values
Data <- Data %>% group_by(patno) %>% filter(n()==1)

Data$foldrise[Data$TREATMNT== "Placebo"] <- 1

Data <- Data %>% 
  mutate(log2foldrise = log2(foldrise)) %>%
  mutate(log2foldrise_sq = log2(foldrise)*log2(foldrise)) %>%
  mutate(log2foldrise_sqrt = sqrt(log2(foldrise))) %>%
  mutate(vaccine = TREATMNT=="Zoster Vaccine") %>%
  mutate(adj_time_event = t_hz_04-rel_dy) # shifting time to week 6 visit (day of foldrise measurement)
Data$vaccine <- as.numeric(Data$vaccine)

Data.control <- Data %>% filter(TREATMNT == "Placebo")
Data.vaccinated <- Data %>% filter(TREATMNT == "Zoster Vaccine")
Data.younger <- Data %>% filter(AGEN<70)
Data.older <- Data %>% filter(AGEN>=70)
Data.younger.vaccinated <- Data.younger %>% filter(TREATMNT == "Zoster Vaccine")
Data.younger.control <- Data.younger %>% filter(TREATMNT == "Placebo")
Data.older.vaccinated <- Data.older %>% filter(TREATMNT == "Zoster Vaccine")
Data.older.control <- Data.older %>% filter(TREATMNT == "Placebo")

sum(Data.younger$f_mitt_04) # number of HZ cases in younger group
sum(Data.older$f_mitt_04) # number of HZ cases in older group

# number of subjects with information needed for this analysis
length(unique(Data$patno)) 
nrow(Data)

length(unique(Data.control$patno)) 
length(unique(Data.vaccinated$patno))

length(unique(Data.younger$patno))
length(unique(Data.older$patno))

# number of subjects experiencing herpes zoster
Data.diseased <- Data %>% filter(f_mitt_04==1)
length(unique(Data.diseased$patno))

# median time-to-HZ, range
summary(Data.diseased$adj_time_event)

# median time-to-event, range
Data.nondiseased <- Data %>% filter(f_mitt_04==0)
summary(Data.nondiseased$adj_time_event)

# age range
summary(Data$AGEN)

# Effect of age group on immunogenicity ---------------------------------------
set.seed(1)
younger.boxplot <- ggplot(data = Data.younger, aes(x=TREATMNT, y=log2foldrise, fill=factor(f_mitt_04))) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitterdodge(), pch=21, aes(alpha = factor(f_mitt_04)))+
  scale_fill_manual(values = c("1" = "firebrick",
                               "0" = "darkgray")) +
  scale_alpha_manual(values = c("1" = 1,
                                "0" = 0.2)) +
  ylim(0,7) +
  ggtitle("Younger") +
  xlab("Study arm") +
  ylab("Log2 Fold rise") +
  theme_bw()

t.test(Data.younger$log2foldrise[Data.younger$TREATMNT=="Zoster Vaccine"], Data.older$log2foldrise[Data.older$TREATMNT=="Zoster Vaccine"])

# Check of t-test assumptions:
# 1. normality
qqnorm(Data.younger$log2foldrise[Data.younger$TREATMNT=="Zoster Vaccine"])
qqline(Data.younger$log2foldrise[Data.younger$TREATMNT=="Zoster Vaccine"])
ks.test(Data.younger$log2foldrise[Data.younger$TREATMNT=="Zoster Vaccine"], "pnorm", mean = mean(Data.younger$log2foldrise[Data.younger$TREATMNT=="Zoster Vaccine"]), sd = sd(Data.younger$log2foldrise[Data.younger$TREATMNT=="Zoster Vaccine"]))

qqnorm(Data.older$log2foldrise[Data.older$TREATMNT=="Zoster Vaccine"])
qqline(Data.older$log2foldrise[Data.older$TREATMNT=="Zoster Vaccine"])
ks.test(Data.older$log2foldrise[Data.younger$TREATMNT=="Zoster Vaccine"], "pnorm", mean = mean(Data.older$log2foldrise[Data.older$TREATMNT=="Zoster Vaccine"]), sd = sd(Data.older$log2foldrise[Data.older$TREATMNT=="Zoster Vaccine"]))

# 2. homogeneity of variances
BDat <- data.frame(
  data = c(Data.younger$log2foldrise[Data.younger$TREATMNT=="Zoster Vaccine"], Data.older$log2foldrise[Data.older$TREATMNT=="Zoster Vaccine"]),
  group = c(rep("younger", length(Data.younger$log2foldrise[Data.younger$TREATMNT=="Zoster Vaccine"])), rep("older", length(Data.older$log2foldrise[Data.older$TREATMNT=="Zoster Vaccine"])))
)
bartlett.test(data ~ group, data = BDat)

## Assumptions violated, thus Man Whitney U test is conducted
wilcox.test(Data.younger$log2foldrise[Data.younger$TREATMNT=="Zoster Vaccine"], Data.older$log2foldrise[Data.older$TREATMNT=="Zoster Vaccine"])


older.boxplot <- ggplot(data = Data.older, aes(x=TREATMNT, y=log2foldrise, fill=factor(f_mitt_04))) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitterdodge(), pch=21, aes(alpha = factor(f_mitt_04)))+
  scale_fill_manual(values = c("1" = "firebrick",
                               "0" = "darkgray")) +
  scale_alpha_manual(values = c("1" = 1,
                                "0" = 0.2)) +
  ylim(0,7) +
  ggtitle("Older") +
  xlab("Study arm") +
  ylab("Log2 Fold rise") +
  theme_bw()

# Correlate of Risk Analysis ---------------------------------------------------
## Cox regression 
# Full covariate model
coxFit.0.a <- coxph(Surv(adj_time_event, f_mitt_04) ~ log2foldrise + AGEN + GENDERN + RACEN, data = Data) 
summary(coxFit.0.a)
cox.zph(coxFit.0.a)
AIC(coxFit.0.a) 

# First find out if "titer" should be a linear or nonlinear predictor
coxFit.1.a <- coxph(Surv(adj_time_event, f_mitt_04) ~ log2foldrise, data = Data)
summary(coxFit.1.a)
cox.zph(coxFit.1.a)
AIC(coxFit.1.a) # winning model because log2foldrise is significant; model (a)

coxFit.1.b <- coxph(Surv(adj_time_event, f_mitt_04) ~ log2foldrise + log2foldrise_sq, data = Data)
summary(coxFit.1.b)
cox.zph(coxFit.1.b)
AIC(coxFit.1.b) # model (b)

# Next, evaluate impact of baseline covariates using model (a)
coxFit.2.a <- coxph(Surv(adj_time_event, f_mitt_04) ~ log2foldrise + age_flag, data = Data) # younger group... age_flag==TRUE
summary(coxFit.2.a)
cox.zph(coxFit.2.a)
AIC(coxFit.2.a) # model (c)

coxFit.2.b <- coxph(Surv(adj_time_event, f_mitt_04) ~ log2foldrise * age_flag, data = Data)
summary(coxFit.2.b)
cox.zph(coxFit.2.b)
AIC(coxFit.2.b) # model (d)

# Correlate of Protection Analysis ---------------------------------------------
# Best-fitting covariate models from CoR analysis
## Cox regression 
coxFit.1.b.CoP <- coxph(Surv(adj_time_event, f_mitt_04) ~ log2foldrise + TREATMNT, data = Data)
summary(coxFit.1.b.CoP)
cox.zph(coxFit.1.b.CoP)
AIC(coxFit.1.b.CoP) # model (e); vaccination status is insignificant

coxFit.1.b.CoP <- coxph(Surv(adj_time_event, f_mitt_04) ~ log2foldrise * TREATMNT, data = Data)
summary(coxFit.1.b.CoP)
cox.zph(coxFit.1.b.CoP)
AIC(coxFit.1.b.CoP) 

# Vaccine Efficacy Prediction --------------------------------------------------
set.seed(1)

## Cox regression 
coxFit.Zost <- coxph(Surv(adj_time_event, f_mitt_04) ~ log2foldrise, data = Data)
summary(coxFit.Zost)
cox.zph(coxFit.Zost)
AIC(coxFit.Zost)

VE.coxFit.Zost <- ve(coxFit.Zost, Data, nboot = 1000)
VE.coxFit.Zost.younger <- ve(coxFit.Zost, Data.younger, nboot = 1000)
VE.coxFit.Zost.older <- ve(coxFit.Zost, Data.older, nboot = 1000)

# Visualizations ---------------------------------------------------------------

# Risk curves for younger and older
## Cox
summary(coxFit.Zost)
log2foldrise <- seq(0,7,0.001)
df <- data.frame(log2foldrise = log2foldrise)

Data$median_obs <- as.numeric(abs(median(Data$log2foldrise)-Data$log2foldrise)==min(abs(median(Data$log2foldrise)-Data$log2foldrise)))

coxFit.Zost$means <- c(log2foldrise = Data$log2foldrise[Data$median_obs==1][1])

lp <- predict(coxFit.Zost, df, type = "lp", se.fit = TRUE)
df$hr <- exp(lp$fit)
df$hr.lb <- exp(lp$fit - 1.96*lp$se.fit)
df$hr.ub <- exp(lp$fit + 1.96*lp$se.fit)

df <- df %>% 
  mutate(hrdiff = abs(hr-1)) %>%
  arrange(hrdiff)
coxFit.Zost$means

set.seed(1)
risk.Zost.younger <- ggplot(df, aes(x = log2foldrise)) +
  labs(x = 'Log2 foldrise', y = 'Hazard ratio') +
  geom_vline(aes(xintercept = 0), linetype=1) +
  geom_hline(aes(yintercept = 1), linetype=3) +
  theme_bw() +
  geom_line(data=df[!is.na(df$hr),], aes(y = hr), color = "darkgray", size=1) +
  geom_ribbon(data = df[!is.na(df$hr),],
              aes(ymin = hr.lb, ymax = hr.ub),
              alpha = 0.3,
              fill = "darkgray")+
  xlim(0,7) +
  ylim(0,1.3) +
  ggtitle("Zost Younger") 
immunogenicity.Zost.younger <- ggplot(Data.younger, aes(x = log2foldrise, y = TREATMNT)) +
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
  xlim(0,7) +
  theme_bw() +
  theme(panel.background = element_blank()) + 
  xlab('Log2 foldrise') +
  ylab("")+
  theme(panel.background = element_blank()) 
grid.arrange(risk.Zost.younger, immunogenicity.Zost.younger, nrow = 2)

set.seed(1)
risk.Zost.older <- ggplot(df, aes(x = log2foldrise)) +
  labs(x = 'Log2 foldrise', y = 'Hazard ratio') +
  geom_vline(aes(xintercept = 0), linetype=1) +
  geom_hline(aes(yintercept = 1), linetype=3) +
  theme_bw() +
  geom_line(data=df[!is.na(df$hr),], aes(y = hr), color = "darkgray", size=1) +
  geom_ribbon(data = df[!is.na(df$hr),],
              aes(ymin = hr.lb, ymax = hr.ub),
              alpha = 0.3,
              fill = "darkgray")+
  xlim(0,7) +
  ylim(0,1.3) +
  ggtitle("Zost Older") 
immunogenicity.Zost.older <- ggplot(Data.older, aes(x = log2foldrise, y = TREATMNT)) +
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
  xlim(0,7) +
  theme_bw() +
  theme(panel.background = element_blank()) + 
  xlab('Log2 foldrise') +
  ylab("")+
  theme(panel.background = element_blank()) 
grid.arrange(risk.Zost.older, immunogenicity.Zost.older, nrow = 2) 

# Tables -----------------------------------------------------------------------

## Table 1 
# Zost Younger
sum(Data.younger.control$adj_time_event)/365.25 # control person-years at risk
sum(Data.younger.control$f_mitt_04) # control cases

sum(Data.younger.vaccinated$adj_time_event)/365.25 # vaccinated person-years at risk
sum(Data.younger.vaccinated$f_mitt_04) # vaccinated cases

lapply(VEwaldCIonDF(A = sum(Data.younger.vaccinated$f_mitt_04==1, na.rm = TRUE),
                    B = sum(Data.younger.vaccinated$adj_time_event)/365.25,
                    C = sum(Data.younger.control$f_mitt_04==1, na.rm = TRUE),
                    D = sum(Data.younger.control$adj_time_event)/365.25),"*",100) # case-count VE
VE.coxFit.Zost.younger  # immunogenicity-based VE

# Zost Older
sum(Data.older.control$adj_time_event)/365.25 # control person-years at risk
sum(Data.older.control$f_mitt_04) # control cases

sum(Data.older.vaccinated$adj_time_event)/365.25 # vaccinated person-years at risk
sum(Data.older.vaccinated$f_mitt_04) # vaccinated cases

lapply(VEwaldCIonDF(A = sum(Data.older.vaccinated$f_mitt_04==1, na.rm = TRUE),
                    B = sum(Data.older.vaccinated$adj_time_event)/365.25,
                    C = sum(Data.older.control$f_mitt_04==1, na.rm = TRUE),
                    D = sum(Data.older.control$adj_time_event)/365.25),"*",100) # case-count VE
VE.coxFit.Zost.older # immunogenicity-based VE

## Case-count overall
lapply(VEwaldCIonDF(A = sum(Data.vaccinated$f_mitt_04==1, na.rm = TRUE),
                    B = sum(Data.vaccinated$adj_time_event)/365.25,
                    C = sum(Data.control$f_mitt_04==1, na.rm = TRUE),
                    D = sum(Data.control$adj_time_event)/365.25),"*",100) # case-count VE
VE.coxFit.Zost

# Numerical results in the text ------------------------------------------------

## Introduction
# Zost case-count in younger and older based on the whole SPS study
EFF.control <- EFF %>% filter(TREAT_CD == 2)
EFF.vaccinated <- EFF %>% filter(TREAT_CD == 1)
EFF.younger <- EFF %>% filter(agen_04<70)
EFF.older <- EFF %>% filter(agen_04>=70)
EFF.younger.vaccinated <- EFF.younger %>% filter(TREAT_CD == 1)
EFF.younger.control <- EFF.younger %>% filter(TREAT_CD == 2)
EFF.older.vaccinated <- EFF.older %>% filter(TREAT_CD == 1)
EFF.older.control <- EFF.older %>% filter(TREAT_CD == 2)

# older
lapply(VEwaldCIonDF(A = sum(EFF.older.vaccinated$f_mitt_04==1, na.rm = TRUE),
                    B = round(sum(EFF.older.vaccinated$t_hz_04)/365.25),
                    C = sum(EFF.older.control$f_mitt_04==1, na.rm = TRUE),
                    D = round(sum(EFF.older.control$t_hz_04)/365.25)),"*",100) # case-count VE

# younger
lapply(VEwaldCIonDF(A = sum(EFF.younger.vaccinated$f_mitt_04==1, na.rm = TRUE),
                    B = sum(EFF.younger.vaccinated$t_hz_04)/365.25,
                    C = sum(EFF.younger.control$f_mitt_04==1, na.rm = TRUE),
                    D = sum(EFF.younger.control$t_hz_04)/365.25),"*",100) # case-count VE

# overall
lapply(VEwaldCIonDF(A = sum(EFF.vaccinated$f_mitt_04==1, na.rm = TRUE),
                    B = round(sum(EFF.vaccinated$t_hz_04)/365.25),
                    C = sum(EFF.control$f_mitt_04==1, na.rm = TRUE),
                    D = round(sum(EFF.control$t_hz_04)/365.25)),"*",100) # case-count VE

# Updated Table 1 ------------------------------------------------------------

VE.RR <- function(model, contrast, model.summary = T, alpha = 0.05){
  coef <- model$coefficients
  if (model.summary == T){print(summary(model))}
  RR <- exp(t(contrast)%*%model$coefficients)
  serr <- sqrt(t(contrast)%*%model$var%*%contrast);
  LL.RR <- exp(t(contrast)%*%model$coefficients - qnorm(1-alpha/2)*serr)
  UL.RR <- exp(t(contrast)%*%model$coefficients + qnorm(1-alpha/2)*serr)
  res.RR <- c(RR, LL.RR, UL.RR)
  VE.RR <- c(1-RR, 1-UL.RR, 1-LL.RR)
  res <- rbind(res.RR, VE.RR)
  rownames(res) <- c('RR', 'VE')
  colnames(res) <- c('Estimate', 'LL','UL')
  return(res)
}

zost.int <- coxph(Surv(adj_time_event, f_mitt_04) ~ TREATMNT*age_flag, data = Data)
summary(zost.int)
cox.zph(zost.int)
AIC(zost.int)

# Estimating RR and VE within the older subjects' age group
VE.RR(model = zost.int, contrast = c(1,0,0))

# Estimating RR and VE within the younger subjects' age group
VE.RR(model = zost.int, contrast = c(1,0,1))
