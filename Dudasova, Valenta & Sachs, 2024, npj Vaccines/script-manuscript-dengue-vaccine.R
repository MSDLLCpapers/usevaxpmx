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

# Data loading
subjectData <- read.csv("subjectData.csv")
symptomaticCases <- read.csv("symptomaticCases.csv")
titerData <- read.csv("titerData.csv")

# Clean data for analysis
Data.all <- ReshapeData(subjectData, symptomaticCases, titerData)

# Number of datapoints for each serotype
df <- symptomaticCases %>% filter(delayPD3>0)
NdiseaseDENV1 <- nrow(filter(df,serotype==1))
NdiseaseDENV2 <- nrow(filter(df,serotype==2))
NdiseaseDENV3 <- nrow(filter(df,serotype==3))
NdiseaseDENV4 <- nrow(filter(df,serotype==4))
df1 <- data.frame("Serotype" = c("DENV1", "DENV2", "DENV3", "DENV4"),
                  "Ndiseased" = c(NdiseaseDENV1, NdiseaseDENV2, NdiseaseDENV3,NdiseaseDENV4)
) 

# Models will combine all serotypes
Data <- Data.all$Data.all

# Adjusted time-to-event to a new time - measurement of day 30 SN titer (specific to each subject)
Data <- Data %>% mutate(time_event_adj = time_event - time_NT)

# Scaled and centered titers
Data$out <- as.factor(Data$out)
m_D1_NT <- mean(Data$D1_NT)
sd_D1_NT <- sd(Data$D1_NT)
m_D2_NT <- mean(Data$D2_NT)
sd_D2_NT <- sd(Data$D2_NT)
m_D3_NT <- mean(Data$D3_NT)
sd_D3_NT <- sd(Data$D3_NT)
m_D4_NT <- mean(Data$D4_NT)
sd_D4_NT <- sd(Data$D4_NT)
Data$sc_D1_NT <- (Data$D1_NT-m_D1_NT)/sd_D1_NT # standard
Data$sc_D2_NT <- (Data$D2_NT-m_D2_NT)/sd_D2_NT
Data$sc_D3_NT <- (Data$D3_NT-m_D3_NT)/sd_D3_NT
Data$sc_D4_NT <- (Data$D4_NT-m_D4_NT)/sd_D4_NT

Data$sc_avg_NT <- (Data$sc_D1_NT+Data$sc_D2_NT+Data$sc_D3_NT+Data$sc_D4_NT)/4
Data$sc_avg_NT_sq <- Data$sc_avg_NT*Data$sc_avg_NT

Data$disease_D1 <- as.numeric(Data$serotype_disease==1)
Data$disease_D1[is.na(Data$disease_D1)] <- 0
Data$disease_D2 <- as.numeric(Data$serotype_disease==2)
Data$disease_D2[is.na(Data$disease_D2)] <- 0
Data$disease_D3 <- as.numeric(Data$serotype_disease==3)
Data$disease_D3[is.na(Data$disease_D3)] <- 0
Data$disease_D4 <- as.numeric(Data$serotype_disease==4)
Data$disease_D4[is.na(Data$disease_D4)] <- 0

Data.control <- Data %>% filter(vaccine == 0)
Data.vaccinated <- Data %>% filter(vaccine == 1)
Data.seropositive <- Data %>% filter(serostatus > 0)
Data.seronegative <- Data %>% filter(serostatus == 0)
Data.seropositive.vaccinated <- Data.vaccinated %>% filter(serostatus > 0)
Data.seropositive.control <- Data.control %>% filter(serostatus > 0)
Data.seronegative.vaccinated <- Data.vaccinated %>% filter(serostatus == 0)
Data.seronegative.control <- Data.control %>% filter(serostatus == 0)

# Decriptive statistics of the data --------------------------------------------
sum(Data$disease_any) # number of DENV-Any cases (first occurrences)
sum(Data$disease_D1) # number of DENV1 cases (first occurrences)
sum(Data$disease_D2) # number of DENV2 cases (first occurrences)
sum(Data$disease_D3) # number of DENV3 cases (first occurrences)
sum(Data$disease_D4) # number of DENV4 cases (first occurrences)

summary(Data$time_event_adj[Data$disease_any==1]) # median time to VCD
summary(Data$time_event_adj) # median time to event
summary(Data$time_NT) # median time of the day 30 visit

nrow(Data) # N subjects
nrow(Data.vaccinated) # N vaccinated
nrow(Data.control) # N control

# Effect of serostatus on immunogenicity ---------------------------------------

Data.seronegative$serotype_disease[is.na(Data.seronegative$serotype_disease)] <- 0
Data.seronegative$serotype_disease[Data.seronegative$serotype_disease==999] <- 0

Data.seropositive$serotype_disease[is.na(Data.seropositive$serotype_disease)] <- 0
Data.seropositive$serotype_disease[Data.seropositive$serotype_disease==999] <- 0

set.seed(1)

## DENV1
DENV1.seronegative.boxplot <- ggplot(data = Data.seronegative, aes(x=treatment, y=D1_NT, fill=factor(disease_D1))) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitterdodge(), pch=21, aes(alpha = factor(disease_D1)))+
  scale_fill_manual(values = c("1" = "firebrick",
                               "0" = "darkgray")) +
  scale_alpha_manual(values = c("1" = 1,
                                "0" = 0.2)) +
  ylim(-0.1,11.5) +
  ggtitle("Seronegative") +
  xlab("Study arm") +
  ylab("Log2 PRNT50") +
  theme_bw()

t.test(Data.seronegative$D1_NT[Data.seronegative$treatment=="vaccine"], 
       Data.seropositive$D1_NT[Data.seropositive$treatment=="vaccine"])

# Check of t-test assumptions:
# 1. normality
qqnorm(Data.seronegative$D1_NT[Data.seronegative$treatment=="vaccine"])
qqline(Data.seronegative$D1_NT[Data.seronegative$treatment=="vaccine"])
ks.test(Data.seronegative$D1_NT[Data.seronegative$treatment=="vaccine"], "pnorm", mean = mean(Data.seronegative$D1_NT[Data.seronegative$treatment=="vaccine"]), sd = sd(Data.seronegative$D1_NT[Data.seronegative$treatment=="vaccine"]))

qqnorm(Data.seropositive$D1_NT[Data.seropositive$treatment=="vaccine"])
qqline(Data.seropositive$D1_NT[Data.seropositive$treatment=="vaccine"])
ks.test(Data.seropositive$D1_NT[Data.seropositive$treatment=="vaccine"], "pnorm", mean = mean(Data.seropositive$D1_NT[Data.seropositive$treatment=="vaccine"]), sd = sd(Data.seropositive$D1_NT[Data.seropositive$treatment=="vaccine"]))

# 2. homogeneity of variances
BDat <- data.frame(
  data = c(Data.seronegative$D1_NT[Data.seronegative$treatment=="vaccine"], Data.seropositive$D1_NT[Data.seropositive$treatment=="vaccine"]),
  group = c(rep("seronegative", length(Data.seronegative$D1_NT[Data.seronegative$treatment=="vaccine"])), rep("seropositive", length(Data.seropositive$D1_NT[Data.seropositive$treatment=="vaccine"])))
)
bartlett.test(data ~ group, data = BDat)

DENV1.seropositive.boxplot <- ggplot(data = Data.seropositive, aes(x=treatment, y=D1_NT, fill=factor(disease_D1))) +  
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitterdodge(), pch=21, aes(alpha = factor(disease_D1)))+
  scale_fill_manual(values = c("1" = "firebrick",
                               "0" = "darkgray")) +
  scale_alpha_manual(values = c("1" = 1,
                                "0" = 0.2)) +
  ylim(-0.1,11.5) +
  ggtitle("Seropositive") +
  xlab("Study arm") +
  ylab("Log2 PRNT50") +
  theme_bw()

## DENV2
DENV2.seronegative.boxplot <- ggplot(data = Data.seronegative, aes(x=treatment, y=D2_NT, fill=factor(disease_D2))) +  
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitterdodge(), pch=21, aes(alpha = factor(disease_D2)))+
  scale_fill_manual(values = c("1" = "darksalmon",
                               "0" = "darkgray")) +
  scale_alpha_manual(values = c("1" = 1,
                                "0" = 0.2)) +
  ylim(-0.1,11.5) +
  ggtitle("Seronegative") +
  xlab("Study arm") +
  ylab("Log2 PRNT50") +
  theme_bw()

t.test(Data.seronegative$D2_NT[Data.seronegative$treatment=="vaccine"], Data.seropositive$D2_NT[Data.seropositive$treatment=="vaccine"])
wilcox.test(Data.seronegative$D2_NT[Data.seronegative$treatment=="vaccine"], Data.seropositive$D2_NT[Data.seropositive$treatment=="vaccine"])

# Check of t-test assumptions:
# 1. normality
qqnorm(Data.seronegative$D2_NT[Data.seronegative$treatment=="vaccine"])
qqline(Data.seronegative$D2_NT[Data.seronegative$treatment=="vaccine"])
ks.test(Data.seronegative$D2_NT[Data.seronegative$treatment=="vaccine"], "pnorm", mean = mean(Data.seronegative$D2_NT[Data.seronegative$treatment=="vaccine"]), sd = sd(Data.seronegative$D2_NT[Data.seronegative$treatment=="vaccine"]))

qqnorm(Data.seropositive$D2_NT[Data.seropositive$treatment=="vaccine"])
qqline(Data.seropositive$D2_NT[Data.seropositive$treatment=="vaccine"])
ks.test(Data.seropositive$D2_NT[Data.seropositive$treatment=="vaccine"], "pnorm", mean = mean(Data.seropositive$D2_NT[Data.seropositive$treatment=="vaccine"]), sd = sd(Data.seropositive$D2_NT[Data.seropositive$treatment=="vaccine"]))

# 2. homogeneity of variances
BDat <- data.frame(
  data = c(Data.seronegative$D2_NT[Data.seronegative$treatment=="vaccine"], Data.seropositive$D2_NT[Data.seropositive$treatment=="vaccine"]),
  group = c(rep("seronegative", length(Data.seronegative$D2_NT[Data.seronegative$treatment=="vaccine"])), rep("seropositive", length(Data.seropositive$D2_NT[Data.seropositive$treatment=="vaccine"])))
)
bartlett.test(data ~ group, data = BDat)

DENV2.seropositive.boxplot <- ggplot(data = Data.seropositive, aes(x=treatment, y=D2_NT, fill=factor(disease_D2))) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitterdodge(), pch=21, aes(alpha = factor(disease_D2)))+
  scale_fill_manual(values = c("1" = "darksalmon",
                               "0" = "darkgray")) +
  scale_alpha_manual(values = c("1" = 1,
                                "0" = 0.2)) +
  ylim(-0.1,11.5) +
  ggtitle("Seropositive") +
  xlab("Study arm") +
  ylab("Log2 PRNT50") +
  theme_bw()

## DENV3
DENV3.seronegative.boxplot <- ggplot(data = Data.seronegative, aes(x=treatment, y=D3_NT, fill=factor(disease_D3))) +  
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitterdodge(), pch=21, aes(alpha = factor(disease_D3)))+
  scale_fill_manual(values = c("1" = "darkred",
                               "0" = "darkgray")) +
  scale_alpha_manual(values = c("1" = 1,
                                "0" = 0.2)) +
  ylim(-0.1,11.5) +
  ggtitle("Seronegative") +
  xlab("Study arm") +
  ylab("Log2 PRNT50") +
  theme_bw()

t.test(Data.seronegative$D3_NT[Data.seronegative$treatment=="vaccine"], Data.seropositive$D3_NT[Data.seropositive$treatment=="vaccine"])

# Check of t-test assumptions:
# 1. normality
qqnorm(Data.seronegative$D3_NT[Data.seronegative$treatment=="vaccine"])
qqline(Data.seronegative$D3_NT[Data.seronegative$treatment=="vaccine"])
ks.test(Data.seronegative$D3_NT[Data.seronegative$treatment=="vaccine"], "pnorm", mean = mean(Data.seronegative$D3_NT[Data.seronegative$treatment=="vaccine"]), sd = sd(Data.seronegative$D1_NT[Data.seronegative$treatment=="vaccine"]))

qqnorm(Data.seropositive$D3_NT[Data.seropositive$treatment=="vaccine"])
qqline(Data.seropositive$D3_NT[Data.seropositive$treatment=="vaccine"])
ks.test(Data.seropositive$D3_NT[Data.seropositive$treatment=="vaccine"], "pnorm", mean = mean(Data.seropositive$D3_NT[Data.seropositive$treatment=="vaccine"]), sd = sd(Data.seropositive$D3_NT[Data.seropositive$treatment=="vaccine"]))

# 2. homogeneity of variances
BDat <- data.frame(
  data = c(Data.seronegative$D3_NT[Data.seronegative$treatment=="vaccine"], Data.seropositive$D3_NT[Data.seropositive$treatment=="vaccine"]),
  group = c(rep("seronegative", length(Data.seronegative$D3_NT[Data.seronegative$treatment=="vaccine"])), rep("seropositive", length(Data.seropositive$D3_NT[Data.seropositive$treatment=="vaccine"])))
)
bartlett.test(data ~ group, data = BDat)

DENV3.seropositive.boxplot <- ggplot(data = Data.seropositive, aes(x=treatment, y=D3_NT, fill=factor(disease_D3))) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitterdodge(), pch=21, aes(alpha = factor(disease_D3)))+
  scale_fill_manual(values = c("1" = "darkred",
                               "0" = "darkgray")) +
  scale_alpha_manual(values = c("1" = 1,
                                "0" = 0.2)) +
  ylim(-0.1,11.5) +
  ggtitle("Seropositive") +
  xlab("Study arm") +
  ylab("Log2 PRNT50") +
  theme_bw()

## DENV4
DENV4.seronegative.boxplot <- ggplot(data = Data.seronegative, aes(x=treatment, y=D4_NT, fill=factor(disease_D4))) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitterdodge(), pch=21, aes(alpha = factor(disease_D4)))+
  scale_fill_manual(values = c("1" = "firebrick1",
                               "0" = "darkgray")) +
  scale_alpha_manual(values = c("1" = 1,
                                "0" = 0.2)) +
  ylim(-0.1,11.5) +
  ggtitle("Seronegative") +
  xlab("Study arm") +
  ylab("Log2 PRNT50") +
  theme_bw()

t.test(Data.seronegative$D4_NT[Data.seronegative$treatment=="vaccine"], Data.seropositive$D4_NT[Data.seropositive$treatment=="vaccine"])
t.test(Data.seronegative$D4_NT[Data.seronegative$treatment=="vaccine"], Data.seropositive$D4_NT[Data.seropositive$treatment=="vaccine"], var.equal = FALSE)

# Check of t-test assumptions:
# 1. normality
qqnorm(Data.seronegative$D4_NT[Data.seronegative$treatment=="vaccine"])
qqline(Data.seronegative$D4_NT[Data.seronegative$treatment=="vaccine"])
ks.test(Data.seronegative$D4_NT[Data.seronegative$treatment=="vaccine"], "pnorm", mean = mean(Data.seronegative$D4_NT[Data.seronegative$treatment=="vaccine"]), sd = sd(Data.seronegative$D4_NT[Data.seronegative$treatment=="vaccine"]))

qqnorm(Data.seropositive$D4_NT[Data.seropositive$treatment=="vaccine"])
qqline(Data.seropositive$D4_NT[Data.seropositive$treatment=="vaccine"])
ks.test(Data.seropositive$D4_NT[Data.seropositive$treatment=="vaccine"], "pnorm", mean = mean(Data.seropositive$D4_NT[Data.seropositive$treatment=="vaccine"]), sd = sd(Data.seropositive$D4_NT[Data.seropositive$treatment=="vaccine"]))

# 2. homogeneity of variances
BDat <- data.frame(
  data = c(Data.seronegative$D4_NT[Data.seronegative$treatment=="vaccine"], Data.seropositive$D4_NT[Data.seropositive$treatment=="vaccine"]),
  group = c(rep("seronegative", length(Data.seronegative$D4_NT[Data.seronegative$treatment=="vaccine"])), rep("seropositive", length(Data.seropositive$D4_NT[Data.seropositive$treatment=="vaccine"])))
)
bartlett.test(data ~ group, data = BDat) # not met, therefore Welsch's t-test conducted above

DENV4.seropositive.boxplot <- ggplot(data = Data.seropositive, aes(x=treatment, y=D4_NT,fill=factor(disease_D4))) + 
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitterdodge(), pch=21, aes(alpha = factor(disease_D4)))+
  scale_fill_manual(values = c("1" = "firebrick1",
                               "0" = "darkgray")) +
  scale_alpha_manual(values = c("1" = 1,
                                "0" = 0.2)) +
  ylim(-0.1,11.5) +
  ggtitle("Seropositive") +
  xlab("Study arm") +
  ylab("Log2 PRNT50") +
  theme_bw()

# Correlate of Risk Analysis ---------------------------------------------------
## 1. Multicollinearity and PCA
NTdata <- Data %>% select(newID, sc_D1_NT, sc_D2_NT, sc_D3_NT, sc_D4_NT)

# VIF
tempFit.a <- glm(disease_any ~ sc_D1_NT + sc_D2_NT + sc_D3_NT + sc_D4_NT + vaccine + serostatus, data = Data, family = binomial()) 
vif(tempFit.a)

# Correlation matrix
r <- cor(NTdata[,2:5], use="complete.obs", method = "pearson")
ggcorrplot(r, 
           hc.order = TRUE,
           lab = TRUE,
           type = "upper",
           outline.color = "white",
           colors = c("#E46726", "white", "#6D9EC1")) 

# PCA
my_pca <- prcomp(NTdata[,2:5], scale = TRUE,
                 center = TRUE, retx = T)
summary(my_pca) # keep first two components
res.pca <- PCA(NTdata[,2:5], ncp = 2, graph = FALSE, scale.unit = TRUE)

Data$PC1 <- my_pca$x[,1]
Data$PC2 <- my_pca$x[,2]
Data$PC3 <- my_pca$x[,3]
Data$PC4 <- my_pca$x[,4]
pcadata <- Data %>% select(newID, PC1, PC2, PC3, PC4)
tempFit.b <- glm(disease_any ~ PC1 + PC2 + PC3 + PC4 + vaccine + serostatus , data = Data, family = binomial()) 
vif(tempFit.b)
r <- cor(pcadata[,2:5], use="complete.obs")
ggcorrplot(r, 
           hc.order = TRUE,
           lab = TRUE)

# Visualization of average titer as the predictor ------------------------------------------------------

set.seed(1)

## DENV-Any
DENVAny.seronegative.boxplot <- ggplot(data = Data.seronegative, aes(x=treatment, y=sc_avg_NT, fill = factor(serotype_disease))) +  
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitterdodge(), pch=21, aes(alpha = factor(serotype_disease)))+
  scale_fill_manual(values = c("0" = "gray",
                               "999" = "gray",
                               "1" = "firebrick",
                               "2" = "darksalmon",
                               "3" = "darkred",
                               "4" = "firebrick1")) +
  scale_alpha_manual(values = c("0" = 0.2,
                                "999" = 0.2,
                                "1" = 1,
                                "2" = 1,
                                "3" = 1,
                                "4" = 1)) +
  ylim(-3.1,2.1) +
  ggtitle("Seronegative") +
  xlab("Study arm") +
  ylab("average titer") +
  theme_bw()

DENVAny.seropositive.boxplot <- ggplot(data = Data.seropositive, aes(x=treatment, y=sc_avg_NT, fill = factor(serotype_disease))) +  
  geom_boxplot(outlier.shape=NA) + 
  geom_point(position=position_jitterdodge(), pch=21, aes(alpha = factor(serotype_disease)))+
  scale_fill_manual(values = c("0" = "gray",
                               "999" = "gray",
                               "1" = "firebrick",
                               "2" = "darksalmon",
                               "3" = "darkred",
                               "4" = "firebrick1")) +
  scale_alpha_manual(values = c("0" = 0.2,
                                "999" = 0.2,
                                "1" = 1,
                                "2" = 1,
                                "3" = 1,
                                "4" = 1)) +
  ylim(-3.1,2.1) +
  ggtitle("Seropositive") +
  xlab("Study arm") +
  ylab("Log2 PRNT50") +
  theme_bw()

## 2. Cox regression with DENV-Any endpoint
############################# PCA 
## First find out if "titer" (represented by principal components) should be a linear or nonlinear predictor
Data$PC1_sq <- Data$PC1*Data$PC1
Data$PC2_sq <- Data$PC2*Data$PC2

coxFit.1.a <- coxph(Surv(time_event_adj, disease_any) ~ PC1 + PC2, data = Data)
summary(coxFit.1.a)
cox.zph(coxFit.1.a)
AIC(coxFit.1.a) 

coxFit.1.b <- coxph(Surv(time_event_adj, disease_any) ~ PC1 + PC1_sq, data = Data)
summary(coxFit.1.b)
cox.zph(coxFit.1.b)
AIC(coxFit.1.b) # best-fitting shape

coxFit.1.c <- coxph(Surv(time_event_adj, disease_any) ~ PC1 + PC2 + PC2_sq, data = Data)
summary(coxFit.1.c)
cox.zph(coxFit.1.c)
AIC(coxFit.1.c)

coxFit.1.d <- coxph(Surv(time_event_adj, disease_any) ~ PC1 + PC1_sq + PC2 + PC2_sq, data = Data)
summary(coxFit.1.d)
cox.zph(coxFit.1.d)
AIC(coxFit.1.d)

## Next, evaluate impact of baseline covariates using full covariate model using model 1.b
coxFit.2.a <- coxph(Surv(time_event_adj, disease_any) ~ PC1 + PC1_sq + serostatus, data = Data)
summary(coxFit.2.a)
cox.zph(coxFit.2.a)
AIC(coxFit.2.a) # winning model, i.e. effect of serostatus is significant, no interactions significant; model 0a-I

coxFit.2.b <- coxph(Surv(time_event_adj, disease_any) ~ PC1 * serostatus + PC1_sq, data = Data)
summary(coxFit.2.b)
cox.zph(coxFit.2.b)
AIC(coxFit.2.b)

coxFit.2.c <- coxph(Surv(time_event_adj, disease_any) ~ PC1 + PC1_sq * serostatus, data = Data)
summary(coxFit.2.c)
cox.zph(coxFit.2.c)
AIC(coxFit.2.c)

############################# average titer 
## First find out if "titer" should be a linear or nonlinear predictor

coxFit.1.A <- coxph(Surv(time_event_adj, disease_any) ~ sc_avg_NT, data = Data)
summary(coxFit.1.A)
cox.zph(coxFit.1.A)
AIC(coxFit.1.A) 

coxFit.1.B <- coxph(Surv(time_event_adj, disease_any) ~ sc_avg_NT + sc_avg_NT_sq, data = Data)
summary(coxFit.1.B)
cox.zph(coxFit.1.B)
AIC(coxFit.1.B) # best-fitting shape

## Next, evaluate impact of baseline covariates using full covariate model using model 1.B
coxFit.2.A <- coxph(Surv(time_event_adj, disease_any) ~ sc_avg_NT + sc_avg_NT_sq + serostatus, data = Data)
summary(coxFit.2.A)
cox.zph(coxFit.2.A)
AIC(coxFit.2.A) # winning model, i.e. effect of serostatus is significant, no interactions significant; model 0b-I

coxFit.2.B <- coxph(Surv(time_event_adj, disease_any) ~ sc_avg_NT * serostatus + sc_avg_NT_sq, data = Data)
summary(coxFit.2.B)
cox.zph(coxFit.2.B)
AIC(coxFit.2.B)

coxFit.2.C <- coxph(Surv(time_event_adj, disease_any) ~ sc_avg_NT + sc_avg_NT_sq * serostatus, data = Data)
summary(coxFit.2.C)
cox.zph(coxFit.2.C)
AIC(coxFit.2.C)

## 3. CSH Cox model with DENV1-4 endpoints
############################# PCA 
Data <- Data %>% 
  mutate(cause = serotype_disease)
Data$cause[is.na(Data$cause)] <- "event-free"
Data$cause[Data$cause==999] <- "event-free"
Data$status <- Data$cause
Data$status[Data$cause=="event-free"] <- 0
Data$cause <- as.factor(Data$cause)
Data$status <- as.factor(Data$status)

table(Data$status, Data$cause, dnn = c("Status", "Cause"))
(tmat <- trans.comprisk(4, names = c("event-free", "1", "2","3","4")))

#### DENV1

## First find out if "titer" should be a linear or nonlinear predictor
CSH.DENV1.1.a <- coxph(Surv(time_event_adj, status == 1) ~ PC1, data = Data); 
summary(CSH.DENV1.1.a)
AIC(CSH.DENV1.1.a) # best-fitting shape

CSH.DENV1.1.b <- coxph(Surv(time_event_adj, status == 1) ~ PC1 + PC1_sq, data = Data); 
summary(CSH.DENV1.1.b)
AIC(CSH.DENV1.1.b)

## Next, evaluate impact of baseline covariates using full covariate model using model 1.a
CSH.DENV1.2.a <- coxph(Surv(time_event_adj, status == 1) ~ PC1 + serostatus, data = Data); 
summary(CSH.DENV1.2.a)
AIC(CSH.DENV1.2.a)
cox.zph(CSH.DENV1.2.a)

CSH.DENV1.2.b <- coxph(Surv(time_event_adj, status == 1) ~ PC1 * serostatus, data = Data); 
summary(CSH.DENV1.2.b)
AIC(CSH.DENV1.2.b) # winning model; model 1a-I
cox.zph(CSH.DENV1.2.b)

#### DENV2

## First find out if "titer" should be a linear or nonlinear predictor
CSH.DENV2.1.a <- coxph(Surv(time_event_adj, status == 2) ~ PC1, data = Data); 
summary(CSH.DENV2.1.a)
AIC(CSH.DENV2.1.a)

CSH.DENV2.1.b <- coxph(Surv(time_event_adj, status == 2) ~ PC1 + PC2, data = Data); 
summary(CSH.DENV2.1.b)
AIC(CSH.DENV2.1.b)

CSH.DENV2.1.c <- coxph(Surv(time_event_adj, status == 2) ~ PC1 + PC1_sq + PC2, data = Data); 
summary(CSH.DENV2.1.c)
AIC(CSH.DENV2.1.c)

CSH.DENV2.1.c <- coxph(Surv(time_event_adj, status == 2) ~ PC1 + PC1_sq , data = Data); 
summary(CSH.DENV2.1.c)
AIC(CSH.DENV2.1.c) # best-fitting shape
cox.zph(CSH.DENV2.1.c)

## Next, evaluate impact of baseline covariates using full covariate model using model 1.c
CSH.DENV2.2.a <- coxph(Surv(time_event_adj, status == 2) ~ PC1 + PC1_sq + serostatus, data = Data); 
summary(CSH.DENV2.2.a)
AIC(CSH.DENV2.2.a) # winning model; 2a-I
cox.zph(CSH.DENV2.2.a)

CSH.DENV2.2.b <- coxph(Surv(time_event_adj, status == 2) ~ PC1 * serostatus + PC1_sq, data = Data); 
summary(CSH.DENV2.2.b)
AIC(CSH.DENV2.2.b)

CSH.DENV2.2.c <- coxph(Surv(time_event_adj, status == 2) ~ PC1  + PC1_sq * serostatus, data = Data); 
summary(CSH.DENV2.2.c)
AIC(CSH.DENV2.2.c)

CSH.DENV2.2.d <- coxph(Surv(time_event_adj, status == 2) ~ PC1  + PC1_sq  + PC2* serostatus, data = Data); 
summary(CSH.DENV2.2.d)
AIC(CSH.DENV2.2.d)

#### DENV3

## First find out if "titer" should be a linear or nonlinear predictor
CSH.DENV3.1.a <- coxph(Surv(time_event_adj, status == 3) ~ PC1, data = Data); 
summary(CSH.DENV3.1.a)
AIC(CSH.DENV3.1.a) # best-fitting shape; winning model; model 3a-I
cox.zph(CSH.DENV3.1.a)

CSH.DENV3.1.b <- coxph(Surv(time_event_adj, status == 3) ~ PC1 + PC2, data = Data); 
summary(CSH.DENV3.1.b)
AIC(CSH.DENV3.1.b)

CSH.DENV3.1.c <- coxph(Surv(time_event_adj, status == 3) ~ PC1 + PC1_sq, data = Data); 
summary(CSH.DENV3.1.c)
AIC(CSH.DENV3.1.c) 

## Next, evaluate impact of baseline covariates using full covariate model using model 1.a
CSH.DENV3.2.a <- coxph(Surv(time_event_adj, status == 3) ~ PC1 + serostatus, data = Data); 
summary(CSH.DENV3.2.a)
AIC(CSH.DENV3.2.a) 

CSH.DENV3.2.b <- coxph(Surv(time_event_adj, status == 3) ~ PC1 * serostatus , data = Data); 
summary(CSH.DENV3.2.b)
AIC(CSH.DENV3.2.b)

#### DENV4

## First find out if "titer" should be a linear or nonlinear predictor
CSH.DENV4.1.a <- coxph(Surv(time_event_adj, status == 4) ~ PC1 + PC2, data = Data); 
summary(CSH.DENV4.1.a)
AIC(CSH.DENV4.1.a)

CSH.DENV4.1.b <- coxph(Surv(time_event_adj, status == 4) ~ PC2, data = Data); 
summary(CSH.DENV4.1.b)
AIC(CSH.DENV4.1.b) # best-fitting shape

CSH.DENV4.1.c <- coxph(Surv(time_event_adj, status == 4) ~ PC2 + PC2_sq, data = Data); 
summary(CSH.DENV4.1.c)
AIC(CSH.DENV4.1.c) 

## Next, evaluate impact of baseline covariates using full covariate model using model 1.b
CSH.DENV4.2.a <- coxph(Surv(time_event_adj, status == 4) ~ PC2 + serostatus, data = Data); 
summary(CSH.DENV4.2.a)
AIC(CSH.DENV4.2.a) # winning model; model 4a-I
cox.zph(CSH.DENV4.2.a)

CSH.DENV4.2.b <- coxph(Surv(time_event_adj, status == 4) ~ PC2* serostatus, data = Data); 
summary(CSH.DENV4.2.b)
AIC(CSH.DENV4.2.b)

############################# average titer 

#### DENV1

## First find out if "titer" should be a linear or nonlinear predictor
CSH.DENV1.1.a <- coxph(Surv(time_event_adj, status == 1) ~ sc_avg_NT, data = Data); 
summary(CSH.DENV1.1.a)
AIC(CSH.DENV1.1.a) # best-fitting shape

CSH.DENV1.1.b <- coxph(Surv(time_event_adj, status == 1) ~ sc_avg_NT + sc_avg_NT_sq, data = Data); 
summary(CSH.DENV1.1.b)
AIC(CSH.DENV1.1.b)

## Next, evaluate impact of baseline covariates using full covariate model using model 1.a
CSH.DENV1.2.a <- coxph(Surv(time_event_adj, status == 1) ~ sc_avg_NT + serostatus, data = Data); 
summary(CSH.DENV1.2.a)
AIC(CSH.DENV1.2.a)
cox.zph(CSH.DENV1.2.a)

CSH.DENV1.2.b <- coxph(Surv(time_event_adj, status == 1) ~ sc_avg_NT * serostatus, data = Data); 
summary(CSH.DENV1.2.b)
AIC(CSH.DENV1.2.b) # winning model; model 1b-I
cox.zph(CSH.DENV1.2.b)

beta <- coef(CSH.DENV1.2.b)
beta.seroneg <- beta[1]
beta.seropos <- beta[1] + beta[3]
cov <- vcov(CSH.DENV1.2.b)
cov.seroneg <- cov[1]
cov.seropos <- cov[1,1]+cov[3,3]+2*cov[1,3]
ci <- 0.95
n <- nrow(Data)
beta.seroneg.CIhigh <- beta.seroneg + qt(ci + (1 - ci)/2, df = n - 1)*sqrt(cov.seroneg)
beta.seroneg.CIlow <- beta.seroneg - qt(ci + (1 - ci)/2, df = n - 1)*sqrt(cov.seroneg)
beta.seropos.CIhigh <- beta.seropos + qt(ci + (1 - ci)/2, df = n - 1)*sqrt(cov.seropos)
beta.seropos.CIlow <- beta.seropos - qt(ci + (1 - ci)/2, df = n - 1)*sqrt(cov.seropos)

exp(beta.seropos)
exp(beta.seropos.CIlow)
exp(beta.seropos.CIhigh)

#### DENV2

## First find out if "titer" should be a linear or nonlinear predictor
CSH.DENV2.1.a <- coxph(Surv(time_event_adj, status == 2) ~ sc_avg_NT, data = Data); 
summary(CSH.DENV2.1.a)
AIC(CSH.DENV2.1.a)

CSH.DENV2.1.b <- coxph(Surv(time_event_adj, status == 2) ~ sc_avg_NT + sc_avg_NT_sq, data = Data); 
summary(CSH.DENV2.1.b)
AIC(CSH.DENV2.1.b) # best-fitting shape

## Next, evaluate impact of baseline covariates using full covariate model using model 1.c
CSH.DENV2.2.a <- coxph(Surv(time_event_adj, status == 2) ~ sc_avg_NT + sc_avg_NT_sq + serostatus, data = Data); 
summary(CSH.DENV2.2.a)
AIC(CSH.DENV2.2.a) # winning model; model 2b-I
cox.zph(CSH.DENV2.2.a)

CSH.DENV2.2.b <- coxph(Surv(time_event_adj, status == 2) ~ sc_avg_NT + sc_avg_NT_sq * serostatus , data = Data); 
summary(CSH.DENV2.2.b)
AIC(CSH.DENV2.2.b)

CSH.DENV2.2.c <- coxph(Surv(time_event_adj, status == 2) ~ sc_avg_NT * serostatus + sc_avg_NT_sq, data = Data); 
summary(CSH.DENV2.2.c)
AIC(CSH.DENV2.2.c)

#### DENV3

## First find out if "titer" should be a linear or nonlinear predictor
CSH.DENV3.1.a <- coxph(Surv(time_event_adj, status == 3) ~ sc_avg_NT, data = Data); 
summary(CSH.DENV3.1.a)
AIC(CSH.DENV3.1.a) # best-fitting shape; winning model; model 3b-I
cox.zph(CSH.DENV3.1.a)

CSH.DENV3.1.b <- coxph(Surv(time_event_adj, status == 3) ~ sc_avg_NT + sc_avg_NT_sq, data = Data); 
summary(CSH.DENV3.1.b)
AIC(CSH.DENV3.1.b)

## Next, evaluate impact of baseline covariates using full covariate model using model 1.a
CSH.DENV3.2.a <- coxph(Surv(time_event_adj, status == 3) ~ sc_avg_NT + serostatus , data = Data); 
summary(CSH.DENV3.2.a)
AIC(CSH.DENV3.2.a) # (not winning, better fit without adjustment)

CSH.DENV3.2.b <- coxph(Surv(time_event_adj, status == 3) ~ sc_avg_NT * serostatus , data = Data); 
summary(CSH.DENV3.2.b)
AIC(CSH.DENV3.2.b)

#### DENV4

## First find out if "titer" should be a linear or nonlinear predictor
CSH.DENV4.1.a <- coxph(Surv(time_event_adj, status == 4) ~ sc_avg_NT, data = Data); 
summary(CSH.DENV4.1.a)
AIC(CSH.DENV4.1.a) # best-fitting shape; winning model; model 4b-I
cox.zph(CSH.DENV4.1.a)

CSH.DENV4.1.b <- coxph(Surv(time_event_adj, status == 4) ~ sc_avg_NT + sc_avg_NT_sq, data = Data); 
summary(CSH.DENV4.1.b)
AIC(CSH.DENV4.1.b) 

## Next, evaluate impact of baseline covariates using full covariate model using model 1.b
CSH.DENV4.2.a <- coxph(Surv(time_event_adj, status == 4) ~ sc_avg_NT + serostatus, data = Data); 
summary(CSH.DENV4.2.a)
AIC(CSH.DENV4.2.a) # not winning model

CSH.DENV4.2.b <- coxph(Surv(time_event_adj, status == 4) ~ sc_avg_NT * serostatus, data = Data); 
summary(CSH.DENV4.2.b)
AIC(CSH.DENV4.2.b)

## 4. Fine-Gray model with DENV1-4 endpoints 
############################# PCA 
nm1 <- c("time_event_adj","status",
         "PC1","PC2",
         "PC1_sq","PC2_sq",
         "sc_avg_NT","sc_avg_NT_sq",
         "vaccine","serostatus")
x1 <- Data[,nm1];

#### DENV1
## Base model
fg.data.1 <- finegray(Surv(time_event_adj, status) ~ .,
                      data=x1, na.action=na.omit, etype = "1")
fg.psurv1.a <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PC1 + PC2,
                     weight=fgwt, data=fg.data.1)
summary(fg.psurv1.a)
AIC(fg.psurv1.a)

fg.psurv1.b <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PC1,
                     weight=fgwt, data=fg.data.1)
summary(fg.psurv1.b)
AIC(fg.psurv1.b) # best-fitting shape

fg.psurv1.c <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PC1 + PC1_sq,
                     weight=fgwt, data=fg.data.1)
summary(fg.psurv1.c)
AIC(fg.psurv1.c)

## Covariate model
fg.psurv1.d <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PC1 + serostatus,
                     weight=fgwt, data=fg.data.1)
summary(fg.psurv1.d)
AIC(fg.psurv1.d)

fg.psurv1.e <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PC1 * serostatus,
                     weight=fgwt, data=fg.data.1)
summary(fg.psurv1.e)
AIC(fg.psurv1.e) # winning model; model 1c-I
cox.zph(fg.psurv1.e)

#### DENV2
## Base model
fg.data.2 <- finegray(Surv(time_event_adj, status) ~ .,
                      data=x1, na.action=na.omit, etype = "2")
fg.psurv2.a <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PC1 + PC2,
                     weight=fgwt, data=fg.data.2)
summary(fg.psurv2.a)
AIC(fg.psurv2.a)

fg.psurv2.b <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PC1 ,
                     weight=fgwt, data=fg.data.2)
summary(fg.psurv2.b)
AIC(fg.psurv2.b)

fg.psurv2.c <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PC1 + PC1_sq + PC2,
                     weight=fgwt, data=fg.data.2)
summary(fg.psurv2.c)
AIC(fg.psurv2.c) 

fg.psurv2.c <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PC1 + PC1_sq,
                     weight=fgwt, data=fg.data.2)
summary(fg.psurv2.c)
AIC(fg.psurv2.c) # best-fitting shape
cox.zph(fg.psurv2.c)

## Covariate model
fg.psurv2.d <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PC1 + PC1_sq+ serostatus,
                     weight=fgwt, data=fg.data.2)
summary(fg.psurv2.d)
AIC(fg.psurv2.d) # winning model; model 2c-I
cox.zph(fg.psurv2.d)

fg.psurv2.e <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PC1 * serostatus+ PC1_sq  ,
                     weight=fgwt, data=fg.data.2)
summary(fg.psurv2.e)
AIC(fg.psurv2.e)

fg.psurv2.f <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PC1 + PC1_sq * serostatus,
                     weight=fgwt, data=fg.data.2)
summary(fg.psurv2.f)
AIC(fg.psurv2.f)

#### DENV3
## Base model
fg.data.3 <- finegray(Surv(time_event_adj, status) ~ .,
                      data=x1, na.action=na.omit, etype = "3")
fg.psurv3.a <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PC1 + PC2 ,
                     weight=fgwt, data=fg.data.3)
summary(fg.psurv3.a)
AIC(fg.psurv3.a)

fg.psurv3.b <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PC1 ,
                     weight=fgwt, data=fg.data.3)
summary(fg.psurv3.b)
AIC(fg.psurv3.b) # best-fitting shape; winning model; model 3c-I
cox.zph(fg.psurv3.b)

fg.psurv3.c <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PC1 + PC1_sq,
                     weight=fgwt, data=fg.data.3)
summary(fg.psurv3.c)
AIC(fg.psurv3.c)

## Covariate model
fg.psurv3.d <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PC1 + serostatus,
                     weight=fgwt, data=fg.data.3)
summary(fg.psurv3.d)
AIC(fg.psurv3.d) # not winning model

fg.psurv3.e <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PC1 * serostatus,
                     weight=fgwt, data=fg.data.3)
summary(fg.psurv3.e)
AIC(fg.psurv3.e)

#### DENV4
## Base model
fg.data.4 <- finegray(Surv(time_event_adj, status) ~ .,
                      data=x1, na.action=na.omit, etype = "4")
fg.psurv4.a <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PC1 + PC2,
                     weight=fgwt, data=fg.data.4)
summary(fg.psurv4.a)
AIC(fg.psurv4.a)

fg.psurv4.b <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PC1,
                     weight=fgwt, data=fg.data.4)
summary(fg.psurv4.b)
AIC(fg.psurv4.b)

fg.psurv4.c <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PC2,
                     weight=fgwt, data=fg.data.4)
summary(fg.psurv4.c)
AIC(fg.psurv4.c) # best-fitting shape

## Covariate model
fg.psurv4.d <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PC2 + serostatus,
                     weight=fgwt, data=fg.data.4)
summary(fg.psurv4.d)
AIC(fg.psurv4.d) # winning model; 4c-I
cox.zph(fg.psurv4.d)

fg.psurv4.e <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PC2 * serostatus,
                     weight=fgwt, data=fg.data.4)
summary(fg.psurv4.e)
AIC(fg.psurv4.e)

############################# average titer
#### DENV1
## Base model
fg.data.1 <- finegray(Surv(time_event_adj, status) ~ .,
                      data=x1, na.action=na.omit, etype = "1")
fg.psurv1.a <- coxph(Surv(fgstart, fgstop, fgstatus) ~ sc_avg_NT,
                     weight=fgwt, data=fg.data.1)
summary(fg.psurv1.a)
AIC(fg.psurv1.a)  # best-fitting shape

fg.psurv1.b <- coxph(Surv(fgstart, fgstop, fgstatus) ~ sc_avg_NT + sc_avg_NT_sq,
                     weight=fgwt, data=fg.data.1)
summary(fg.psurv1.b)
AIC(fg.psurv1.b)

## Covariate model
fg.psurv1.d <- coxph(Surv(fgstart, fgstop, fgstatus) ~  sc_avg_NT + serostatus,
                     weight=fgwt, data=fg.data.1)
summary(fg.psurv1.d)
AIC(fg.psurv1.d)

fg.psurv1.e <- coxph(Surv(fgstart, fgstop, fgstatus) ~ sc_avg_NT * serostatus,
                     weight=fgwt, data=fg.data.1)
summary(fg.psurv1.e)
AIC(fg.psurv1.e) # winning model; model 1d-I
cox.zph(fg.psurv1.e)

beta <- coef(fg.psurv1.e)
beta.seroneg <- beta[1]
beta.seropos <- beta[1] + beta[3]
cov <- vcov(fg.psurv1.e)
cov.seroneg <- cov[1]
cov.seropos <- cov[1,1]+cov[3,3]+2*cov[1,3]
ci <- 0.95
n <- nrow(Data)
beta.seroneg.CIhigh <- beta.seroneg + qt(ci + (1 - ci)/2, df = n - 1)*sqrt(cov.seroneg)
beta.seroneg.CIlow <- beta.seroneg - qt(ci + (1 - ci)/2, df = n - 1)*sqrt(cov.seroneg)
beta.seropos.CIhigh <- beta.seropos + qt(ci + (1 - ci)/2, df = n - 1)*sqrt(cov.seropos)
beta.seropos.CIlow <- beta.seropos - qt(ci + (1 - ci)/2, df = n - 1)*sqrt(cov.seropos)

exp(beta.seropos)
exp(beta.seropos.CIlow)
exp(beta.seropos.CIhigh)

#### DENV2
## Base model
fg.data.2 <- finegray(Surv(time_event_adj, status) ~ .,
                      data=x1, na.action=na.omit, etype = "2")
fg.psurv2.a <- coxph(Surv(fgstart, fgstop, fgstatus) ~ sc_avg_NT + sc_avg_NT_sq,
                     weight=fgwt, data=fg.data.2)
summary(fg.psurv2.a)
AIC(fg.psurv2.a) # best-fitting shape

fg.psurv2.b <- coxph(Surv(fgstart, fgstop, fgstatus) ~ sc_avg_NT ,
                     weight=fgwt, data=fg.data.2)
summary(fg.psurv2.b)
AIC(fg.psurv2.b)

## Covariate model
fg.psurv2.d <- coxph(Surv(fgstart, fgstop, fgstatus) ~ sc_avg_NT + sc_avg_NT_sq + serostatus,
                     weight=fgwt, data=fg.data.2)
summary(fg.psurv2.d)
AIC(fg.psurv2.d) # winning model; model 2d-I
cox.zph(fg.psurv2.d)

fg.psurv2.e <- coxph(Surv(fgstart, fgstop, fgstatus) ~ sc_avg_NT * serostatus+ sc_avg_NT_sq ,
                     weight=fgwt, data=fg.data.2)
summary(fg.psurv2.e)
AIC(fg.psurv2.e)

fg.psurv2.f <- coxph(Surv(fgstart, fgstop, fgstatus) ~ sc_avg_NT + sc_avg_NT_sq * serostatus ,
                     weight=fgwt, data=fg.data.2)
summary(fg.psurv2.f)
AIC(fg.psurv2.f)

#### DENV3
## Base model
fg.data.3 <- finegray(Surv(time_event_adj, status) ~ .,
                      data=x1, na.action=na.omit, etype = "3")
fg.psurv3.a <- coxph(Surv(fgstart, fgstop, fgstatus) ~ sc_avg_NT +sc_avg_NT_sq,
                     weight=fgwt, data=fg.data.3)
summary(fg.psurv3.a)
AIC(fg.psurv3.a)

fg.psurv3.b <- coxph(Surv(fgstart, fgstop, fgstatus) ~ sc_avg_NT ,
                     weight=fgwt, data=fg.data.3)
summary(fg.psurv3.b)
AIC(fg.psurv3.b) # best-fitting shape; winning model; model 3d-I
cox.zph(fg.psurv3.b)

## Covariate model
fg.psurv3.d <- coxph(Surv(fgstart, fgstop, fgstatus) ~ sc_avg_NT + serostatus,
                     weight=fgwt, data=fg.data.3)
summary(fg.psurv3.d)
AIC(fg.psurv3.d) # not winning model

fg.psurv3.e <- coxph(Surv(fgstart, fgstop, fgstatus) ~ sc_avg_NT * serostatus,
                     weight=fgwt, data=fg.data.3)
summary(fg.psurv3.e)
AIC(fg.psurv3.e)

#### DENV4
## Base model
fg.data.4 <- finegray(Surv(time_event_adj, status) ~ .,
                      data=x1, na.action=na.omit, etype = "4")
fg.psurv4.a <- coxph(Surv(fgstart, fgstop, fgstatus) ~ sc_avg_NT + sc_avg_NT_sq,
                     weight=fgwt, data=fg.data.4)
summary(fg.psurv4.a)
AIC(fg.psurv4.a)

fg.psurv4.b <- coxph(Surv(fgstart, fgstop, fgstatus) ~ sc_avg_NT,
                     weight=fgwt, data=fg.data.4)
summary(fg.psurv4.b)
AIC(fg.psurv4.b) # best-fitting shape; winning model; model 4d-I
cox.zph(fg.psurv4.b)

## Covariate model
fg.psurv4.d <- coxph(Surv(fgstart, fgstop, fgstatus) ~ sc_avg_NT + serostatus,
                     weight=fgwt, data=fg.data.4)
summary(fg.psurv4.d)
AIC(fg.psurv4.d) # not winning model

fg.psurv4.e <- coxph(Surv(fgstart, fgstop, fgstatus) ~ sc_avg_NT * serostatus,
                     weight=fgwt, data=fg.data.4)
summary(fg.psurv4.e)
AIC(fg.psurv4.e)

# Correlate of Protection Analysis ---------------------------------------------
# Best-fitting covariate models from CoR analysis
## 2. Cox regression with DENV-Any endpoint
############################# PCA 
coxFit.2.a.CoP <- coxph(Surv(time_event_adj, disease_any) ~ PC1 + PC1_sq + serostatus + vaccine, data = Data)
summary(coxFit.2.a.CoP)
cox.zph(coxFit.2.a.CoP)
AIC(coxFit.2.a.CoP) # vaccination status is insignificant; model 0a-II

############################# average titer
coxFit.2.A.CoP <- coxph(Surv(time_event_adj, disease_any) ~ sc_avg_NT + sc_avg_NT_sq + serostatus + vaccine, data = Data)
summary(coxFit.2.A.CoP)
cox.zph(coxFit.2.A.CoP)
AIC(coxFit.2.A.CoP) # vaccination status is insignificant; model 0b-II

## 3. CSH Cox model with DENV1-4 endpoints 
############################# PCA
#### DENV1
CSH.DENV1.2.b.CoP <- coxph(Surv(time_event_adj, status == 1) ~ PC1 * serostatus + vaccine, data = Data); 
summary(CSH.DENV1.2.b.CoP)
AIC(CSH.DENV1.2.b.CoP) # vaccine insignificant; model 1a-II
cox.zph(CSH.DENV1.2.b.CoP)

#### DENV2
CSH.DENV2.2.a.CoP <- coxph(Surv(time_event_adj, status == 2) ~ PC1 + PC1_sq + serostatus + vaccine, data = Data); 
summary(CSH.DENV2.2.a.CoP)
AIC(CSH.DENV2.2.a.CoP) # vaccine insignificant; model 2a-II
cox.zph(CSH.DENV2.2.a.CoP)

#### DENV3
CSH.DENV3.2.a.CoP <- coxph(Surv(time_event_adj, status == 3) ~ PC1 + vaccine , data = Data); 
summary(CSH.DENV3.2.a.CoP)
AIC(CSH.DENV3.2.a.CoP) # vaccine insignificant; model 3a-II
cox.zph(CSH.DENV3.2.a.CoP)

#### DENV4
CSH.DENV4.2.a.CoP <- coxph(Surv(time_event_adj, status == 4) ~ PC2 + serostatus + vaccine, data = Data); 
summary(CSH.DENV4.2.a.CoP)
AIC(CSH.DENV4.2.a.CoP) # vaccine marginally significant; model 4a-II
cox.zph(CSH.DENV4.2.a.CoP)

############################# average titer 
#### DENV1
CSH.DENV1.2.B.CoP <- coxph(Surv(time_event_adj, status == 1) ~ sc_avg_NT * serostatus + vaccine, data = Data); 
summary(CSH.DENV1.2.B.CoP)
AIC(CSH.DENV1.2.B.CoP) # vaccine insignificant; model 1b-II
cox.zph(CSH.DENV1.2.B.CoP) 

#### DENV2
CSH.DENV2.2.A.CoP <- coxph(Surv(time_event_adj, status == 2) ~ sc_avg_NT + sc_avg_NT_sq + serostatus + vaccine, data = Data); 
summary(CSH.DENV2.2.A.CoP)
AIC(CSH.DENV2.2.A.CoP) # vaccine insignificant; model 2b-II
cox.zph(CSH.DENV2.2.A.CoP)

#### DENV3
CSH.DENV3.2.A.CoP <- coxph(Surv(time_event_adj, status == 3) ~ sc_avg_NT + vaccine, data = Data); 
summary(CSH.DENV3.2.A.CoP)
AIC(CSH.DENV3.2.A.CoP) # vaccine insignificant; model 3b-II
cox.zph(CSH.DENV3.2.A.CoP)

#### DENV4
CSH.DENV4.2.A.CoP <- coxph(Surv(time_event_adj, status == 4) ~ sc_avg_NT + vaccine, data = Data); 
summary(CSH.DENV4.2.A.CoP)
AIC(CSH.DENV4.2.A.CoP) # vaccine marginally significant; model 4b-II
cox.zph(CSH.DENV4.2.A.CoP)

## 4. Fine-Gray model with DENV1-4 endpoints 
############################# PCA 
#### DENV1
fg.psurv1.e.CoP <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PC1 * serostatus + vaccine,
                         weight=fgwt, data=fg.data.1)
summary(fg.psurv1.e.CoP)
AIC(fg.psurv1.e.CoP) # vaccine is insignificant; model 1c-II
cox.zph(fg.psurv1.e.CoP)

#### DENV2
fg.psurv2.d.CoP <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PC1 + PC1_sq + serostatus + vaccine,
                         weight=fgwt, data=fg.data.2)
summary(fg.psurv2.d.CoP)
AIC(fg.psurv2.d.CoP) # vaccine is insignificant; model 2c-II
cox.zph(fg.psurv2.d.CoP)

#### DENV3
fg.psurv3.d.CoP <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PC1 + vaccine,
                         weight=fgwt, data=fg.data.3)
summary(fg.psurv3.d.CoP)
AIC(fg.psurv3.d.CoP) # vaccine is insignificant; model 3c-II
cox.zph(fg.psurv3.d.CoP)

#### DENV4
fg.psurv4.d.CoP <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PC2 + serostatus + vaccine,
                         weight=fgwt, data=fg.data.4)
summary(fg.psurv4.d.CoP)
AIC(fg.psurv4.d.CoP) # vaccine is significant; model 4c-II
cox.zph(fg.psurv4.d.CoP)

############################# average titer
#### DENV1
fg.psurv1.E.CoP <- coxph(Surv(fgstart, fgstop, fgstatus) ~ sc_avg_NT * serostatus + vaccine,
                         weight=fgwt, data=fg.data.1)
summary(fg.psurv1.E.CoP)
AIC(fg.psurv1.E.CoP) # vaccine is insignificant; model 1d-II
cox.zph(fg.psurv1.E.CoP)

#### DENV2
fg.psurv2.D.CoP <- coxph(Surv(fgstart, fgstop, fgstatus) ~ sc_avg_NT + sc_avg_NT_sq + serostatus + vaccine,
                         weight=fgwt, data=fg.data.2)
summary(fg.psurv2.D.CoP)
AIC(fg.psurv2.D.CoP) # vaccine is insignificant; model 2d-II
cox.zph(fg.psurv2.D.CoP)

#### DENV3
fg.psurv3.D.CoP <- coxph(Surv(fgstart, fgstop, fgstatus) ~ sc_avg_NT + vaccine,
                         weight=fgwt, data=fg.data.3)
summary(fg.psurv3.D.CoP)
AIC(fg.psurv3.D.CoP) # vaccine is insignificant; model 3d-II
cox.zph(fg.psurv3.D.CoP)

#### DENV4
fg.psurv4.D.CoP <- coxph(Surv(fgstart, fgstop, fgstatus) ~ sc_avg_NT + vaccine,
                         weight=fgwt, data=fg.data.4)
summary(fg.psurv4.D.CoP)
AIC(fg.psurv4.D.CoP) # vaccine is significant; model 4d-II
cox.zph(fg.psurv4.D.CoP)

# Vaccine Efficacy Prediction --------------------------------------------------
set.seed(1)
## Data 
Data.control <- Data %>% filter(vaccine == 0)
Data.vaccinated <- Data %>% filter(vaccine == 1)
Data.seropositive <- Data %>% filter(serostatus > 0)
Data.seronegative <- Data %>% filter(serostatus == 0)
Data.seropositive.vaccinated <- Data.vaccinated %>% filter(serostatus > 0)
Data.seropositive.control <- Data.control %>% filter(serostatus > 0)
Data.seronegative.vaccinated <- Data.vaccinated %>% filter(serostatus == 0)
Data.seronegative.control <- Data.control %>% filter(serostatus == 0)

## 2. Cox regression with DENV-Any endpoint 
############################# PCA 
coxFit.DENVAny.pca <- coxph(Surv(time_event_adj, disease_any) ~ PC1 + PC1_sq + serostatus, data = Data)
summary(coxFit.DENVAny.pca)
cox.zph(coxFit.DENVAny.pca)
AIC(coxFit.DENVAny.pca)

VE.coxFit.DENVAny.pca <- ve(coxFit.DENVAny.pca, Data, nboot = 1000)
VE.coxFit.DENVAny.pca.seronegative <- ve(coxFit.DENVAny.pca, Data.seronegative, nboot = 1000)
VE.coxFit.DENVAny.pca.seropositive <- ve(coxFit.DENVAny.pca, Data.seropositive, nboot = 1000)

############################# average titer 
coxFit.DENVAny.avg <- coxph(Surv(time_event_adj, disease_any) ~ sc_avg_NT + sc_avg_NT_sq + serostatus, data = Data)
summary(coxFit.DENVAny.avg)
cox.zph(coxFit.DENVAny.avg)
AIC(coxFit.DENVAny.avg)

VE.coxFit.DENVAny.avg <- ve(coxFit.DENVAny.avg, Data, nboot = 1000)
VE.coxFit.DENVAny.avg.seronegative <- ve(coxFit.DENVAny.avg, Data.seronegative, nboot = 1000)
VE.coxFit.DENVAny.avg.seropositive <- ve(coxFit.DENVAny.avg, Data.seropositive, nboot = 1000)

## 3. CSH Cox model with DENV1-4 endpoints 
############################# PCA
#### DENV1
CSH.DENV1.pca <- coxph(Surv(time_event_adj, status == 1) ~ PC1 * serostatus, data = Data); 
summary(CSH.DENV1.pca)
AIC(CSH.DENV1.pca)

VE.CSH.DENV1.pca <- ve(CSH.DENV1.pca, Data, nboot = 1000)
VE.CSH.DENV1.pca.seronegative <- ve(CSH.DENV1.pca, Data.seronegative, nboot = 1000)
VE.CSH.DENV1.pca.seropositive <- ve(CSH.DENV1.pca, Data.seropositive, nboot = 1000)

#### DENV2
CSH.DENV2.pca <- coxph(Surv(time_event_adj, status == 2) ~ PC1 + PC1_sq + serostatus, data = Data); 
summary(CSH.DENV2.pca)
AIC(CSH.DENV2.pca)

VE.CSH.DENV2.pca <- ve(CSH.DENV2.pca, Data, nboot = 1000)
VE.CSH.DENV2.pca.seronegative <- ve(CSH.DENV2.pca, Data.seronegative, nboot = 1000)
VE.CSH.DENV2.pca.seropositive <- ve(CSH.DENV2.pca, Data.seropositive, nboot = 1000)

#### DENV3
CSH.DENV3.pca <- coxph(Surv(time_event_adj, status == 3) ~ PC1, data = Data); 
summary(CSH.DENV3.pca)
AIC(CSH.DENV3.pca)

VE.CSH.DENV3.pca <- ve(CSH.DENV3.pca, Data, nboot = 1000)
VE.CSH.DENV3.pca.seronegative <- ve(CSH.DENV3.pca, Data.seronegative, nboot = 1000)
VE.CSH.DENV3.pca.seropositive <- ve(CSH.DENV3.pca, Data.seropositive, nboot = 1000)

#### DENV4
CSH.DENV4.pca <- coxph(Surv(time_event_adj, status == 4) ~ PC2 + serostatus + vaccine, data = Data); 
summary(CSH.DENV4.pca)
AIC(CSH.DENV4.pca)

VE.CSH.DENV4.pca.CoP <- ve(CSH.DENV4.pca, Data, nboot = 1000) 
VE.CSH.DENV4.pca.CoP.seronegative <- ve(CSH.DENV4.pca, Data.seronegative, nboot = 1000)
VE.CSH.DENV4.pca.CoP.seropositive <- ve(CSH.DENV4.pca, Data.seropositive, nboot = 1000)

############################# average titer 
#### DENV1
CSH.DENV1.avg <- coxph(Surv(time_event_adj, status == 1) ~ sc_avg_NT * serostatus, data = Data); 
summary(CSH.DENV1.avg)
AIC(CSH.DENV1.avg)

VE.CSH.DENV1.avg <- ve(CSH.DENV1.avg, Data, nboot = 1000)  
VE.CSH.DENV1.avg.seronegative <- ve(CSH.DENV1.avg, Data.seronegative, nboot = 1000) 
VE.CSH.DENV1.avg.seropositive <- ve(CSH.DENV1.avg, Data.seropositive, nboot = 1000) 

#### DENV2
CSH.DENV2.avg <- coxph(Surv(time_event_adj, status == 2) ~ sc_avg_NT + sc_avg_NT_sq + serostatus, data = Data); 
summary(CSH.DENV2.avg)
AIC(CSH.DENV2.avg)

VE.CSH.DENV2.avg <- ve(CSH.DENV2.avg, Data, nboot = 1000)  
VE.CSH.DENV2.avg.seronegative <- ve(CSH.DENV2.avg, Data.seronegative, nboot = 1000) 
VE.CSH.DENV2.avg.seropositive <- ve(CSH.DENV2.avg, Data.seropositive, nboot = 1000) 

#### DENV3
CSH.DENV3.avg <- coxph(Surv(time_event_adj, status == 3) ~ sc_avg_NT, data = Data); 
summary(CSH.DENV3.avg)
AIC(CSH.DENV3.avg)

VE.CSH.DENV3.avg <- ve(CSH.DENV3.avg, Data, nboot = 1000) 
VE.CSH.DENV3.avg.seronegative <- ve(CSH.DENV3.avg, Data.seronegative, nboot = 1000)
VE.CSH.DENV3.avg.seropositive <- ve(CSH.DENV3.avg, Data.seropositive, nboot = 1000)

#### DENV4
CSH.DENV4.avg <- coxph(Surv(time_event_adj, status == 4) ~ sc_avg_NT, data = Data); 
summary(CSH.DENV4.avg)
AIC(CSH.DENV4.avg)

VE.CSH.DENV4.avg <- ve(CSH.DENV4.avg, Data, nboot = 1000) 
VE.CSH.DENV4.avg.seronegative <- ve(CSH.DENV4.avg, Data.seronegative, nboot = 1000) 
VE.CSH.DENV4.avg.seropositive <- ve(CSH.DENV4.avg, Data.seropositive, nboot = 1000) 

## 4. Fine-Gray model with DENV1-4 endpoints
############################# PCA
#### DENV1
FG.DENV1.pca <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PC1 * serostatus,
                      weight=fgwt, data=fg.data.1)
summary(FG.DENV1.pca)
AIC(FG.DENV1.pca)

VE.FG.DENV1.pca <- ve(FG.DENV1.pca, Data, nboot = 1000) 
VE.FG.DENV1.pca.seronegative <- ve(FG.DENV1.pca, Data.seronegative, nboot = 1000) 
VE.FG.DENV1.pca.seropositive <- ve(FG.DENV1.pca, Data.seropositive, nboot = 1000) 

#### DENV2
FG.DENV2.pca <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PC1 + PC1_sq + serostatus,
                      weight=fgwt, data=fg.data.2)
summary(FG.DENV2.pca)
AIC(FG.DENV2.pca)

VE.FG.DENV2.pca <- ve(FG.DENV2.pca, Data, nboot = 1000) 
VE.FG.DENV2.pca.seronegative <- ve(FG.DENV2.pca, Data.seronegative, nboot = 1000) 
VE.FG.DENV2.pca.seropositive <- ve(FG.DENV2.pca, Data.seropositive, nboot = 1000) 

#### DENV3
FG.DENV3.pca <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PC1,
                      weight=fgwt, data=fg.data.3)
summary(FG.DENV3.pca)
AIC(FG.DENV3.pca)

VE.FG.DENV3.pca <- ve(FG.DENV3.pca, Data, nboot = 1000)
VE.FG.DENV3.pca.seronegative <- ve(FG.DENV3.pca, Data.seronegative, nboot = 1000)
VE.FG.DENV3.pca.seropositive <- ve(FG.DENV3.pca, Data.seropositive, nboot = 1000)

#### DENV4
FG.DENV4.pca <- coxph(Surv(fgstart, fgstop, fgstatus) ~ PC2 + serostatus + vaccine,
                      weight=fgwt, data=fg.data.4)
summary(FG.DENV4.pca)
AIC(FG.DENV4.pca)

VE.FG.DENV4.pca <- ve(FG.DENV4.pca, Data, nboot = 1000)
VE.FG.DENV4.pca.seronegative <- ve(FG.DENV4.pca, Data.seronegative, nboot = 1000)
VE.FG.DENV4.pca.seropositive <- ve(FG.DENV4.pca, Data.seropositive, nboot = 1000)

############################# average titer 
#### DENV1
FG.DENV1.avg <- coxph(Surv(fgstart, fgstop, fgstatus) ~ sc_avg_NT * serostatus,
                      weight=fgwt, data=fg.data.1)
summary(FG.DENV1.avg)
AIC(FG.DENV1.avg)

VE.FG.DENV1.avg <- ve(FG.DENV1.avg, Data, nboot = 1000)
VE.FG.DENV1.avg.seronegative <- ve(FG.DENV1.avg, Data.seronegative, nboot = 1000)
VE.FG.DENV1.avg.seropositive <- ve(FG.DENV1.avg, Data.seropositive, nboot = 1000)

#### DENV2
FG.DENV2.avg <- coxph(Surv(fgstart, fgstop, fgstatus) ~ sc_avg_NT + sc_avg_NT_sq + serostatus,
                      weight=fgwt, data=fg.data.2)
summary(FG.DENV2.avg)
AIC(FG.DENV2.avg)

VE.FG.DENV2.avg <- ve(FG.DENV2.avg, Data, nboot = 1000)
VE.FG.DENV2.avg.seronegative <- ve(FG.DENV2.avg, Data.seronegative, nboot = 1000)
VE.FG.DENV2.avg.seropositive <- ve(FG.DENV2.avg, Data.seropositive, nboot = 1000)

#### DENV3
FG.DENV3.avg <- coxph(Surv(fgstart, fgstop, fgstatus) ~ sc_avg_NT,
                      weight=fgwt, data=fg.data.3)
summary(FG.DENV3.avg)
AIC(FG.DENV3.avg)

VE.FG.DENV3.avg <- ve(FG.DENV3.avg, Data, nboot = 1000)
VE.FG.DENV3.avg.seronegative <- ve(FG.DENV3.avg, Data.seronegative, nboot = 1000)
VE.FG.DENV3.avg.seropositive <- ve(FG.DENV3.avg, Data.seropositive, nboot = 1000)

#### DENV4
FG.DENV4.avg <- coxph(Surv(fgstart, fgstop, fgstatus) ~ sc_avg_NT,
                      weight=fgwt, data=fg.data.4)
summary(FG.DENV4.avg)
AIC(FG.DENV4.avg)

VE.FG.DENV4.avg <- ve(FG.DENV4.avg, Data, nboot = 1000)
VE.FG.DENV4.avg.seronegative <- ve(FG.DENV4.avg, Data.seronegative, nboot = 1000)
VE.FG.DENV4.avg.seropositive <- ve(FG.DENV4.avg, Data.seropositive, nboot = 1000)

# Visualizations ---------------------------------------------------------------

# Risk curves for DENV-Any and DENV1-4 using average titer

# DENV-Any
set.seed(1)
## Cox, DENV-Any, average titer
summary(coxFit.DENVAny.avg)
sc_avg_nt <- seq(-3,3,0.001)
sc_avg_nt_sq <- sc_avg_nt*sc_avg_nt
serostatus <- rbinom(length(sc_avg_nt),1,1/2)
df <- data.frame(sc_avg_NT = sc_avg_nt,
                 sc_avg_NT_sq = sc_avg_nt_sq,
                 serostatus = serostatus)

Data$median_obs <- as.numeric(abs(median(Data$sc_avg_NT)-Data$sc_avg_NT)==min(abs(median(Data$sc_avg_NT)-Data$sc_avg_NT)))

coxFit.DENVAny.avg$means <- c(sc_avg_NT = Data$sc_avg_NT[Data$median_obs==1][1],
                              sc_avg_NT_sq = Data$sc_avg_NT_sq[Data$median_obs==1][1],
                              serostatus = Data$serostatus[Data$median_obs==1][1])

lp <- predict(coxFit.DENVAny.avg, df, type = "lp", se.fit = TRUE)
df$hr <- exp(lp$fit)
df$hr.lb <- exp(lp$fit - 1.96*lp$se.fit)
df$hr.ub <- exp(lp$fit + 1.96*lp$se.fit)

df <- df %>% 
  mutate(hrdiff = abs(hr-1)) %>%
  arrange(hrdiff)
coxFit.DENVAny.avg$means

df$hr.seropositive <- df$hr
df$hr.seropositive[df$serostatus==0] <- NA
df$hr.lb.seropositive <- df$hr.lb
df$hr.lb.seropositive[df$serostatus==0] <- NA
df$hr.ub.seropositive <- df$hr.ub
df$hr.ub.seropositive[df$serostatus==0] <- NA
df$hr.seronegative <- df$hr
df$hr.seronegative[df$serostatus==1] <- NA
df$hr.lb.seronegative <- df$hr.lb
df$hr.lb.seronegative[df$serostatus==1] <- NA
df$hr.ub.seronegative <- df$hr.ub
df$hr.ub.seronegative[df$serostatus==1] <- NA

line <- data.frame(vline = coxFit.DENVAny.avg$means[1])

risk.DENVany.seropositive <- ggplot(df, aes(x = sc_avg_NT)) +
  labs(x = 'Average SN titer', y = 'Relative hazard') +
  theme_bw() +
  geom_vline(data = line, 
             aes(xintercept = vline), linetype=1) +
  scale_y_continuous(limits = c(0, 13),
                     breaks = c(0,1,5,10),
                     minor_breaks = c(2.5,7.5,12.5),
                     labels = c('0', '1', '5', '10')) +
  geom_hline(aes(yintercept = 1), linetype=3) +
  geom_line(data=df[!is.na(df$hr.seropositive),], aes(y = hr.seropositive), color = "darkgray", size=1) +
  geom_ribbon(data = df[!is.na(df$hr.seropositive),],
              aes(ymin = hr.lb.seropositive, ymax = hr.ub.seropositive),
              alpha = 0.3,
              fill = "darkgray")+
  xlim(-3.2,2.2) +
  ggtitle("DENV-Any Seropositive") +
  theme(axis.text.x=element_blank(), 
        axis.text.y=element_blank()  
  )
immunogenicity.DENVany.seropositive <- ggplot(Data.seropositive, aes(x = sc_avg_NT, y = treatment)) +
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
  xlim(-3.2,2.2) +
  theme_bw() +
  theme(panel.background = element_blank()) + 
  xlab('Average SN titer') +
  ylab("")+
  theme(panel.background = element_blank()) 
grid.arrange(risk.DENVany.seropositive, immunogenicity.DENVany.seropositive, nrow = 2) 

risk.DENVany.seronegative <- ggplot(df, aes(x = sc_avg_NT)) +
  labs(x = 'Average SN titer', y = 'Relative hazard') +
  theme_bw() +
  geom_line(data=df[!is.na(df$hr.seronegative),], aes(y = hr.seronegative), color = "darkgray", size=1) +
  geom_ribbon(data = df[!is.na(df$hr.seronegative),],
              aes(ymin = hr.lb.seronegative, ymax = hr.ub.seronegative),
              alpha = 0.3,
              fill = "darkgray")+
  xlim(-3.2,2.2) +
  ylim(0,13) +
  ggtitle("DENV-Any Seronegative") +
  theme(panel.background = element_blank()) +
  theme(axis.text.x=element_blank(), 
        axis.text.y=element_blank() 
  )
immunogenicity.DENVany.seronegative <- ggplot(Data.seronegative, aes(x = sc_avg_NT, y = treatment)) +
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
  xlim(-3.2,2.2) +
  theme_bw() +
  theme(panel.background = element_blank()) + 
  xlab('Average SN titer') +
  ylab("")+
  theme(panel.background = element_blank()) 
grid.arrange(risk.DENVany.seronegative, immunogenicity.DENVany.seronegative, nrow = 2)

# DENV 1
set.seed(1)
## Cox CSH, DENV1, average titer
summary(CSH.DENV1.avg)
sc_avg_nt <- seq(-3,3,0.001)
serostatus <- rbinom(length(sc_avg_nt),1,1/2)
df <- data.frame(sc_avg_NT = sc_avg_nt,
                 serostatus = serostatus)

Data$sc_avg_NT_serostatus <- Data$sc_avg_NT*Data$serostatus

CSH.DENV1.avg$means <- c(sc_avg_NT = Data$sc_avg_NT[Data$median_obs==1][1],
                         serostatus = Data$serostatus[Data$median_obs==1][1],
                         sc_avg_NT_serostatus = Data$sc_avg_NT_serostatus[Data$median_obs==1][1])

lp <- predict(CSH.DENV1.avg, df, type = "lp", se.fit = TRUE)
df$hr <- exp(lp$fit)
df$hr.lb <- exp(lp$fit - 1.96*lp$se.fit)
df$hr.ub <- exp(lp$fit + 1.96*lp$se.fit)

df <- df %>% 
  mutate(hrdiff = abs(hr-1)) %>%
  arrange(hrdiff)
CSH.DENV1.avg$means

df$hr.seropositive <- df$hr
df$hr.seropositive[df$serostatus==0] <- NA
df$hr.lb.seropositive <- df$hr.lb
df$hr.lb.seropositive[df$serostatus==0] <- NA
df$hr.ub.seropositive <- df$hr.ub
df$hr.ub.seropositive[df$serostatus==0] <- NA
df$hr.seronegative <- df$hr
df$hr.seronegative[df$serostatus==1] <- NA
df$hr.lb.seronegative <- df$hr.lb
df$hr.lb.seronegative[df$serostatus==1] <- NA
df$hr.ub.seronegative <- df$hr.ub
df$hr.ub.seronegative[df$serostatus==1] <- NA

risk.DENV1.seropositive <- ggplot(df, aes(x = sc_avg_NT)) +
  labs(x = 'Average SN titer', y = 'Relative hazard') +
  theme_bw() +
  geom_vline(data = line, 
             aes(xintercept = vline), linetype=1) +
  geom_hline(aes(yintercept = 1), linetype=3) +
  geom_ribbon(data = df[!is.na(df$hr.seropositive),],
              aes(ymin = hr.lb.seropositive, ymax = hr.ub.seropositive),
              fill = "cadetblue3",
              alpha = 0.3)+
  geom_line(data=df[!is.na(df$hr.seropositive),], aes(y = hr.seropositive), 
            color = "cadetblue3", size=1) +
  ggtitle("DENV1 Seropositive") +
  theme(panel.background = element_blank())
immunogenicity.DENV1.seropositive <- ggplot(Data.seropositive, aes(x = sc_avg_NT, y = treatment)) +
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
  theme_bw() +
  theme(panel.background = element_blank()) + 
  xlab('Average SN titer') +
  ylab("")+
  theme(panel.background = element_blank())
grid.arrange(risk.DENV1.seropositive, immunogenicity.DENV1.seropositive, nrow = 2) 

risk.DENV1.seronegative <- ggplot(df, aes(x = sc_avg_NT)) +
  labs(x = 'Average SN titer', y = 'Relative hazard') +
  theme_bw() +
  geom_line(data=df[!is.na(df$hr.seronegative),], aes(y = hr.seronegative), color = "cadetblue3", size=1) +
  geom_ribbon(data = df[!is.na(df$hr.seronegative),],
              aes(ymin = hr.lb.seronegative, ymax = hr.ub.seronegative),
              alpha = 0.3,
              fill = "cadetblue3")+
  ggtitle("DENV1 Seronegative") +
  theme(panel.background = element_blank()) 
immunogenicity.DENV1.seronegative <- ggplot(Data.seronegative, aes(x = sc_avg_NT, y = treatment)) +
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
  theme_bw() +
  theme(panel.background = element_blank()) + 
  xlab('Average SN titer') +
  ylab("")+
  theme(panel.background = element_blank())
grid.arrange(risk.DENV1.seronegative, immunogenicity.DENV1.seronegative, nrow = 2)

# DENV 2
set.seed(1)
## Cox CSH, DENV2, average titer
summary(CSH.DENV2.avg)
sc_avg_nt <- seq(-3,3,0.001)
sc_avg_nt_sq <- sc_avg_nt*sc_avg_nt
serostatus <- rbinom(length(sc_avg_nt),1,1/2)
df <- data.frame(sc_avg_NT = sc_avg_nt,
                 sc_avg_NT_sq = sc_avg_nt_sq,
                 serostatus = serostatus)

CSH.DENV2.avg$means <- c(sc_avg_NT = Data$sc_avg_NT[Data$median_obs==1][1],
                         sc_avg_NT_sq = Data$sc_avg_NT_sq[Data$median_obs==1][1],
                         serostatus = Data$serostatus[Data$median_obs==1][1])

lp <- predict(CSH.DENV2.avg, df, type = "lp", se.fit = TRUE)
df$hr <- exp(lp$fit)
df$hr.lb <- exp(lp$fit - 1.96*lp$se.fit)
df$hr.ub <- exp(lp$fit + 1.96*lp$se.fit)

df <- df %>% 
  mutate(hrdiff = abs(hr-1)) %>%
  arrange(hrdiff)
CSH.DENV2.avg$means

df$hr.seropositive <- df$hr
df$hr.seropositive[df$serostatus==0] <- NA
df$hr.lb.seropositive <- df$hr.lb
df$hr.lb.seropositive[df$serostatus==0] <- NA
df$hr.ub.seropositive <- df$hr.ub
df$hr.ub.seropositive[df$serostatus==0] <- NA
df$hr.seronegative <- df$hr
df$hr.seronegative[df$serostatus==1] <- NA
df$hr.lb.seronegative <- df$hr.lb
df$hr.lb.seronegative[df$serostatus==1] <- NA
df$hr.ub.seronegative <- df$hr.ub
df$hr.ub.seronegative[df$serostatus==1] <- NA

risk.DENV2.seropositive <- ggplot(df, aes(x = sc_avg_NT)) +
  labs(x = 'Average SN titer', y = 'Relative hazard') +
  theme_bw() +
  geom_vline(data = line, 
             aes(xintercept = vline), linetype=1) +
  geom_hline(aes(yintercept = 1), linetype=3) +
  geom_line(data=df[!is.na(df$hr.seropositive),], aes(y = hr.seropositive), color = "cadetblue4", size=1) +
  geom_ribbon(data = df[!is.na(df$hr.seropositive),],
              aes(ymin = hr.lb.seropositive, ymax = hr.ub.seropositive),
              alpha = 0.3,
              fill = "cadetblue4")+
  ggtitle("DENV2 Seropositive") +
  theme(panel.background = element_blank()) 
immunogenicity.DENV2.seropositive <- ggplot(Data.seropositive, aes(x = sc_avg_NT, y = treatment)) +
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
  theme_bw() +
  theme(panel.background = element_blank()) + 
  xlab('Average SN titer') +
  ylab("")+
  theme(panel.background = element_blank()) 
grid.arrange(risk.DENV2.seropositive, immunogenicity.DENV2.seropositive, nrow = 2) 

risk.DENV2.seronegative <- ggplot(df, aes(x = sc_avg_NT)) +
  labs(x = 'Average SN titer', y = 'Relative hazard') +
  theme_bw() +
  geom_line(data=df[!is.na(df$hr.seronegative),], aes(y = hr.seronegative), color = "cadetblue4", size=1) +
  geom_ribbon(data = df[!is.na(df$hr.seronegative),],
              aes(ymin = hr.lb.seronegative, ymax = hr.ub.seronegative),
              alpha = 0.3,
              fill = "cadetblue4")+
  ggtitle("DENV2 Seronegative") +
  theme(panel.background = element_blank())
immunogenicity.DENV2.seronegative <- ggplot(Data.seronegative, aes(x = sc_avg_NT, y = treatment)) +
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
  theme_bw() +
  theme(panel.background = element_blank()) + 
  xlab('Average SN titer') +
  ylab("")+
  theme(panel.background = element_blank())
grid.arrange(risk.DENV2.seronegative, immunogenicity.DENV2.seronegative, nrow = 2)

# DENV 3
set.seed(1)
## Cox CSH, DENV3, average titer
summary(CSH.DENV3.avg)
sc_avg_nt <- seq(-3,3,0.001)
sc_avg_nt_sq <- sc_avg_nt*sc_avg_nt
serostatus <- rbinom(length(sc_avg_nt),1,1/2)
df <- data.frame(sc_avg_NT = sc_avg_nt,
                 sc_avg_NT_sq = sc_avg_nt_sq,
                 serostatus = serostatus)

CSH.DENV3.avg$means <- c(sc_avg_NT = Data$sc_avg_NT[Data$median_obs==1][1])

lp <- predict(CSH.DENV3.avg, df, type = "lp", se.fit = TRUE)
df$hr <- exp(lp$fit)
df$hr.lb <- exp(lp$fit - 1.96*lp$se.fit)
df$hr.ub <- exp(lp$fit + 1.96*lp$se.fit)

df <- df %>% 
  mutate(hrdiff = abs(hr-1)) %>%
  arrange(hrdiff)
CSH.DENV3.avg$means

df$hr.seropositive <- df$hr
df$hr.seropositive[df$serostatus==0] <- NA
df$hr.lb.seropositive <- df$hr.lb
df$hr.lb.seropositive[df$serostatus==0] <- NA
df$hr.ub.seropositive <- df$hr.ub
df$hr.ub.seropositive[df$serostatus==0] <- NA
df$hr.seronegative <- df$hr
df$hr.seronegative[df$serostatus==1] <- NA
df$hr.lb.seronegative <- df$hr.lb
df$hr.lb.seronegative[df$serostatus==1] <- NA
df$hr.ub.seronegative <- df$hr.ub
df$hr.ub.seronegative[df$serostatus==1] <- NA

risk.DENV3.seropositive <- ggplot(df, aes(x = sc_avg_NT)) +
  labs(x = 'Average SN titer', y = 'Relative hazard') +
  theme_bw() +
  geom_vline(data = line, 
             aes(xintercept = vline), linetype=1) +
  geom_hline(aes(yintercept = 1), linetype=3) +
  geom_line(data=df[!is.na(df$hr.seropositive),], aes(y = hr.seropositive), color = "cornflowerblue", size=1) +
  geom_ribbon(data = df[!is.na(df$hr.seropositive),],
              aes(ymin = hr.lb.seropositive, ymax = hr.ub.seropositive),
              alpha = 0.3,
              fill = "cornflowerblue")+
  ggtitle("DENV3 Seropositive") +
  theme(panel.background = element_blank())
immunogenicity.DENV3.seropositive <- ggplot(Data.seropositive, aes(x = sc_avg_NT, y = treatment)) +
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
  theme_bw() +
  theme(panel.background = element_blank()) + 
  xlab('Average SN titer') +
  ylab("")+
  theme(panel.background = element_blank())
grid.arrange(risk.DENV3.seropositive, immunogenicity.DENV3.seropositive, nrow = 2) 

risk.DENV3.seronegative <- ggplot(df, aes(x = sc_avg_NT)) +
  labs(x = 'Average SN titer', y = 'Relative hazard') +
  theme_bw() +
  geom_line(data=df[!is.na(df$hr.seronegative),], aes(y = hr.seronegative), color = "cornflowerblue", size=1) +
  geom_ribbon(data = df[!is.na(df$hr.seronegative),],
              aes(ymin = hr.lb.seronegative, ymax = hr.ub.seronegative),
              alpha = 0.3,
              fill = "cornflowerblue")+
  ggtitle("DENV3 Seronegative") +
  theme(panel.background = element_blank())
immunogenicity.DENV3.seronegative <- ggplot(Data.seronegative, aes(x = sc_avg_NT, y = treatment)) +
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
  theme_bw() +
  theme(panel.background = element_blank()) + 
  xlab('Average SN titer') +
  ylab("")+
  theme(panel.background = element_blank()) 
grid.arrange(risk.DENV3.seronegative, immunogenicity.DENV3.seronegative, nrow = 2)

# DENV 4
set.seed(1)
## Cox CSH, DENV4, average titer
summary(CSH.DENV4.avg)
sc_avg_nt <- seq(-3,3,0.001)
sc_avg_nt_sq <- sc_avg_nt*sc_avg_nt
serostatus <- rbinom(length(sc_avg_nt),1,1/2)
df <- data.frame(sc_avg_NT = sc_avg_nt,
                 sc_avg_NT_sq = sc_avg_nt_sq,
                 serostatus = serostatus)

CSH.DENV4.avg$means <- c(sc_avg_NT = Data$sc_avg_NT[Data$median_obs==1][1])

lp <- predict(CSH.DENV4.avg, df, type = "lp", se.fit = TRUE)
df$hr <- exp(lp$fit)
df$hr.lb <- exp(lp$fit - 1.96*lp$se.fit)
df$hr.ub <- exp(lp$fit + 1.96*lp$se.fit)

df <- df %>% 
  mutate(hrdiff = abs(hr-1)) %>%
  arrange(hrdiff)
CSH.DENV4.avg$means

df$hr.seropositive <- df$hr
df$hr.seropositive[df$serostatus==0] <- NA
df$hr.lb.seropositive <- df$hr.lb
df$hr.lb.seropositive[df$serostatus==0] <- NA
df$hr.ub.seropositive <- df$hr.ub
df$hr.ub.seropositive[df$serostatus==0] <- NA
df$hr.seronegative <- df$hr
df$hr.seronegative[df$serostatus==1] <- NA
df$hr.lb.seronegative <- df$hr.lb
df$hr.lb.seronegative[df$serostatus==1] <- NA
df$hr.ub.seronegative <- df$hr.ub
df$hr.ub.seronegative[df$serostatus==1] <- NA

risk.DENV4.seropositive <- ggplot(df, aes(x = sc_avg_NT)) +
  labs(x = 'Average SN titer', y = 'Relative hazard') +
  theme_bw() +
  geom_vline(data = line, 
             aes(xintercept = vline), linetype=1) +
  geom_hline(aes(yintercept = 1), linetype=3) +
  geom_line(data=df[!is.na(df$hr.seropositive),], aes(y = hr.seropositive), color = "deepskyblue4", size=1) +
  geom_ribbon(data = df[!is.na(df$hr.seropositive),],
              aes(ymin = hr.lb.seropositive, ymax = hr.ub.seropositive),
              alpha = 0.3,
              fill = "deepskyblue4")+
  ggtitle("DENV4 Seropositive") +
  theme(panel.background = element_blank()) 
immunogenicity.DENV4.seropositive <- ggplot(Data.seropositive, aes(x = sc_avg_NT, y = treatment)) +
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
  theme_bw() +
  theme(panel.background = element_blank()) + 
  xlab('Average SN titer') +
  ylab("")+
  theme(panel.background = element_blank()) 
grid.arrange(risk.DENV4.seropositive, immunogenicity.DENV4.seropositive, nrow = 2) 

risk.DENV4.seronegative <- ggplot(df, aes(x = sc_avg_NT)) +
  labs(x = 'Average SN titer', y = 'Relative hazard') +
  theme_bw() +
  geom_line(data=df[!is.na(df$hr.seronegative),], aes(y = hr.seronegative), color = "deepskyblue4", size=1) +
  geom_ribbon(data = df[!is.na(df$hr.seronegative),],
              aes(ymin = hr.lb.seronegative, ymax = hr.ub.seronegative),
              alpha = 0.3,
              fill = "deepskyblue4")+
  ggtitle("DENV4 Seronegative") +
  theme(panel.background = element_blank()) 
immunogenicity.DENV4.seronegative <- ggplot(Data.seronegative, aes(x = sc_avg_NT, y = treatment)) +
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
  theme_bw() +
  theme(panel.background = element_blank()) + 
  xlab('Average SN titer') +
  ylab("")+
  theme(panel.background = element_blank()) 
grid.arrange(risk.DENV4.seronegative, immunogenicity.DENV4.seronegative, nrow = 2)

# Tables -----------------------------------------------------------------------

## Table 2 
# DENV-Any Seropositive
sum(Data.seropositive.control$time_event)/365.25 # control person-years at risk
sum(Data.seropositive.control$disease_any) # control cases

sum(Data.seropositive.vaccinated$time_event)/365.25 # vaccinated person-years at risk
sum(Data.seropositive.vaccinated$disease_any) # vaccinated cases

lapply(VEwaldCIonDF(A = sum(Data.seropositive.vaccinated$disease_any),
                    B = sum(Data.seropositive.vaccinated$time_event)/365.25,
                    C = sum(Data.seropositive.control$disease_any),
                    D = sum(Data.seropositive.control$time_event)/365.25),"*",100) # case-count VE
VE.coxFit.DENVAny.avg.seropositive # immunogenicity-based VE

# DENV-Any Seronegative
sum(Data.seronegative.control$time_event)/365.25 # control person-years at risk
sum(Data.seronegative.control$disease_any) # control cases

sum(Data.seronegative.vaccinated$time_event)/365.25 # vaccinated person-years at risk
sum(Data.seronegative.vaccinated$disease_any) # vaccinated cases

lapply(VEwaldCIonDF(A = sum(Data.seronegative.vaccinated$disease_any),
                    B = sum(Data.seronegative.vaccinated$time_event)/365.25,
                    C = sum(Data.seronegative.control$disease_any),
                    D = sum(Data.seronegative.control$time_event)/365.25),"*",100) # case-count VE
VE.coxFit.DENVAny.avg.seronegative # immunogenicity-based VE

# DENV1 Seropositive
sum(Data.seropositive.control$serotype_disease==1, na.rm = TRUE) # control cases
sum(Data.seropositive.vaccinated$serotype_disease==1, na.rm = TRUE) # vaccinated cases
lapply(VEwaldCIonDF(A = sum(Data.seropositive.vaccinated$serotype_disease==1, na.rm = TRUE),
                    B = sum(Data.seropositive.vaccinated$time_event)/365.25,
                    C = sum(Data.seropositive.control$serotype_disease==1, na.rm = TRUE),
                    D = sum(Data.seropositive.control$time_event)/365.25),"*",100) # case-count VE
VE.CSH.DENV1.avg.seropositive # immunogenicity-based VE

# DENV1 Seronegative
sum(Data.seronegative.control$serotype_disease==1, na.rm = TRUE) # control cases
sum(Data.seronegative.vaccinated$serotype_disease==1, na.rm = TRUE) # vaccinated cases
lapply(VEwaldCIonDF(A = sum(Data.seronegative.vaccinated$serotype_disease==1, na.rm = TRUE),
                    B = sum(Data.seronegative.vaccinated$time_event)/365.25,
                    C = sum(Data.seronegative.control$serotype_disease==1, na.rm = TRUE),
                    D = sum(Data.seronegative.control$time_event)/365.25),"*",100) # case-count VE
VE.CSH.DENV1.avg.seronegative # immunogenicity-based VE

# DENV2 Seropositive
sum(Data.seropositive.control$serotype_disease==2, na.rm = TRUE) # control cases
sum(Data.seropositive.vaccinated$serotype_disease==2, na.rm = TRUE) # vaccinated cases
lapply(VEwaldCIonDF(A = sum(Data.seropositive.vaccinated$serotype_disease==2, na.rm = TRUE),
                    B = sum(Data.seropositive.vaccinated$time_event)/365.25,
                    C = sum(Data.seropositive.control$serotype_disease==2, na.rm = TRUE),
                    D = sum(Data.seropositive.control$time_event)/365.25),"*",100) # case-count VE
VE.CSH.DENV2.avg.seropositive # immunogenicity-based VE

# DENV2 Seronegative
sum(Data.seronegative.control$serotype_disease==2, na.rm = TRUE) # control cases
sum(Data.seronegative.vaccinated$serotype_disease==2, na.rm = TRUE) # vaccinated cases
lapply(VEwaldCIonDF(A = sum(Data.seronegative.vaccinated$serotype_disease==2, na.rm = TRUE),
                    B = sum(Data.seronegative.vaccinated$time_event)/365.25,
                    C = sum(Data.seronegative.control$serotype_disease==2, na.rm = TRUE),
                    D = sum(Data.seronegative.control$time_event)/365.25),"*",100) # case-count VE
VE.CSH.DENV2.avg.seronegative # immunogenicity-based VE

# DENV3 Seropositive
sum(Data.seropositive.control$serotype_disease==3, na.rm = TRUE) # control cases
sum(Data.seropositive.vaccinated$serotype_disease==3, na.rm = TRUE) # vaccinated cases
lapply(VEwaldCIonDF(A = sum(Data.seropositive.vaccinated$serotype_disease==3, na.rm = TRUE),
                    B = sum(Data.seropositive.vaccinated$time_event)/365.25,
                    C = sum(Data.seropositive.control$serotype_disease==3, na.rm = TRUE),
                    D = sum(Data.seropositive.control$time_event)/365.25),"*",100) # case-count VE
VE.CSH.DENV3.avg.seropositive # immunogenicity-based VE

# DENV3 Seronegative
sum(Data.seronegative.control$serotype_disease==3, na.rm = TRUE) # control cases
sum(Data.seronegative.vaccinated$serotype_disease==3, na.rm = TRUE) # vaccinated cases
lapply(VEwaldCIonDF(A = sum(Data.seronegative.vaccinated$serotype_disease==3, na.rm = TRUE),
                    B = sum(Data.seronegative.vaccinated$time_event)/365.25,
                    C = sum(Data.seronegative.control$serotype_disease==3, na.rm = TRUE),
                    D = sum(Data.seronegative.control$time_event)/365.25),"*",100) # case-count VE
VE.CSH.DENV3.avg.seronegative # immunogenicity-based VE

# DENV4 Seropositive
sum(Data.seropositive.control$serotype_disease==4, na.rm = TRUE) # control cases
sum(Data.seropositive.vaccinated$serotype_disease==4, na.rm = TRUE) # vaccinated cases
lapply(VEwaldCIonDF(A = sum(Data.seropositive.vaccinated$serotype_disease==4, na.rm = TRUE),
                    B = sum(Data.seropositive.vaccinated$time_event)/365.25,
                    C = sum(Data.seropositive.control$serotype_disease==4, na.rm = TRUE),
                    D = sum(Data.seropositive.control$time_event)/365.25),"*",100) # case-count VE
VE.CSH.DENV4.avg.seropositive # immunogenicity-based VE

# DENV4 Seronegative
sum(Data.seronegative.control$serotype_disease==4, na.rm = TRUE) # control cases
sum(Data.seronegative.vaccinated$serotype_disease==4, na.rm = TRUE) # vaccinated cases
lapply(VEwaldCIonDF(A = sum(Data.seronegative.vaccinated$serotype_disease==4, na.rm = TRUE),
                    B = sum(Data.seronegative.vaccinated$time_event)/365.25,
                    C = sum(Data.seronegative.control$serotype_disease==4, na.rm = TRUE),
                    D = sum(Data.seronegative.control$time_event)/365.25),"*",100) # case-count VE
VE.CSH.DENV4.avg.seronegative # immunogenicity-based VE

# Updated Table 2 --------------------------------------------------------------

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

## DENV-Any
denv.any.int <- coxph(Surv(time_event_adj, disease_any) ~ vaccine*serostatus, data = Data)
summary(denv.any.int)
cox.zph(denv.any.int)
AIC(denv.any.int)

# Estimating RR and VE within the seronegative group
VE.RR(model = denv.any.int, contrast = c(1,0,0))

# Estimating RR and VE within the seropositive group
VE.RR(model = denv.any.int, contrast = c(1,0,1))

## DENV1
denv1.int <- coxph(Surv(time_event_adj, status == 1) ~ vaccine*serostatus, data = Data)
summary(denv1.int)
cox.zph(denv1.int)
AIC(denv1.int)

# Estimating RR and VE within the seronegative group
VE.RR(model = denv1.int, contrast = c(1,0,0))

# Estimating RR and VE within the seropositive group
VE.RR(model = denv1.int, contrast = c(1,0,1))

## DENV2
denv2.int <- coxph(Surv(time_event_adj, status == 2) ~ vaccine*serostatus, data = Data)
summary(denv2.int)
cox.zph(denv2.int)
AIC(denv2.int)

# Estimating RR and VE within the seronegative group
VE.RR(model = denv2.int, contrast = c(1,0,0))

# Estimating RR and VE within the seropositive group
VE.RR(model = denv2.int, contrast = c(1,0,1))

## DENV3
denv3.int <- coxph(Surv(time_event_adj, status == 3) ~ vaccine*serostatus, data = Data)
summary(denv3.int)
cox.zph(denv3.int)
AIC(denv3.int)

# Estimating RR and VE within the seronegative group
VE.RR(model = denv3.int, contrast = c(1,0,0))

# Estimating RR and VE within the seropositive group
VE.RR(model = denv3.int, contrast = c(1,0,1))

## DENV4
denv4.int <- coxph(Surv(time_event_adj, status == 4) ~ vaccine*serostatus, data = Data)
summary(denv4.int)
cox.zph(denv4.int)
AIC(denv4.int)

# Estimating RR and VE within the seronegative group
VE.RR(model = denv4.int, contrast = c(1,0,0))

# Estimating RR and VE within the seropositive group
VE.RR(model = denv4.int, contrast = c(1,0,1))
