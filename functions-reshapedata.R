# Reference: Julie Dudášová, Zdeněk Valenta, Jeffrey R. Sachs (2024).
# "Improving precision of vaccine efficacy evaluation using immune correlate data in time-to-event models"
# This program is released under the GNU GPL3 license. Copyright © 2024 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.

ReshapeDataSPS <- function(IMM, EFF){
  ## IMM
  # patno...subject identifier
  # cmi_flag...lines that contain string "HZ Sample" in this column should be discarded
  # result...value of the assay
  # assaynam...specifies whether the assay is gpELISA, RCF, or ELISPOT
  # AGEN...age of subject at entry
  # TREATMNT...string placebo or vaccine group
  # visit...string that specifies when the sample was taken
  # b_val...titer value at "visit 1" - baseline titer
  # foldrise...foldrise in titers - result/b_val
  # GENDERN...gender
  # RACEN...race
  # exam_dt...date of visit
  # RDate...date of visit 1 (baseline, vaccination)
  # AGESTRTN...age group (1 for younger, 2 for older)
  # rel_dy...days post vaccination
  
  # remove invalid data
  IMM <- filter(IMM)
  
  # remove invalid data + lines that contain string "HZ Sample" + remove irrelevant columns
  IMM <- IMM %>% 
    filter(cmi_flag=="CMI Sample", valid == 1) %>% 
    select(patno, result, assaynam, AGEN, TREATMNT, visit, rel_dy, foldrise, b_val, GENDERN, RACEN, exam_dt, RDate)
  
  tbl_df(IMM)    
  
  ## EFF
  # patno...subject identifier
  # agen_04...age at which subject entered study SPS
  # mitt_04...1 if subject was included in SPS study
  # f_mitt_04...1 if subject developed HZ in SPS study
  # phn90_04...if subject had post herpetic neuralgia in SPS study
  # t_phn90_04...time till subject developed post herpetic neuralgia in SPS study
  # mitt_04_y1...1 if subject was included in the first year of SPS study
  # f_mitt_04_y1...1 if subject developed HZ the in first year of SPS study
  # t_hz_04_y1...number of days in year 1 when subject was healthy
  # 0405 identifies variables from STPS study
  # 013 identifies variables from LTPS study
  # DOBDATE...birth date

  # remove irrelevant columns
  EFF <- select(EFF, patno, f_mitt_04, t_hz_04) 
  
  tbl_df(EFF)
  
  # Merge information for IMM patients with EFF information
  tables <- merge(x = IMM, y = EFF, by = "patno") 
  
  # We use gpELISA assay (antibody titers) in our analysis
  # For Cox analysis we use constant biomarker
  data <- tables %>% 
    filter(assaynam == "gpELISA", visit == "week 6")
  
  data$result <- as.numeric(data$result)
  data <- data %>% mutate(age_flag = AGEN<70)
  data <- data %>% mutate(b_val_flag = b_val>400)
  data <- data %>% filter(!is.na(foldrise))
  
  # foldrise lower than 1 indicates a measurement error
  data$foldrise[data$foldrise < 1] <- 1
  data$result[data$foldrise < 1] <- data$b_val[data$foldrise < 1]
  
  tbl_df(data)
  
  data$exam_dt <- as.Date(data$exam_dt, format='%m/%d/%Y')
  data$RDate <- as.Date(data$RDate, format='%m/%d/%Y')

  return(data)
}

ReshapeData <- function(subjectData, symptomaticCases, titerData){
  ## subjectData
  # X...row identifier
  # newID...subject identifier
  # sero...serostatus pre vaccination (1 = seropositive, 0 = seronegative); seropositive subject has developed antibodies against at least one serotype
  # vaccine...vaccination status (1 = vaccinated, 0 = placebo)
 
  as_tibble(subjectData)    
  
  ## symptomaticCases
  # X...row identifier
  # newID...subject identifier
  # delayPD3...time to symptomatic infection (i.e. disease) since dose 3, in days
  # serotype...serotype of confirmed DENV infection (not clear what 999 means - maybe unidentified?)
  
  as_tibble(symptomaticCases)

  ## titerData
  # X...row identifier
  # newID...subject identifier
  # delayPD3...time to symptomatic infection (i.e. disease) since dose 3, in days
  # D1_NT...level of netralizing titer (PRNT50 assay) against DENV1
  # D2_NT...level of netralizing titer (PRNT50 assay) against DENV2
  # D3_NT...level of netralizing titer (PRNT50 assay) against DENV3
  # D4_NT...level of netralizing titer (PRNT50 assay) against DENV4
  
  as_tibble(titerData)

  # Select relevant columns and edit their names
  subjectData <- subjectData %>% select("newID","sero","vaccine")
  names(subjectData) [2] <- "serostatus" # use "serostatus" instead of "sero"
  
  symptomaticCases <- symptomaticCases %>% select("newID","delayPD3","serotype")
  names(symptomaticCases) [c(2,3)] <- c("time_disease","serotype_disease") # use "time_disease" instead of "delayPD3", "serotype_disease" instead of "serotype"
  symptomaticCases <- symptomaticCases %>% mutate("disease" = 1) # add a column for disease status
  
  titerData <- titerData %>% select("newID","delayPD3","D1_NT","D2_NT","D3_NT","D4_NT")
  names(titerData) [2] <- "time_NT" # use "time_NT" instead of "delayPD3"
  
  # Merge all information from original published data into one data.frame
  Data.raw <- titerData
  Data.raw <- merge(Data.raw, subjectData, by="newID", all=T)
  Data.raw <- merge(Data.raw, symptomaticCases, by="newID", all=T)
  
  ## Derived data using only one measurement of NT (peak titer, ~30 days post dose)
  # Add derived columns to the original published data
  Data <- Data.raw
  Data <- Data %>%
    group_by(newID) %>% 
    mutate("time_event" = max(time_NT)) # assuming that each subject was observed until her/his last NT measurement
  
  # For analyses DENV2-specific
  Data.DENV2 <- Data
  # All symptomatic infections with other than DENV2 serotype will be considered healthy
  Data.DENV2$disease[Data.DENV2$serotype_disease!=2] <- 0 # diseased with DENV1, DENV3, DENV4, and unidentified serotype are considered healthy (or non-diseased with DENV2)
  Data.DENV2$disease[is.na(Data.DENV2$disease)]<-0
  Data.DENV2$time_event[Data.DENV2$disease==1] <- Data.DENV2$time_disease[Data.DENV2$disease==1] # censoring of observations after the first symptomatic disease

  Data.DENV2 <- Data.DENV2 %>%
    group_by(newID) %>%
    filter(time_NT>0) %>% #  Remove subjects where disease occurred prior to dose 3
    mutate("day30PD3"=min(time_NT)==time_NT) %>% # for selection of NT measurement ~30 days post dose 3
    filter(day30PD3==TRUE) %>% # use only NT measurements ~30 days post dose 3 
    select(-one_of("day30PD3")) %>%
    filter(time_event>0) %>% # Remove subjects where disease occurred prior to dose 3
    filter(time_NT<300) # Remove subjects with ~30 days visit missing
  
  # Remove rows with other infections for subjects who got sick with DENV2
  c <- count(Data.DENV2,newID)
  out <- c$newID[c$n>1]
  for (i in 1:length(out)){
    if (sum(Data.DENV2$newID==out[i] & Data.DENV2$serotype_disease==2)>0){
      Data.DENV2 <- Data.DENV2 %>% filter(newID!=out[i]|serotype_disease==2)
    }
  }
  
  # Remove multiple rows for non-diseased (with DENV2) subjects
  Data.DENV2 <- distinct(Data.DENV2, newID, .keep_all = TRUE)
  
  Data.DENV2$titer <- Data.DENV2$D2_NT
  
  Data.DENV2$serostatus[Data.DENV2$serostatus>0] <- 1 # seropositive is defined as the subject with previous exposure to at least one DENV serotype
  Data.DENV2 <- Data.DENV2 %>% 
    mutate(treatment = 'vaccine') %>%
    mutate(sero = 'seropositive') 
  Data.DENV2$treatment[Data.DENV2$vaccine==0] <- 'placebo'
  Data.DENV2$sero[Data.DENV2$serostatus==0] <- 'seronegative'
  
  Data.durability <- NA
  
  # For analyses DENV1-specific
  Data.DENV1 <- Data
  Data.DENV1$disease[Data.DENV1$serotype_disease!=1] <- 0 # diseased with DENV1, DENV3, DENV4, and unidentified serotype are considered healthy (or non-diseased with DENV2)
  Data.DENV1$disease[is.na(Data.DENV1$disease)]<-0
  Data.DENV1$time_event[Data.DENV1$disease==1] <- Data.DENV1$time_disease[Data.DENV1$disease==1] # censoring of observations after the first symptomatic disease
  Data.DENV1 <- Data.DENV1 %>%
    group_by(newID) %>%
    filter(time_NT>0) %>% #  Remove subjects where disease occurred prior to dose 3
    mutate("day30PD3"=min(time_NT)==time_NT) %>% # for selection of NT measurement ~30 days post dose 3
    filter(day30PD3==TRUE) %>% # use only NT measurements ~30 days post dose 3 
    select(-one_of("day30PD3")) %>%
    filter(time_event>0) %>% # Remove subjects where disease occurred prior to dose 3
    filter(time_NT<300) # Remove subjects with ~30 days visit missing
  c <- count(Data.DENV1,newID)
  out <- c$newID[c$n>1]
  for (i in 1:length(out)){
    if (sum(Data.DENV1$newID==out[i] & Data.DENV1$serotype_disease==1)>0){
      Data.DENV1 <- Data.DENV1 %>% filter(newID!=out[i]|serotype_disease==1)
    }
  }
  Data.DENV1 <- distinct(Data.DENV1, newID, .keep_all = TRUE)
  Data.DENV1$titer <- Data.DENV1$D1_NT
  Data.DENV1$serostatus[Data.DENV1$serostatus>0] <- 1 # seropositive is defined as the subject with previous exposure to at least one DENV serotype
  Data.DENV1 <- Data.DENV1 %>% 
    mutate(treatment = 'vaccine') %>%
    mutate(sero = 'seropositive') 
  Data.DENV1$treatment[Data.DENV1$vaccine==0] <- 'placebo'
  Data.DENV1$sero[Data.DENV1$serostatus==0] <- 'seronegative'
  
  # For analyses DENV3-specific
  Data.DENV3 <- Data
  Data.DENV3$disease[Data.DENV3$serotype_disease!=3] <- 0 # diseased with DENV1, DENV3, DENV4, and unidentified serotype are considered healthy (or non-diseased with DENV2)
  Data.DENV3$disease[is.na(Data.DENV3$disease)]<-0
  Data.DENV3$time_event[Data.DENV3$disease==1] <- Data.DENV3$time_disease[Data.DENV3$disease==1] # censoring of observations after the first symptomatic disease
  Data.DENV3 <- Data.DENV3 %>%
    group_by(newID) %>%
    filter(time_NT>0) %>% #  Remove subjects where disease occurred prior to dose 3
    mutate("day30PD3"=min(time_NT)==time_NT) %>% # for selection of NT measurement ~30 days post dose 3
    filter(day30PD3==TRUE) %>% # use only NT measurements ~30 days post dose 3 
    select(-one_of("day30PD3")) %>%
    filter(time_event>0) %>% # Remove subjects where disease occurred prior to dose 3
    filter(time_NT<300) # Remove subjects with ~30 days visit missing
  c <- count(Data.DENV3,newID)
  out <- c$newID[c$n>1]
  for (i in 1:length(out)){
    if (sum(Data.DENV3$newID==out[i] & Data.DENV3$serotype_disease==3)>0){
      Data.DENV3 <- Data.DENV3 %>% filter(newID!=out[i]|serotype_disease==3)
    }
  }
  Data.DENV3 <- distinct(Data.DENV3, newID, .keep_all = TRUE)
  Data.DENV3$titer <- Data.DENV3$D3_NT
  Data.DENV3$serostatus[Data.DENV3$serostatus>0] <- 1 # seropositive is defined as the subject with previous exposure to at least one DENV serotype
  Data.DENV3 <- Data.DENV3 %>% 
    mutate(treatment = 'vaccine') %>%
    mutate(sero = 'seropositive') 
  Data.DENV3$treatment[Data.DENV3$vaccine==0] <- 'placebo'
  Data.DENV3$sero[Data.DENV3$serostatus==0] <- 'seronegative'
  
  # For analyses DENV4-specific
  Data.DENV4 <- Data
  Data.DENV4$disease[Data.DENV4$serotype_disease!=4] <- 0 # diseased with DENV1, DENV3, DENV4, and unidentified serotype are considered healthy (or non-diseased with DENV2)
  Data.DENV4$disease[is.na(Data.DENV4$disease)]<-0
  Data.DENV4$time_event[Data.DENV4$disease==1] <- Data.DENV4$time_disease[Data.DENV4$disease==1] # censoring of observations after the first symptomatic disease
  Data.DENV4 <- Data.DENV4 %>%
    group_by(newID) %>%
    filter(time_NT>0) %>% #  Remove subjects where disease occurred prior to dose 3
    mutate("day30PD3"=min(time_NT)==time_NT) %>% # for selection of NT measurement ~30 days post dose 3
    filter(day30PD3==TRUE) %>% # use only NT measurements ~30 days post dose 3 
    select(-one_of("day30PD3")) %>%
    filter(time_event>0) %>% # Remove subjects where disease occurred prior to dose 3
    filter(time_NT<300) # Remove subjects with ~30 days visit missing
  c <- count(Data.DENV4,newID)
  out <- c$newID[c$n>1]
  for (i in 1:length(out)){
    if (sum(Data.DENV4$newID==out[i] & Data.DENV4$serotype_disease==4)>0){
      Data.DENV4 <- Data.DENV4 %>% filter(newID!=out[i]|serotype_disease==4)
    }
  }
  Data.DENV4 <- distinct(Data.DENV4, newID, .keep_all = TRUE)
  Data.DENV4$titer <- Data.DENV4$D4_NT
  Data.DENV4$serostatus[Data.DENV4$serostatus>0] <- 1 # seropositive is defined as the subject with previous exposure to at least one DENV serotype
  Data.DENV4 <- Data.DENV4 %>% 
    mutate(treatment = 'vaccine') %>%
    mutate(sero = 'seropositive') 
  Data.DENV4$treatment[Data.DENV4$vaccine==0] <- 'placebo'
  Data.DENV4$sero[Data.DENV4$serostatus==0] <- 'seronegative'
  
  # For analyses with all serotypes
  Data.all <- Data
  Data.all$disease[is.na(Data.all$disease)]<-0
  Data.all$time_event[Data.all$disease==1] <- Data.all$time_disease[Data.all$disease==1] 
  Data.all <- Data.all %>%
    group_by(newID) %>%
    filter(time_NT>0) %>% #  Remove subjects where disease occurred prior to dose 3
    mutate("day30PD3"=min(time_NT)==time_NT) %>% # for selection of NT measurement ~30 days post dose 3
    filter(day30PD3==TRUE) %>% # use only NT measurements ~30 days post dose 3 
    select(-one_of("day30PD3")) %>%
    filter(time_event>0) %>% # Remove subjects where disease occurred prior to dose 3
    filter(time_NT<300) # Remove subjects with ~30 days visit missing
  Data.all$serostatus[Data.all$serostatus>0] <- 1 # seropositive is defined as the subject with previous exposure to at least one DENV serotype
  Data.all <- Data.all %>% 
    mutate(treatment = 'vaccine') %>%
    mutate(sero = 'seropositive') %>%
    mutate(out = serotype_disease) %>%
    mutate(disease_any = out>0)
  Data.all$treatment[Data.all$vaccine==0] <- 'placebo'
  Data.all$sero[Data.all$serostatus==0] <- 'seronegative'
  Data.all$out[is.na(Data.all$out)] <- 0
  Data.all$out[Data.all$out==999] <- 0 # diseased with unspecified serotype are considered healthy
  Data.all$disease_any[Data.all$disease_any==TRUE] <- 1
  Data.all$disease_any[is.na(Data.all$disease_any)] <- 0
  # censoring due to competing risks (observations are censored after the first VCD)
  Data.all <- Data.all %>% 
    arrange(newID,time_event) %>%
   distinct(newID, disease, .keep_all = TRUE)
  
  return(list(Data.raw = Data.raw,
              Data.DENV2 = Data.DENV2,
              Data.DENV1 = Data.DENV1,
              Data.DENV3 = Data.DENV3,
              Data.DENV4 = Data.DENV4,
              Data.all = Data.all))
}

VEwaldCIonDF <- function (A, B, C, D, confLevel = 0.95) 
{
  N. <- 1 - ((1 - confLevel)/2)
  z <- stats::qnorm(N., mean = 0, sd = 1)
  N_VACC <- B
  N_CONTROL <- D
  RR_point_t <- (A/N_VACC)/(C/N_CONTROL)
  ln_RR_point <- log(RR_point_t)
  ln_RR_point_se <- sqrt((1/A) - (1/N_VACC) + (1/C) - (1/N_CONTROL))
  RR_point_se <- exp(ln_RR_point_se)
  RR_lower <- 1 - exp(ln_RR_point + (z * ln_RR_point_se))
  RR_upper <- 1 - exp(ln_RR_point - (z * ln_RR_point_se))
  RR_point <- 1 - RR_point_t
  wald_ci_point <- RR_point
  wald_ci_lower <- RR_lower
  wald_ci_upper <- RR_upper
  return(list(VE = RR_point,
              lowerBound = wald_ci_lower, 
              upperBound = wald_ci_upper))
}

ConvertNormalToGeometric <- function(mean,
                                     sd,
                                     n,
                                     ci = 0.95){
  se <- sd/sqrt(n) 
  GMT.upperBound <- 2^(mean+(se * qt(ci + (1 - ci)/2, df = n - 1)))
  GMT.lowerBound <- 2^(mean-(se * qt(ci + (1 - ci)/2, df = n - 1)))
  GMT <- 2^mean
  
  return(list(GMT = GMT,
              GMT.upperBound = GMT.upperBound,
              GMT.lowerBound = GMT.lowerBound))
}