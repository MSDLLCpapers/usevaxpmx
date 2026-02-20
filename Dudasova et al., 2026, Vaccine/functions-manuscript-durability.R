# Reference: Julie Dudášová et al. (2026).
# "Annual immunogenicity measurements in randomized clinical trials help reveal heterogeneity in long-term vaccine efficacy: a pooled logistic regression modeling study"
# This program is released under the GNU GPL3 license. Copyright © 2026 Merck & Co., Inc., Rahway, NJ, USA and its affiliates. All rights reserved.

ReshapeDataDurabilitySPS <- function(IMM, EFF){
  # remove invalid data
  IMM <- filter(IMM)
  
  # remove invalid data + lines that contain string "HZ Sample" + remove irrelevant columns
  IMM <- IMM %>% 
    filter(cmi_flag=="CMI Sample", valid == 1) %>% 
    select(patno, result, assaynam, AGEN, TREATMNT, visit, foldrise, b_val, GENDERN, RACEN, rel_dy)
  tbl_df(IMM)
  
  # remove irrelevant columns
  EFF <- select(EFF, patno, f_mitt_04, t_hz_04, 
                f_mitt_04_y1, f_mitt_04_y2, f_mitt_04_y3, f_mitt_04_y4,
                t_hz_04_y1, t_hz_04_y2, t_hz_04_y3, t_hz_04_y4) 
  tbl_df(EFF)
  
  # Merge information for IMM patients with EFF information
  tables <- merge(x = IMM, y = EFF, by = "patno") 
  
  # We use gpELISA assay (antibody titers) in our analysis
  # For Cox analysis we use constant biomarker
  data <- tables %>% 
    filter(assaynam == "gpELISA", rel_dy>=0)
  data <- data %>%
    group_by(patno) %>%
    arrange(rel_dy, .by_group = TRUE) %>%
    mutate(time_1 = rel_dy) %>%
    mutate(time_2 = rel_dy) %>%
    mutate(cens = 0)
  # time_2
  data$time_2 <- lead(data$time_1,1)
  data$time_2[is.na(data$time_2)==TRUE] <- 0
  data$time_2[data$time_2==0] <- data$t_hz_04[data$time_2==0] # last time 2 for each subject is the time of censoring
  # zoster cases:
  data$cens[data$f_mitt_04==1 & data$time_1<=data$t_hz_04 & data$time_2>=data$t_hz_04] <- 1
  # noncases:
  data$time_2[data$f_mitt_04==0 & data$time_2>data$t_hz_04] <- data$t_hz_04[data$f_mitt_04==0 & data$time_2>data$t_hz_04]
  data <- data %>% 
    filter(time_2 > time_1) # if time of censoring (disease) happened before the ELISA measurement, remove the ELISA measurement record
  
  # check: number of events
  length(unique(data$patno[data$f_mitt_04==1]))
  sum(data$cens)
  data.durability <- data
  return(data.durability)
}

ReshapeDataDurability <- function(subjectData, symptomaticCases, titerData){
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
  
  # For analyses with all serotypes
  Data.all <- Data
  Data.all$disease[is.na(Data.all$disease)]<-0
  Data.all$time_event[Data.all$disease==1] <- Data.all$time_disease[Data.all$disease==1] # for disease cases time-to-event is time-to-disease
  Data.all <- Data.all %>%
    group_by(newID) %>%
    filter(time_disease>0 | is.na(time_disease)) #  Remove subjects where disease occurred prior to dose 3
  Data.all$serostatus[Data.all$serostatus>0] <- 1 # seropositive is defined as the subject with previous exposure to at least one DENV serotype
  Data.all <- Data.all %>% 
    mutate(treatment = 'vaccine') %>%
    mutate(sero = 'seropositive') %>%
    mutate(out = serotype_disease) %>%
    mutate(disease_any = out>0)
  Data.all$treatment[Data.all$vaccine==0] <- 'placebo'
  Data.all$sero[Data.all$serostatus==0] <- 'seronegative'
  Data.all$out[is.na(Data.all$out)] <- 0
  Data.all$disease_any[Data.all$disease_any==TRUE] <- 1
  Data.all$disease_any[is.na(Data.all$disease_any)] <- 0
  # censoring due to competing risks (observations are censored after the first VCD)
  Data.all <- Data.all %>% 
    arrange(newID,time_event) %>%
    group_by(newID) %>%
    mutate("first_disease"=min(time_disease)==time_disease) %>% # for selection of first disease 
    filter(first_disease==TRUE | is.na(first_disease)) %>% # use only first disease and non-diseased 
    arrange(time_NT, .by_group = TRUE) %>%
    mutate(time_1 = time_NT) %>%
    mutate(time_2 = time_NT) %>%
    mutate(cens = 0)
  # time_2
  Data.all$time_2 <- lead(Data.all$time_1,1)
  Data.all$time_2[Data.all$time_2<Data.all$time_1 | is.na(Data.all$time_2)] <- Data.all$time_1[Data.all$time_2<Data.all$time_1| is.na(Data.all$time_2)] # correct time_2 for the last observation for each subject 
  # cens: value of 1 in the cens column indicates if the time_event is within the <time_1;time_2) interval
  Data.all$cens[Data.all$time_1<=Data.all$time_disease & Data.all$time_2>Data.all$time_disease] <- 1
  Data.all <- Data.all %>% 
    mutate(last_obs = max(time_disease, time_NT))
  # for diseased time_2 equals to time_disease
  Data.all$time_2[Data.all$cens==1] <- Data.all$time_disease[Data.all$cens==1]
  # for diseased, data after the first infection are no longer used, i.e. event=first symptomatic infection
  Data.all <- group_by(Data.all, newID) %>%
    mutate(first2 = min(which(cens == 1 | row_number() == n()))) %>%
    filter(row_number() <= first2) %>%
    select(-first2) %>%
    mutate("day30PD3"=min(time_NT))
  # adding additional record for diseased with disease occurring post last titer measurement (e.g., newID=76)
  Data.all <- Data.all %>%
    group_by(newID) %>%
    mutate("lastNT"=max(time_NT)==time_NT) # for selection of the last NT measurement 
  
  subjID <- Data.all$newID[Data.all$lastNT==TRUE & Data.all$disease==1 & Data.all$cens==0]
  for (i in subjID){
    df1 <- Data.all[(Data.all$newID==i & Data.all$lastNT==TRUE),]
    df1.new <- df1
    df1.new4time_1 <- df1$time_2
    df1.new$time_2 <- df1$time_disease
    df1.new$cens <- 1
    Data.all <- rbind(Data.all,df1.new)
  }
  Data.all <- arrange(Data.all, newID, time_NT) %>%
    select(-first_disease, -last_obs, -lastNT)
  return(list(Data.raw = Data.raw,
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

caseCountVE <- function(Data){
  Data.control <- Data %>% filter(vaccine == 0)
  Data.vaccinated <- Data %>% filter(vaccine == 1)
  Data.seropositive <- Data %>% filter(serostatus > 0)
  Data.seronegative <- Data %>% filter(serostatus == 0)
  Data.seropositive.vaccinated <- Data.vaccinated %>% filter(serostatus > 0)
  Data.seropositive.control <- Data.control %>% filter(serostatus > 0)
  Data.seronegative.vaccinated <- Data.vaccinated %>% filter(serostatus == 0)
  Data.seronegative.control <- Data.control %>% filter(serostatus == 0)
  
  # case-count 
  # seropositive:
  case.Data.seropositive.control <- Data.seropositive.control %>% filter(disease_D2==1)
  noncase.Data.seropositive.control <- Data.seropositive.control %>% filter(disease_D2==0)
  case.Data.seropositive.vaccinated <- Data.seropositive.vaccinated %>% filter(disease_D2==1)
  noncase.Data.seropositive.vaccinated <- Data.seropositive.vaccinated %>% filter(disease_D2==0)
  seropositive.control.PY <- sum(noncase.Data.seropositive.control$person_years_D2) # control person-years at risk
  seropositive.control.cases <- sum(case.Data.seropositive.control$disease_D2) # control cases
  seropositive.vaccinated.PY <- sum(noncase.Data.seropositive.vaccinated$person_years_D2) # vaccinated person-years at risk
  seropositive.vaccinated.cases <- sum(case.Data.seropositive.vaccinated$disease_D2) # vaccinated cases  
  seropositive.caseCount <- lapply(VEwaldCIonDF(A = seropositive.vaccinated.cases,
                                                B = seropositive.vaccinated.PY,
                                                C = seropositive.control.cases,
                                                D = seropositive.control.PY),"*",100) # case-count VE
  # seronegative:
  case.Data.seronegative.control <- Data.seronegative.control %>% filter(disease_D2==1)
  noncase.Data.seronegative.control <- Data.seronegative.control %>% filter(disease_D2==0)
  case.Data.seronegative.vaccinated <- Data.seronegative.vaccinated %>% filter(disease_D2==1)
  noncase.Data.seronegative.vaccinated <- Data.seronegative.vaccinated %>% filter(disease_D2==0)
  seronegative.control.PY <- sum(noncase.Data.seronegative.control$person_years_D2) # control person-years at risk
  seronegative.control.cases <- sum(case.Data.seronegative.control$disease_D2) # control cases
  seronegative.vaccinated.PY <- sum(noncase.Data.seronegative.vaccinated$person_years_D2) # vaccinated person-years at risk
  seronegative.vaccinated.cases <- sum(case.Data.seronegative.vaccinated$disease_D2) # vaccinated cases  
  seronegative.caseCount <- lapply(VEwaldCIonDF(A = seronegative.vaccinated.cases,
                                                B = seronegative.vaccinated.PY,
                                                C = seronegative.control.cases,
                                                D = seronegative.control.PY),"*",100) # case-count VE
  return(list(seronegative.VE = seronegative.caseCount, 
              seropositive.VE = seropositive.caseCount,
              seronegative.control.PY = seronegative.control.PY, 
              seronegative.control.cases = seronegative.control.cases, 
              seronegative.vaccinated.PY = seronegative.vaccinated.PY, 
              seronegative.vaccinated.cases = seronegative.vaccinated.cases,
              seropositive.control.PY = seropositive.control.PY, 
              seropositive.control.cases = seropositive.control.cases, 
              seropositive.vaccinated.PY = seropositive.vaccinated.PY, 
              seropositive.vaccinated.cases = seropositive.vaccinated.cases))
}

caseCountVEZostavax <- function(Data){
  Data.control <- Data %>% filter(vaccine == 0)
  Data.vaccinated <- Data %>% filter(vaccine == 1)
  Data.younger <- Data %>% filter(AGEN<70)
  Data.older <- Data %>% filter(AGEN>=70)
  Data.younger.vaccinated <- Data.vaccinated %>% filter(AGEN<70)
  Data.younger.control <- Data.control %>% filter(AGEN<70)
  Data.older.vaccinated <- Data.vaccinated %>% filter(AGEN>=70)
  Data.older.control <- Data.control %>% filter(AGEN>=70)
  
  # case-count 
  # younger:
  case.Data.younger.control <- Data.younger.control %>% filter(cens==1)
  noncase.Data.younger.control <- Data.younger.control %>% filter(cens==0)
  case.Data.younger.vaccinated <- Data.younger.vaccinated %>% filter(cens==1)
  noncase.Data.younger.vaccinated <- Data.younger.vaccinated %>% filter(cens==0)
  younger.control.PY <- sum(noncase.Data.younger.control$person_years) # control person-years at risk
  younger.control.cases <- sum(case.Data.younger.control$cens) # control cases
  younger.vaccinated.PY <- sum(noncase.Data.younger.vaccinated$person_years) # vaccinated person-years at risk
  younger.vaccinated.cases <- sum(case.Data.younger.vaccinated$cens) # vaccinated cases  
  younger.caseCount <- lapply(VEwaldCIonDF(A = younger.vaccinated.cases,
                                           B = younger.vaccinated.PY,
                                           C = younger.control.cases,
                                           D = younger.control.PY),"*",100) # case-count VE
  # older:
  case.Data.older.control <- Data.older.control %>% filter(cens==1)
  noncase.Data.older.control <- Data.older.control %>% filter(cens==0)
  case.Data.older.vaccinated <- Data.older.vaccinated %>% filter(cens==1)
  noncase.Data.older.vaccinated <- Data.older.vaccinated %>% filter(cens==0)
  older.control.PY <- sum(noncase.Data.older.control$person_years) # control person-years at risk
  older.control.cases <- sum(case.Data.older.control$cens) # control cases
  older.vaccinated.PY <- sum(noncase.Data.older.vaccinated$person_years) # vaccinated person-years at risk
  older.vaccinated.cases <- sum(case.Data.older.vaccinated$cens) # vaccinated cases  
  older.caseCount <- lapply(VEwaldCIonDF(A = older.vaccinated.cases,
                                         B = older.vaccinated.PY,
                                         C = older.control.cases,
                                         D = older.control.PY),"*",100) # case-count VE
  return(list(older.VE = older.caseCount, 
              younger.VE = younger.caseCount,
              older.control.PY = older.control.PY, 
              older.control.cases = older.control.cases, 
              older.vaccinated.PY = older.vaccinated.PY, 
              older.vaccinated.cases = older.vaccinated.cases,
              younger.control.PY = younger.control.PY, 
              younger.control.cases = younger.control.cases, 
              younger.vaccinated.PY = younger.vaccinated.PY, 
              younger.vaccinated.cases = younger.vaccinated.cases))
}

OverlayObservedProportions <- function(plot_obj, data, predictor = "D2_NT", N = 10) {
  data$predictor <- data[[predictor]]
  
  # Create intervals based on the log-transformed predictor
  data$interval <- cut(data$predictor, breaks = N, include.lowest = TRUE)
  
  # Calculate observed proportions within each interval
  observed_data <- data %>%
    group_by(interval) %>%
    reframe(
      mid_point = {
        bounds <- str_extract_all(as.character(interval), "-?\\d+\\.?\\d*")[[1]]
        if (length(bounds) == 2) {
          ((as.numeric(bounds[1]) + as.numeric(bounds[2])) / 2)
        } else {
          NA  # If bounds are not as expected, set midpoint to NA
        }
      },
      # observed_prop = sum(disease_D2, na.rm = TRUE) / sum(!is.na(disease_D2)),
      weighted_sum = sum(person_years_D2[disease_D2 == 0], na.rm = TRUE),
      weighted_prop = sum(disease_D2, na.rm = TRUE) / weighted_sum,
      # count = n(),
      count = round(weighted_sum),
      successes = sum(disease_D2, na.rm = TRUE),
      prop_ci_lower = if (count > 0) binom.test(successes, count)$conf.int[1] else NA,
      prop_ci_upper = if (count > 0) binom.test(successes, count)$conf.int[2] else NA
    ) %>%
    mutate(mid_point = mid_point) %>%
    filter(count > 0)
  
  # Add observed proportions to the plot as points
  plot_obj + 
    #geom_point(data = observed_data, aes(x = mid_point, y = observed_prop), color = "blue", size = 2) +
    geom_point(data = observed_data, aes(x = mid_point, y = weighted_prop), color = "black", size = 2) +
    geom_errorbar(data = observed_data, aes(x = mid_point, ymin = prop_ci_lower, ymax = prop_ci_upper), width = 0.1, color = "black", inherit.aes = FALSE) +
    # geom_text(data = observed_data, aes(x = mid_point, y = observed_prop, label = count), 
    #           vjust = -0.7, color = "blue", size = 3) +  # Adjust position and size as needed
    labs(subtitle = paste("Observed proportions (", N, " intervals)", sep = ""))
  
}

OverlayObservedProportionsZostavax <- function(plot_obj, data, predictor = "log2result", N = 10) {
  data$predictor <- data[[predictor]]
  
  # Create intervals based on the log-transformed predictor
  data$interval <- cut(data$predictor, breaks = N, include.lowest = TRUE)
  
  # Calculate observed proportions within each interval
  observed_data <- data %>%
    group_by(interval, vaccine) %>%
    reframe(
      mid_point = {
        bounds <- str_extract_all(as.character(interval), "-?\\d+\\.?\\d*")[[1]]
        if (length(bounds) == 2) {
          ((as.numeric(bounds[1]) + as.numeric(bounds[2])) / 2)
        } else {
          NA  # If bounds are not as expected, set midpoint to NA
        }
      },
      # observed_prop = sum(disease_D2, na.rm = TRUE) / sum(!is.na(disease_D2)),
      weighted_sum = sum(person_years[cens == 0], na.rm = TRUE),
      weighted_prop = sum(cens, na.rm = TRUE) / weighted_sum,
      # count = n(),
      count = round(weighted_sum),
      successes = sum(cens, na.rm = TRUE),
      prop_ci_lower = if (count > 0) binom.test(successes, count)$conf.int[1] else NA,
      prop_ci_upper = if (count > 0) binom.test(successes, count)$conf.int[2] else NA
    ) %>%
    mutate(mid_point = mid_point) %>%
    filter(count > 0)
  
  # Add observed proportions to the plot as points
  plot_obj + 
    geom_point(data = observed_data[observed_data$vaccine == 1,], 
               aes(x = mid_point + 0.05, y = weighted_prop), 
               color = "gray20", size = 2) +
    geom_errorbar(data = observed_data[observed_data$vaccine == 1,], 
                  aes(x = mid_point + 0.05, ymin = prop_ci_lower, ymax = prop_ci_upper), 
                  width = 0.1, color = "gray20", inherit.aes = FALSE) +
    geom_point(data = observed_data[observed_data$vaccine == 0,], 
               aes(x = mid_point - 0.05, y = weighted_prop), 
               color = "darkslategray4", size = 2) +
    geom_errorbar(data = observed_data[observed_data$vaccine == 0,], 
                  aes(x = mid_point - 0.05, ymin = prop_ci_lower, ymax = prop_ci_upper), 
                  width = 0.1, color = "darkslategray4", inherit.aes = FALSE) +
    labs(subtitle = paste("Observed proportions by group (", N, " intervals)", sep = ""))
}

robust_summary <- function(model, 
                           cluster = "newID") {
  
  # Ensure the cluster variable is passed correctly
  if (missing(cluster)) {
    stop("Please provide a cluster variable.")
  }
  
  # Calculate clustered variance-covariance matrix
  vcov_clustered <- vcovCL(model, cluster = as.formula(paste("~", cluster)))
  
  
  # Get coefficients
  coef_est <- coef(model)
  
  # Standard errors (square root of diagonal)
  robust_se <- sqrt(diag(vcov_clustered))
  
  # z-values and p-values
  z_vals <- coef_est / robust_se
  p_vals <- 2 * (1 - pnorm(abs(z_vals)))
  
  # Combine in a table
  robust_results <- data.frame(
    Estimate = coef_est,
    "Robust SE" = robust_se,
    "z value" = z_vals,
    "Pr(>|z|)" = p_vals
  )
  
  return(robust_results)
}