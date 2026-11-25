source("06_fever_analysis.R")
#==========================TREATMENT OUTCOME ANALYSIS===========================

#----------------Microscopy availability during early follow-up-----------------

early_microscopy_visits <- c(
  "H0", "H8", "H24", "H36",
  "H48", "H60", "H72"
)

SMART_analysis |>
  filter(
    study_visit %in% early_microscopy_visits
  ) |>
  group_by(
    treatment_arm,
    study_visit
  ) |>
  summarise(
    n_records = n(),
    n_with_microscopy = sum(!is.na(microscopy_sample_id)),
    n_with_parasite_density = sum(!is.na(parasite_density)),
    n_positive = sum(parasite_density > 0, na.rm = TRUE),
    n_zero = sum(parasite_density == 0, na.rm = TRUE),
    n_missing = sum(is.na(parasite_density)),
    .groups = "drop"
  ) |>
  arrange(
    treatment_arm,
    match(study_visit, early_microscopy_visits)
  ) |>
  print(n = Inf)

#-------By site
SMART_analysis |>
  filter(
    study_visit %in% early_microscopy_visits
  ) |>
  group_by(
    site,
    study_visit
  ) |>
  summarise(
    n_records = n(),
    n_with_parasite_density = sum(!is.na(parasite_density)),
    n_positive = sum(parasite_density > 0, na.rm = TRUE),
    n_zero = sum(parasite_density == 0, na.rm = TRUE),
    n_missing = sum(is.na(parasite_density)),
    .groups = "drop"
  ) |>
  arrange(
    site,
    match(study_visit, early_microscopy_visits)
  ) |>
  print(n = Inf)

#========Create data for only participants with parasite density >0 for H0======
eligible_participants <- smart_analysis_final |>
  filter(
    study_visit == "H0",
    !is.na(parasite_density),
    parasite_density > 0
  ) |>
  distinct(
    participant_id,
    treatment_arm,
    site
  )

# Number eligible
eligible_participants |>
  count(treatment_arm)

#----Check microscopy available for eligible participants-----------------------
smart_analysis_final |>
  semi_join(
    eligible_participants,
    by = "participant_id"
  ) |>
  filter(
    study_visit %in% c(
      "H0", "H8", "H24", "H36",
      "H48", "H60", "H72",
      "D7", "D14", "D21", "D28",
      "D35", "D42"
    )
  ) |>
  group_by(
    treatment_arm,
    study_visit
  ) |>
  summarise(
    n_participants = n_distinct(participant_id),
    n_density_available = sum(!is.na(parasite_density)),
    n_positive = sum(parasite_density > 0, na.rm = TRUE),
    n_zero = sum(parasite_density == 0, na.rm = TRUE),
    n_missing = sum(is.na(parasite_density)),
    .groups = "drop"
  ) |>
  print(n = Inf)

#===================Early Treatment Failure(ETF)================================
etf_data <- smart_analysis_final |> 
  semi_join(
    eligible_participants,
    by = "participant_id"
  ) |> 
  filter(
    study_visit %in% c("H0", "H48", "H72")
  ) |> 
  select(
    participant_id,
    treatment_arm,
    site,
    study_visit,
    parasite_density,
    temperature
  ) |> 
  pivot_wider(
    names_from = study_visit,
    values_from = c(
      parasite_density,
      temperature
    ),
    names_glue = "{.value}_{study_visit}"
  )

#check
etf_data |> 
  summarise(
    n = n(),
    missing_h0 = sum(is.na(parasite_density_H0)),
    missing_h48 = sum(is.na(parasite_density_H48)),
    missing_h72 = sum(is.na(parasite_density_H72)),
    missing_temp_h72 = sum(is.na(temperature_H72))
  )

#--------------Derive the 3 ETF criteria-----------
etf_data <- etf_data |>
  mutate(
    
    # Criterion 1:
    # Day 2 parasitaemia greater than Day 0
    etf_day2_higher = case_when(
      !is.na(parasite_density_H48) &
        parasite_density_H48 > parasite_density_H0 ~ TRUE,
      
      !is.na(parasite_density_H48) ~ FALSE,
      
      TRUE ~ NA
    ),
    
    # Criterion 2:
    # Day 3 parasitaemia with temperature >= 37.5 C
    etf_day3_fever = case_when(
      !is.na(parasite_density_H72) &
        parasite_density_H72 > 0 &
        !is.na(temperature_H72) &
        temperature_H72 >= 37.5 ~ TRUE,
      
      !is.na(parasite_density_H72) &
        !is.na(temperature_H72) ~ FALSE,
      
      TRUE ~ NA
    ),
    
    # Criterion 3:
    # Day 3 parasitaemia >= 25% of Day 0 count
    etf_day3_25pct = case_when(
      !is.na(parasite_density_H72) &
        parasite_density_H72 >=
        0.25 * parasite_density_H0 ~ TRUE,
      
      !is.na(parasite_density_H72) ~ FALSE,
      
      TRUE ~ NA
    )
  )

#------------count them
etf_data |>
  summarise(
    day2_higher = sum(etf_day2_higher %in% TRUE),
    day3_fever = sum(etf_day3_fever %in% TRUE),
    day3_25pct = sum(etf_day3_25pct %in% TRUE)
  )

etf_data |>
  group_by(treatment_arm) |>
  summarise(
    n = n(),
    day2_higher = sum(etf_day2_higher %in% TRUE),
    day3_fever = sum(etf_day3_fever %in% TRUE),
    day3_25pct = sum(etf_day3_25pct %in% TRUE),
    .groups = "drop"
  )

#---------------check danger sign or severe malaria----------
treatment_failure_check |>
  semi_join(
    eligible_participants,
    by = "participant_id"
  ) |>
  count(
    #tf_danger_sign,
    tf_severe_manifestation,
    #tf_malaria_diagnosis,
    .drop = FALSE
  )

#-------check if there was treatment failure on day 1 -3 ----------------------
treatment_failure_check |>
  semi_join(
    eligible_participants,
    by = "participant_id"
  ) |>
  filter(
    tf_visit %in% c(
      "H24", "H48", "H72"
    )
  ) |>
  select(
    participant_id,
    tf_visit,
    tf_microscopy_result,
    tf_danger_sign,
    tf_severe_manifestation,
    tf_malaria_diagnosis
  ) |>
  print(n = Inf)

#---------------------Final ETF classification----------------------------------

etf_data <- etf_data |>
  mutate(
    etf_status = case_when(
      
      # Meets at least one ETF criterion
      etf_day2_higher == TRUE |
        etf_day3_fever == TRUE |
        etf_day3_25pct == TRUE ~ "ETF",
      
      # All criteria assessed and none met
      etf_day2_higher == FALSE &
        etf_day3_fever == FALSE &
        etf_day3_25pct == FALSE ~ "No ETF",
      
      # At least one required measurement unavailable
      TRUE ~ "Indeterminate"
    )
  )

etf_data |>
  count(
    treatment_arm,
    etf_status
  )

etf_data |>
  count(etf_status)


#------------Participants with ETF indeterminate--------------------------------
etf_data |> 
  filter(etf_status == "Indeterminate") |> 
  select(
    participant_id,
    starts_with("parasite_density"),
    starts_with("temperature"),
    etf_status
  )

#=========================LATE CLINICAL FAILURE (LCF)===========================

no_etf_participants <- etf_data |>
  filter(etf_status == "No ETF") |>
  select(
    participant_id,
    treatment_arm,
    site
  )

#---------------------Create Day 4-42 lcf dataset----------------------------------

lcf_data <- smart_analysis_final |>
  semi_join(
    no_etf_participants,
    by = "participant_id"
  ) |>
  mutate(
    analysis_day = as.numeric(
      str_extract(study_day, "\\d+")
    )
  ) |>
  filter(
    analysis_day >= 4,
    analysis_day <= 42
  )

#---------------------Presence of parasitaemia----------------------------------

lcf_data <- lcf_data |>
  mutate(
    parasitaemia_present = case_when(
      
      # Parasite count confirms parasitaemia
      !is.na(parasite_density) &
        parasite_density > 0 ~ TRUE,
      
      # Microscopy confirms parasitaemia even if count unavailable
      microscopy_positive == 1 ~ TRUE,
      
      # Parasite count confirms absence
      !is.na(parasite_density) &
        parasite_density == 0 ~ FALSE,
      
      # Microscopy confirms absence
      microscopy_positive == 0 ~ FALSE,
      
      TRUE ~ NA
    )
  )

#----------------LCF fever + parasitaemia criteria
lcf_data <- lcf_data |>
  mutate(
    lcf_fever = case_when(
      
      parasitaemia_present == TRUE &
        !is.na(temperature) &
        temperature >= 37.5 ~ TRUE,
      
      !is.na(parasitaemia_present) &
        !is.na(temperature) ~ FALSE,
      
      TRUE ~ NA
    )
  )

#------------------Check severe malaria in treatment-failure visits-------------
#Check actual day of severe treatment failures

severe_tf_lcf <- treatment_failure_check |>
  semi_join(
    no_etf_participants,
    by = "participant_id"
  ) |>
  filter(
    tf_severe_manifestation == 1,
    tf_microscopy_result == 1
  ) |>
  select(
    participant_id,
    tf_visit,
    tf_visit_date,
    tf_microscopy_result,
    tf_danger_sign,
    tf_severe_manifestation,
    tf_malaria_diagnosis
  )

#----------------Check fever + parasitaemia LCF  criteria 2--------------------------

lcf_data |>
  filter(lcf_fever == TRUE) |>
  select(
    participant_id,
    treatment_arm,
    analysis_day,
    study_visit,
    parasite_density,
    temperature
  ) |>
  arrange(
    participant_id,
    analysis_day
  ) |>
  print(n = Inf)

lcf_data |>
  ungroup() |> 
  summarise(
    n_lcf_fever = n_distinct(
      participant_id[lcf_fever == TRUE]
    )
  )

#----------------LCF danger sign/severe malaria + parasitaemia criteria 1------------------

lcf_danger_severe <- lcf_data |>
  filter(
    
    # Parasitaemia confirmed by treatment-failure microscopy
    tf_microscopy_result == 1,
    
    # Danger sign OR severe malaria
    (
      tf_danger_sign == 1 |
        tf_severe_manifestation == 1
    ),
    
    # Keep only the actual treatment-failure visit
    (
      study_visit == tf_visit |
        (
          tf_visit == "Unscheduled visit" &
            followup_date == tf_visit_date
        )
    )
  ) |>
  select(
    participant_id,
    treatment_arm,
    analysis_day,
    study_visit,
    tf_visit,
    tf_visit_date,
    tf_microscopy_result,
    tf_danger_sign,
    tf_severe_manifestation
  ) |>
  arrange(participant_id, analysis_day)

lcf_danger_severe |>
  print(n = Inf)

#------------------------Final LCF participant----------------------------------

lcf_participants <- bind_rows(
  
  # Criterion 1: danger sign/severe malaria + parasitaemia
  lcf_danger_severe |>
    select(
      participant_id,
      treatment_arm
    ),
  
  # Criterion 2: fever + parasitaemia
  lcf_data |>
    filter(lcf_fever == TRUE) |>
    select(
      participant_id,
      treatment_arm
    )
  
) |>
  distinct(
    participant_id,
    treatment_arm
  )

# Confirm one row per participant
lcf_participants |>
  count(participant_id) |>
  filter(n > 1)

# Confirm all LCF participants came from the No-ETF population
lcf_participants |>
  anti_join(
    no_etf_participants,
    by = "participant_id"
  )

# Total LCF participants
lcf_participants |>
  ungroup() |> 
  summarise(
    n_lcf = n_distinct(participant_id)
  )

# LCF by treatment arm
lcf_participants |>
  count(treatment_arm)

#====================LATE PARASITOLOGICAL FAILURE (LPF)=========================

# Participants who had neither ETF nor LCF
no_etf_lcf_participants <- no_etf_participants |>
  anti_join(
    lcf_participants,
    by = "participant_id"
  )

# Check
no_etf_lcf_participants |>
  count(treatment_arm)

# Confirm no LCF participant remains
no_etf_lcf_participants |>
  semi_join(
    lcf_participants,
    by = "participant_id"
  )

#--------------create day 7 - 42 lpf dataset------------------------------------
lpf_data <- lcf_data |>
  anti_join(
    lcf_participants,
    by = "participant_id"
  ) |>
  filter(
    analysis_day >= 7,
    analysis_day <= 42
  )

lpf_data |>
  ungroup() |> 
  summarise(
    n_participants = n_distinct(participant_id)
  )

#----------------Check  no fever + parasitaemia LPF criteria -------------------
lpf_data <- lpf_data |>
  mutate(
    lpf = case_when(
      
      parasitaemia_present == TRUE &
        !is.na(temperature) &
        temperature < 37.5 ~ TRUE,
      
      !is.na(parasitaemia_present) &
        !is.na(temperature) ~ FALSE,
      
      TRUE ~ NA
    )
  )

#---------------------------check LPF cases-------------------------------------
lpf_data |>
  filter(lpf == TRUE) |>
  select(
    participant_id,
    treatment_arm,
    analysis_day,
    study_visit,
    parasite_density,
    temperature,
    microscopy_positive
  ) |>
  arrange(
    participant_id,
    analysis_day
  ) |>
  print(n = Inf)

#count
lpf_data |>
  ungroup() |> 
  summarise(
    n_lpf = n_distinct(
      participant_id[lpf == TRUE]
    )
  )

#---------------------------Final LPF participants------------------------------
lpf_participants <- lpf_data |>
  filter(lpf == TRUE) |>
  distinct(
    participant_id,
    treatment_arm
  )

# LPF by treatment arm
lpf_participants |>
  count(treatment_arm)

# Total LPF participants
lpf_participants |>
  ungroup() |> 
  summarise(
    n_lpf = n_distinct(participant_id)
  )

#========== ADEQUATE CLINICAL AND PARASITOLOGICAL RESPONSE (ACPR)===============
no_etf_lcf_lpf_participants <- no_etf_lcf_participants |> 
  anti_join(
    lpf_participants,
    by = "participant_id"
  )

#check remaining participant by arm
no_etf_lcf_lpf_participants |> 
  count(treatment_arm)

# confirm no lpf participant exist 
no_etf_lcf_lpf_participants |> 
  semi_join(
    lcf_participants,
     by = "participant_id"
  )

#-----------------create day 42 ACPR dataset------------------------------------
acpr_data <- lpf_data |>
  anti_join(
    lpf_participants,
    by = "participant_id"
  ) |>
  filter(
    analysis_day == 42
  )

#Inspect microscopy at day 42
acpr_data |>
  select(
    participant_id,
    treatment_arm,
    study_visit,
    parasite_density,
    microscopy_positive
  ) |>
  print(n = Inf)

#----------check absence of parasitaemia acpr criteria--------------------------
acpr_data <- acpr_data |>
  mutate(
    parasitaemia_absent = case_when(
      
      # Numerical parasite density confirms absence
      !is.na(parasite_density) &
        parasite_density == 0 ~ TRUE,
      
      # Microscopy confirms absence
      microscopy_positive == 0 ~ TRUE,
      
      # Parasitaemia is present
      !is.na(parasite_density) &
        parasite_density > 0 ~ FALSE,
      
      microscopy_positive == 1 ~ FALSE,
      
      # Cannot determine
      TRUE ~ NA
    )
  )

#count
acpr_data |>
  count(
    treatment_arm,
    parasitaemia_absent,
    .drop = FALSE
  )

# Indeterminate outcome
acpr_data <- acpr_data |>
  mutate(
    acpr_status = case_when(
      parasitaemia_absent == TRUE ~ "ACPR",
      parasitaemia_absent == FALSE ~ "Not ACPR",
      TRUE ~ "Indeterminate"
    )
  )

#check
acpr_data |>
  count(
    treatment_arm,
    acpr_status
  )

#-------------------------Final ACPR participant--------------------------------
acpr_participants <- acpr_data |>
  filter(
    parasitaemia_absent == TRUE
  ) |>
  distinct(
    participant_id,
    treatment_arm
  )

# ACPR by treatment arm
acpr_participants |>
  count(treatment_arm)

# Total ACPR participants
acpr_participants |>
  ungroup() |> 
  summarise(
    n_acpr = n_distinct(participant_id)
  )

#==================FINAL TREATMENT OUTCOME CLASSIFICATION=======================

treatment_outcome <- eligible_participants |>
  select(
    participant_id,
    treatment_arm,
    site
  ) |>
  mutate(
    treatment_outcome = case_when(
      
      # ETF
      participant_id %in%
        etf_data$participant_id[
          etf_data$etf_status == "ETF"
        ] ~ "ETF",
      
      # LCF
      participant_id %in%
        lcf_participants$participant_id ~ "LCF",
      
      # LPF
      participant_id %in%
        lpf_participants$participant_id ~ "LPF",
      
      # ACPR
      participant_id %in%
        acpr_participants$participant_id ~ "ACPR",
      
      # Could not be classified
      TRUE ~ "Indeterminate"
    )
  )

#-------------------------------check final outcome-----------------------------
treatment_outcome |>
  count(treatment_outcome)

treatment_outcome |>
  count(
    treatment_arm,
    treatment_outcome
  )

treatment_outcome |>
  ungroup() |> 
  summarise(
    n_participants = n(),
    n_unique = n_distinct(participant_id)
  )

#check everything add up
treatment_outcome |>
  count(treatment_outcome) |>
  summarise(
    total = sum(n)
  )

#=================DESCRIPTIVE SUMMARY OF ETF/LCF/LPF/ACPR=======================
#------------Add age group to treatment outcome---------------------------------

age_reference <- smart_analysis_final |>
  ungroup() |>
  filter(study_visit == "H0") |>
  select(
    participant_id,
    age
  ) |>
  distinct(participant_id, .keep_all = TRUE)

age_reference <- age_reference |> 
  mutate(
    age_group = case_when(
      age < 60  ~ "<5 years",
      age >= 60 ~ ">=5 years",
      TRUE ~ NA_character_
    )
  )

treatment_outcome <- treatment_outcome |>
  left_join(
    age_reference,
    by = "participant_id"
  ) 

#confirm
treatment_outcome |>
  count(age_group, .drop = FALSE)

#----------------Treatment-outcome proportions by arm---------------------------
outcome_by_arm <-treatment_outcome |> 
  count(
    treatment_arm,
    treatment_outcome
  ) |> 
  group_by(treatment_arm) |> 
  mutate(
    total = sum(n),
    proportion = n / total,
    percentage = round(100* proportion, 1)
  ) |> 
  ungroup()
outcome_by_arm

#----------------Treatment-outcome proportions by age subgroup------------------
outcome_by_age <- treatment_outcome |> 
  filter(!is.na(age_group)) |> 
  count(
    age_group,
    treatment_outcome
  ) |> 
  group_by(age_group) |> 
  mutate(
    total = sum(n),
    proportion = n/total,
    percentage = round(100* proportion, 1)
  ) |> 
  ungroup()
outcome_by_age

#------Treatment-outcome proportions by age subgroup within treatment arm-------
outcome_by_age_arm <- treatment_outcome |> 
  filter(!is.na(age_group)) |> 
  count(
    treatment_arm,
    age_group,
    treatment_outcome
  ) |> 
  group_by(
    treatment_arm,
    age_group
  ) |> 
  mutate(
    total = sum(n),
    proportion = n/total,
    percentage = round(100* proportion, 1)
  ) |> 
  ungroup()
outcome_by_age_arm

#------Treatment-outcome proportions by treatment arm within  site--------------

outcome_by_arm_site <- treatment_outcome |> 
  #filter(!is.na(site)) |> 
  count(
    site,
    treatment_arm,
    treatment_outcome
  ) |> 
  group_by(
    site,
    treatment_arm
  ) |> 
  mutate(
    total = sum(n),
    proportion = n/total,
    percentage = round(100* proportion, 1)
  ) |> 
  ungroup()
outcome_by_arm_site

#-----------------------Label variables-----------------------------------------
treatment_outcome <- treatment_outcome |>
  mutate(
    treatment_outcome = factor(
      treatment_outcome,
      levels = c(
        "ETF",
        "LCF",
        "LPF",
        "ACPR",
        "Indeterminate"
      )
    ),
    
    age_group = factor(
      age_group,
      levels = c("<5 years", ">=5 years")
    ),
    
    treatment_arm = factor(
      treatment_arm,
      levels = c("AL", "AP"),
      labels = c("Artemether-Lumefantrine", "Artesunate-Pyranoridine")
    ),
    
    site = factor(
      site,
      levels = c(1, 2),
      labels = c("Assin Fosu", "Agogo")
    )
  )

#==========================TREATMENT OUTCOME TABLES==============================
#-------------------1. Outcome by arm descriptive table----------------------------

table_outcome_arm <- treatment_outcome |> 
  select(
    treatment_arm,
    treatment_outcome
  ) |> 
  tbl_summary(
    by = treatment_arm,
    
    statistic = list(
      all_categorical() ~ "{n} ({p}%)"
    ),
    
    missing = "no",
    
    label = list(
      treatment_outcome ~"Treatment outcome"
    )
  ) |> 
  add_overall() |> 
  bold_labels()

table_outcome_arm_gt <- table_outcome_arm |> 
  as_gt() |> 
  gt:: tab_header(
    title = gt::md("**PCR-uncorrected treatment outcomes by treatment arm**")
  )
table_outcome_arm_gt

#save the table
gt::gtsave(
  table_outcome_arm_gt, "figures/treatment_outcomes_arm.png"
)

#-------------------2. Outcome by age group descriptive table-------------------------
table_outcome_age <- treatment_outcome |>
  select(
    age_group,
    treatment_outcome
  ) |>
  tbl_summary(
    by = age_group,
    
    statistic = list(
      all_categorical() ~ "{n} ({p}%)"
    ),
    
    missing = "no",
    
    label = list(
      treatment_outcome ~ "Treatment outcome"
    )
  ) |>
  add_overall() |>
  bold_labels()

table_outcome_age_gt <- table_outcome_age |>
  as_gt() |>
  gt::tab_header(
    title = gt::md("**PCR-uncorrected treatment outcomes by age group**")
  )
table_outcome_age_gt

#save the table
gt::gtsave(
  table_outcome_age_gt, "figures/treatment_outcomes_age.png"
)

#-------------------3. Outcome by age group with arm descriptive table-------------------------
table_outcome_age_arm <- treatment_outcome |>
  select(
    treatment_arm,
    age_group,
    treatment_outcome
  ) |>
  tbl_strata(
    strata = treatment_arm,

    .tbl_fun = ~ .x |>
      tbl_summary(
        by = age_group,

        statistic = list(
          all_categorical() ~ "{n} ({p}%)"
        ),

        missing = "no",

        label = list(
          treatment_outcome ~ "Treatment outcome"
        )
      ) |> 
      bold_labels(),
    .combine_with = "tbl_merge"
    
  )

table_outcome_age_arm_gt <- table_outcome_age_arm |> 
  as_gt() |> 
  gt::tab_header(
    title = gt::md("**PCR-uncorrected treatment outcomes by treatment arm and age group**")
  )
table_outcome_age_arm_gt

#save
gt::gtsave(
  table_outcome_age_arm_gt,
  "figures/treatment_outcomes_age_arm.png"
)

#-------------------4. Outcome by arm within site descriptive table-------------------------
table_outcome_arm_site <- treatment_outcome |>
  select(
    site,
    treatment_arm,
    treatment_outcome
  ) |>
  tbl_strata(
    strata = site,
    
    .tbl_fun = ~ .x |>
      tbl_summary(
        by = treatment_arm,
        
        statistic = list(
          all_categorical() ~ "{n} ({p}%)"
        ),
        
        missing = "no",
        
        label = list(
          treatment_outcome ~ "Treatment outcome"
        )
      ) |> 
      bold_labels(),
    .combine_with = "tbl_merge"
    
  )

table_outcome_arm_site_gt <- table_outcome_arm_site |> 
  as_gt() |> 
  gt::tab_header(
    title = gt::md("**PCR-uncorrected treatment outcomes by site and treatment arm**")
  )
table_outcome_arm_site_gt

#save
gt::gtsave(
  table_outcome_arm_site_gt,
  "figures/treatment_outcomes_arm_site.png"
)

