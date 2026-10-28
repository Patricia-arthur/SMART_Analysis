source("Data/SMARTUncomplicatedMa-Datamain2_R_2026-10-27_0948.r")

library(tidyverse)
library(lubridate)
library(survival)
library(survminer)
library(janitor)
library(naniar)

#------------------STUDY DATASET-----------------------------------------------------#
#       
#-----------------------------------------------------------------------#
SMARTUncomplicated <- data

#View(SMARTUncomplicated)
dim(SMARTUncomplicated)

SMARTUncomplicated |>
  summarise(
    rows = n(),
    columns = ncol(SMARTUncomplicated),
    total_records = n(),
    unique_participants = n_distinct(id_um, na.rm = TRUE),
    missing_id = sum(is.na(id_um) | id_um == "")
  )

#Checking for missing ID
SMARTUncomplicated |>
  filter(is.na(id_um) | id_um =="") |>
  select(
    redcap_event_name,
  )

#Participants with more than 1 event
SMARTUncomplicated |> 
  count(id_um, redcap_event_name, sort =TRUE) |>
  filter(n > 1)

SMARTUncomplicated |> 
  count(id_um, redcap_event_name) |>
  filter(n > 1) |>
  count(redcap_event_name, sort=TRUE)

#Number per each event
SMARTUncomplicated |>
  count(redcap_event_name, sort=TRUE)

#check duplicate records
nrow(SMARTUncomplicated) - nrow(distinct(SMARTUncomplicated))  


#------ Participant in each treatment group ------#
SMARTUncomplicated |>
  filter(redcap_event_name =="hour_0_day_0_arm_1") |>
  count(
    j01_which_treatment_has_be,
    j01_which_treatment_has_be.factor,
    .drop = FALSE
  )

SMARTUncomplicated |>
  filter(redcap_event_name =="hour_0_day_0_arm_1") |>
  count(j01_which_treatment_has_be)

#Participants with missing treatment
SMARTUncomplicated |>
  filter (
    redcap_event_name == "hour_0_day_0_arm_1",
    is.na(j01_which_treatment_has_be)
  ) |>
  select(
    id_um,
    h01_3um,
    j01_which_treatment_has_be,
    j01_which_treatment_has_be.factor
  )

#check treatment existence
missing_treatment_id <- SMARTUncomplicated |>
  filter(
    redcap_event_name == "hour_0_day_0_arm_1",
    is.na(j01_which_treatment_has_be)
  ) |>
  pull(id_um)

SMARTUncomplicated |>
  filter(id_um %in% missing_treatment_id) |>
  select(
    id_um,
    h01_3um,
    redcap_event_name,
    starts_with("da_al_"),
    starts_with("da_ap_")
  )

SMARTUncomplicated |>
  filter(id_um %in% missing_treatment_id) |>
  select(
    id_um,
    h01_3um,
    redcap_event_name,
    j01_which_treatment_has_be,
    j01_which_treatment_has_be.factor,
    e04um,
    e04um.factor,
    e04um_2,
    e04um_2.factor,
    da_ap_pid,
    dadmin_ap,
    da_ap_b04,
    da_ap_b05,
    da_al_pid,
    dadmin1,
    da_al_b04,
    da_al_b05
  )

SMARTUncomplicated |>
  filter(id_um %in% missing_treatment_id) |>
  count(redcap_event_name)

#checking for other participants
SMARTUncomplicated |>
  filter(redcap_event_name == "hour_0_day_0_arm_1") |>
  filter(!is.na(j01_which_treatment_has_be)) |>
  group_by(id_um) |>
  summarise(
    n_treatments = n_distinct(j01_which_treatment_has_be),
    .groups = "drop"
  ) |>
  filter(n_treatments > 1)

#check if a participant has more than one baseline record
SMARTUncomplicated |>
  filter(redcap_event_name == "hour_0_day_0_arm_1") |>
  count(id_um) |>
  filter(n > 1)

#creating treatment for participant across all records
treatment_lookup <- SMARTUncomplicated |>
  filter(redcap_event_name == "hour_0_day_0_arm_1") |>
  transmute(
    id_um,
    treatment_arm = case_when(
      j01_which_treatment_has_be == 1 ~ "AL",
      j01_which_treatment_has_be == 2 ~ "AP",
      TRUE ~ NA_character_
    )
  )

SMARTUncomplicated <- SMARTUncomplicated |>
  left_join(
    treatment_lookup,
    by = "id_um"
  )

SMARTUncomplicated |>
  distinct(id_um, treatment_arm) |>
  count(treatment_arm, .drop = FALSE)

SMARTUncomplicated |> select(
  id_um,
  h01_3um,
  redcap_event_name,
  j01_which_treatment_has_be.factor,
  treatment_arm
) #|>
  #View()

#--------------------------Mapping visit schedules------------------------------
SMARTUncomplicated <- SMARTUncomplicated |>
  mutate(
    study_day = case_when(
      redcap_event_name %in% c(
        "hour_0_day_0_arm_1", 
        "hour_8_day_0_arm_1"
        ) ~ "Day 0",
      redcap_event_name %in% c(
        "hour_24_day_1_arm_1", 
        "hour_36_day_1_arm_1"
        ) ~ "Day 1",
      redcap_event_name %in% c(
        "hour_48_day_2_arm_1", 
        "hour_60_day_2_arm_1"
        ) ~"Day 2",
      
      redcap_event_name == "hour_72_day_3_arm_1" ~ "Day 3",
      redcap_event_name == "day_7_arm_1" ~ "Day 7",
      redcap_event_name == "day_14_arm_1" ~ "Day 14",
      redcap_event_name == "day_21_arm_1" ~ "Day 21",
      redcap_event_name == "day_28_arm_1" ~ "Day 28",
      redcap_event_name == "day_35_arm_1" ~ "Day 35",
      redcap_event_name == "day_42_arm_1" ~ "Day 42",
       
      redcap_event_name == "unsheduled_visit_arm_1" ~ "Unscheduled visit",
      redcap_event_name == "treatment_failure_arm_1" ~ "Treatment failure",
      redcap_event_name == "concomitant_medica_arm_1" ~ "Concomitant medication",
      redcap_event_name == "schedule_visit_arm_1" ~ "Scheduled visit",
      
      TRUE ~NA_character_
    ),
    
    study_visit = case_when(
      redcap_event_name == "hour_0_day_0_arm_1" ~ "H0",
      redcap_event_name == "hour_8_day_0_arm_1" ~ "H8",
      redcap_event_name == "hour_24_day_1_arm_1" ~ "H24",
      redcap_event_name == "hour_36_day_1_arm_1" ~ "H36",
      redcap_event_name == "hour_48_day_2_arm_1" ~ "H48",
      redcap_event_name == "hour_60_day_2_arm_1" ~ "H60",
      redcap_event_name == "hour_72_day_3_arm_1" ~ "H72",
      redcap_event_name == "day_7_arm_1" ~ "D7",
      redcap_event_name == "day_14_arm_1" ~ "D14",
      redcap_event_name == "day_21_arm_1" ~ "D21",
      redcap_event_name == "day_28_arm_1" ~ "D28",
      redcap_event_name == "day_35_arm_1" ~ "D35",
      redcap_event_name == "day_42_arm_1" ~ "D42",
      
      redcap_event_name == "unsheduled_visit_arm_1" ~ "Unscheduled visit",
      redcap_event_name == "treatment_failure_arm_1" ~ "Treatment failure",
      redcap_event_name == "concomitant_medica_arm_1" ~ "Concomitant medication",
      redcap_event_name == "schedule_visit_arm_1" ~ "Scheduled visit",
      
      TRUE ~NA_character_
    ) 
  )


SMARTUncomplicated |>
  count(study_day, study_visit, sort = TRUE)

# c("study_day", "study_visit", "sheduled_visit") %in% names(SMARTUncomplicated)
# SMARTUncomplicated |>
#   distinct(redcap_event_name) |>
#   filter(
#     str_detect(redcap_event_name, "^visit")
#   )

# SMARTUncomplicated <- SMARTUncomplicated |>
#   mutate(
#     study_day = factor(
#       study_day,
#       levels = c(
#         "Day 0", "Day 1", "Day 2", "Day 3", "Day 7","Day 14", "Day 21", "Day 28", "Day 35", 
#         "Day 42", "Unscheduled visit", "Treatment failure", "Concomitant medication", "Scheduled visit"
#       )
#     ),
#     study_visit = factor(
#       study_visit,
#       levels = c("H0","H8","H24","H36","H48","H60","H72","D7","D14","D21","D28","D35","D42",
#                  "Unscheduled visit", "Treatment failure", "Concomitant medication", "Scheduled visit")
#     )
#   )


#------------------Temperature variables within data---------------------------
temperature_mapping <- SMARTUncomplicated |>
  select(
    redcap_event_name,
    study_day,
    study_visit,
    
    #baseline temp
    c03_1um,
    c03_2um_2,
    
    #followup H8-H60
    a04_1_fuv,
    a04_3_2_fuv,
    
    #followup D3-42
    a03_1_d_um,
    a03_2_d_um,
    
    #treatment failure
    b01_1_tf_um_sm,
    b01_2_tf_um_sm,
    
    #unscheduled visit
    b03_1_uv,
    b03_2_uv
  )

temperature_mapping |>
  summarise(
    across(
      c( c03_1um, c03_2um_2, a04_1_fuv, a04_3_2_fuv, a03_1_d_um,  
         a03_2_d_um, b01_1_tf_um_sm, b01_2_tf_um_sm, b03_1_uv, b03_2_uv
      ),
      ~ sum(!is.na(.))
    ),
    .by = redcap_event_name
  )

#------------------Hb variables---------------------
names(SMARTUncomplicated) |>
  str_subset(
    regex(
      "d05um|f02um|a07_d_um|c03_haemoglobin|haemoglobin|hgb_edt|haemoglobin|hgb",
      ignore_case = TRUE
    )
  )

SMARTUncomplicated |>
  select(
    redcap_event_name,
    d05um, a07_d_um, c03_haemoglobin,
    
  )

SMARTUncomplicated |>
  select(
    redcap_event_name, d05um, a07_d_um,
    c03_haemoglobin, hgb_edta
  ) |>
  summarise(
    across(
      everything(),
      ~sum(!is.na(.))
    ),
    .by = redcap_event_name
  )

#-------------------------Biochemistry variables-----------
names(SMARTUncomplicated) |>
  str_subset(
    regex(
      paste(
        "c01_date_and_time_of_bioch",
        "c03_1_bilt", "c03_2_if_done",
        "c04_1_bild", "c04_1",
        "alb_heparin", "c05_2_if_done",
        "alt_heparin", "c05_2_if_done_2",
        "c07_heparin", "c07_2_if_done",
        "alp_heparin", "c08_2_if_done",
        "ggt_heparin", "c09_2_if_done",
        "pro_heparin", "c10_2_if_done",
        sep = "|"
      ),
      ignore_case = TRUE
    )
  )

SMARTUncomplicated |>
  select(
    redcap_event_name,
    c01_date_and_time_of_bioch,
    
    # Actual biochemistry results
    c03_2_if_done,       # Total bilirubin
    c04_1,               # Direct bilirubin
    c05_2_if_done,       # Albumin
    c05_2_if_done_2,     # ALT
    c07_2_if_done,       # AST
    c08_2_if_done,       # ALP
    c09_2_if_done,       # GGT
    c10_2_if_done        # Total protein
  ) |>
  summarise(
    across(
      everything(),
      ~sum(!is.na(.))
    ),
    .by = redcap_event_name
  )

#--------------Drug administration variables-------------------------------
names(SMARTUncomplicated) |>
  str_subset(
    regex(
      "da_al_|^da_ap_|^dadmin",
      ignore_case = TRUE
    )
  )

SMARTUncomplicated|>
  select(
    redcap_event_name,
    
    # AP
    da_ap_b01_1,
    da_ap_b01_2,
    da_ap_b02,
    dadmin_ap,
    da_ap_b03_2,
    da_ap_b04,
    da_ap_b05,
    da_ap_b07_1,
    da_ap_b07_2,
    da_ap_b07_3,
    
    # AL
    da_al_b01_1,
    da_al_b01_2,
    da_al_b02,
    dadmin1,
    b03_da,
    da_al_b04,
    da_al_b05,
    da_al_b07_1,
    da_al_b07_2,
    da_al_b07_3
    
  ) |>
  summarise(
    across(
      everything(),
      ~sum(!is.na(.))
    ),
    .by = redcap_event_name
  )


#-------------------Drug administered for AL and AP--------------
SMARTUncomplicated |>
  count(
    redcap_event_name,
    da_ap_b01_1,
    sort = TRUE
  )

SMARTUncomplicated |>
  count(
    redcap_event_name,
    da_al_b01_1,
    sort = TRUE
  )

SMARTUncomplicated |>
  count(
    redcap_event_name,
    da_ap_b01_1,
    da_ap_b01_2,
    sort = TRUE
  )

SMARTUncomplicated |>
  count(
    redcap_event_name,
    da_al_b01_1,
    da_al_b01_2,
    sort = TRUE
  )

SMARTUncomplicated |>
  count(
    redcap_event_name,
    da_ap_b01_1,
    dadmin_ap,
    sort =TRUE
  )

SMARTUncomplicated |>
  count(
    redcap_event_name,
    da_al_b01_1,
    da_al_b02,
    dadmin1,
    sort = TRUE
  )

# vomitting occured
SMARTUncomplicated |>
  filter(da_ap_b01_1 == 1) |>
  count(
    redcap_event_name,
    da_ap_b03_2,
    da_ap_b07_1,
    da_ap_b07_3,
    sort = TRUE
  )

SMARTUncomplicated |>
  filter(da_al_b01_1 == 1) |>
  count(
    redcap_event_name,
    b03_da,
    da_al_b07_1,
    da_al_b07_3,
    sort = TRUE
  )

#Drug re-administered (second dose)
# AP re-administration
SMARTUncomplicated |>
  filter(da_ap_c01 == 1) |>
  count(
    redcap_event_name,
    rdadmin,
    da_ap_c03_2,
    da_ap_c07_1,
    sort = TRUE
  )

# AL re-administration
SMARTUncomplicated |>
  filter(da_al_c01 == 1) |>
  count(
    redcap_event_name,
    dadmin,
    da_al_c03,
    da_al_c07_1,
    sort = TRUE
  )


#-------------------Adverse events (AE)-----

SMARTUncomplicated |>
  count(
    redcap_event_name,
    a05_1_fuv,
    a04_1_d_um,
    b04_any_new_complain_since,
    sort = TRUE
  )

SMARTUncomplicated |>
  filter(a05_1_fuv == 1) |>
  select(
    redcap_event_name,
    starts_with("a05_")
  ) |>
  glimpse()

SMARTUncomplicated |>
  filter(a04_1_d_um == 1) |>
  select(
    redcap_event_name,
    a04_1_d_um,
    starts_with("a04_2_d_um"),
    if_nsd_specify,
    if_gid_specify,
    if_rd_specify,
    if_md_specify,
    if_hlsd_specify,
    if_sd_specify,
    if_id_specify,
    if_gud_specify,
    if_ed_specify,
    if_enod_specify,
    if_end_specify,
    if_cd_specify,
    if_others_specify_complain
  ) |>
  glimpse()

SMARTUncomplicated |>
  filter(b04_any_new_complain_since == 1) |>
  select(
    redcap_event_name,
    starts_with("b04_")
  ) |>
  glimpse()

#-------------------Concomitant medications----
names(SMARTUncomplicated) |>
  str_subset(
    regex(
      "cm_|medication|medname|route|dose|concom",
      ignore_case = TRUE
    )
  )

SMARTUncomplicated |>
  filter(redcap_event_name == "concomitant_medica_arm_1") |>
  select(
    cm_pid_rm,
    cmno_rm,
    medication_name_rm,
    medname_rm,
    cm_ind_rm,
    cm_oth_rm,
    route_rm,
    route_oth_rm,
    sdate_rm,
    ts_rm,
    edate_rm
  ) |>
  summarise(
    across(
      everything(),
      ~sum(!is.na(.) & .!="")
    )
  )

SMARTUncomplicated |>
  filter(redcap_event_name == "concomitant_medica_arm_1") |>
  count(cm_ind_rm, cm_ind_rm.factor, sort = TRUE)

SMARTUncomplicated |>
  filter(redcap_event_name == "concomitant_medica_arm_1") |>
  count(medname_rm, medname_rm.factor, sort = TRUE)

### Microscopy/ Parasitaemia ###
names(SMARTUncomplicated) |>
  str_subset(
    regex(
      "paras|micros|smear|slide|density|troph|gamet|malaria|blood.?film",
      ignore_case = TRUE
    )
  )

names(SMARTUncomplicated) |>
  str_subset(
    regex(
      "^d02_1um_2$|^d02_2um_2$|^d02_4um$|malaria_micros",
      ignore_case = TRUE
    )
  )

SMARTUncomplicated |>
  select(
    redcap_event_name,
    d02_1um_2,
    d02_2um_2,
    d02_4um,
    b03_malaria_microscopy_pos_sm,
    b03_malaria_microscopy_pos_2_sm,
    b03_malaria_microscopy_pos_3_sm,
    c02_malaria_microscopy_don
  )|>
  summarise(
    across(
      everything(),
      ~sum(!is.na(.))
    ),
    .by = redcap_event_name
  )

#-------------------Pharmacokinetic variables------

names(SMARTUncomplicated)|>
  str_subset(
    regex(
      "pk|pharmac|plasma|lumefantrine|pyronaridine|aliquot",
      ignore_case = TRUE
    )
  )

#-------------------Check PK variables identified from the codebook----------

names(SMARTUncomplicated) |>
  str_subset(
    regex(
      "sptaken_heparin_fo|blab_dat_fo|pk_fo|b02_2_if_yes|bfx2_fo|b04_dtstorage_fo|
      hep_bio_3|blood_heparin|b05_biochemistry_performed|b05_1_if_yes|
      b06_biochemistry_performed|b06_1_if_yes",
      ignore_case = TRUE
    )
  )


#================DATA DICTIONARY======================================

SMART_dictionary <- tribble(
  ~domain, ~variable, ~description, ~visit_source, ~analysis_role,
  
  #-----------------------------------------------------------------------
  #   Identifiers and study structure
  #-----------------------------------------------------------------------
  
  "Identifiers & structure", "h01_3um", "Participant ID", "Participant-level identifier", 
  "Primary identifier and dataset linkage key",
  
  "Identifiers & structure", "id_um", "Screening ID", "All records", 
  "Secondary identifier",
  
  "Identifiers & structure", "redcap_event_name", "Original REDCap event name", 
  "All records",
  "Original visit/event structure",
  
  "Identifiers & structure", "study_day", "Derived study day", "Day 0 to Day 42, 
  treatment failure, unscheduled and other study events",
  "Main study-day variable for analysis",
  
  "Identifiers & structure", "study_visit", "Derived study visit", "H0, H8, H24,
  H36, H48, H60, H72,
  D7, D14, D28,D35, D42", "Detailed visit/time-point variable",
  
  "Indentifiers & structure", "treatment arm", "Derived treatment arm: 
  Artemether-Lumefantrine or Artesunate- Pyronaridine",
  "Participant level", "Main treatment exposure",
  
  #-----------------------------------------------------------------------
  #   Temperature
  #-----------------------------------------------------------------------
  
  "Temperature", "c03_1um", "Temperature", "H0", "Baseline temperature",
  
  "Temperature", "c03_2um_2", "Temperature route", "H0", 
  "Temperature measurement information",
  
  "Temperature", "a04_1_fuv", "Temperature", "H8-H60", 
  "Early fever-resolution assessment",
  
  "Temperature", "a04_3_2_fuv", "Temperature route", "H8-H60", 
  "Temperature measurement information",
  
  "Temperature", "a03_1_d_um", "Temperature", "H72 and D7-D47", 
  "Follow-up fever/clinical assessment",
  
  "Temperature", "a03_2_d_um", "Temperature route", "H72 and D7-D42", 
  "Temperature measurement information",
  
  "Temperature", "b01_1_tf_um_sm", "Temperature", "Treatment failure", 
  "Clinical assessment at treatment failure",
  
  "Temperature", "b01_2_tf_um_sm", "Temperature route", "Treatment failure", 
  "Temperature measurement information",
  
  "Temperature", "b03_1_uv", "Temperature", "Unscheduled visit", 
  "Clinical assessment at unschduled visit",
  
  "Temperature", "b03_2_uv", "Temperature route", "Unscheduled visit", 
  "Temperature measurement information",
  
  # ---------------------------------------------------------------------------
  #    Haemoglobin
  # ---------------------------------------------------------------------------
  "Haemoglobin", "d05um", "Haemoglobin (g/dL)", "H0", "Baseline haemoglobin",
  
  "Haemoglobin", "a07_d_um", "Follow-up haemoglobin (g/dL)", "H72 and D7-D42", 
  "Follow-up haematological assessment",
  
  "Haemoglobin", "c03_haemoglobin", "Haemoglobin (g/dL)", "Unscheduled visit", 
  "Unscheduled haematological assessment",
  
  "Haemoglobin", "hgb_edta", "Laboratory HGB from EDTA sample (g/dL)", 
  "Laboratory/EDTA records", "Laboratory haemoglobin assessment",
  
  # ---------------------------------------------------------------------------
  #    Microscopy
  # ---------------------------------------------------------------------------
  "Microscopy", "d02_1um_2", "Malaria microscopy done?", "H0", 
  "Baseline microscopy assessment",
  
  "Microscopy", "d02_2um_2", "Malaria microscopy result", "H0", 
  "Baseline malaria confirmation",
  
  "Microscopy", "d02_4um", "P. falciparum asexual parasite only?", "H0", 
  "Baseline parasite characterization",
  
  "Microscopy", "b03_malaria_microscopy_pos_sm", "Malaria microscopy positive?",
  "Treatment failure", "Treatment-failure confirmation",
  
  "Microscopy", "b03_malaria_microscopy_pos_2_sm", "RDT positive?", "Treatment failure", 
  "Supporting malaria diagnosis",
  
  "Microscopy", "b03_malaria_microscopy_pos_3_sm", "Presence of severe manifestation
  of P. falciparum malaria?", "Treatment failure", "Severity/treatment-failure assessment",
  
  "Microscopy", "c02_malaria_microscopy_don", "Malaria microscopy done?", 
  "Unscheduled visit", "Assessment for possible recurrence",
  
  # ---------------------------------------------------------------------------
  # Parasitaemia
  # ---------------------------------------------------------------------------
  "Parasitaemia", NA_character_, "Quantitative parasite density/count", 
  "Expected H0-H72 and later relevant follow-up visits",
  "Pending separate parasitaemia dataset; required for parasite clearance and recurrent parasitaemia",
  
  # ---------------------------------------------------------------------------
  # AP drug administration
  # ---------------------------------------------------------------------------
  "Drug administration - AP", "da_ap_b01_1", "Study drug given?", "AP dosing visits",
  "Confirms scheduled treatment receipt",
  
  "Drug administration - AP", "dadmin_ap", "Was Artesunate-Pyronaridine administered?",
  "AP dosing visits", "Confirms actual AP administration",
  
  "Drug administration - AP", "da_ap_b03_2", "AP dose administered",
  "AP dosing visits", "Treatment dose/exposure",
  
  "Drug administration - AP", "da_ap_b04", "Date of AP administration",
  "AP dosing visits", "Treatment timing",
  
  "Drug administration - AP", "da_ap_b05",
  "Time of AP administration/start of observation","AP dosing visits", 
  "Treatment timing",
  
  "Drug administration - AP", "da_ap_b07_1", "Vomited after AP administration?",
  "AP dosing visits", "Treatment tolerance/exposure",
  
  "Drug administration - AP", "da_ap_b07_2", "Time of vomiting after AP",
  "AP dosing visits", "Treatment exposure assessment",
  
  "Drug administration - AP", "da_ap_b07_3", "Was AP re-administered after vomiting?",
  "AP dosing visits", "Determines effective treatment exposure",
  
  "Drug administration - AP", "da_ap_c01","Is a second AP dose being administered?", 
  "AP re-administration", "Re-administration indicator",
  
  "Drug administration - AP", "rdadmin","Was AP re-administered?",
  "AP re-administration", "Confirms replacement administration",
  
  "Drug administration - AP", "da_ap_c03_2", "Replacement AP dose",
  "AP re-administration", "Replacement treatment exposure",
  
  "Drug administration - AP", "da_ap_c04", "Date of AP re-administration",
  "AP re-administration", "Replacement treatment timing",
  
  "Drug administration - AP", "da_ap_c05", "Time of AP re-administration",
  "AP re-administration", "Replacement treatment timing",
  
  "Drug administration - AP", "da_ap_c07_1", "Vomited after AP re-administration?",
  "AP re-administration", "Tolerance of replacement dose",
  
  # ---------------------------------------------------------------------------
  # AL drug administration
  # ---------------------------------------------------------------------------
  "Drug administration - AL", "da_al_b01_1", "Study drug given?",
  "AL dosing visits", "Confirms scheduled treatment receipt",
  
  "Drug administration - AL", "dadmin1", "Was Artemether-Lumefantrine administered?",
  "AL dosing visits", "Confirms actual AL administration",
  
  "Drug administration - AL", "b03_da", "AL dose administered",
  "AL dosing visits", "Treatment dose",
  
  "Drug administration - AL", "da_al_b04", "Date of AL administration",
  "AL dosing visits", "Treatment timing",
  
  "Drug administration - AL", "da_al_b05",
  "Time of AL administration/start of observation",
  "AL dosing visits", "Treatment timing",
  
  "Drug administration - AL", "da_al_b07_1", "Vomited after AL administration?",
  "AL dosing visits", "Treatment tolerance",
  
  "Drug administration - AL", "da_al_b07_2", "Time of vomiting after AL",
  "AL dosing visits", "Treatment exposure assessment",
  
  "Drug administration - AL", "da_al_b07_3", "Was AL re-administered after vomiting?",
  "AL dosing visits", "Determines effective treatment exposure",
  
  "Drug administration - AL", "da_al_c01", 
  "Is a second AL dose being administered?", "AL re-administration",
  "Re-administration indicator",
  
  "Drug administration - AL", "dadmin", "Was AL re-administered?",
  "AL re-administration", "Confirms replacement administration",
  
  "Drug administration - AL", "da_al_c03", "Replacement AL dose",
  "AL re-administration", "Replacement treatment exposure",
  
  "Drug administration - AL", "da_al_c04", "Date of AL re-administration",
  "AL re-administration", "Replacement treatment timing",
  
  "Drug administration - AL", "da_al_c05", "Time of AL re-administration",
  "AL re-administration", "Replacement treatment timing",
  
  "Drug administration - AL", "da_al_c07_1", "Vomited after AL re-administration?",
  "AL re-administration","Tolerance of replacement dose",
  
  # ---------------------------------------------------------------------------
  # New clinical complaints / adverse-event indicators
  # ---------------------------------------------------------------------------
  "Clinical complaints / AE indicators", "a05_1_fuv",
  "Any new complaint since last visit?",
  "H8-H60", "Early follow-up safety",
  
  "Clinical complaints / AE indicators", "a05_2_fuv___1 to ___13",
  "Affected system(s)", "H8-H60", "Classification of new complaint",
  
  "Clinical complaints / AE indicators", "a05_*_specify",
  "Specific complaint description",
  "H8-H60", "Detailed clinical complaint",
  
  "Clinical complaints / AE indicators", "a04_1_d_um",
  "Any new complaint since last visit?",
  "H72 and D7-D42", "Follow-up safety",
  
  "Clinical complaints / AE indicators", "a04_2_d_um___1 to ___13",
  "Affected system(s)", "H72 and D7-D42", "Classification of new complaint",
  
  "Clinical complaints / AE indicators", "if_*_specify",
  "Specific complaint description",
  "H72 and D7-D42", "Detailed clinical complaint",
  
  "Clinical complaints / AE indicators", "b04_any_new_complain_since",
  "Any new complaint since last visit?",
  "Unscheduled visit","Unscheduled safety",
  
  "Clinical complaints / AE indicators", "b04_2_if_yes_tick_all_syst___1 to ___13",
  "Affected system(s)", "Unscheduled visit", "Classification of new complaint",
  
  "Clinical complaints / AE indicators", "b04_3_if_*_specify",
  "Specific complaint description", "Unscheduled visit",
  "Detailed clinical complaint",
  
  # ---------------------------------------------------------------------------
  # Concomitant medications
  # ---------------------------------------------------------------------------
  "Concomitant medication", "cm_pid_rm", "Participant ID",
  "Concomitant medication form", "Record linkage",
  
  "Concomitant medication", "cmno_rm", "Concomitant medication number",
  "Concomitant medication form", "Medication record identifier",
  
  "Concomitant medication", "medication_name_rm", "Medication generic name",
  "Concomitant medication form", "Concomitant medication exposure",
  
  "Concomitant medication", "medname_rm", "Forbidden medication?",
  "Concomitant medication form", "Protocol compliance assessment",
  
  "Concomitant medication", "cm_ind_rm", "Medication indication",
  "Concomitant medication form", "Reason for medication use",
  
  "Concomitant medication", "cm_oth_rm", "Other indication specified",
  "Concomitant medication form", "Additional indication information",
  
  "Concomitant medication", "route_rm",
  "Medication route", "Concomitant medication form",
  "Medication administration information",
  
  "Concomitant medication", "route_oth_rm", "Other medication route specified",
  "Concomitant medication form", "Additional administration information",
  
  "Concomitant medication", "sdate_rm", "Medication start date",
  "Concomitant medication form", "Medication exposure timing",
  
  "Concomitant medication", "ts_rm", "Treatment stopped?",
  "Concomitant medication form", "Medication status",
  
  "Concomitant medication", "edate_rm", "Medication end date",
  "Concomitant medication form", "Medication exposure duration"
)

write_csv(SMART_dictionary, "Discrepancy data/SMART_dictionary.csv")

SMART_supporting_dictionary <- tribble(
  ~domain, ~variable, ~description, ~visit_source, ~analysis_role,
  
  # ---------------------------------------------------------------------------
  # Biochemistry
  # ---------------------------------------------------------------------------
  "Biochemistry", "c01_date_and_time_of_bioch",
  "Date and time biochemistry was performed",
  "Laboratory records", "Laboratory timing",
  
  "Biochemistry", "c03_2_if_done",
  "Total bilirubin",
  "H0, H72, D7, D14, D28, D42 and selected other visits",
  "Potential safety/liver-function assessment",
  
  "Biochemistry", "c04_1",
  "Direct bilirubin",
  "H0, H72, D7, D14, D28, D42 and selected other visits",
  "Potential safety/liver-function assessment",
  
  "Biochemistry", "c05_2_if_done",
  "Albumin",
  "H0, H72, D7, D14, D28, D42 and selected other visits",
  "Potential safety/liver-function assessment",
  
  "Biochemistry", "c05_2_if_done_2",
  "Alanine aminotransferase (ALT)",
  "H0, H72, D7, D14, D28, D42 and selected other visits",
  "Potential liver-safety assessment",
  
  "Biochemistry", "c07_2_if_done",
  "Aspartate aminotransferase (AST)",
  "H0, H72, D7, D14, D28, D42 and selected other visits",
  "Potential liver-safety assessment",
  
  "Biochemistry", "c08_2_if_done",
  "Alkaline phosphatase (ALP)",
  "H0, H72, D7, D14, D28, D42 and selected other visits",
  "Potential liver-safety assessment",
  
  "Biochemistry", "c09_2_if_done",
  "Gamma-glutamyl transferase (GGT)",
  "H0, H72, D7, D14, D28, D42 and selected other visits",
  "Potential liver-safety assessment",
  
  "Biochemistry", "c10_2_if_done",
  "Total protein",
  "H0, H72, D7, D14, D28, D42 and selected other visits",
  "Potential biochemical/safety assessment",
  
  # ---------------------------------------------------------------------------
  # Pharmacokinetic sample collection and processing - AP
  # ---------------------------------------------------------------------------
  "PK sample processing - AP", "sptaken_heparin_fo",
  "Fluoride oxalate sample taken for AP pharmacokinetics",
  "PK sampling records", "Specimen availability",
  
  "PK sample processing - AP", "blab_dat_fo",
  "Date and time AP PK sample received in laboratory",
  "PK laboratory records", "Specimen timing",
  
  "PK sample processing - AP", "pk_fo", "AP PK sample aliquoted (BFX1)",
  "PK laboratory records", "Specimen processing indicator",
  
  "PK sample processing - AP", "b02_2_if_yes",
  "Volume of AP PK aliquot", "PK laboratory records",
  "Specimen processing information",
  
  "PK sample processing - AP", "bfx2_fo",
  "AP backup PK sample prepared (BFX2)",
  "PK laboratory records",
  "Backup specimen indicator",
  
  "PK sample processing - AP", "b02_2_if_yes_2",
  "Volume of AP backup PK aliquot",
  "PK laboratory records",
  "Specimen processing information",
  
  "PK sample processing - AP", "b04_dtstorage_fo",
  "Date and time AP PK samples stored at -20°C",
  "PK laboratory records",
  "Specimen storage information",
  
  # ---------------------------------------------------------------------------
  # Pharmacokinetic sample collection and processing - AL
  # ---------------------------------------------------------------------------
  "PK sample processing - AL", "hep_bio_3",
  "Heparin blood sample collected for AL pharmacokinetics",
  "PK sampling records",
  "Specimen availability",
  
  "PK sample processing - AL", "blood_heparin",
  "Blood centrifugation performed for biochemistry and AL PK",
  "Laboratory records",
  "Specimen processing indicator",
  
  "PK sample processing - AL", "b05_biochemistry_performed",
  "AL PK sample aliquoted (BHP2)",
  "PK laboratory records",
  "Specimen processing indicator",
  
  "PK sample processing - AL", "b05_1_if_yes",
  "Volume of AL PK aliquot",
  "PK laboratory records",
  "Specimen processing information",
  
  "PK sample processing - AL", "b06_biochemistry_performed",
  "AL backup PK sample prepared (BHP3)",
  "PK laboratory records",
  "Backup specimen indicator",
  
  "PK sample processing - AL", "b06_1_if_yes",
  "Volume of AL backup PK aliquot",
  "PK laboratory records",
  "Specimen processing information"
)

write_csv(SMART_supporting_dictionary, "Discrepancy data/SMART_supporting_dictionary.csv")

### Complete Dictionary 
SMART_full_dictionary <- bind_rows(
  SMART_dictionary,
  SMART_supporting_dictionary
)
write_csv(SMART_full_dictionary, "Discrepancy data/SMART_complete_data_dictionary.csv")


#========================ANALYSIS DATASET=======================================

SMART_analysis <- SMARTUncomplicated |>
  select(
    h01_3um,
    id_um,
    redcap_event_name,
    treatment_arm,
    study_day,
    study_visit,
    site,
    age = b02_1um,
    weight = b04um,
    sex = b05um,
    
    #------ Adding temperature ------#
    c03_1um,
    c03_2um_2,
    a04_1_fuv,
    a04_3_2_fuv,
    a03_1_d_um,
    a03_2_d_um,
    b01_1_tf_um_sm,
    b01_2_tf_um_sm,
    b03_1_uv,
    b03_2_uv,
    
    #------ Adding haemoglobin ------#
    d05um,
    a07_d_um,
    c03_haemoglobin,
    hgb_edta,
    
    #------- Adding RDT-------------#
    d02_1um,
    d02_2um,
    c01_malaria_rdt_done,
    c01_2_uv,
    
    #------ Adding microscopy ------#
    d02_1um_2,
    d02_2um_2,
    d02_4um,
    b03_malaria_microscopy_pos_sm,
    b03_malaria_microscopy_pos_2_sm,
    b03_malaria_microscopy_pos_3_sm,
    c02_malaria_microscopy_don,
    
    #------ Adding drug administration ------#
    # AP
    da_ap_b01_1,
    dadmin_ap,
    da_ap_b03_2,
    da_ap_b04,
    da_ap_b05,
    da_ap_b07_1,
    da_ap_b07_3,
    da_ap_c01,
    rdadmin,
    da_ap_c03_2,
    
    # AL
    da_al_b01_1,
    dadmin1,
    b03_da,
    da_al_b04,
    da_al_b05,
    da_al_b07_1,
    da_al_b07_3,
    da_al_c01,
    dadmin,
    da_al_c03,
    
    #------ Adding adverse event ------#
    a05_1_fuv,
    a04_1_d_um,
    b04_any_new_complain_since,
    
    visitp_fup,
    visitp_fupd,
    
    #---- A dding Biochemistry variables-----#
    # Biochemistry sample / processing
    samplec_heparin,
    vdate_heparin,
    hep_bio,
    blab_dat,
    blood_heparin,
    b03_biochemistry_performed,
    c01_date_and_time_of_bioch,
    
    # Biochemistry test status + results
    c03_1_bilt,
    c03_2_if_done,
    
    c04_1_bild,
    c04_1,
    
    alb_heparin,
    c05_2_if_done,
    
    alt_heparin,
    c05_2_if_done_2,
    
    c07_heparin,
    c07_2_if_done,
    
    alp_heparin,
    c08_2_if_done,
    
    ggt_heparin,
    c09_2_if_done,
    
    pro_heparin,
    c10_2_if_done,
    
    #----Adding adverse event selection fields----
    a05_2_fuv___1,
    a05_2_fuv___2,
    a05_2_fuv___3,
    a05_2_fuv___4,
    a05_2_fuv___5,
    a05_2_fuv___6,
    a05_2_fuv___7,
    a05_2_fuv___8,
    a05_2_fuv___9,
    a05_2_fuv___10,
    a05_2_fuv___11,
    a05_2_fuv___12,
    a05_2_fuv___13,
    
    # Later follow-up complaint details: H72-D42
    a04_2_d_um___1,
    a04_2_d_um___2,
    a04_2_d_um___3,
    a04_2_d_um___4,
    a04_2_d_um___5,
    a04_2_d_um___6,
    a04_2_d_um___7,
    a04_2_d_um___8,
    a04_2_d_um___9,
    a04_2_d_um___10,
    a04_2_d_um___11,
    a04_2_d_um___12,
    a04_2_d_um___13,
    
    # Unscheduled visit complaint details
    b04_2_if_yes_tick_all_syst___1,
    b04_2_if_yes_tick_all_syst___2,
    b04_2_if_yes_tick_all_syst___3,
    b04_2_if_yes_tick_all_syst___4,
    b04_2_if_yes_tick_all_syst___5,
    b04_2_if_yes_tick_all_syst___6,
    b04_2_if_yes_tick_all_syst___7,
    b04_2_if_yes_tick_all_syst___8,
    b04_2_if_yes_tick_all_syst___9,
    b04_2_if_yes_tick_all_syst___10,
    b04_2_if_yes_tick_all_syst___11,
    b04_2_if_yes_tick_all_syst___12,
    b04_2_if_yes_tick_all_syst___13,
    
    #unscheduled visit date
    a02_uv,
    a01_uv,
    a02_uv,
    
    a01_fuv,
    a01_d_um,
    
    #scheduled visit
    d01_date_of_study_end,
    d3_1_was_the_participant_fol,
    d3_2_if_no_select_reason_f,
    sremarks
  )|>
  
  mutate(
    h01_3um = na_if(as.character(h01_3um), "")
  )|>
  group_by(id_um)|>
  mutate(
    h01_3um = first(na.omit(h01_3um))
  ) |>
  ungroup()|>
  
  rename(
    participant_id = h01_3um,
    screening_id = id_um
  )|>
  
  mutate(
    across(
      c(
        age,
        weight,
        sex,
        site,
        c03_1um,
        c03_2um_2,
        a04_1_fuv,
        a04_3_2_fuv,
        a03_1_d_um,
        a03_2_d_um,
        b01_1_tf_um_sm,
        b01_2_tf_um_sm,
        b03_1_uv,
        b03_2_uv,
        d05um,
        a07_d_um,
        c03_haemoglobin,
        hgb_edta,
        d02_1um,
        d02_2um,
        d02_1um_2,
        d02_2um_2,
        d02_4um,
        b03_malaria_microscopy_pos_sm,
        b03_malaria_microscopy_pos_2_sm,
        b03_malaria_microscopy_pos_3_sm,
        c01_malaria_rdt_done,
        c01_2_uv,
        c02_malaria_microscopy_don,
        da_ap_b01_1,
        dadmin_ap,
        da_ap_b03_2,
        da_ap_b07_1,
        da_ap_b07_3,
        da_ap_c01,
        rdadmin,
        da_ap_c03_2,
        da_al_b01_1,
        dadmin1,
        b03_da,
        da_al_b07_1,
        da_al_b07_3,
        da_al_c01,
        dadmin,
        da_al_c03,
        a05_1_fuv,
        a04_1_d_um,
        b04_any_new_complain_since,
        visitp_fup,
        visitp_fupd,
        c03_2_if_done,
        c04_1,
        c05_2_if_done,
        c05_2_if_done_2,
        c07_2_if_done,
        c08_2_if_done,
        c09_2_if_done,
        c10_2_if_done,
        a01_uv,
        d3_1_was_the_participant_fol,
        d3_2_if_no_select_reason_f
      ),
      as.numeric
    ),
    
    across(
      c(
        da_ap_b04,
        da_ap_b05,
        da_al_b04,
        da_al_b05,
        a01_fuv,
        a01_d_um,
        a02_uv,
        d01_date_of_study_end
      ),
      as.character
    )
  )|>
  mutate(
    temperature = dplyr::case_when(
      study_visit == "H0" ~ c03_1um,
      study_visit %in% c("H8", "H24", "H36", "H48", "H60") ~ a04_1_fuv,
      study_visit %in% c("H72","D7","D14","D21","D28","D35","D42")~ a03_1_d_um,
      study_visit == "Treatment failure" ~ b01_1_tf_um_sm,
      study_visit == "Unscheduled visit" ~ b03_1_uv,
      TRUE ~ NA_real_
    ),
    temperature_route = dplyr::case_when(
      study_visit == "H0" ~ c03_2um_2,
      study_visit %in% c("H8", "H24", "H36", "H48", "H60") ~ a04_3_2_fuv,
      study_visit %in% c("H72","D7","D14","D21","D28","D35","D42")~ a03_2_d_um,
      study_visit == "Treatment failure" ~ b01_2_tf_um_sm,
      study_visit == "Unscheduled visit" ~ b03_2_uv,
      TRUE ~ NA_real_
    ),
    haemoglobin = dplyr::case_when(
      study_visit == "H0" ~ d05um,
      study_visit %in% c("H72","D7","D14","D21","D28","D35","D42")~ a07_d_um,
      study_visit == "Unscheduled visit" ~ c03_haemoglobin,
      TRUE ~ NA_real_
    ),
    microscopy_done = dplyr::case_when(
      study_visit == "H0" ~ d02_1um_2,
      study_visit == "Unscheduled visit" ~ c02_malaria_microscopy_don,
      TRUE ~ NA_real_
    ),
    microscopy_result = dplyr::case_when(
      study_visit == "H0" ~ d02_2um_2,
      study_visit == "Treatment failure" ~ b03_malaria_microscopy_pos_sm,
      TRUE ~ NA_real_
    ),
    rdt_done = case_when(
      study_visit == "H0" ~ d02_1um,
      study_visit == "Unscheduled visit" ~ c01_malaria_rdt_done,
      TRUE ~ NA_real_
    ),
    
    rdt_result = case_when(
      study_visit == "H0" ~ d02_2um,
      study_visit == "Unscheduled visit" ~ c01_2_uv,
      TRUE ~ NA_real_
    ),
    pf_asexual_only = dplyr::case_when(
      study_visit == "H0" ~ d02_4um,
      TRUE ~ NA_real_
    ),
    drug_given = dplyr:: case_when(
    treatment_arm == "AP" ~ da_ap_b01_1,
    treatment_arm == "AL" ~da_al_b01_1,
    TRUE ~ NA_real_
    ),
    drug_administered = dplyr::case_when(
      treatment_arm == "AP" ~dadmin_ap,
      treatment_arm == "AL" ~ dadmin1,
      TRUE ~ NA_real_
    ),
    drug_dose = dplyr::case_when(
      treatment_arm == "AP" ~ da_ap_b03_2,
      treatment_arm == "AL" ~ b03_da,
      TRUE ~ NA_real_
    ),
    drug_admin_date = dplyr::case_when(
      treatment_arm == "AP" ~ da_ap_b04,
      treatment_arm == "AL" ~ da_al_b04,
      TRUE ~ NA_character_
    ),
    drug_admin_time = dplyr::case_when(
      treatment_arm == "AP" ~ da_ap_b05,
      treatment_arm == "AL" ~ da_al_b05,
      TRUE ~ NA_character_
    ),
    vomited = dplyr::case_when(
      treatment_arm == "AP" ~ da_ap_b07_1,
      treatment_arm == "AL" ~ da_al_b07_1,
      TRUE ~ NA_real_
    ),  
    drug_readministered = dplyr::case_when(
      treatment_arm == "AP" ~ da_ap_b07_3,
      treatment_arm == "AL" ~ da_al_b07_3,
      TRUE ~ NA_real_
    ),
    new_complaint = dplyr::case_when(
      study_visit %in% c("H8", "H24", "H36", "H48", "H60") ~ a05_1_fuv,
      study_visit %in% c("H72","D7","D14","D21","D28","D35","D42")~ a04_1_d_um,
      study_visit == "Unscheduled visit" ~ b04_any_new_complain_since,
      TRUE ~ NA_real_
    ),
    unscheduled_visit = dplyr::case_when(
      !is.na(a01_uv) ~ paste0(
        "U",
        stringr::str_pad(a01_uv, width = 2, pad = "0")
      ),
      TRUE ~ NA_character_
    )
  ) |>
  mutate(
    drug_admin_date = na_if(drug_admin_date, ""),
    drug_admin_date = as.Date(drug_admin_date),
    drug_admin_time = na_if(drug_admin_time, ""),
    a02_uv = as.Date(a02_uv),
    a01_fuv = as.Date(a01_fuv),
    a01_d_um = as.Date(a01_d_um),
    a02_uv = as.Date(a02_uv),
    d01_date_of_study_end = as.Date(d01_date_of_study_end)
  ) |>
  # Fill baseline participant characteristics across all rows
  group_by(participant_id) |>
    fill(
      site,
      age,
      weight,
      sex,
      .direction = "downup"
    ) |>
    ungroup()
  
#Remove participant with missing treatment S24-20057-AG(S24-S0276-AG--1)
SMART_analysis <- SMART_analysis |>
  filter(participant_id != "S24-20057-AG")

# #Check available temperature and haemoglobin 
# SMART_analysis |>
#   group_by(study_visit) |>
#   summarise(
#     n = n(),
#     temp_available = sum(!is.na(temperature)),
#     route_available = sum(!is.na(temperature_route)),
#     hb_available = sum(!is.na(haemoglobin)),
#     hgb_edta_available = sum(!is.na(hgb_edta))
#   )

# # Check available microscopy data
# SMART_analysis |>
#   group_by(study_visit) |>
#   summarise(
#     n = n(),
#     microscopy_done_available = sum(!is.na(microscopy_done)),
#     microscopy_result_available = sum(!is.na(microscopy_result)),
#     pf_asexual_available = sum(!is.na(pf_asexual_only)),
#     tf_micro_2_available = sum(!is.na(b03_malaria_microscopy_pos_2_sm)),
#     tf_micro_3_available = sum(!is.na(b03_malaria_microscopy_pos_3_sm))
#   )

# # Check drug admin
# SMART_analysis |>
#   group_by(study_visit) |>
#   summarise(
#     n = n(),
#     ap_drug_given = sum(!is.na(da_ap_b01_1)),
#     ap_administered = sum(!is.na(dadmin_ap)),
#     al_drug_given = sum(!is.na(da_al_b01_1)),
#     al_administered = sum(!is.na(dadmin1)),
#     .groups = "drop" 
#   ) |>
#   filter(
#     ap_drug_given > 0 |
#     ap_administered > 0 |
#     al_drug_given > 0 |
#     al_administered > 0
#   )
# SMART_analysis |>
#   group_by(treatment_arm, study_visit) |>
#   summarise(
#     n = n(),
#     drug_given_available = sum(!is.na(drug_given)),
#     administered_available = sum(!is.na(drug_administered)),
#     dose_available = sum(!is.na(drug_dose)),
#     vomiting_recorded = sum(!is.na(vomited)),
#     readministration_recorded = sum(!is.na(drug_readministered)),
#     .groups = "drop"
#   )|>
#   filter(drug_given_available > 0)
# 
# SMART_analysis |>
#   filter (!is.na(vomited)) |>
#   count(treatment_arm, study_visit, vomited)
# 
# SMART_analysis |>
#   filter (!is.na(drug_readministered)) |>
#   count(treatment_arm, study_visit, drug_readministered)
# 
# SMART_analysis |>
#   filter(!is.na(new_complaint)) |>
#   count(study_visit, new_complaint)

#=======================CREATING CONCOMITANT DATASET============================

SMART_conmed <- SMARTUncomplicated |>
  mutate(
    h01_3um = na_if(as.character(h01_3um), "")
  )|>
  group_by(id_um)|>
  mutate(
    h01_3um = first(na.omit(h01_3um))
  )|>
  ungroup()|>
  
  filter(redcap_event_name == "concomitant_medica_arm_1") |>
  select(
    participant_id = h01_3um,
    screening_id = id_um,
    cm_pid_rm,
    cmno_rm,
    medication_name_rm,
    medname_rm,
    cm_ind_rm,
    cm_oth_rm,
    route_rm,
    route_oth_rm,
    sdate_rm,
    ts_rm,
    edate_rm
  )
glimpse(SMART_conmed)

SMART_conmed |>
  summarise(
    n_records = n(),
    n_participants = n_distinct(participant_id),
    participant_id_missing = sum(is.na(participant_id))
  )

SMART_conmed|>
  select(
    participant_id,
    screening_id,
    cmno_rm,
    medication_name_rm
  )|>
  slice_head(n=20)

SMART_analysis |>
  select(
    screening_id,
    treatment_arm,
    study_visit,
    drug_admin_date,
    drug_admin_time
  ) |>
  filter(!is.na(drug_admin_date)) |>
  slice_head(n = 21)


#=============================MICROSCOPY DATASET================================
microscopy_data <- readxl::read_excel(
  "Data/Malaria_microscopy_Updated_original.xlsx",
  sheet = "MalariaSlideReadingD-AllSMARTS2"
) |>
  clean_names()

dim(microscopy_data)
names(microscopy_data)

microscopy_data |>
  select(
    sample_id,
    contains("reading_date"),
    contains("detection"),
    contains("count"),
    contains("microscope_factor"),
    falciparum,
    malariae,
    ovale,
    vivax,
    contains("final_parasitemia")
  ) |>
  glimpse()

microscopy_data |>
  count(
    final_parasitemia_over_all_readings_in_parasites_m_l,
    sort = TRUE
  )

microscopy_data |>
  select(
    detection_6,
    detection_15,
    detection_29,
    detection_38,
    detection_53,
    detection_62
  )|>
  summarise(
    across(
      everything(),
      ~paste(
        sort(unique(na.omit(.))),
        collapse = ", "
      )
    )
  )
  
microscopy_data |>
  count(sample_id, sort = TRUE) |>
  slice_head(n = 20)

# extracting sample_id and visit day/hour
microscopy_data <- microscopy_data|>
  mutate(
    participant_id = str_remove(sample_id, "-[^-]+$"),
    microscopy_visit = str_extract(sample_id, "[^-]+$")
  )

microscopy_data |>
  count(microscopy_visit, sort = TRUE)
microscopy_data |>
  summarise(
    n_samples = n(),
    n_participants = n_distinct(participant_id),
    n_visits = n_distinct(microscopy_visit)
  )

SMART_microscopy <- microscopy_data |>
  select(
    sample_id,
    participant_id,
    microscopy_visit,
    reading_date_2,
    reading_date_3,
    reading_date_4,
    detection_6,
    detection_15,
    detection_29,
    detection_38,
    x3rd_reading_done,
    falciparum,
    malariae,
    ovale,
    vivax,
    final_parasitemia_over_all_readings_in_parasites_m_l
  )
glimpse(SMART_microscopy)

write_csv(SMART_microscopy, "Discrepancy data/SMART_microscopy.csv")

SMART_microscopy |>
  count(sample_id) |>
  filter(n > 1)

Microscopy_with_discrepancy <- microscopy_data |>
  filter(final_parasitemia_over_all_readings_in_parasites_m_l == -1,
         malariae == -1,
         falciparum == -1,
         ovale == -1,
         vivax == -1) |>
  select(
    sample_id,
    reading_date_2,
    reading_date_3,
    reading_date_4,
    slide_reader_5,
    detection_6,
    count_7,
    count_on_8,
    detection_15,
    counted_16,
    count_17,
    count_on_18,
    data_entry_clerk_26,
    slide_reader_28,
    detection_29,
    count_30,
    count_on_31,
    detection_38,
    counted_39,
    count_40,
    count_on_41,
    data_entry_clerk_49,
    x3rd_reading_done,
    slide_reader_52,
    detection_53,
    count_54,
    count_on_55,
    detection_62,
    counted_63,
    count_64,
    count_on_65,
    data_entry_clerk_72,
    falciparum,
    malariae,
    ovale,
    vivax,
    final_parasitemia_over_all_readings_in_parasites_m_l
  )
write_csv(Microscopy_with_discrepancy, "Discrepancy data/Microscopy_with_discrepancy.csv")

# SMART_analysis %>% group_by(redcap_event_name="hour_0_day_0_arm_1") %>% miss_var_table()
#gg_miss_var(SMART_analysis)
# gg_miss_var(SMART_analysis, facet = redcap_event_name)


