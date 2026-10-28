#===================DATA CLEANING===============================================

source("02_data_preparation.R")
library(writexl)


#=====================1.PARTICIPANT IDENTIFICATION==========================================================

# basic id check
SMART_analysis |> 
  summarise(
    n_records = n(),
    n_particpant = n_distinct(participant_id),
    n_screening_ids = n_distinct(screening_id),
    missing_participant_id = sum(is.na(participant_id)),
    missing_screening_id = sum(is.na(screening_id)),
  )

#check whether one participant ID is linked to multiple screening IDs
SMART_analysis |>
  distinct(participant_id, screening_id) |>
  count(participant_id, name = "n_screening_ids") |>
  filter(n_screening_ids > 1)

#check whether one screening ID is linked to multiple participant IDs
SMART_analysis |>
  distinct(participant_id, screening_id) |>
  count(screening_id, name = "n_participant_ids") |>
  filter(n_participant_ids >1)

#==================================2.TREATMENT ALLOCATION=======================
# Distribution of treatment arm
SMART_analysis |>
  distinct(participant_id, treatment_arm) |>
  count(treatment_arm)

#check missing treatment arm
SMART_analysis|>
  summarise(
    missing_treatment_arm = sum(is.na(treatment_arm)),
    blank_treatment_arm = sum(
      trimws(treatment_arm) == "",
      na.rm = TRUE
    )
  )

# check whether a participant has more than one treatment arm
SMART_analysis |>
  filter(!is.na(treatment_arm)) |>
  distinct(participant_id, treatment_arm) |>
  count(participant_id, name = "n_treatment_arms") |>
  filter(n_treatment_arms >1)

# Inspect all records for participant with missing treatment allocation
SMART_analysis |>
  filter(is.na(treatment_arm)) |>
  select(
    participant_id,
    screening_id,
    redcap_event_name,
    study_day,
    study_visit,
    treatment_arm,
    drug_given,
    drug_administered,
    drug_dose,
    drug_admin_date,
    drug_admin_time
  )
# Identify participant with missing treatment arm
SMART_analysis |>
  group_by(participant_id, screening_id) |>
  summarise(
    treatment_values = paste(
      sort(unique(na.omit(treatment_arm))),
      collapse = ", "
    ),
    n_missing_treatment = sum(is.na(treatment_arm)),
    .groups = "drop"
  ) |>
  filter(treatment_values == "")

#=================3.  STUDY VISITS==============================================================
SMART_analysis |>
  count(study_day, study_visit, sort = TRUE)

#check missing study day or visit
SMART_analysis |>
  summarise(
    missing_study_day = sum(is.na(study_day)),
    missing_study_visit = sum(is.na(study_visit)),
    blank_study_day = sum(trimws(study_day)=="", na.rm = TRUE),
    blank_study_vist = sum(trimws(study_visit)=="", na.rm = TRUE)
  )

#Inspect actual event mapping
SMART_analysis |>
  distinct(redcap_event_name, study_day, study_visit) |>
  arrange(study_day, study_visit)

#Check participant with repeated visit records
SMART_analysis |>
  count(
    participant_id,
    study_visit,
    name = "n_records"
  ) |>
  filter(n_records > 1) |>
  count(
    study_visit,
    name = "n_paticipant_visits"
  )|>
  arrange(desc(n_paticipant_visits)) 

#check whether each REDCap event maps to only one study day/visit
SMART_analysis |>
  distinct(
    redcap_event_name,
    study_day,
    study_visit
  ) |>
  count(
    redcap_event_name,
    name = "n_mappings"
  ) |>
  filter(n_mappings > 1)

#Visit attendance H8 to H60
visit_status_early <- SMART_analysis|>
  filter(study_visit %in% c("H8", "H24","H36","H48", "H60")) |>
  count(study_visit, visitp_fup, name = "n") |>
  arrange(study_visit, visitp_fup)
visit_status_early

#Visit attendance H72 to D42
visit_status_late <- SMART_analysis|>
  filter(
    study_visit %in% c("H72","D3","D7","D14","D21", "D28", "D35", "D42")
  ) |>
  count(study_visit, visitp_fupd, name = "n") |>
  arrange(study_visit, visitp_fupd)
visit_status_late |>
  print(n=Inf)

#Number of participant represented at each scheduled visit
SMART_analysis |>
  filter(
    study_visit %in% c(
      "H0","H8","H24","H36","H48","H60","H72","D7","D14","D21",
      "D28","D35", "D42"
    )
  ) |>
  group_by(study_visit)|>
  summarise(
    n_records = n(),
    n_participant = n_distinct(participant_id),
    .groups = "drop"
  )|>
  arrange(match(
    study_visit,
    c("H0","H8","H24","H36","H48","H60","H72","D7","D14","D21",
      "D28","D35", "D42")
  ))

# Participant with a treatment failure record
treatment_failure_participant <- SMART_analysis|>
  filter(study_visit == "Treatment failure") |>
  distinct(participant_id, screening_id)

treatment_failure_participant |>
  summarise(
    n_treatment_failure = n_distinct(participant_id)
  )

#Check participant without DAY 42   record
# Participants at baseline
baseline_participants <- SMART_analysis |>
  filter(study_visit == "H0") |>
  distinct(participant_id)

# Participants with a D42 record
d42_participants <- SMART_analysis |>
  filter(study_visit == "D42") |>
  distinct(participant_id)

# Participants without a D42 record
no_d42 <- baseline_participants |>
  anti_join(d42_participants, by = "participant_id")

no_d42 |>
  summarise(n_without_D42 = n())

#How many participants without D42 had a treatment failure?
no_d42 |>
  left_join(
    treatment_failure_participant |>
      mutate(treatment_failure = "Yes"),
    by = "participant_id"
  ) |>
    mutate(
      treatment_failure = replace_na(treatment_failure, "No")
    ) |>
   count(treatment_failure)

#Treatment arm for day 42 participant
SMART_analysis |>
  filter(study_visit == "D42") |>
  count(treatment_arm, useNA = "ifany")

#======================4.TEMPERATURE============================================

SMART_analysis |>
  summarise(
    n_available = sum(!is.na(temperature)),
    n_missing = sum(is.na(temperature)),
    minimum = min(temperature, na.rm = TRUE),
    q1 = quantile( temperature, 0.25, na.rm = TRUE),
    median = median(temperature, na.rm = TRUE),
    q3 = quantile(temperature, 0.75, na.rm =TRUE),
    maximum = max(temperature, na.rm = TRUE)
  )

#Temperature availability by study visit
SMART_analysis |>
  group_by(study_visit) |>
  summarise(
    n_records = n(),
    n_available = sum(!is.na(temperature)),
    n_missing = sum(is.na(temperature)),
    percent_missing = round(100 * n_missing / n_records, 1),
    .groups = "drop"
  )

#Inspect lowest temperature values
SMART_analysis |>
  filter(!is.na(temperature)) |>
  arrange(temperature) |>
  select(
    participant_id,
    study_visit,
    temperature,
    temperature_route
  ) |>
  slice_head(n = 10)

#Inspect highest temperature values
SMART_analysis |>
  filter(!is.na(temperature)) |>
  arrange(desc(temperature)) |>
  select(
    participant_id,
    study_visit,
    temperature,
    temperature_route
  ) |>
  slice_head(n = 10)

# Check completeness of temperature and temperature route
SMART_analysis |>
  summarise(
    temp_present_route_missing =
      sum(!is.na(temperature) & is.na(temperature_route)),
    
    temp_missing_route_present =
      sum(is.na(temperature) & !is.na(temperature_route))
  )

# Identify temperature with missing route

missing_temperature_route <- SMART_analysis |>
  filter(
    !is.na(temperature),
    is.na(temperature_route)
  ) |>
  select(
    participant_id,
    screening_id,
    study_visit,
    temperature,
    temperature_route
  )
missing_temperature_route

# Examine temperature route codes
SMART_analysis|>
  filter(!is.na(temperature_route)) |>
  count(temperature_route, sort = TRUE)

#Check missing temperature at visits where temperature is expected
missing_temperature <- SMART_analysis |>
  filter(
    study_visit %in% c(
      "H0","H8","H24","H36","H48","H60","H72","D7","D14","D21",
      "D28","D35","D42", "Treatment failure", "Unscheduled visit"
    ),
    is.na(temperature),
    is.na(temperature_route)
  ) 
missing_temperature |>
  count(study_visit, name = "n_missing_temperature")

missing_temperature |>
  count(study_visit, "n_missing_temperature") |>
  arrange(match(
    study_visit,
    c("H8","H24","H36","H48","H60","H72","D7","D14","D21",
      "D28","D35","D42", "Treatment failure", "Unscheduled visit")
  ))


# Check reason for missing temperature
missing_temperature_status <- missing_temperature |>
  mutate(
    followup_status = case_when(
      
      study_visit %in% c("H8", "H24", "H36", "H48", "H60") &
        visitp_fup == 1 ~ "Follow-up performed",
      
      study_visit %in% c("H8", "H24", "H36", "H48", "H60") &
        visitp_fup == 0 ~ "Follow-up not performed",
      
      study_visit %in% c("H8", "H24", "H36", "H48", "H60") &
        is.na(visitp_fup) ~ "Follow-up status missing",
      
      study_visit %in% c("H72", "D7", "D14", "D21",
                         "D28", "D35", "D42") &
        visitp_fupd == 1 ~ "Follow-up performed",
      
      study_visit %in% c("H72", "D7", "D14", "D21",
                         "D28", "D35", "D42") &
        visitp_fupd == 0 ~ "Follow-up not performed",
      
      study_visit %in% c("H72", "D7", "D14", "D21",
                         "D28", "D35", "D42") &
        is.na(visitp_fupd) ~ "Follow-up status missing",
      
      TRUE ~ "Not applicable"
    )
  )

missing_temperature_status |>
  count(study_visit, followup_status, name = "n") |>
  arrange(
    match(
      study_visit,
      c("H8", "H24", "H36", "H48", "H60", "H72",
        "D7", "D14", "D21", "D28", "D35", "D42",
        "Treatment failure", "Unscheduled visit")
    ),
    followup_status
  ) |>
  print(n = Inf)

#===============================================================================
#             MISSING TEMPERATURE REPORT
#===============================================================================
temperature_missingness_summary <- SMART_analysis |>
  filter(
    study_visit %in% c(
      "H8", "H24", "H36", "H48", "H60", "H72",
      "D7", "D14", "D21", "D28", "D35", "D42",
      "Treatment failure", "Unscheduled visit"
    )
  ) |>
  group_by(study_visit) |>
  summarise(
    n_records = n(),
    n_missing = sum(is.na(temperature)),
    percent_missing = round(100 * n_missing / n_records, 1),
    .groups = "drop"
  ) |>
  arrange(match(
    study_visit,
    c("H8", "H24", "H36", "H48", "H60", "H72",
      "D7", "D14", "D21", "D28", "D35", "D42",
      "Treatment failure", "Unscheduled visit")
  ))

temperature_missingness_summary

#Temperature missing despite follow-up being performed
temperature_queries <- missing_temperature_status |>
  filter(followup_status == "Follow-up performed")|>
  select(
    participant_id,
    screening_id,
    study_visit,
    temperature
  )
temperature_queries

#-----------------------PARTICIPANT WITH MISSING TEMPERATURE--------------------
participants_missing_temperature <- SMART_analysis|>
  filter(
    study_visit %in% c(
      "H8", "H24", "H36", "H48", "H60", "H72",
      "D7", "D14", "D21", "D28", "D35", "D42",
      "Treatment failure", "Unscheduled visit"
    ),
    is.na(temperature)
  )|>
  select(
    participant_id,
    screening_id,
    #treatment_arm,
    study_visit,
    temperature
  ) |>
  arrange(
    match(
      study_visit,
      c("H8", "H24", "H36", "H48", "H60", "H72",
        "D7", "D14", "D21", "D28", "D35", "D42",
        "Treatment failure", "Unscheduled visit"
      )
    ),
    participant_id
  )
participants_missing_temperature 

write_csv(participants_missing_temperature, "Discrepancy data/Participant_with_missing_temperature.csv")


# # Temperature  and follow-up status missing 
# missing_temperature_status |>
#   filter(followup_status == "Follow-up status missing") |>
#   select(
#     participant_id,
#     screening_id,
#     study_visit,
#     temperature,
#     visitp_fup,
#     visitp_fupd
#   ) |>
#   arrange(study_visit, participant_id)

#===============================================================================
# missingness_summary <- SMART_analysis |>
#   summarise(
#     across(
#       c(
#         participant_id,
#         treatment_arm,
#         study_day,
#         study_visit,
#         temperature,
#         temperature_route,
#         haemoglobin,
#         hgb_edta,
#         microscopy_done,
#         microscopy_result,
#         pf_asexual_only,
#         drug_given,
#         drug_administered,
#         drug_dose,
#         drug_admin_date,
#         drug_admin_time,
#         vomited,
#         drug_readministered,
#         new_complaint
#       ),
#       ~ sum(is.na(.))
#     )
#   )
# 
# missingness_summary <- SMART_analysis |>
#   summarise(
#     across(
#       c(
#         participant_id,
#         treatment_arm,
#         study_day,
#         study_visit,
#         temperature,
#         temperature_route,
#         haemoglobin,
#         hgb_edta,
#         microscopy_done,
#         microscopy_result,
#         pf_asexual_only,
#         drug_given,
#         drug_administered,
#         drug_dose,
#         drug_admin_date,
#         drug_admin_time,
#         vomited,
#         drug_readministered,
#         new_complaint
#       ),
#       list(
#         n_missing = ~ sum(is.na(.)),
#         n_available = ~ sum(!is.na(.))
#       )
#     )
#   ) |>
#   pivot_longer(
#     everything(),
#     names_to = c("variable", ".value"),
#     names_pattern = "(.*)_(n_missing|n_available)"
#   ) |>
#   mutate(
#     total = n_missing + n_available,
#     percent_missing = round(100 * n_missing / total, 1)
#   ) |>
#   arrange(desc(percent_missing))
# missingness_summary
#===============================================================================

#======================HAEMOGLOBIN MISSINGNESS==================================

#Check total haemoglobin available for each visit
SMART_analysis|>
  group_by(study_visit) |>
  summarise(
    n_records = n(),
    n_available = sum(!is.na(haemoglobin)),
    n_missing = sum(is.na(haemoglobin)),
    .groups = "drop"
  ) |>
  arrange(
    match(
      study_visit,
      c("H0","H8", "H24", "H36", "H48", "H60", "H72", "D7", "D14",
         "D21", "D28", "D35", "D42","Treatment failure", "Scheduled visit",
         "Unscheduled visit", "Concomitant medication"
      )
    )
  )

# Check where Haemoglobin is populated
SMART_analysis |>
  select(
    participant_id,
    study_visit,
    d05um,
    a07_d_um,
    c03_haemoglobin,
    hgb_edta,
    haemoglobin
  ) |>
  group_by(study_visit) |>
  summarise(
    n_records = n(),
    d05um_n = sum(!is.na(d05um)),
    a07_d_um_n = sum(!is.na(a07_d_um)),
    c03_haemoglobin_n = sum(!is.na(c03_haemoglobin)),
    hgb_edta_n = sum(!is.na(hgb_edta)),
    haemoglobin_n = sum(!is.na(haemoglobin)),
    .groups = "drop"
  ) |>
  arrange(
    match(
      study_visit,
      c("H0","H8", "H24", "H36", "H48", "H60", "H72", "D7", "D14",
        "D21", "D28", "D35", "D42","Treatment failure", "Scheduled visit",
        "Unscheduled visit", "Concomitant medication"
      )
    )
  ) |> 
  print(n=Inf)

#-----------Document missing haemoglobin----------------------------------------
haemoglobin_missing_summary <- SMART_analysis|>
  filter(
    study_visit %in% c("H0", "H72", "D7", "D14",
        "D21", "D28", "D35", "D42", "Unscheduled visit"
    )
  )|>
  group_by(study_visit) |>
  summarise(
    n_records = n(),
    n_missing = sum(is.na(haemoglobin)),
    percentage_missing = round(100* n_missing / n_records,1),
    .groups = "drop"
  )|>
  arrange(
    match(
      study_visit,
      c("H0", "H72", "D7", "D14",
        "D21", "D28", "D35", "D42", "Unscheduled visit"
      )
    )
  )
haemoglobin_missing_summary

#-------PATICIPANTS WITH MISSING HAEMOGLOBIN-----------------------------------
participant_missing_haemoglobin <- SMART_analysis|>
  filter(
    study_visit %in% c("H0", "H72", "D7", "D14",
                       "D28", "D42", "Unscheduled visit"
    ),
    is.na(haemoglobin)
  )|>
  select(
    participant_id,
    screening_id,
    study_visit,
    haemoglobin,
  ) |>
  arrange(
    match(
      study_visit,
      c("H0", "H72", "D7", "D14","D21", 
        "D28", "D35", "D42", "Unscheduled visit")
    ),
    participant_id
  )
participant_missing_haemoglobin

write_csv(participant_missing_haemoglobin, "Discrepancy data/Participants_with_missing_haemoglobin.csv")


#=============MALARIA MICROSCOPY MISSINGNESS=====================================

#Check where miscroscopy results are populated
SMART_analysis |>
  group_by(study_visit)|>
  summarise(
    n_records = n(),
    microscopy_done_available = sum(!is.na(microscopy_done)),
    microscopy_results_available = sum(!is.na(microscopy_result)),
    pf_asexual_available = sum(!is.na(pf_asexual_only)),
    rdt_done_available = sum(!is.na(rdt_done)),
    rdt_result_available = sum(!is.na(rdt_result)),
    .groups = "drop"
  ) |>
  arrange(
    match(
      study_visit,
      c(
        "H0", "H8", "H24", "H36", "H48", "H60", "H72",
        "D7", "D14", "D21", "D28", "D35", "D42",
        "Treatment failure", "Scheduled visit",
        "Unscheduled visit", "Concomitant medication"
      )
    )
  ) |>
  print( n = Inf)

#Missing microscopy results at treatment failure
missing_tf_microscopy_result <- SMART_analysis |>
  filter(
    study_visit == "Treatment failure",
    is.na(microscopy_result)
  ) |>
  select (
    participant_id,
    screening_id,
    study_visit,
    microscopy_result
  )
missing_tf_microscopy_result

# Missing microscopy-done information at unscheduled visit
missing_unscheduled_microscopy <- SMART_analysis |>
  filter(
    study_visit == "Unscheduled visit",
    is.na(microscopy_done)
  )|>
  select(
    participant_id,
    screening_id,
    study_visit,
    microscopy_done
  )
missing_unscheduled_microscopy

# Check RDT count
SMART_analysis |>
  filter(study_visit == "Unscheduled visit") |>
  count(
    rdt_done,
    rdt_result,
    .drop = FALSE
  ) 

# Missing RDT done at Unscheduled visits
#RDT not done and NA
participant_missing_rdt <- SMART_analysis |>
  filter(
    #(study_visit ==  "Unscheduled visit" &
    #rdt_done == 0) |
    (study_visit ==  "Unscheduled visit"&
     is.na(rdt_done)&
     is.na(rdt_result))
    
  ) |>
  select(
    participant_id,
    screening_id,
    study_visit,
    rdt_done,
    rdt_result
  )
participant_missing_rdt


#------------MICROSCOPY AND RDT MISSINGNESS SUMMARY----------------------------

microscopy_missingness_summary <- tibble(
  study_visit = c(
    "Treatment failure",
    "Unscheduled visit",
    "Unscheduled visit"
  ),
  variable = c(
    "microscopy_result",
    "microscopy_done",
    "rdt_done"
  ),
  n_records = c(124, 25, 25),
  n_missing = c(2, 7, 7)
) |>
  mutate(
    percent_missing = round(100 * n_missing / n_records, 1)
  )

microscopy_missingness_summary


#==================== DRUG ADMINISTRATION MISSINGNESS===========================

missing_drug_admin <- SMART_analysis |>
  filter(
    (treatment_arm == "AP" &
       study_visit %in% c("H0", "H24", "H48")
     ) |
    (treatment_arm == "AL" &
       study_visit %in% c("H0","H8", "H24","H36", "H48","H60"))
  )
missing_drug_admin |>
  group_by(treatment_arm, study_visit) |>
  summarise(
    n_records = n(),
    missing_drug_given        = sum(is.na(drug_given)),
    missing_drug_administered = sum(is.na(drug_administered)),
    missing_dose              = sum(is.na(drug_dose)),
    missing_date              = sum(is.na(drug_admin_date)),
    missing_time              = sum(is.na(drug_admin_time)),
    missing_vomitting         = sum(is.na(vomited)),
    #drug_readministered = sum(is.na(drug_readministered)),
    .groups = "drop"
  ) |>
  arrange(
    match(
      study_visit,
      c(
        "H0", "H8", "H24", "H36", "H48", "H60"
      )
    )
  ) |>
  print(n= Inf)

#--------Participants with missing drug admin-----
participant_missing_drug_admin <- missing_drug_admin |>
  filter(
    is.na(drug_given) |
    is.na(drug_administered) |
    is.na(drug_dose) |
    is.na(drug_admin_date) |
    is.na(drug_admin_time) |
    is.na(vomited)
  ) |>
  select(
    participant_id,
    screening_id,
    study_visit,
    drug_given,
    drug_administered,
    drug_dose,
    drug_admin_date,
    drug_admin_time,
    vomited
  )|>
  arrange(
    match(study_visit,
      c("H0", "H8", "H24", "H36", "H48", "H60")),
    participant_id
  )
participant_missing_drug_admin

# Check re-administration information where vomiting was recorded
SMART_analysis |>
  filter(
    study_visit %in% c("H0", "H8", "H24", "H36", "H48", "H60"),
    vomited == 1
  ) |>
  summarise(
    n_vomited = n(),
    n_readmin_available = sum(!is.na(drug_readministered)),
    n_readmin_missing = sum(is.na(drug_readministered))
  )

#===============BIOCHEMISTRY MISSINGNESS=======================================

SMART_analysis |>
  group_by(study_visit) |>
  summarise(
    n_records = n(),
    c03_2_n = sum(!is.na(c03_2_if_done)),
    c04_1_n = sum(!is.na(c04_1)),
    c05_2_n = sum(!is.na(c05_2_if_done)),
    c05_2_2_n = sum(!is.na(c05_2_if_done_2)),
    c07_2_n = sum(!is.na(c07_2_if_done)),
    c08_2_n = sum(!is.na(c08_2_if_done)),
    c09_2_n = sum(!is.na(c09_2_if_done)),
    c10_2_n = sum(!is.na(c10_2_if_done)),
    .groups = "drop"
  ) |>
  arrange(
    match(
      study_visit,
      c("H0", "H8", "H24", "H36", "H48", "H60", "H72",
        "D7", "D14", "D21", "D28", "D35", "D42",
        "Treatment failure", "Unscheduled visit")
    )
  ) |>
  print(n = Inf)

SMART_analysis |>
  filter(
    study_visit %in% c("H0","H72","D7","D14","D14","D28","D42","Unscheduled visit")
  ) |>
  summarise(
    n_records = n(),
    
    total_bili_done = sum(c03_1_bilt == 1, na.rm = TRUE),
    total_bili_result = sum(!is.na(c03_2_if_done)),
    
    direct_bili_done = sum(c04_1_bild == 1, na.rm = TRUE),
    direct_bili_result = sum(!is.na(c04_1)),
    
    albumin_done = sum(alb_heparin == 1, na.rm = TRUE),
    albumin_result = sum(!is.na(c05_2_if_done)),
    
    alt_done = sum(alt_heparin == 1, na.rm = TRUE),
    alt_result = sum(!is.na(c05_2_if_done_2)),
    
    ast_done = sum(c07_heparin == 1, na.rm = TRUE),
    ast_result = sum(!is.na(c07_2_if_done)),
    
    alp_done = sum(alp_heparin == 1, na.rm = TRUE),
    alp_result = sum(!is.na(c08_2_if_done)),
    
    ggt_done = sum(ggt_heparin == 1, na.rm = TRUE),
    ggt_result = sum(!is.na(c09_2_if_done)),
    
    protein_done = sum(pro_heparin == 1, na.rm = TRUE),
    protein_result = sum(!is.na(c10_2_if_done))
  ) 

#Verify available bilirubin result
# Direct bilirubin result available but test not recorded as Done
direct_bili_inconsistency <- SMART_analysis |>
  filter(
    study_visit %in% c("H0", "H72", "D7", "D14", "D28", "D42", "Unscheduled visit"),
    !is.na(c04_1),
    is.na(c04_1_bild) | c04_1_bild != 1
  ) |>
  select(
    participant_id,
    screening_id,
    study_visit,
    c04_1_bild,
    c04_1
  )
direct_bili_inconsistency

# Check biochemistry tests with missing result
biochem_missing_results <- SMART_analysis |>
  filter(
    (c03_1_bilt == 1 & is.na(c03_2_if_done)) |
      (c04_1_bild == 1 & is.na(c04_1)) |
      (alb_heparin == 1 & is.na(c05_2_if_done)) |
      (alt_heparin == 1 & is.na(c05_2_if_done_2)) |
      (c07_heparin == 1 & is.na(c07_2_if_done)) |
      (alp_heparin == 1 & is.na(c08_2_if_done)) |
      (ggt_heparin == 1 & is.na(c09_2_if_done)) |
      (pro_heparin == 1 & is.na(c10_2_if_done))
  ) |>
  select(
    participant_id,
    screening_id,
    study_visit,
    c03_1_bilt, c03_2_if_done,
    c04_1_bild, c04_1,
    alb_heparin, c05_2_if_done,
    alt_heparin, c05_2_if_done_2,
    c07_heparin, c07_2_if_done,
    alp_heparin, c08_2_if_done,
    ggt_heparin, c09_2_if_done,
    pro_heparin, c10_2_if_done
  )

biochem_missing_results

#=================ADVERSE EVENTS MISSINGNESS======================================


# Check availability of the main complaint indicator by visit
SMART_analysis |>
  group_by(study_visit) |>
  summarise(
    n_records = n(),
    
    early_complaint_n = sum(!is.na(a05_1_fuv)),
    late_complaint_n = sum(!is.na(a04_1_d_um)),
    unscheduled_complaint_n =
      sum(!is.na(b04_any_new_complain_since)),
    
    .groups = "drop"
  ) |>
  arrange(
    match(
      study_visit,
      c(
        "H0", "H8", "H24", "H36", "H48", "H60",
        "H72", "D7", "D14", "D21", "D28", "D35", "D42",
        "Treatment failure", "Unscheduled visit"
      )
    )
  ) |>
  print(n = Inf)

#Identify missing clinical complain indicators
missing_complaint <- SMART_analysis |>
  filter(
    # Early follow-up
    (study_visit %in% c("H8", "H24", "H36", "H48", "H60") &
       is.na(a05_1_fuv)) |
      
      # Later follow-up
      (study_visit %in% c("H72", "D7", "D14", "D21",
                          "D28", "D35", "D42") &
         is.na(a04_1_d_um)) |
      
      # Unscheduled visit
      (study_visit == "Unscheduled visit" &
         is.na(b04_any_new_complain_since))
  ) |>
  select(
    participant_id,
    screening_id,
    treatment_arm,
    study_visit
  ) |>
  arrange(
    match(
      study_visit,
      c("H8", "H24", "H36", "H48", "H60",
        "H72", "D7", "D14", "D21", "D28",
        "D35", "D42", "Unscheduled visit")
    ),
    participant_id
  )

missing_complaint

complaint_missingness_summary <- missing_complaint |>
  count(study_visit, name = "n_missing") |>
  arrange(
    match(
      study_visit,
      c("H8", "H24", "H36", "H48", "H60",
        "H72", "D7", "D14", "D21", "D28",
        "D35", "D42", "Unscheduled visit")
    )
  )

complaint_missingness_summary

#-----------------Checking NA
# Early follow-up: H8-H60
SMART_analysis |>
  filter(study_visit %in% c("H8", "H24", "H36", "H48", "H60")) |>
  count(study_visit, a05_1_fuv, name = "n") |>
  arrange(
    match(study_visit, c("H8", "H24", "H36", "H48", "H60")),
    a05_1_fuv
  ) |>
  print(n = Inf)

# Later follow-up: H72-D42
SMART_analysis |>
  filter(
    study_visit %in% c("H72", "D7", "D14", "D21", "D28", "D35", "D42")
  ) |>
  count(study_visit, a04_1_d_um, name = "n") |>
  arrange(
    match(
      study_visit,
      c("H72", "D7", "D14", "D21", "D28", "D35", "D42")
    ),
    a04_1_d_um
  ) |>
  print(n = Inf)

# Unscheduled visits
SMART_analysis |>
  filter(study_visit == "Unscheduled visit") |>
  count(b04_any_new_complain_since, name = "n") |>
  print(n = Inf)


#       CHECK DETAILS WHEN  COMPLAINT = YES
#    H8-H60
early_complaint_details <- SMART_analysis |>
  filter(
    study_visit %in% c("H8", "H24", "H36", "H48", "H60"),
    a05_1_fuv == 1
  ) |>
  mutate(
    n_complaints_recorded = rowSums(
      across(
        c(
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
          a05_2_fuv___13
        ),
        ~ . == 1
      ),
      na.rm = TRUE
    )
  )

early_complaint_details |>
  summarise(
    complaint_yes = n(),
    details_available = sum(n_complaints_recorded > 0),
    details_missing = sum(n_complaints_recorded == 0)
  )

#     D3-D42
late_complaint_details <- SMART_analysis |>
  filter(
    study_visit %in% c("H72", "D7", "D14", "D21", "D28", "D35", "D42"),
    a04_1_d_um == 1
  ) |>
  mutate(
    n_complaints_recorded = rowSums(
      across(
        starts_with("a04_2_d_um___"),
        ~ . == 1
      ),
      na.rm = TRUE
    )
  )

late_complaint_details |>
  summarise(
    complaint_yes = n(),
    details_available = sum(n_complaints_recorded > 0),
    details_missing = sum(n_complaints_recorded == 0)
  )

# Unscheduled visit
unscheduled_complaint_details <- SMART_analysis |>
  filter(
    study_visit == "Unscheduled visit",
    b04_any_new_complain_since == 1
  ) |>
  mutate(
    n_complaints_recorded = rowSums(
      across(
        starts_with("b04_2_if_yes_tick_all_syst___"),
        ~ . == 1
      ),
      na.rm = TRUE
    )
  )


unscheduled_complaint_details |>
  summarise(
    complaint_yes = n(),
    details_available = sum(n_complaints_recorded > 0),
    details_missing = sum(n_complaints_recorded == 0)
  )

# Identify participants with complaint reported but complaint details missing
missing_complaint_details <- late_complaint_details |>
  filter(n_complaints_recorded == 0) |>
  select(
    participant_id,
    screening_id,
    study_visit,
    a04_1_d_um,
    starts_with("a04_2_d_um___")
  
  )

missing_complaint_details

missing_complaint_status <- SMART_analysis |>
  filter(
    (study_visit %in% c("H8", "H24", "H36", "H48", "H60") &
       is.na(a05_1_fuv)) |
      
      (study_visit %in% c("H72", "D7", "D14", "D21",
                          "D28", "D35", "D42") &
         is.na(a04_1_d_um)) |
      
      (study_visit == "Unscheduled visit" &
         is.na(b04_any_new_complain_since))
  ) |>
  select(
    participant_id,
    screening_id,
    study_visit,
    a05_1_fuv,
    a04_1_d_um,
    b04_any_new_complain_since
  )
missing_complaint_status


#=============== Concomitant Medication Missingness=============================

# Check structure of the concomitant-medication dataset
names(SMART_conmed)

# Number of records and participants
SMART_conmed |>
  summarise(
    n_records = n(),
    n_participants = n_distinct(participant_id)
  )

SMART_conmed |>
  summarise(
    n_records = n(),
    
    medication_name_missing =
      sum(is.na(medication_name_rm) | medication_name_rm == ""),
    
    medication_status_missing =
      sum(is.na(medname_rm)),
    
    indication_missing =
      sum(is.na(cm_ind_rm)),
    
    route_missing =
      sum(is.na(route_rm)),
    
    start_date_missing =
      sum(is.na(sdate_rm)),
    
    start_time_missing =
      sum(is.na(ts_rm)),
    
    end_date_missing =
      sum(is.na(edate_rm))
  )

# Identify concomitant medication records with missing indication
missing_conmed_indication <- SMART_conmed |>
  filter(is.na(cm_ind_rm)) |>
  select(
    participant_id,
    screening_id,
    cmno_rm,
    medication_name_rm,
    #medname_rm,
     cm_ind_rm,
    # cm_oth_rm,
    # route_rm,
    # sdate_rm,
    # ts_rm,
    # edate_rm
  )

missing_conmed_indication |>
  print(n = Inf)

# Participant scheduled visit
visit_order <- c(
  "H0", "H8", "H24", "H36", "H48", "H60", "H72",
  "D7", "D14", "D21", "D28", "D35", "D42"
)

participant_followup <- SMART_analysis |>
  filter(study_visit %in% visit_order) |>
  distinct(participant_id, screening_id, study_visit) |>
  mutate(
    visit_number = match(study_visit, visit_order)
  ) |>
  arrange(participant_id, visit_number)
#View(participant_followup)

#Participants last follow-up day  
last_followup <- participant_followup |>
  group_by( participant_id, screening_id) |>
  slice_max(visit_number, n =1, with_ties = FALSE)|>
  ungroup() |>
  transmute(
    participant_id,
    screening_id,
    last_recorded_visit = study_visit,
    last_visit_number = visit_number
  )

# Total Number of participants with each visit day as their last recorded follow-up
last_followup |>
  count(last_recorded_visit)|>
  mutate(
    last_recorded_visit = factor(
      last_recorded_visit,
      levels = visit_order
    )
  ) |>
  arrange(last_recorded_visit)

SMART_missingness <- SMART_analysis |>
  left_join(
    last_followup |>
      select(
        participant_id,
        last_visit_number
      ),
    by = "participant_id"
  ) |>
  mutate(
    current_visit_number = match(study_visit, visit_order)
  ) 

SMART_missingness <- SMART_missingness |>
  mutate(
    within_followup = case_when(
      study_visit %in% visit_order &
        current_visit_number <= last_visit_number ~ TRUE,
      
      study_visit %in% visit_order &
        current_visit_number > last_visit_number ~ FALSE,
      
      TRUE ~ NA
    )
  )

SMART_missingness |>
  select(
    participant_id,
    study_visit,
    current_visit_number,
    last_visit_number,
    within_followup
  ) #|>
  #View()

#======================CHECKING REAL MISSINGNESS==============================
#------Removed visit days after participants stopped followup
# 1. Temperature
participants_missing_temperature_true <- participants_missing_temperature |>
  left_join(
    last_followup |>
      select(participant_id, last_visit_number),
    by = "participant_id"
  ) |>
  mutate(
    missing_visit_number = match( study_visit, visit_order)
  ) |>
  filter(
    missing_visit_number <= last_visit_number
  ) |>
  select(-last_visit_number, -missing_visit_number)

#Compare with old missing temperature
tibble(
  old_missing = nrow(participants_missing_temperature),
  true_missing = nrow(participants_missing_temperature_true),
  removed_after_followup = old_missing - true_missing
)

#2. Haemoglobin
participant_missing_haemoglobin_true <- participant_missing_haemoglobin |>
  left_join(
    last_followup |>
      select(participant_id, last_visit_number),
    by = "participant_id"
  ) |>
  mutate(
    missing_visit_number = match(study_visit, visit_order)
  ) |>
  filter(
    missing_visit_number <= last_visit_number
  ) |>
  select(-last_visit_number, -missing_visit_number)

#Compare with old missing haemoglobin
tibble(
  old_missing = nrow(participant_missing_haemoglobin),
  true_missing = nrow(participant_missing_haemoglobin_true),
  removed_after_followup = old_missing - true_missing
)

# Function to identify whether a scheduled visit was explicitly NOT performed
not_performed <- function(study_visit, visitp_fup, visitp_fupd) {
  case_when(
    study_visit %in% c("H8", "H24", "H36", "H48", "H60") &
      visitp_fup == 0 ~ TRUE,
    
    study_visit %in% c("H72", "D7", "D14", "D21", "D28", "D35", "D42") &
      visitp_fupd == 0 ~ TRUE,
    
    TRUE ~ FALSE
  )
}

#--------Removing visits that wasn't performed------------------------------
 #Temperature again
participants_missing_temperature_final <-
  participants_missing_temperature_true |>
  left_join(
    SMART_analysis |>
      select(
        participant_id,
        study_visit,
        visitp_fup,
        visitp_fupd
      ) |>
      distinct(),
    by = c("participant_id", "study_visit")
  ) |>
  mutate(
    visit_not_performed = not_performed(
      study_visit,
      visitp_fup,
      visitp_fupd
    )
  ) |>
  filter(!visit_not_performed) |>
  select(
    -visitp_fup,
    -visitp_fupd,
    -visit_not_performed
  )

#Haemoglobin again
participant_missing_haemoglobin_final <-
  participant_missing_haemoglobin_true |>
  left_join(
    SMART_analysis |>
      select(
        participant_id,
        study_visit,
        visitp_fup,
        visitp_fupd
      ) |>
      distinct(),
    by = c("participant_id", "study_visit")
  ) |>
  mutate(
    visit_not_performed = not_performed(
      study_visit,
      visitp_fup,
      visitp_fupd
    )
  ) |>
  filter(!visit_not_performed) |>
  select(
    -visitp_fup,
    -visitp_fupd,
    -visit_not_performed
  )

#Compare both old and new missingness of Temperature and Haemoglobin
tibble(
  variable = c("Temperature", "Haemoglobin"),
  
  originally_flagged = c(
    nrow(participants_missing_temperature),
    nrow(participant_missing_haemoglobin)
  ),
  
  after_last_followup = c(
    nrow(participants_missing_temperature_true),
    nrow(participant_missing_haemoglobin_true)
  ),
  
  final_missing = c(
    nrow(participants_missing_temperature_final),
    nrow(participant_missing_haemoglobin_final)
  )
)

#===================COMPLAINT STATUS: FINAL MISSINGNESS=======================

missing_complaint_status_final <- missing_complaint_status |>
  left_join(
    last_followup |>
      select(participant_id, last_visit_number),
    by = "participant_id"
  ) |>
  mutate(
    missing_visit_number = match(study_visit, visit_order)
  ) |>
  filter(
    missing_visit_number <= last_visit_number
  ) |>
  left_join(
    SMART_analysis |>
      select(
        participant_id,
        study_visit,
        visitp_fup,
        visitp_fupd
      ) |>
      distinct(),
    by = c("participant_id", "study_visit")
  ) |>
  mutate(
    visit_not_performed = not_performed(
      study_visit,
      visitp_fup,
      visitp_fupd
    )
  ) |>
  filter(!visit_not_performed) |>
  select(
    -last_visit_number,
    -missing_visit_number,
    -visitp_fup,
    -visitp_fupd,
    -visit_not_performed
  )


#=====================COMPLAINT DETAILS: FINAL MISSINGNESS======================

missing_complaint_details_final <- missing_complaint_details |>
  left_join(
    last_followup |>
      select(participant_id, last_visit_number),
    by = "participant_id"
  ) |>
  mutate(
    missing_visit_number = match(study_visit, visit_order)
  ) |>
  filter(
    missing_visit_number <= last_visit_number
  ) |>
  left_join(
    SMART_analysis |>
      select(
        participant_id,
        study_visit,
        visitp_fup,
        visitp_fupd
      ) |>
      distinct(),
    by = c("participant_id", "study_visit")
  ) |>
  mutate(
    visit_not_performed = not_performed(
      study_visit,
      visitp_fup,
      visitp_fupd
    )
  ) |>
  filter(!visit_not_performed) |>
  select(
    -last_visit_number,
    -missing_visit_number,
    -visitp_fup,
    -visitp_fupd,
    -visit_not_performed
  )

# final count of all the missingness/issues identified
missingness_final_counts <- tibble(
  variable = c(
    "Temperature",
    "Temperature route",
    "Haemoglobin",
    "Treatment failure microscopy",
    "Unscheduled microscopy",
    "RDT",
    "Drug administration",
    "Direct bilirubin inconsistency",
    "Complaint status",
    "Complaint details",
    "Concomitant medication indication"
  ),
  
  n_flagged = c(
    nrow(participants_missing_temperature_final),
    nrow(missing_temperature_route),
    nrow(participant_missing_haemoglobin_final),
    nrow(missing_tf_microscopy_result),
    nrow(missing_unscheduled_microscopy),
    nrow(participant_missing_rdt),
    nrow(participant_missing_drug_admin),
    nrow(direct_bili_inconsistency),
    nrow(missing_complaint_status_final),
    nrow(missing_complaint_details_final),
    nrow(missing_conmed_indication)
  )
)

missingness_final_counts

#=========================CALCULATING DENOMINATORS AND PERCENTAGES==============

eligible_followup <- SMART_missingness |>
  mutate(
    visit_not_performed = not_performed(
      study_visit,
      visitp_fup,
      visitp_fupd
    )
  ) |>
  filter(
    within_followup == TRUE,
    !visit_not_performed
  )

# Temperature
temperature_missingness <- eligible_followup |>
  filter(
    study_visit %in% c(
      "H8", "H24", "H36", "H48", "H60", "H72",
      "D7", "D14", "D21", "D28", "D35", "D42"
    )
  ) |>
  summarise(
    expected_n = n(),
    missing_n = sum(is.na(temperature)),
    missing_percent = round(100 * missing_n / expected_n, 1),
    .by = study_visit
  ) |>
  mutate(
    study_visit = factor(study_visit, levels = visit_order)
  ) |>
  arrange(study_visit)

temperature_missingness

#Haemoglobin
haemoglobin_missingness <- eligible_followup |>
  filter(
    study_visit %in% c(
      "H72", "D7", "D14", "D21",
      "D28", "D35", "D42"
    )
  ) |>
  summarise(
    expected_n = n(),
    missing_n = sum(is.na(haemoglobin)),
    missing_percent = round(100 * missing_n / expected_n, 1),
    .by = study_visit
  ) |>
  mutate(
    study_visit = factor(study_visit, levels = visit_order)
  ) |>
  arrange(study_visit)

haemoglobin_missingness

# Complains
complaint_missingness <- eligible_followup |>
  filter(
    study_visit %in% c(
      "H8", "H24", "H36", "H48", "H60",
      "H72", "D7", "D14", "D21",
      "D28", "D35", "D42"
    )
  ) |>
  mutate(
    complaint_status = case_when(
      study_visit %in% c("H8", "H24", "H36", "H48", "H60") ~
        as.character(a05_1_fuv),
      
      study_visit %in% c(
        "H72", "D7", "D14", "D21", "D28", "D35", "D42"
      ) ~
        as.character(a04_1_d_um),
      
      TRUE ~ NA_character_
    )
  ) |>
  summarise(
    expected_n = n(),
    missing_n = sum(is.na(complaint_status)),
    missing_percent = round(100 * missing_n / expected_n, 1),
    .by = study_visit
  ) |>
  mutate(
    study_visit = factor(study_visit, levels = visit_order)
  ) |>
  arrange(study_visit)

complaint_missingness

tibble(
  variable = c(
    "Temperature",
    "Haemoglobin",
    "Complaint status"
  ),
  
  from_final_lists = c(
    nrow(participants_missing_temperature_final),
    nrow(participant_missing_haemoglobin_final),
    nrow(missing_complaint_status_final)
  ),
  
  from_summary = c(
    sum(temperature_missingness$missing_n),
    sum(haemoglobin_missingness$missing_n),
    sum(complaint_missingness$missing_n)
  )
)

#==============Overall Missingness===========================================
main_missingness_summary <- bind_rows(
  
  temperature_missingness |>
    summarise(
      variable = "Temperature",
      expected_n = sum(expected_n),
      missing_n = sum(missing_n)
    ),
  
  haemoglobin_missingness |>
    summarise(
      variable = "Haemoglobin",
      expected_n = sum(expected_n),
      missing_n = sum(missing_n)
    ),
  
  complaint_missingness |>
    summarise(
      variable = "Complaint status",
      expected_n = sum(expected_n),
      missing_n = sum(missing_n)
    )
  
) |>
  mutate(
    missing_percent = round(100 * missing_n / expected_n, 1)
  )

main_missingness_summary

#====================Retain by visit==========================
drug_given_missingness <- missing_drug_admin |>
  group_by(study_visit) |>
  summarise(
    expected_n = n(),
    missing_n = sum(is.na(drug_given)),
    missing_percent = round(100 * missing_n / expected_n, 1),
    .groups = "drop"
  )

drug_dose_missingness <- missing_drug_admin |>
  filter(drug_given == 1) |>
  group_by(study_visit) |>
  summarise(
    expected_n = n(),
    missing_n = sum(is.na(drug_dose)),
    missing_percent = round(100 * missing_n / expected_n, 1),
    .groups = "drop"
  )

#===========================MISSINGNESS OF KEY VARIABLES BY VISIT===============

missingness_by_visit <- bind_rows(
  
  temperature_missingness |>
    mutate(variable = "Temperature"),
  
  haemoglobin_missingness |>
    mutate(variable = "Haemoglobin"),
  
  drug_given_missingness |>
    mutate(variable = "Study drug given"),
  
  drug_dose_missingness |>
    mutate(variable = "Drug dose")
  
) |>
  select(
    variable,
    study_visit,
    expected_n,
    missing_n,
    missing_percent
  ) |>
  mutate(
    study_visit = factor(
      study_visit,
      levels = c(
        "H0", "H8", "H24", "H36", "H48", "H60",
        "H72", "D7", "D14", "D21", "D28", "D35", "D42"
      )
    )
  ) |>
  arrange(variable, study_visit) |>
  mutate(study_visit = as.character(study_visit))

#View(missingness_by_visit)

# Adding the other variables
#=========================OTHER MISSINGNESS SUMMARY=============================

other_missingness_summary <- bind_rows(
  
  #-----------------------------------------------------------------------------
  # Temperature route
  # Expected whenever temperature was actually recorded
  #-----------------------------------------------------------------------------
  
  SMART_analysis |>
    filter(!is.na(temperature)) |>
    summarise(
      variable = "Temperature route",
      expected_n = n(),
      missing_n = sum(is.na(temperature_route))
    ),
  
  
  #-----------------------------------------------------------------------------
  # Treatment-failure microscopy
  #-----------------------------------------------------------------------------
  
  SMART_analysis |>
    filter(study_visit == "Treatment failure") |>
    summarise(
      variable = "Treatment failure microscopy",
      expected_n = n(),
      missing_n = sum(is.na(microscopy_result))
    ),
  
  
  #-----------------------------------------------------------------------------
  # Unscheduled microscopy
  #-----------------------------------------------------------------------------
  
  SMART_analysis |>
    filter(study_visit == "Unscheduled visit") |>
    summarise(
      variable = "Unscheduled microscopy",
      expected_n = n(),
      missing_n = sum(is.na(microscopy_done))
    ),
  
  
  #=============================================================================
  # DRUG ADMINISTRATION
  #=============================================================================
  
  #-----------------------------------------------------------------------------
  # Study drug given
  # Expected at every scheduled dosing visit
  #-----------------------------------------------------------------------------
  
  missing_drug_admin |>
    summarise(
      variable = "Study drug given",
      expected_n = n(),
      missing_n = sum(is.na(drug_given))
    ),
  
  
  # #-----------------------------------------------------------------------------
  # # Drug administered
  # # Expected only when study drug was recorded as given
  # #-----------------------------------------------------------------------------
  # 
  # missing_drug_admin |>
  #   filter(drug_given == 1) |>
  #   summarise(
  #     variable = "Drug administered",
  #     expected_n = n(),
  #     missing_n = sum(is.na(drug_administered))
  #   ),
  # 
  
  #-----------------------------------------------------------------------------
  # Drug dose
  # Expected only when study drug was recorded as given
  #-----------------------------------------------------------------------------
  
  missing_drug_admin |>
    filter(drug_given == 1) |>
    summarise(
      variable = "Drug dose",
      expected_n = n(),
      missing_n = sum(is.na(drug_dose))
    ),
  
  
  # #-----------------------------------------------------------------------------
  # # Drug administration date
  # # Expected only when study drug was recorded as given
  # #-----------------------------------------------------------------------------
  # 
  # missing_drug_admin |>
  #   filter(drug_given == 1) |>
  #   summarise(
  #     variable = "Drug administration date",
  #     expected_n = n(),
  #     missing_n = sum(is.na(drug_admin_date))
  #   ),
  # 
  # 
  # #-----------------------------------------------------------------------------
  # # Drug administration time
  # # Expected only when study drug was recorded as given
  # #-----------------------------------------------------------------------------
  # 
  # missing_drug_admin |>
  #   filter(drug_given == 1) |>
  #   summarise(
  #     variable = "Drug administration time",
  #     expected_n = n(),
  #     missing_n = sum(is.na(drug_admin_time))
  #   ),
  # 
  # 
  # #-----------------------------------------------------------------------------
  # # Vomiting after drug administration
  # # Expected only when study drug was recorded as given
  # #-----------------------------------------------------------------------------
  # 
  # missing_drug_admin |>
  #   filter(drug_given == 1) |>
  #   summarise(
  #     variable = "Vomiting after drug administration",
  #     expected_n = n(),
  #     missing_n = sum(is.na(vomited))
  #   ),
  # 
  # 
  # #-----------------------------------------------------------------------------
  # # Drug re-administration
  # # Expected only when vomiting was recorded
  # #-----------------------------------------------------------------------------
  # 
  # missing_drug_admin |>
  #   filter(
  #     drug_given == 1,
  #     vomited == 1
  #   ) |>
  #   summarise(
  #     variable = "Drug re-administration after vomiting",
  #     expected_n = n(),
  #     missing_n = sum(is.na(drug_readministered))
  #   ),
  # 
  
  #-----------------------------------------------------------------------------
  # Complaint details
  # Expected only when a complaint was reported
  #-----------------------------------------------------------------------------
  
  eligible_followup |>
    filter(
      study_visit %in% c(
        "H72", "D7", "D14", "D21", "D28", "D35", "D42"
      ),
      a04_1_d_um == 1
    ) |>
    mutate(
      n_complaints_recorded = rowSums(
        across(
          starts_with("a04_2_d_um___"),
          ~ . == 1
        ),
        na.rm = TRUE
      )
    ) |>
    summarise(
      variable = "Complaint details",
      expected_n = n(),
      missing_n = sum(n_complaints_recorded == 0)
    ),
  
  
  #-----------------------------------------------------------------------------
  # Concomitant medication indication
  #-----------------------------------------------------------------------------
  
  SMART_conmed |>
    summarise(
      variable = "Concomitant medication indication",
      expected_n = n(),
      missing_n = sum(is.na(cm_ind_rm))
    )
  
) |>
  mutate(
    missing_percent = round(
      100 * missing_n / expected_n,
      1
    )
  )

other_missingness_summary

final_missingness_summary <- bind_rows(
  main_missingness_summary,
  other_missingness_summary
) |>
  mutate(
    missing_percent = case_when(
      missing_n > 0 & missing_percent == 0 ~
        round(100 * missing_n / expected_n, 2),
      TRUE ~ missing_percent
    )
  ) |>
  rename(
    `Variable` = variable,
    `Expected records (N)` = expected_n,
    `Missing frequency (n)` = missing_n,
    `Missing (%)` = missing_percent
  )

final_missingness_summary

print(missingness_by_visit, n = Inf)

#-----------RDT
SMART_analysis |>
  filter(study_visit =="Unscheduled visit") |>
  count(rdt_done,rdt_result)

rdt_summary <- SMART_analysis |>
  filter(study_visit == "Unscheduled visit") |>
  summarise(
    variable = "RDT status",
    expected_n = n(),
    missing_n = sum( is.na(rdt_done) & is.na(rdt_result)),
    missing_percent = round(100* missing_n / expected_n, 1)
  )
rdt_summary

#-- Bilirubin
direct_bilirubin_check <- SMART_analysis |>
  filter(
    !is.na(c04_1) & is.na(c04_1_bild)
  ) |>
  select(
    participant_id,
    screening_id,
    study_visit,
    c04_1_bild,
    c04_1
  )

direct_bilirubin_check

SMART_analysis |>
  filter(
    c04_1_bild == 1,
    is.na(c04_1)
  ) |>
  select(
    participant_id,
    screening_id,
    study_visit,
    c04_1_bild,
    c04_1
  )

# Adding RDT summary to the final missingness summary
final_missingness_summary <- bind_rows(
  final_missingness_summary,
  rdt_summary |>
    rename(
      `Variable` = variable,
      `Expected records (N)` = expected_n,
      `Missing frequency (n)` = missing_n,
      `Missing (%)` = missing_percent
    )
)

final_missingness_summary


#===============================================================================
#        EXCEL FOR REAL MISSINGNESS
#===============================================================================

excel_na <- function(data) {
  data |>
    mutate(
      across(
        everything(),
        ~ replace_na(as.character(.), "NA")
      )
    )
}

write_xlsx(
  list(
    
    # Corrected scheduled-visit missingness
    Temperature =
      excel_na(participants_missing_temperature_final),
    
    Temp_Route =
      excel_na(missing_temperature_route),
    
    Haemoglobin =
      excel_na(participant_missing_haemoglobin_final),
    
    # Event-specific records
    TF_Microscopy =
      excel_na(missing_tf_microscopy_result),
    
    Unscheduled_Microscopy =
      excel_na(missing_unscheduled_microscopy),
    
    RDT =
      excel_na(participant_missing_rdt),
    
    # Treatment administration
    Drug_Admin =
      excel_na(participant_missing_drug_admin),
    
    # Data consistency issue
    Direct_Bilirubin =
      excel_na(direct_bili_inconsistency),
    
    # Corrected scheduled follow-up missingness
    Complaint_Status =
      excel_na(missing_complaint_status_final),
    
    Complaint_Details =
      excel_na(missing_complaint_details_final),
    
    # Concomitant medication
    Conmed_Indication =
      excel_na(missing_conmed_indication)
  ),
  
  path = "Discrepancy data/SMART_Final_Missing_Data.xlsx"
)


#===============================================================================
#        EXCEL FOR REAL MISSINGNESS - REDCap LABELS
#===============================================================================

excel_na <- function(data) {
  data |>
    mutate(
      across(
        everything(),
        ~ replace_na(as.character(.), "NA")
      )
    )
}

# 1. TEMPERATURE
temperature_excel <- participants_missing_temperature_final |>
  mutate(
    `REDCap Field Label` = case_when(
      study_visit %in% c("H8", "H24", "H36", "H48", "H60") ~
        "A04.1. Temperature",
      
      study_visit %in% c("H72", "D7", "D14", "D21", "D28", "D35", "D42") ~
        "A03.1. Temperature",
      
      study_visit == "Treatment failure" ~
        "B01.1. Temperature",
      
      study_visit == "Unscheduled visit" ~
        "B03.1. Temperature",
      
      TRUE ~ NA_character_
    )
  ) |>
  select(
    `Participant ID` = participant_id,
    `Screening ID` = screening_id,
    `Study Visit` = study_visit,
    `REDCap Field Label`,
    `Recorded Value` = temperature
  )
# 2. TEMPERATURE ROUTE
temperature_route_excel <- missing_temperature_route |>
  mutate(
    `REDCap Field Label` = case_when(
      study_visit %in% c("H8", "H24", "H36", "H48", "H60") ~
        "A04.2. Temperature route",
      
      study_visit %in% c("H72", "D7", "D14", "D21", "D28", "D35", "D42") ~
        "A03.2. Temperature route:",
      
      study_visit == "Treatment failure" ~
        "B01.2. Temperature route",
      
      study_visit == "Unscheduled visit" ~
        "B03.2. Temperature route",
      
      TRUE ~ NA_character_
    )
  ) |>
  select(
    `Participant ID` = participant_id,
    `Screening ID` = screening_id,
    `Study Visit` = study_visit,
    `Temperature` = temperature,
    `REDCap Field Label`,
    `Recorded Value` = temperature_route
  )

# 3. HAEMOGLOBIN
haemoglobin_excel <- participant_missing_haemoglobin_final |>
  mutate(
    `REDCap Field Label` = case_when(
      study_visit == "H0" ~
        "C03. Haemoglobin",
      
      study_visit %in% c("H72", "D7", "D14", "D28", "D42") ~
        "A08. Haemoglobin:",
      
      study_visit == "Unscheduled visit" ~
        "C03. Haemoglobin",
      
      TRUE ~ NA_character_
    )
  ) |>
  select(
    `Participant ID` = participant_id,
    `Screening ID` = screening_id,
    `Study Visit` = study_visit,
    `REDCap Field Label`,
    `Recorded Value` = haemoglobin
  )

# 4. TREATMENT FAILURE MICROSCOPY
tf_microscopy_excel <- missing_tf_microscopy_result |>
  transmute(
    `Participant ID` = participant_id,
    `Screening ID` = screening_id,
    `Study Visit` = study_visit,
    `REDCap Field Label` = "B03.2. Malaria microscopy positive?",
    `Recorded Value` = microscopy_result
  )

# 5. UNSCHEDULED-VISIT MICROSCOPY
unscheduled_microscopy_excel <- missing_unscheduled_microscopy |>
  transmute(
    `Participant ID` = participant_id,
    `Screening ID` = screening_id,
    `Study Visit` = study_visit,
    `REDCap Field Label` = "C02.1. Malaria microscopy done?",
    `Recorded Value` = microscopy_done
  )

# 6. RDT
rdt_excel <- participant_missing_rdt |>
  transmute(
    `Participant ID` = participant_id,
    `Screening ID` = screening_id,
    `Study Visit` = study_visit,
    `C01.1. Malaria RDT done?` = rdt_done,
    `C01.2. →If yes, results` = rdt_result
  )

# 7. DRUG ADMINISTRATION
drug_admin_excel <- participant_missing_drug_admin |>
  transmute(
    `Participant ID` = participant_id,
    `Screening ID` = screening_id,
    `Study Visit` = study_visit,
    `B01.1. Study drug given?` = drug_given,
    `B03.1. Study drug administered?` = drug_administered,
    `B03.2. Dose to be administered` = drug_dose,
    `B04. Date of drug administration:` = drug_admin_date,
    `B05. Time of drug administration/Start of observation:` = drug_admin_time,
    `Vomiting after drug administration` = vomited
  )

# 8. DIRECT BILIRUBIN CONSISTENCY

direct_bilirubin_excel <- direct_bili_inconsistency |>
  transmute(
    `Participant ID` = participant_id,
    `Screening ID` = screening_id,
    `Study Visit` = study_visit,
    `C04.1. Direct bilirubin (BIL-D), Value:` = c04_1_bild,
    `C04.2. If done:` = c04_1
  )

# 9. COMPLAINT STATUS
complaint_status_excel <- missing_complaint_status_final |>
  transmute(
    `Participant ID` = participant_id,
    `Screening ID` = screening_id,
    `Study Visit` = study_visit,
    `A05.1. Any new complain since last visit` = a05_1_fuv,
    `A04.1. Any new complain since last visit` = a04_1_d_um,
    `B04.1. Any new complain since last visit` =
      b04_any_new_complain_since
  )

# 10. COMPLAINT DETAILS
complaint_details_excel <- missing_complaint_details_final |>
  rename(
    `Participant ID` = participant_id,
    `Screening ID` = screening_id,
    `Study Visit` = study_visit,
    `A04.1. Any new complain since last visit` = a04_1_d_um
  )

# 11. CONCOMITANT MEDICATION INDICATION
conmed_indication_excel <- missing_conmed_indication |>
  transmute(
    `Participant ID` = participant_id,
    `Screening ID` = screening_id,
    `Medication Number` = cmno_rm,
    `Medication Name` = medication_name_rm,
    `Indication:` = cm_ind_rm
  )

#===============================================================================
# WRITE FINAL EXCEL WORKBOOK
#===============================================================================

write_xlsx(
  list(
    Temperature =
      excel_na(temperature_excel),
    
    # Temp_Route =
    #   excel_na(temperature_route_excel),
    # 
    Haemoglobin =
      excel_na(haemoglobin_excel),
    
    TF_Microscopy =
      excel_na(tf_microscopy_excel),
    
    Unscheduled_Microscopy =
      excel_na(unscheduled_microscopy_excel),
    
    RDT =
      excel_na(rdt_excel),
    
    Drug_Admin =
      excel_na(drug_admin_excel),
    
    Direct_Bilirubin =
      excel_na(direct_bilirubin_excel),
    
    Complaint_Status =
      excel_na(complaint_status_excel),
    
    Complaint_Details =
      excel_na(complaint_details_excel),
    
    Conmed_Indication =
      excel_na(conmed_indication_excel)
  ),
  
  path = "Discrepancy data/SMART_Final_Missing_Data1.xlsx"
)



# #====================DIVED INTO LOST TO FOLLOW UP=================XXXXX
# #--------Check participants without day 42 record for treatment failure
# #IDs of participants who have a treatment-failure(tf) record
# tf_participants <- SMART_analysis |>
#   filter(study_visit == "Treatment failure") |>
#   distinct(participant_id) |>
#   pull(participant_id)
# 
# #participant whose last schedule visit was beofre D42
# early_end_followup <- last_followup |>
#   mutate(
#     treatment_failure = if_else(
#       participant_id %in% tf_participants,
#       "Yes",
#       "No"
#     )
#   )
# 
# #overall counts
# early_end_followup |>
#   count(treatment_failure)
# 
# #Check if the total number of participant with each visit day as their
# #last follow-up had treatment failure
# early_end_followup |>
#   count(last_recorded_visit, treatment_failure) |>
#   mutate(
#     last_recorded_visit = factor(
#       last_recorded_visit,
#       levels = visit_order
#     )
#   ) |>
#   arrange(last_recorded_visit, treatment_failure)
# 
# #Check participant with no D42 and treatment failure record (lost to follow-up (ltfu))
# potential_ltfu <- last_followup |>
#    filter(
#      last_recorded_visit != "D42",
#      !participant_id %in% tf_participants,
#    )
# View(potential_ltfu)
# 
# # Records for ltfu
# potential_ltfu_records <- SMART_analysis |>
#   filter(participant_id %in% potential_ltfu $participant_id) |>
#   select(
#     participant_id,
#     screening_id,
#     study_visit,
#     redcap_event_name,
#     everything()
#   ) |>
#   arrange(participant_id)
# View(potential_ltfu_records)
# 
# 
# #----Checking for any form of record or reason for ltfu
# 
# SMART_analysis |>
#   select(visitp_fup, visitp_fupd) |>
#   distinct()
# 
# SMART_analysis |>
#   count(study_visit, visitp_fup, visitp_fupd) |>
#   arrange(study_visit)
# 
# SMART_analysis |>
#   filter(participant_id %in% potential_ltfu$participant_id) |>
#   select(
#     participant_id,
#     screening_id,
#     study_visit,
#     visitp_fup,
#     visitp_fupd,
#     #temperature
#   ) |>
#   arrange(
#     match(participant_id, potential_ltfu$participant_id),
#     match(study_visit, visit_order)
#   )
# 
# 
# SMARTUncomplicated |>
#   select(
#     id_um,
#     redcap_event_name,
#     visitp_fup,
#     vperformed_fup,
#     others_fupv,
#     visitp_fupd,
#     notperform_fupd,
#     sreason
#   )|>
#   filter(id_um %in% potential_ltfu$screening_id)|>
#   View()
#   
# SMARTUncomplicated|>
#   select(
#     visitp_fup,
#     vperformed_fup,
#     others_fupv,
#     visitp_fup,
#     notperform_fupd,
#     sreason,
#   ) |>
#   summarise(
#     across(
#       everything(),
#       ~paste(sort(unique(na.omit(.))), collapse = "")
#     )
#   )
# 
# #Checking visit by visit for 19 participants without D42 and treatment failure record
# ltfu_check <- SMARTUncomplicated |>
#   filter(id_um %in% potential_ltfu$screening_id) |>
#   filter(
#     redcap_event_name %in% c(
#       "hour_0_day_0_arm_1",
#       "hour_8_day_0_arm_1",
#       "hour_24_day_1_arm_1",
#       "hour_36_day_1_arm_1",
#       "hour_48_day_2_arm_1",
#       "hour_60_day_2_arm_1",
#       "hour_72_day_3_arm_1",
#       "day_7_arm_1",
#       "day_14_arm_1",
#       "day_21_arm_1",
#       "day_28_arm_1",
#       "day_35_arm_1",
#       "day_42_arm_1"
#     )
#   ) |>
#   select(
#     id_um,
#     redcap_event_name,
#     visitp_fup,
#     vperformed_fup,
#     others_fupv,
#     visitp_fupd,
#     notperform_fupd,
#     sreason
#   )
# 
# ltfu_check <- ltfu_check |>
#   mutate(
#     study_visit = case_when(
#       redcap_event_name == "hour_0_day_0_arm_1"  ~ "H0",
#       redcap_event_name == "hour_8_day_0_arm_1"  ~ "H8",
#       redcap_event_name == "hour_24_day_1_arm_1" ~ "H24",
#       redcap_event_name == "hour_36_day_1_arm_1" ~ "H36",
#       redcap_event_name == "hour_48_day_2_arm_1" ~ "H48",
#       redcap_event_name == "hour_60_day_2_arm_1" ~ "H60",
#       redcap_event_name == "hour_72_day_3_arm_1" ~ "H72",
#       redcap_event_name == "day_7_arm_1"  ~ "D7",
#       redcap_event_name == "day_14_arm_1" ~ "D14",
#       redcap_event_name == "day_21_arm_1" ~ "D21",
#       redcap_event_name == "day_28_arm_1" ~ "D28",
#       redcap_event_name == "day_35_arm_1" ~ "D35",
#       redcap_event_name == "day_42_arm_1" ~ "D42"
#     )
#   ) |>
#   select(
#     id_um,
#     study_visit,
#     visitp_fup,
#     vperformed_fup,
#     others_fupv,
#     visitp_fupd,
#     notperform_fupd,
#     sreason
#   ) |>
#   arrange(
#     id_um,
#     match(study_visit, visit_order)
#   )
# 
# View(ltfu_check)
# 
# #Check participant last followup and reason
# potential_ltfu_final <- ltfu_check |>
#   inner_join(
#     potential_ltfu |>
#       select(screening_id,last_recorded_visit),
#     by = c("id_um" = "screening_id")
#   )|>
#   filter(study_visit == last_recorded_visit) |>
#   select(
#     id_um,
#     last_recorded_visit,
#     visitp_fup,
#     vperformed_fup,
#     others_fupv,
#     visitp_fupd,
#     notperform_fupd,
#     sreason
#   )
# 
# View(potential_ltfu_final)

# View(SMARTUncomplicated)
# SMARTUncomplicated |> 
#   filter(
#     redcap_event_name == "unsheduled_visit_arm_1",
#     is.na(c02_malaria_microscopy_don)
#     ) |> 
#   select(
#     id_um,
#     c02_malaria_microscopy_don
#   )
