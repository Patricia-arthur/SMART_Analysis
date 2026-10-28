#===============================================================================
#==========SMART ANALYSIS-READY DATASET==========================================

source("03_data_cleaning.R")
library(readxl)
#===================1. MICROSCOPY DATA==========================================


#-------------------------------------------------------------------------------
#       Prepare microscopy data
#-------------------------------------------------------------------------------
a <- SMART_microscopy |>
  filter(
    microscopy_visit %in% c(
      "H8", "H08", "D3","D03","D7","D07"
    )
  ) |>
  group_by(microscopy_visit) |>
  summarise(
    n_records = n(),
    
    final_parasite_available =
      sum(!is.na(final_parasitemia_over_all_readings_in_parasites_m_l)),
    
    final_parasite_zero =
      sum((final_parasitemia_over_all_readings_in_parasites_m_l == 0),
          na.rm = TRUE),
    
    final_parasite_discrepancy = sum(
      final_parasitemia_over_all_readings_in_parasites_m_l== -1,
      na.rm = TRUE
    ),
    
    final_parasite_valid =
      sum(
        final_parasitemia_over_all_readings_in_parasites_m_l >= 0,
        na.rm = TRUE
      ),
    .groups = "drop"
  )

#------------Check participant overlap between duplicate visit formats----------
visit_overlap <- bind_rows(
  
  SMART_microscopy |>
    filter(microscopy_visit %in% c("H8", "H08")) |>
    distinct(participant_id, microscopy_visit) |>
    count(participant_id) |>
    summarise(
      visit = "H8/H08",
      participants = n(),
      in_both = sum(n>1)
      
    ),
  
  SMART_microscopy |>
    filter(microscopy_visit %in% c("D3", "D03")) |>
    distinct(participant_id, microscopy_visit) |>
    count(participant_id) |>
    summarise(
      visit = "D3/D03",
      participants = n(),
      in_both = sum(n>1)
    ),
  
  SMART_microscopy |>
    filter(microscopy_visit %in% c("D7", "D07")) |>
    distinct(participant_id, microscopy_visit) |>
    count(participant_id) |>
    summarise(
      visit = "D7/D07",
      participants = n(),
      in_both = sum(n>1)
    )
)
visit_overlap

#-------------Check if participants have duplicate record for a visit
SMART_microscopy |>
  filter(
    microscopy_visit %in% c(
      "H8", "H08",
      "D3", "D03",
      "D7", "D07"
    )
  ) |>
  count(participant_id, microscopy_visit) |>
  filter(n > 1) |>
  arrange(participant_id, microscopy_visit)

#=================HARMONISE MICROSCOPY VISITS:==================================
#         
# choose visit with valid final microscopy in the case of H8/H08,D3/D03,D7/D07,
#===============================================================================
# Map the days to align with the study days

SMART_microscopy_clean <- SMART_microscopy |>
  mutate(
    microscopy_visit_original = microscopy_visit,
    
    microscopy_visit = case_when(
      microscopy_visit %in% c("H8", "H08") ~ "H8",
      microscopy_visit %in% c("D3", "D03") ~ "H72",
      microscopy_visit %in% c("D7", "D07") ~ "D7",
      TRUE ~ microscopy_visit
      
    )
  )

# Give each record a final-parasitaemia quality status
SMART_microscopy_clean <- SMART_microscopy_clean |>
  mutate(
    final_parasite_status = case_when(
      final_parasitemia_over_all_readings_in_parasites_m_l >= 0 ~
        "Valid",
      
      final_parasitemia_over_all_readings_in_parasites_m_l %in% c(-1, -2) ~
        "Discrepancy",
      
      is.na(final_parasitemia_over_all_readings_in_parasites_m_l) ~
        "Missing",
      
      TRUE ~ "Other"
    ),
    
    final_parasite_priority = case_when(
      final_parasite_status == "Valid" ~ 1,
      final_parasite_status == "Discrepancy" ~ 2,
      final_parasite_status == "Missing" ~ 3,
      TRUE ~ 4
    )
  )

#Check number of participants with record in both visits (H8/H08,D3/D03,D7/D07)
SMART_microscopy_clean |>
  count(participant_id, microscopy_visit) |>
  filter(n > 1) |>
  count(microscopy_visit, name = "participants_with_multiple_records")

# Check conflicting records for final parasite on visit days
SMART_microscopy_clean |>
  filter(microscopy_visit %in% c("H8", "H72", "D7")) |>
  group_by(participant_id, microscopy_visit) |>
  filter(n() > 1) |>
  summarise(
    statuses = paste(
      microscopy_visit_original,
      final_parasite_status,
      sep = ": ",
      collapse = " | "
    ),
    .groups = "drop"
  ) |>
  count(microscopy_visit, statuses, sort = TRUE) |>
  print(n = Inf)

#----------------Check final parasite for valid duplicate records--------------

# check if participants with duplicate records for a visit have the same final 
#parasite or different value in both visits

valid_duplicate_check <- SMART_microscopy_clean |>
  filter(
    microscopy_visit %in% c("H72", "D7"),
    final_parasite_status == "Valid"
  ) |>
  group_by(participant_id, microscopy_visit) |>
  filter(n() > 1) |>
  summarise(
    n_records = n(),
    n_different_values =
      n_distinct(final_parasitemia_over_all_readings_in_parasites_m_l),
    
    values = paste(
      microscopy_visit_original,
      final_parasitemia_over_all_readings_in_parasites_m_l,
      sep = " = ",
      collapse = " | "
    ),
    
    .groups = "drop"
  ) |>
  mutate(
    agreement = if_else(
      n_different_values == 1,
      "Same",
      "Different"
    )
  )

valid_duplicate_check |>
  count(microscopy_visit, agreement)


#----------------Resolve duplicate visit format---------------------------------


SMART_microscopy_clean <- SMART_microscopy_clean |>
  arrange(
    participant_id,
    microscopy_visit,
    final_parasite_priority
  ) |>
  group_by(participant_id, microscopy_visit) |>
  slice_head(n = 1) |>
  ungroup()

SMART_microscopy_clean |>
  count(participant_id, microscopy_visit) |>
  filter(n > 1)

#-------------------------------------------------------------------------------
#-------------Create analysis variable for final parasite-----------------------
#-------------------------------------------------------------------------------

SMART_microscopy_clean <- SMART_microscopy_clean |>
  mutate(
    
    parasite_density = case_when(
      final_parasitemia_over_all_readings_in_parasites_m_l >= 0 ~
        final_parasitemia_over_all_readings_in_parasites_m_l,
      
      TRUE ~ NA_real_
    )
  )

SMART_microscopy_clean |>
  count(microscopy_visit, sort = TRUE) |>
  print(n = Inf)

SMART_microscopy_clean |>
  summarise(
    n_records = n(),
    n_participants = n_distinct(participant_id),
    valid_parasite_density = sum(!is.na(parasite_density)),
    zero_parasitaemia = sum(parasite_density == 0, na.rm = TRUE),
    positive_parasitaemia = sum(parasite_density > 0, na.rm = TRUE),
    unresolved_discrepancies =
      sum(final_parasite_status == "Discrepancy")
  )


#=================== UNUSUAL SAMPLE IDS=========================================
#            
#===============================================================================

unusual_id <- SMART_microscopy_clean|>
  filter((microscopy_visit == "AG")|
           (microscopy_visit == "H18")|
           (microscopy_visit == "08")|
           (microscopy_visit == "A0")|
           (microscopy_visit == "AFD42")|
           (microscopy_visit == "HO8"))|>
  select(
    sample_id,
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
  ) |>
  print(N = Inf)
  #View()
#write_csv(unusual_id,"unsual_microscopy_sample_id.csv")
#===============================================================================


#-------------------------------------------------------------------------------
#           Inspect detection fields
#-------------------------------------------------------------------------------
SMART_microscopy_clean |>
  summarise(
    across(
      c(
        detection_6,
        detection_15,
        detection_29,
        detection_38
      ),
      ~paste(sort(unique(na.omit(.))), collapse = ", ")  
      
    )
  )

SMART_microscopy_clean |>
  select(
    detection_6,
    detection_15,
    detection_29,
    detection_38,
    final_parasitemia_over_all_readings_in_parasites_m_l,
    parasite_density
  ) |>
  filter(
    !is.na(detection_6)|
      !is.na(detection_15)|
      !is.na(detection_29)|
      !is.na(detection_38)
  )|>
  head(30)


#---------Create microscopy analysis variable-----------------------------------
SMART_microscopy_clean <- SMART_microscopy_clean |>
  mutate(
    # 0 = no parasites detected
    # 1 = parasites detected
    # NA = final parasite density currently unresolved/unavailable
    
    microscopy_positive = case_when(
      parasite_density == 0 ~ 0,
      parasite_density > 0 ~ 1,
      TRUE ~ NA_real_
    ) 
  )

#=========================MICROSCOPY VARIABLES FOR ANALYSIS=====================

SMART_microscopy_clean <- SMART_microscopy_clean |>
  mutate(
    # 0 = no parasites detected
    # 1 = parasites detected
    # NA = final parasite density unresolved/unavailable
    microscopy_positive = case_when(
      parasite_density == 0 ~ 0,
      parasite_density > 0  ~ 1,
      TRUE ~ NA_real_
    )
  )


#========CHECK U01 MICROSCOPY AGAINST UNSCHEDULED CLINICAL VISITS===============

u01_microscopy <- SMART_microscopy_clean |>
  filter(microscopy_visit == "U01") |>
  select(
    participant_id,
    sample_id,
    microscopy_visit,
    parasite_density,
    microscopy_positive,
    falciparum,
    malariae,
    ovale,
    vivax
  )

u01_microscopy

SMART_analysis |>
  filter(
    participant_id %in% u01_microscopy$participant_id,
    study_visit == "Unscheduled visit"
  ) |>
  select(
    participant_id,
    screening_id,
    study_visit,
    a02_uv
  ) |>
  arrange(participant_id, a02_uv) |>
  print(n = Inf)

# Standard visits that can currently be linked to SMART_analysis
standard_microscopy_visits <- c(
  "H0", "H8", "H24", "H36", "H48", "H60", "H72",
  "D7", "D14", "D21", "D28", "D35", "D42"
)

# create two microscopy dataset
SMART_microscopy_standard <- SMART_microscopy_clean |>
  filter(
    microscopy_visit %in% standard_microscopy_visits
  ) |>
  select(
    microscopy_sample_id = sample_id,
    participant_id,
    study_visit = microscopy_visit,
    parasite_density,
    microscopy_positive,
    falciparum,
    malariae,
    ovale,
    vivax
  )

SMART_microscopy_u01 <- SMART_microscopy_clean |>
  filter(microscopy_visit == "U01") |>
  transmute(
    microscopy_sample_id = sample_id,
    participant_id,
    unscheduled_visit = microscopy_visit,
    parasite_density,
    microscopy_positive,
    falciparum,
    malariae,
    ovale,
    vivax
  )

#check if it matches
SMART_microscopy_u01 |>
  left_join(
    SMART_analysis |>
      filter(study_visit == "Unscheduled visit") |>
      select(
        participant_id,
        study_visit,
        unscheduled_visit,
        a02_uv
      ),
    by = c("participant_id", "unscheduled_visit")
  ) |>
  select(
    participant_id,
    unscheduled_visit,
    study_visit,
    a02_uv,
    microscopy_sample_id,
    parasite_density
  ) |>
  print(n = Inf)




#----------------- Create microscopy dataset for analysis-----------------------
SMART_microscopy_analysis <- SMART_microscopy_clean |>
  filter(
    microscopy_visit %in% c(
      standard_microscopy_visits,
      "U01"
    )
  ) |>
  transmute(
    microscopy_sample_id = sample_id,
    participant_id,
    
    study_visit = case_when(
      microscopy_visit == "U01" ~ "Unscheduled visit",
      TRUE ~ microscopy_visit
    ),
    
    unscheduled_visit = case_when(
      microscopy_visit == "U01" ~ "U01",
      TRUE ~ NA_character_
    ),
    
    parasite_density,
    microscopy_positive,
    falciparum,
    malariae,
    ovale,
    vivax
  )

#create a temporary join key in both dataset so we can join unsheduled visit in 
#microscopy dataset
SMART_analysis <- SMART_analysis |>
  mutate(
    microscopy_join_key = case_when(
      study_visit == "Unscheduled visit" &
        !is.na(unscheduled_visit) ~
        paste(participant_id, unscheduled_visit, sep = "_"),
      
      TRUE ~
        paste(participant_id, study_visit, sep = "_")
    )
  )

SMART_microscopy_analysis <- SMART_microscopy_analysis |>
  mutate(
    microscopy_join_key = case_when(
      study_visit == "Unscheduled visit" ~
        paste(participant_id, unscheduled_visit, sep = "_"),
      
      TRUE ~
        paste(participant_id, study_visit, sep = "_")
    )
  )

# safety check
# Check participant-visit uniqueness
SMART_microscopy_analysis |>
  count(participant_id, study_visit) |>
  filter(n > 1)

# Check microscopy records that do not match SMART_analysis
SMART_microscopy_analysis |>
  anti_join(
    SMART_analysis |>
      distinct(participant_id, study_visit),
    by = c("participant_id", "study_visit")
  ) |>
  count(study_visit, sort = TRUE) |>
  print(n = Inf)


#===============================================================================
# Check whether unmatched microscopy records belong to participants
# outside the SMART analysis population
#===============================================================================
SMART_microscopy_analysis |>
  anti_join(
    SMART_analysis |>
      distinct(participant_id),
    by = "participant_id"
  ) |>
  summarise(
    n_records = n(),
    n_participants = n_distinct(participant_id)
  )

microscopy_from_DA <- SMART_microscopy_analysis |>
  anti_join(
    SMART_analysis |>
      distinct(participant_id),
    by = "participant_id"
  ) |>
  distinct(participant_id) |>
  arrange(participant_id) |>
  print(n = Inf)
  #View()

#=====Microscopy records for participants who don't have corresponding visit in
#SMART_analysis
#
check_ids <- SMART_microscopy_analysis |>
  semi_join(
    SMART_analysis |>
      distinct(participant_id),
    by = "participant_id"
  ) |>
  anti_join(
    SMART_analysis |>
      distinct(participant_id, study_visit),
    by = c("participant_id", "study_visit")
  ) |>
  select(
    participant_id,
    study_visit,
    microscopy_sample_id,
    parasite_density,
    microscopy_positive
  ) |>
  arrange(participant_id, study_visit) |>
  print(n = Inf)

# Participants whose microscopy visit does not match a SMART_analysis visit
 SMART_microscopy_analysis |>
  semi_join(
    SMART_analysis |> distinct(participant_id),
    by = "participant_id"
  ) |>
  anti_join(
    SMART_analysis |> distinct(participant_id, study_visit),
    by = c("participant_id", "study_visit")
  ) |>
  distinct(participant_id)

# Show all clinical visits available for these participants
SMART_analysis |>
  filter(participant_id %in% check_ids$participant_id) |>
  select(participant_id, study_visit) |>
  distinct() |>
  arrange(participant_id, study_visit) |>
  print(n = Inf)

SMART_analysis |>
  filter(
    participant_id %in% check_ids$participant_id,
    study_visit == "Treatment failure"
  ) |>
  select(
    participant_id,
    study_day,
    study_visit,
    microscopy_done,
    microscopy_result
  ) |>
  arrange(participant_id) |>
  print(n = Inf)

#===============================================================================
#        CHECK ACTUAL VISIT OF TREATMENT FAILURE
#===============================================================================

treatment_failure_check <- SMARTUncomplicated |>
  filter(redcap_event_name == "treatment_failure_arm_1") |>
  select(
    participant_id_tf_um_sm,
    a01_visit_sm,
    a01_2_if_unscheduled_state_sm,
    a01_tf_um_sm,
    b03_malaria_microscopy_pos_sm
  ) |>
  arrange(participant_id_tf_um_sm)

#============= PREPARE TREATMENT FAILURE VISITS=================================

treatment_failure_check <- SMARTUncomplicated |>
  filter(redcap_event_name == "treatment_failure_arm_1") |>
  transmute(
    screening_id = id_um,
    
    tf_visit = case_when(
      as.numeric(a01_visit_sm) == 1  ~ "H0",
      as.numeric(a01_visit_sm) == 2  ~ "H24",
      as.numeric(a01_visit_sm) == 3  ~ "H48",
      as.numeric(a01_visit_sm) == 4  ~ "H72",
      as.numeric(a01_visit_sm) == 5  ~ "D7",
      as.numeric(a01_visit_sm) == 6  ~ "D14",
      as.numeric(a01_visit_sm) == 7  ~ "D21",
      as.numeric(a01_visit_sm) == 8  ~ "D28",
      as.numeric(a01_visit_sm) == 9  ~ "D35",
      as.numeric(a01_visit_sm) == 10 ~ "D42",
      as.numeric(a01_visit_sm) == 11 ~ "Unscheduled visit",
      TRUE ~ NA_character_
    ),
    
    tf_unscheduled_visit = if_else(
      as.numeric(a01_visit_sm) == 11,
      as.character(a01_2_if_unscheduled_state_sm),
      NA_character_
    ),
    
    tf_visit_date = as.Date(a01_tf_um_sm),
    tf_microscopy_result = b03_malaria_microscopy_pos_sm
  ) |>
  left_join(
    SMART_analysis |>
      distinct(screening_id, participant_id),
    by = "screening_id"
  )

treatment_failure_check |>
  summarise(
    n_tf_records = n(),
    n_with_participant_id = sum(!is.na(participant_id)),
    n_without_participant_id = sum(is.na(participant_id))
  )


#=======COMPARE UNMATCHED MICROSCOPY WITH TREATMENT FAILURE VISITS===============


a <-SMART_microscopy_analysis |>
  semi_join(
    SMART_analysis |> distinct(participant_id),
    by = "participant_id"
  ) |>
  anti_join(
    SMART_analysis |> distinct(participant_id, study_visit),
    by = c("participant_id", "study_visit")
  ) |>
  left_join(
    treatment_failure_check,
    by = "participant_id"
  ) |>
  select(
    participant_id,
    microscopy_visit = study_visit,
    parasite_density,
    tf_visit,
    tf_unscheduled_visit,
    tf_visit_date,
    tf_microscopy_result
  ) |>
  arrange(participant_id, microscopy_visit) |>
  #View()
  print(n = Inf)


# ---------Confirm the remaining unmatched participant-visits-------------------
SMART_microscopy_analysis |>
  anti_join(
    SMART_analysis |>
      distinct(participant_id, study_visit),
    by = c("participant_id", "study_visit")
  ) |>
  count(study_visit, sort = TRUE) |>
  print(n = Inf)


#==========JOIN MICROSCOPY DATA TO SMART ANALYSIS===============================

SMART_analysis <- SMART_analysis |>
  left_join(
    SMART_microscopy_analysis |>
      select(
        -participant_id,
        -study_visit,
        -unscheduled_visit
      ),
    by = "microscopy_join_key"
  ) |>
  select(-microscopy_join_key)


# Check total records in the analysis data now
SMART_analysis |>
  summarise(
    n_microscopy_records = sum(!is.na(microscopy_sample_id)),
    n_parasite_density = sum(!is.na(parasite_density)),
    n_positive = sum(microscopy_positive == 1, na.rm = TRUE),
    n_negative = sum(microscopy_positive == 0, na.rm = TRUE)
  )

SMART_analysis |>
  filter(!is.na(microscopy_sample_id)) |>
  count(study_visit, sort = FALSE) |>
  print(n = Inf)
#View(SMART_analysis)

SMART_analysis |>
  filter(
    study_visit == "Unscheduled visit",
    !is.na(microscopy_sample_id)
  ) |>
  select(
    participant_id,
    unscheduled_visit,
    microscopy_sample_id,
    parasite_density
  )



#-------FINAL MICROSCOPY DATA-QUALITY SUMMARY-------------

# 1. Exact breakdown of usable, missing and discrepant final parasite density
SMART_microscopy_clean |>
  count(final_parasite_status, sort = TRUE) |>
  print(n = Inf)


# 2. Final microscopy information successfully linked to SMART_analysis
SMART_analysis |>
  summarise(
    n_microscopy_records = sum(!is.na(microscopy_sample_id)),
    n_parasite_density   = sum(!is.na(parasite_density)),
    n_positive           = sum(microscopy_positive == 1, na.rm = TRUE),
    n_negative           = sum(microscopy_positive == 0, na.rm = TRUE)
  )


# 3. Recreate the FINAL unmatched microscopy dataset after adding U01
unmatched_microscopy <- SMART_microscopy_analysis |>
  anti_join(
    SMART_analysis |>
      distinct(
        participant_id,
        study_visit,
        unscheduled_visit
      ),
    by = c(
      "participant_id",
      "study_visit",
      "unscheduled_visit"
    )
  )

unmatched_microscopy |>
  summarise(
    n_records = n(),
    n_participants = n_distinct(participant_id)
  )

unmatched_microscopy |>
  count(
    study_visit,
    unscheduled_visit,
    sort = TRUE
  ) |>
  print(n = Inf)
#View(unmatched_microscopy)

# 4. Distinguish unmatched participants from unmatched visits
unmatched_microscopy |>
  mutate(
    unmatched_reason = if_else(
      participant_id %in% SMART_analysis$participant_id,
      "Participant present but visit not matched",
      "Participant not in SMART analysis population"
    )
  ) |>
  count(unmatched_reason)

SMART_microscopy_clean |>
  filter(final_parasite_status == "Other") |>
  select(
    sample_id,
    participant_id,
    microscopy_visit,
    final_parasitemia_over_all_readings_in_parasites_m_l
  ) |>
  print(n = Inf)

SMART_microscopy_clean |>
  count(final_parasite_status, sort = TRUE)

#Count parasite counts that are less than 0
SMART_microscopy_clean |>
  filter(
    final_parasitemia_over_all_readings_in_parasites_m_l < 0
  ) |>
  count(
    final_parasitemia_over_all_readings_in_parasites_m_l
  )


#======CREATING PARTICIPANT LAST FOLLOW-UP VISIT DAY & DATE=====================

#---------------------Create follow-up date-------------------------------------

SMART_analysis <- SMART_analysis |>
  mutate(
    followup_date = case_when(
      
      study_visit %in% c(
        "H8", "H24", "H36", "H48", "H60"
      ) ~ a01_fuv,
      
      study_visit %in% c(
        "H72", "D7", "D14", "D21",
        "D28", "D35", "D42"
      ) ~ a01_d_um,
      
      study_visit == "Unscheduled visit" ~ a02_uv,
      
      TRUE ~ as.Date(NA)
    ),
    
    followup_visit = case_when(
      study_visit == "Unscheduled visit" ~ unscheduled_visit,
      !is.na(followup_date) ~ study_visit,
      TRUE ~ NA_character_
    )
  )

#---------------------Create  last follow-up day--------------------------------
last_followup <- SMART_analysis |>
  filter(
    !is.na(followup_date)
  ) |>
  group_by(participant_id) |>
  filter(
    followup_date == max(followup_date)
  ) |>
  summarise(
    last_followup_date = last(followup_date),
    last_followup_visit = last(followup_visit),
    .groups = "drop"
  )
#View(last_followup)

#--------------join to smart_analysis--------------
SMART_analysis <- SMART_analysis |>
  left_join(
    last_followup,
    by = "participant_id"
  )

#check
SMART_analysis |>
  select(
    participant_id,
    last_followup_visit,
    last_followup_date
  ) |>
  distinct() |>
  arrange(participant_id) |>
  print(n = Inf)

#-----------Check for duplicate last follow-up visit
SMART_analysis |>
  select(
    participant_id,
    last_followup_visit,
    last_followup_date
  ) |>
  distinct() |>
  summarise(
    n_participants = n(),
    n_with_last_followup = sum(!is.na(last_followup_date)),
    n_without_last_followup = sum(is.na(last_followup_date))
  )

#--------------check participants with no followup------------------------------
SMART_analysis |>
  select(
    participant_id,
    last_followup_visit,
    last_followup_date
  ) |>
  distinct() |>
  filter(is.na(last_followup_date)) |>
  print(n = Inf)

#------------check participants with more that one visit day--------------------
SMART_analysis |>
  select(
    participant_id,
    last_followup_visit,
    last_followup_date
  ) |>
  distinct() |>
  filter(str_detect(last_followup_visit, " / ")) |>
  print(n = Inf)

#----------CHECK PARTICIPANTS WITH NO IDENTIFIED FOLLOW-UP----------------------

SMART_analysis |>
  filter(
    participant_id %in% c(
      "S24-10062-AG",
      "S24-10182-AF"
    )
  ) |>
  select(
    participant_id,
    redcap_event_name,
    study_visit,
    a01_fuv,
    a01_d_um,
    a01_uv,
    a02_uv,
    followup_visit,
    followup_date
  ) |>
  arrange(participant_id) |>
  print(n = Inf, width = Inf)

# check for other follow-up details for these two participants
no_last_followup <- SMART_analysis |>
  filter(
    participant_id %in% c(
      "S24-10062-AG",
      "S24-10182-AF"
    )
  ) |>
  select(
    participant_id,
    redcap_event_name,
    study_visit,
    visitp_fup,
    visitp_fupd,
    a01_fuv,
    a01_d_um,
    a01_uv,
    a02_uv,
    d01_date_of_study_end,
    d3_1_was_the_participant_fol,
    d3_2_if_no_select_reason_f,
    sremarks
  ) |>
  arrange(participant_id) |>
  print(n = Inf, width = Inf)

#------check scheduled form for their study-end information--------------
SMART_analysis |>
  filter(
    participant_id %in% c(
      "S24-10062-AG",
      "S24-10182-AF"
    ),
    study_visit == "Scheduled visit"
  ) |>
  select(
    participant_id,
    d01_date_of_study_end,
    d3_1_was_the_participant_fol,
    d3_2_if_no_select_reason_f,
    sremarks
  ) |>
  print(n = Inf, width = Inf)
 
#--------------check treatment failure form for these participants
treatment_failure_check |>
  filter(
    participant_id %in% c(
      "S24-10062-AG",
      "S24-10182-AF"
    )
  ) |>
  print(n = Inf, width = Inf)

#==========================PARTICIPANT STUDY-END INFORMATION====================

study_end_info <- SMART_analysis |>
  filter(study_visit == "Scheduled visit") |>
  transmute(
    participant_id,
    
    study_end_date = d01_date_of_study_end,
    
    study_end_reason = case_when(
      
      # Participant completed normal Day 42 follow-up
      d3_1_was_the_participant_fol == 1 ~
        "Regular study end",
      
      # Participant did not complete normal Day 42 follow-up
      d3_1_was_the_participant_fol == 0 &
        d3_2_if_no_select_reason_f == 1 ~
        "Treatment failure",
      
      d3_1_was_the_participant_fol == 0 &
        d3_2_if_no_select_reason_f == 2 ~
        "Voluntary withdrawal by participant",
      
      d3_1_was_the_participant_fol == 0 &
        d3_2_if_no_select_reason_f == 3 ~
        "Withdrawal based on investigator decision",
      
      d3_1_was_the_participant_fol == 0 &
        d3_2_if_no_select_reason_f == 4 ~
        "Loss to follow-up",
      
      d3_1_was_the_participant_fol == 0 &
        d3_2_if_no_select_reason_f == 5 ~
        "Participant died",
      
      TRUE ~ NA_character_
    )
  )

study_end_info |>
  count(study_end_reason, sort = TRUE) |>
  print(n = Inf)

study_end_info |>
  count(participant_id) |>
  filter(n > 1)

# cross checking 124 treatment failure participants against the 108 selected as 
#a reason
tf_participants <- treatment_failure_check |>
  distinct(participant_id)

tf_participants |>
  left_join(
    study_end_info,
    by = "participant_id"
  ) |>
  count(study_end_reason, sort = TRUE) 

#check those not classified as treatment failure at the end of the study
tf_end_mismatch <- tf_participants |>
  left_join(
    study_end_info,
    by = "participant_id"
  ) |>
  filter(
    is.na(study_end_reason) |
      study_end_reason != "Treatment failure"
  )

tf_end_mismatch 
tf_end_mismatch |>
  count(study_end_reason, sort = TRUE)

#----------check across record for these participants---------------
tf_end_investigation <- treatment_failure_check |>
  select(
    participant_id,
    tf_visit,
    tf_unscheduled_visit,
    tf_visit_date
  ) |>
  left_join(
    study_end_info,
    by = "participant_id"
  ) |>
  left_join(
    SMART_analysis |>
      select(
        participant_id,
        last_followup_visit,
        last_followup_date
      ) |>
      distinct(),
    by = "participant_id"
  ) |>
  filter(
    is.na(study_end_reason) |
      study_end_reason != "Treatment failure"
  ) |>
  mutate(
    tf_matches_study_end_date =
      tf_visit_date == study_end_date,
    
    tf_matches_last_followup_date =
      tf_visit_date == last_followup_date
  ) |>
  arrange(participant_id)

tf_end_investigation |>
  select(
    participant_id,
    tf_visit,
    tf_unscheduled_visit,
    tf_visit_date,
    last_followup_visit,
    last_followup_date,
    study_end_date,
    study_end_reason,
    tf_matches_study_end_date,
    tf_matches_last_followup_date
  ) 

#==========================JOIN STUDY-END INFORMATION===========================

SMART_analysis <- SMART_analysis |>
  left_join(
    study_end_info,
    by = "participant_id"
  )

SMART_analysis |>
  summarise(
    n_rows = n(),
    n_participants = n_distinct(participant_id),
    n_with_study_end_date = n_distinct(
      participant_id[!is.na(study_end_date)]
    )
  )

#===================MICROSCOPY VISIT MISMATCH INVESTIGATION====================

microscopy_visit_mismatch <- unmatched_microscopy |>
  semi_join(
    SMART_analysis |>
      distinct(participant_id),
    by = "participant_id"
  )

#check actual visit information
microscopy_visit_mismatch_check <- microscopy_visit_mismatch |>
  left_join(
    SMART_analysis |>
      select(
        participant_id,
        last_followup_visit,
        last_followup_date
      ) |>
      distinct(),
    by = "participant_id"
  ) |>
  left_join(
    treatment_failure_check |>
      select(
        participant_id,
        tf_visit,
        tf_unscheduled_visit,
        tf_visit_date
      ) |>
      distinct(),
    by = "participant_id"
  )

microscopy_visit_mismatch_check |>
  select(
    participant_id,
    microscopy_sample_id,
    study_visit,
    unscheduled_visit,
    last_followup_visit,
    last_followup_date,
    tf_visit,
    tf_unscheduled_visit,
    tf_visit_date
  ) |>
  arrange(participant_id) |>
  print(n = Inf, width = Inf)
 
# CHECK day of unscheduled failure for this participant if it maybe D28 or not
SMART_analysis |> 
  filter(participant_id =="S24-10033-AG") |> 
  select(
    participant_id,
    redcap_event_name,
    study_visit
  )

SMART_analysis |>
  filter(
    participant_id == "S24-10033-AG",
    study_visit %in% c("D21", "Unscheduled visit")
  ) |>
  select(
    participant_id,
    study_visit,
    unscheduled_visit,
    a01_d_um,
    a02_uv
  )

as.Date("2024-06-22") - as.Date("2024-06-18")


#===========================2.PCR DATA==========================================

PCR_raw <- read_csv(
  "Data/pcr_correction_results.csv",
  show_col_types = FALSE
)

dim(PCR_raw)
names(PCR_raw)

PCR_raw |>
  slice_head(n = 10)

#   Clean PCR object
SMART_pcr <- PCR_raw |>
  select(-`...1`) |>
  rename(
    participant_id = Screening,
    pcr_failure_visit = `Failure Date`,
    pcr_failure_sample_id = Failure,
    pcr_msp1 = msp1,
    pcr_msp2 = msp2,
    pcr_ta60 = ta60,
    pcr_conclusion_who = `Conclusion WHO 3/3`,
    pcr_conclusion_2of3 = `Conclusion 2/3`,
    pcr_execution_time = `Execution Time`,
    pcr_delta_value = `Delta Value`,
    pcr_peak_height = `Peak Height`
  )

SMART_pcr |>
  summarise(
    n_records = n(),
    n_participants = n_distinct(participant_id),
    n_failure_samples = n_distinct(pcr_failure_sample_id)
  )

# Check if all participants in PCR data exist in SMART_analysis
SMART_pcr |>
  anti_join(
    SMART_analysis |>
      distinct(participant_id),
    by = "participant_id"
  ) |>
  select(
    participant_id,
    pcr_failure_visit,
    pcr_failure_sample_id
  ) |>
  print(n = Inf)

# Exclude PCR records that do not belong to the SMART analysis population
SMART_pcr <- SMART_pcr |>
  filter(
    !str_detect(participant_id, "-TP$"),
    !str_detect(participant_id, "^S24-20")
  )

SMART_pcr |>
  summarise(
    n_records = n(),
    n_participants = n_distinct(participant_id),
    n_failure_samples = n_distinct(pcr_failure_sample_id)
  )

# Check again if all participants belong to SMART_analysis
SMART_pcr |>
  anti_join(
    SMART_analysis |>
      distinct(participant_id),
    by = "participant_id"
  ) |>
  select(
    participant_id,
    pcr_failure_visit,
    pcr_failure_sample_id
  ) |>
  print(n = Inf)

# Compare PCR failure visit with the treatment failure visit in SMART data
pcr_tf_check <- SMART_pcr |>
  left_join(
    treatment_failure_check |>
      select(
        participant_id,
        tf_visit,
        tf_unscheduled_visit,
        tf_visit_date,
        tf_microscopy_result
      ),
    by = "participant_id"
  )

pcr_tf_check |>
  count(
    pcr_failure_visit,
    tf_visit,
    tf_unscheduled_visit,
    sort = TRUE
  ) |>
  print(n = Inf)

#-------------- Check any disagreement--------------------
pcr_tf_check |>
  filter(
    is.na(tf_visit) |
      (
        tf_visit != "Unscheduled visit" &
          pcr_failure_visit != tf_visit
      ) |
      (
        tf_visit == "Unscheduled visit" &
          pcr_failure_visit != tf_unscheduled_visit
      )
  ) |>
  select(
    participant_id,
    pcr_failure_visit,
    pcr_failure_sample_id,
    tf_visit,
    tf_unscheduled_visit,
    tf_visit_date,
    tf_microscopy_result
  ) |>
  arrange(participant_id) |>
  print(n = Inf)
#----------------------------------------------------------
  
names(SMART_pcr)
glimpse(SMART_pcr)

#---------------Inspect PCR marker results----------------------------

SMART_pcr |>
  count(pcr_msp1, sort = TRUE) |>
  print(n = Inf)

SMART_pcr |>
  count(pcr_msp2, sort = TRUE) |>
  print(n = Inf)

SMART_pcr |>
  count(pcr_ta60, sort = TRUE) |>
  print(n = Inf)

#----------------Inspect the two PCR conclusions--------------------------------
SMART_pcr |>
  count(pcr_conclusion_who, sort = TRUE) |>
  print(n = Inf)

SMART_pcr |>
  count(pcr_conclusion_2of3, sort = TRUE) |>
  print(n = Inf)

#-----------Check disagreement between the two conclusions
SMART_pcr |>
  count(
    pcr_conclusion_who,
    pcr_conclusion_2of3,
    sort = TRUE
  ) |>
  print(n = Inf)


# Create the visit that should correspond to SMART_analysis
SMART_pcr <- SMART_pcr |>
  mutate(
    study_visit = case_when(
      pcr_failure_visit == "U01" ~ "Unscheduled visit",
      TRUE ~ pcr_failure_visit
    )
  )

# Check participant + visit matching
SMART_pcr |>
  anti_join(
    SMART_analysis |>
      distinct(participant_id, study_visit),
    by = c("participant_id", "study_visit")
  ) |>
  select(
    participant_id,
    pcr_failure_visit,
    study_visit,
    pcr_failure_sample_id,
    pcr_conclusion_who,
    pcr_conclusion_2of3
  ) |>
  print(n = Inf)


#====================Create PCR dataset for analysis============================

# Create PCR analysis dataset while retaining visit discrepancies
SMART_pcr_analysis <- SMART_pcr |>
  left_join(
    treatment_failure_check |>
      select(
        participant_id,
        tf_visit,
        tf_unscheduled_visit,
        tf_visit_date,
        tf_microscopy_result
      ),
    by = "participant_id"
  ) |>
  mutate(
    pcr_tf_visit_match = case_when(
      tf_visit == pcr_failure_visit ~ TRUE,
      
      tf_visit == "Unscheduled visit" &
        pcr_failure_visit == tf_unscheduled_visit ~ TRUE,
      
      TRUE ~ FALSE
    )
  ) |>
  select(
    participant_id,
    pcr_failure_visit,
    pcr_failure_sample_id,
    pcr_msp1,
    pcr_msp2,
    pcr_ta60,
    pcr_conclusion_who,
    pcr_conclusion_2of3,
    tf_visit,
    tf_unscheduled_visit,
    tf_visit_date,
    tf_microscopy_result,
    pcr_tf_visit_match
  )

SMART_pcr_analysis |>
  summarise(
    n_records = n(),
    n_participants = n_distinct(participant_id),
    n_samples = n_distinct(pcr_failure_sample_id),
    n_visit_matches = sum(pcr_tf_visit_match),
    n_visit_discrepancies = sum(!pcr_tf_visit_match)
  )

#---------------------DOCUMENT UNMATCHED PCR RECORD-----------------------------

SMART_pcr_analysis |>
  filter(!pcr_tf_visit_match) |>
  select(
    participant_id,
    pcr_failure_visit,
    tf_visit,
    tf_unscheduled_visit,
    tf_visit_date,
    pcr_conclusion_who,
    pcr_conclusion_2of3
  )

#-------------------Select PCR variables needed for analysis---------------------
SMART_pcr_analysis <- SMART_pcr_analysis |>
  select(
    participant_id,
    pcr_failure_visit,
    pcr_failure_sample_id,
    
    # PCR marker results
    pcr_msp1,
    pcr_msp2,
    pcr_ta60,
    
    # PCR classifications
    pcr_conclusion_who,
    pcr_conclusion_2of3,
    
    # Treatment-failure information
    tf_visit,
    tf_unscheduled_visit,
    tf_visit_date,
    tf_microscopy_result,
    
    # PCR/clinical visit consistency check
    pcr_tf_visit_match
  )

glimpse(SMART_pcr_analysis)

#==============================JOIN PCR TO ANALYSIS DATA========================

#-------------MATCH U01 PCR FAILURES TO UNSCHEDULED VISIT DATES-----------------
SMART_pcr_analysis |>
  filter(pcr_failure_visit == "U01") |>
  select(
    participant_id,
    pcr_failure_visit,
    tf_visit_date
  ) |>
  left_join(
    SMART_analysis |>
      filter(study_visit == "Unscheduled visit") |>
      select(
        participant_id,
        study_visit,
        a02_uv
      ),
    by = "participant_id"
  ) |>
  mutate(
    date_match = tf_visit_date == a02_uv
  ) |>
  arrange(participant_id, a02_uv) |>
  print(n = Inf)


#                  DEFINE PCR JOIN TARGETS


SMART_pcr_join <- SMART_pcr_analysis |>
  filter(pcr_tf_visit_match == TRUE) |>
  mutate(
    pcr_target_visit = case_when(
      participant_id == "S24-10002-AG" ~ "Treatment failure",
      pcr_failure_visit == "U01" ~ "Unscheduled visit",
      TRUE ~ pcr_failure_visit
    ),
    
    pcr_target_date = case_when(
      pcr_failure_visit == "U01" ~ as.Date(tf_visit_date),
      TRUE ~ as.Date(NA)
    )
  )
SMART_pcr_join |>
  count(pcr_target_visit)

#check match of these pcr records
SMART_pcr_join |>
  rowwise() |>
  mutate(
    n_matching_rows = sum(
      SMART_analysis$participant_id == participant_id &
        SMART_analysis$study_visit == pcr_target_visit &
        (
          pcr_target_visit != "Unscheduled visit" |
            (!is.na(SMART_analysis$a02_uv) &
               SMART_analysis$a02_uv == pcr_target_date)
        ),
      na.rm = TRUE
    )
  ) |>
  ungroup() |>
  count(n_matching_rows)


# Create a unique key for the existing SMART_analysis rows
SMART_analysis <- SMART_analysis |>
  mutate(
    pcr_join_key = case_when(
      study_visit == "Unscheduled visit" & !is.na(a02_uv) ~
        paste(participant_id, study_visit, a02_uv, sep = "_"),
      
      TRUE ~
        paste(participant_id, study_visit, sep = "_")
    )
  )

# Create the corresponding key for the 116 confirmed PCR records
SMART_pcr_join <- SMART_pcr_join |>
  mutate(
    pcr_join_key = case_when(
      pcr_target_visit == "Unscheduled visit" ~
        paste(participant_id, pcr_target_visit, pcr_target_date, sep = "_"),
      
      TRUE ~
        paste(participant_id, pcr_target_visit, sep = "_")
    )
  )

# SELECT PCR VARIABLES TO JOIN  TO SMART_analysis
SMART_pcr_for_join <- SMART_pcr_join |>
  select(
    pcr_join_key,
    pcr_failure_visit,
    pcr_failure_sample_id,
    pcr_msp1,
    pcr_msp2,
    pcr_ta60,
    pcr_conclusion_who,
    pcr_conclusion_2of3,
    pcr_tf_visit_match
  )


# JOIN TO SMART_analysis
SMART_analysis <- SMART_analysis |>
  left_join(
    SMART_pcr_for_join,
    by = "pcr_join_key"
  ) |>
  select(-pcr_join_key)

# Confirm join
SMART_analysis |>
  summarise(
    n_rows = n(),
    n_pcr_records = sum(!is.na(pcr_failure_sample_id)),
    n_pcr_participants = n_distinct(
      participant_id[!is.na(pcr_failure_sample_id)]
    )
  )

SMART_analysis |>
  filter(!is.na(pcr_failure_sample_id)) |>
  count(study_visit, pcr_failure_visit)


#------------------THE TWO PCR VISIT DISCREPANCIES--------------------------

pcr_visit_mismatch_check <- SMART_pcr_analysis |>
  filter(!pcr_tf_visit_match) |>
  select(
    participant_id,
    pcr_failure_visit,
    pcr_failure_sample_id,
    pcr_conclusion_who,
    pcr_conclusion_2of3
  ) |>
  left_join(
    SMART_analysis |>
      select(
        participant_id,
        last_followup_visit,
        last_followup_date,
        study_end_date,
        study_end_reason
      ) |>
      distinct(),
    by = "participant_id"
  ) |>
  left_join(
    treatment_failure_check |>
      select(
        participant_id,
        tf_visit,
        tf_unscheduled_visit,
        tf_visit_date
      ),
    by = "participant_id"
  )

pcr_visit_mismatch_check |>
  select(
    participant_id,
    pcr_failure_visit,
    last_followup_visit,
    last_followup_date,
    tf_visit,
    tf_unscheduled_visit,
    tf_visit_date,
    study_end_date,
    study_end_reason
  ) |>
  print(n = Inf, width = Inf)

#check unscheduled visit date and the last scheduled visit date
SMART_analysis |>
  filter(
    participant_id == "S24-10033-AG",
    study_visit %in% c("D21", "Unscheduled visit")
  ) |>
  select(
    participant_id,
    study_visit,
    unscheduled_visit,
    a01_d_um,
    a02_uv
  ) |>
  print(n = Inf)

#calculate their interval
SMART_analysis |>
  filter(
    participant_id == "S24-10033-AG",
    study_visit %in% c("D21", "Unscheduled visit")
  ) |>
  summarise(
    d21_date = a01_d_um[study_visit == "D21"],
    u01_date = a02_uv[
      study_visit == "Unscheduled visit" &
        unscheduled_visit == "U01"
    ],
    days_between = as.numeric(u01_date - d21_date)
  )
  

#==========================3.PK DATA============================================

pk_data <- read_excel(
  "Data/PK_Lumefantrine_Pyronaridine_levels.xlsx",
  sheet = "Uncomplicated_AL_AP_combined",
  col_names = FALSE
)

dim(pk_data)

pk_data |>
  slice(1:10) |>
  print(width = Inf)

#-------------------------RESHAPE PK DATA---------------------------------------

# Day 07
pk_day07 <- pk_data |>
  slice(-(1:3)) |>
  select(
    pk_date = ...1,
    participant_id = ...2,
    pk_concentration = ...3
  ) |>
  mutate(
    pk_day = "D7"
  )

# Day 03
pk_day03 <- pk_data |>
  slice(-(1:3)) |>
  select(
    pk_date = ...5,
    participant_id = ...6,
    pk_concentration = ...7
  ) |>
  mutate(
    pk_day = "D3"
  )

# H00
pk_h00 <- pk_data |>
  slice(-(1:3)) |>
  select(
    pk_date = ...9,
    participant_id = ...10,
    pk_concentration = ...11
  ) |>
  mutate(
    pk_day = "H0"
  )

#------Stack them
SMART_pk <- bind_rows(
  pk_day07,
  pk_day03,
  pk_h00
) |>
  filter(!is.na(participant_id)) |>
  mutate(
    participant_id = str_trim(participant_id),
    pk_concentration = as.numeric(pk_concentration),
    pk_date = as.Date(as.numeric(pk_date), origin = "1899-12-30")
  ) |>
  select(
    participant_id,
    pk_day,
    pk_date,
    pk_concentration
  )

dim(SMART_pk)

SMART_pk |>
  count(pk_day)

SMART_pk |>
  group_by(pk_day) |>
  summarise(
    n_records = n(),
    n_participants = n_distinct(participant_id),
    n_missing_concentration = sum(is.na(pk_concentration)),
    min_concentration = min(pk_concentration, na.rm = TRUE),
    max_concentration = max(pk_concentration, na.rm = TRUE)
  )

# Check participants with more than one PK measurement on the same PK day
SMART_pk |>
  count(participant_id, pk_day) |>
  filter(n > 1) |>
  arrange(participant_id, pk_day)

# Inspect records for them
SMART_pk |>
  semi_join(
    SMART_pk |>
      count(participant_id, pk_day) |>
      filter(n > 1),
    by = c("participant_id", "pk_day")
  ) |>
  arrange(participant_id, pk_day, pk_date) |>
  print(n = Inf)


#-----------------------------CHECK CLEAN PK DATA-------------------------------
SMART_pk |>
  summarise(
    n_records = n(),
    n_participants = n_distinct(participant_id),
    n_missing_concentration = sum(is.na(pk_concentration))
  )

SMART_pk |>
  count(pk_day)

# Check PK participants not found in SMART_analysis
SMART_pk |>
  anti_join(
    SMART_analysis |>
      select(participant_id) |>
      distinct(),
    by = "participant_id"
  ) |>
  distinct(participant_id)


# Check PK measurements by treatment arm and sampling day
SMART_pk |>
  left_join(
    SMART_analysis |>
      select(participant_id, treatment_arm) |>
      distinct(),
    by = "participant_id"
  ) |>
  count(treatment_arm, pk_day)

#Final PK count

SMART_pk |>
  summarise(
    n_records = n(),
    n_participants = n_distinct(participant_id),
    n_missing_concentration = sum(is.na(pk_concentration))
  )

SMART_pk |>
  count(pk_day)

SMART_pk |>
  count(participant_id, pk_day) |>
  filter(n > 1)

# Creating a day match to smart_analysis
SMART_pk <- SMART_pk |>
  mutate(
    pk_target_visit = case_when(
      pk_day == "H0" ~ "H0",
      pk_day == "D3" ~ "H72",
      pk_day == "D7" ~ "D7",
      TRUE ~ NA_character_
    )
  )

# check PK data has exactly one row in smart_analysis
SMART_pk |>
  rowwise() |>
  mutate(
    n_matching_rows = sum(
      SMART_analysis$participant_id == participant_id &
        SMART_analysis$study_visit == pk_target_visit,
      na.rm = TRUE
    )
  ) |>
  ungroup() |>
  count(n_matching_rows)

#=============================JOIN PK TO SMART_ANALYSIS=========================

SMART_pk_join <- SMART_pk |>
  transmute(
    participant_id,
    study_visit = pk_target_visit,
    pk_day,
    pk_date,
    pk_concentration
  )

SMART_analysis <- SMART_analysis |>
  left_join(
    SMART_pk_join,
    by = c("participant_id", "study_visit")
  )

SMART_analysis |>
  summarise(
    n_rows = n(),
    n_pk_records = sum(!is.na(pk_concentration)),
    n_pk_participants = n_distinct(
      participant_id[!is.na(pk_concentration)]
    )
  )

SMART_analysis |>
  filter(!is.na(pk_concentration)) |>
  count(study_visit, pk_day)


#===========================ADD TREATMENT FAILURE INFORMATION===================

SMART_tf <- treatment_failure_check |>
  select(
    participant_id,
    tf_visit,
    tf_unscheduled_visit,
    tf_visit_date,
    tf_microscopy_result
  )

SMART_tf |>
  count(participant_id) |>
  filter(n > 1)

#-------------------Join treatment failure to smart_analysis--------------------
SMART_analysis <- SMART_analysis |>
  left_join(
    SMART_tf,
    by = "participant_id"
  )
SMART_analysis |>
  summarise(
    n_rows = n(),
    n_participants = n_distinct(participant_id),
    n_tf_participants = n_distinct(
      participant_id[!is.na(tf_visit)]
    )
  )


#=======================SELECTING FINAL VARIABLES FOR ANALYSIS================
dim(SMART_analysis)

names(SMART_analysis)
glimpse(SMART_analysis)

SMART_analysis |>
  summarise(
    n_rows = n(),
    n_participants = n_distinct(participant_id),
    missing_site = sum(is.na(site)),
    missing_age = sum(is.na(age)),
    missing_weight = sum(is.na(weight)),
    missing_sex = sum(is.na(sex))
  )

#-------------Checking Unscheduled visits info for participants-----------------
  # Participants with unscheduled visit rows
  SMART_analysis |>
  filter(study_visit == "Unscheduled visit") |>
  count(participant_id, sort = TRUE) |>
  print(n = Inf)

multiple_unscheduled <- SMART_analysis |>
  filter(study_visit == "Unscheduled visit") |>
  count(participant_id) |>
  filter(n > 1)

multiple_unscheduled

#check information recorded for such participants
SMART_analysis |>
  semi_join(
    multiple_unscheduled,
    by = "participant_id"
  ) |>
  filter(study_visit == "Unscheduled visit") |>
  select(
    participant_id,
    unscheduled_visit,
    a01_uv,
    a02_uv,
    temperature,
    temperature_route,
    haemoglobin,
    new_complaint,
    microscopy_sample_id,
    parasite_density,
    microscopy_positive
  ) |>
  arrange(participant_id, unscheduled_visit) |>
  print(n = Inf, width = Inf)

#=====REPLACING TREATMENT FAILURE AND UNSCHEDULED VISITS WITH EXACT DAYS========

#-------Reference for deriving the study day of unscheduled visits--------------

scheduled_reference <- SMART_analysis |>
  filter(
    study_visit %in% c(
      "H0", "H8", "H24", "H36", "H48", "H60", "H72",
      "D7", "D14", "D21", "D28", "D35", "D42"
    ),
    !is.na(followup_date)
  ) |>
  mutate(
    protocol_day = case_when(
      study_visit %in% c("H0", "H8")   ~ 0,
      study_visit %in% c("H24", "H36") ~ 1,
      study_visit %in% c("H48", "H60") ~ 2,
      study_visit == "H72"              ~ 3,
      study_visit == "D7"               ~ 7,
      study_visit == "D14"              ~ 14,
      study_visit == "D21"              ~ 21,
      study_visit == "D28"              ~ 28,
      study_visit == "D35"              ~ 35,
      study_visit == "D42"              ~ 42
    )
  ) |>
  select(
    participant_id,
    scheduled_visit = study_visit,
    scheduled_date = followup_date,
    protocol_day
  )  

#------------Create a function for an unscheduled visit-------------------------
derive_unscheduled_day <- function(id, event_date) {
  
  if (is.na(event_date)) {
    return(NA_integer_)
  }
  
  previous_visit <- scheduled_reference |>
    filter(
      participant_id == id,
      scheduled_date <= event_date
    ) |>
    arrange(
      desc(scheduled_date),
      desc(protocol_day)
    ) |>
    slice(1)
  
  if (nrow(previous_visit) == 0) {
    return(NA_integer_)
  }
  
  previous_visit$protocol_day +
    as.integer(event_date - previous_visit$scheduled_date)
}

#=====================FINAL ANALYSIS-READY DATASET==============================
#===============================================================================

smart_analysis_final <- SMART_analysis |>
  
  
  filter(
    #remove rows that are not needed for analysis
    !study_visit %in% c(
      "Scheduled visit",
      "Concomitant medication"
    ),
    
    # Remove empty unscheduled-visit placeholder rows
    !(
      study_visit == "Unscheduled visit" &
        is.na(a01_uv) &
        is.na(a02_uv)
    )
  ) |>
  
  mutate(
    original_study_visit = study_visit
  ) |> 
  
  rowwise() |>
  
  mutate(
    
    # Derive the day only where necessary
    derived_day = case_when(
      
      # Ordinary unscheduled visit
      study_visit == "Unscheduled visit" ~
        derive_unscheduled_day(
          participant_id,
          a02_uv
        ),
      
      # Treatment failure occurring at an unscheduled visit
      study_visit == "Treatment failure" &
        tf_visit == "Unscheduled visit" ~
        derive_unscheduled_day(
          participant_id,
          tf_visit_date
        ),
      
      TRUE ~ NA_integer_
    ),
    
    #---------------------------------------------------------------------------
    # REPLACE STUDY DAY ONLY IN THE FINAL DATASET
    #---------------------------------------------------------------------------
    
    study_day = case_when(
      
      # Ordinary unscheduled visit
      study_visit == "Unscheduled visit" &
        !is.na(derived_day) ~
        paste0("Day ", derived_day),
      
      # TF occurring during an unscheduled visit
      study_visit == "Treatment failure" &
        tf_visit == "Unscheduled visit" &
        !is.na(derived_day) ~
        paste0("Day ", derived_day),
      
      # TF occurring at scheduled visits
      study_visit == "Treatment failure" &
        tf_visit == "H0" ~ "Day 0",
      
      study_visit == "Treatment failure" &
        tf_visit == "H24" ~ "Day 1",
      
      study_visit == "Treatment failure" &
        tf_visit == "H48" ~ "Day 2",
      
      study_visit == "Treatment failure" &
        tf_visit == "H72" ~ "Day 3",
      
      study_visit == "Treatment failure" &
        str_detect(tf_visit, "^D\\d+$") ~
        paste0(
          "Day ",
          str_remove(tf_visit, "D")
        ),
      
      # Normal scheduled rows stay unchanged
      TRUE ~ study_day
    ),
    
    #---------------------------------------------------------------------------
    # REPLACE STUDY VISIT ONLY IN THE FINAL DATASET
    #---------------------------------------------------------------------------
    
    study_visit = case_when(
      
      # Ordinary unscheduled visit
      study_visit == "Unscheduled visit" &
        !is.na(derived_day) ~
        paste0("D", derived_day),
      
      # TF during unscheduled visit
      study_visit == "Treatment failure" &
        tf_visit == "Unscheduled visit" &
        !is.na(derived_day) ~
        paste0("D", derived_day),
      
      # TF already linked to a scheduled visit
      study_visit == "Treatment failure" &
        !is.na(tf_visit) &
        tf_visit != "Unscheduled visit" ~
        tf_visit,
      
      # Normal scheduled visits remain unchanged
      TRUE ~ study_visit
    )
  ) |>
  
  select(
    participant_id,
    screening_id,
    study_day,
    study_visit,
    original_study_visit,
    treatment_arm,
    site,
    age,
    sex,
    weight,
    
    
    #unscheduled_visit,
    
    temperature,
    temperature_route,
    
    haemoglobin,
    
    drug_administered,
    drug_dose,
    drug_admin_date,
    drug_admin_time,
    #vomited,
    #drug_readministered,
    
    followup_date,
    last_followup_visit,
    last_followup_date,
    study_end_date,
    study_end_reason,
    
    tf_visit,
    #tf_unscheduled_visit,
    tf_visit_date,
    #tf_microscopy_result,
    
    #microscopy_sample_id,
    parasite_density,
    #microscopy_positive,
    
    pcr_failure_visit,
    #pcr_failure_sample_id,
    #pcr_msp1,
    #pcr_msp2,
    #pcr_ta60,
    pcr_conclusion_who,
    pcr_conclusion_2of3,
    
    pk_day,
    #pk_date,
    pk_concentration,
    
    new_complaint,
    
    total_bilirubin = c03_2_if_done,
    direct_bilirubin = c04_1,
    albumin = c05_2_if_done,
    alt = c05_2_if_done_2,
    ast = c07_2_if_done,
    alp = c08_2_if_done,
    ggt = c09_2_if_done,
    total_protein = c10_2_if_done,
  )

# Identify duplicated treatment-failure + unscheduled-visit pairs
tf_unscheduled_pairs <- smart_analysis_final |>
  group_by(
    participant_id,
    study_day,
    study_visit
  ) |>
  filter(
    any(original_study_visit == "Unscheduled visit") &
      any(original_study_visit == "Treatment failure")
  ) |>
  ungroup()

#check
tf_unscheduled_pairs |>
  select(
    participant_id,
    study_day,
    study_visit,
    original_study_visit,
    temperature,
    haemoglobin,
    parasite_density,
    pcr_failure_visit,
    pcr_conclusion_who
  ) |>
  arrange(participant_id, study_visit) |>
  print(n = Inf, width = Inf)

tf_unscheduled_pairs |>
  count(
    participant_id,
    study_day,
    study_visit
  )

#------merge rows with duplicte visit as a result of treatment failure and unscheduled visit
# Function to take the first non-missing value
first_non_missing <- function(x) {
  
  i <- which(!is.na(x))
  
  if (length(i) == 0) {
    return(x[NA_integer_][1])
  }
  
  x[i[1]]
}

#now merge
merged_tf_unscheduled <- tf_unscheduled_pairs |>
  
  # Treatment-failure row comes first
  arrange(
    participant_id,
    study_day,
    study_visit,
    desc(original_study_visit == "Treatment failure")
  ) |>
  
  group_by(
    participant_id,
    study_day,
    study_visit
  ) |>
  
  summarise(
    across(
      -original_study_visit,
      first_non_missing
    ),
    .groups = "drop"
  )

#check the merged results
merged_tf_unscheduled |>
  select(
    participant_id,
    study_day,
    study_visit,
    temperature,
    haemoglobin,
    parasite_density,
    pcr_failure_visit,
    pcr_conclusion_who
  ) |>
  arrange(participant_id) |>
  print(n = Inf, width = Inf)

#---Replace the duplicate in the data with the merged reult
smart_analysis_final <- smart_analysis_final |>
  
  anti_join(
    tf_unscheduled_pairs |>
      distinct(
        participant_id,
        study_day,
        study_visit
      ),
    by = c(
      "participant_id",
      "study_day",
      "study_visit"
    )
  ) |>
  
  bind_rows(
    merged_tf_unscheduled
  ) 


#---check if unscheduled visit and treatment failure still remain in the dataset
smart_analysis_final |>
  filter(
    study_visit %in% c(
      "Unscheduled visit",
      "Treatment failure"
    )
  ) |>
  select(
    participant_id,
    study_day,
    study_visit,
    tf_visit,
    tf_visit_date,
    followup_date
  ) |>
  print(n = Inf, width = Inf)

#---check derived unscheduled visit days
smart_analysis_final |>
  filter(
    str_detect(study_visit, "^D\\d+$")
  ) |>
  count(study_day, study_visit) |>
  arrange(
    as.numeric(str_remove(study_visit, "D"))
  ) |>
  print(n = Inf)

#---check there is no duplicate visits
smart_analysis_final |>
  count(
    participant_id,
    study_day,
    study_visit
  ) |>
  filter(n > 1) |>
  print(n = Inf)

#--- check the main dataset of these observed duplicates
duplicate_visits <- smart_analysis_final |>
  count(
    participant_id,
    study_day,
    study_visit
  ) |>
  filter(n > 1) |>
  select(
    participant_id,
    study_visit
  )

SMART_analysis |>
  semi_join(
    duplicate_visits,
    by = "participant_id"
  ) |>
  filter(
    study_visit == "Treatment failure" |
      study_visit %in% duplicate_visits$study_visit
  ) |>
  select(
    participant_id,
    study_visit,
    tf_visit,
    temperature,
    haemoglobin,
    parasite_density
  ) |>
  arrange(participant_id) |>
  print(n = Inf, width = Inf)


#===MERGE SCHEDULED VISIT + TREATMENT FAILURE ROWS ON THE SAME STUDY DAY========

# Identify the duplicate scheduled + treatment-failure pairs
tf_scheduled_pairs <- smart_analysis_final |>
  group_by(
    participant_id,
    study_day,
    study_visit
  ) |>
  filter(
    any(original_study_visit == "Treatment failure") &
      any(original_study_visit != "Treatment failure")
  ) |>
  ungroup()

# define variable to choose highest value
measurement_vars <- c(
  "age",
  "weight",
  "temperature",
  "haemoglobin",
  "drug_dose",
  "parasite_density",
  "pk_concentration",
  "total_bilirubin",
  "direct_bilirubin",
  "albumin",
  "alt",
  "ast",
  "alp",
  "ggt",
  "total_protein"
)

#now merge
merged_tf_scheduled <- tf_scheduled_pairs |>
  
  # Put the normal scheduled row first for non-measurement variables
  arrange(
    participant_id,
    study_day,
    study_visit,
    original_study_visit == "Treatment failure"
  ) |>
  
  group_by(
    participant_id,
    study_day,
    study_visit
  ) |>
  
  summarise(
    
    # For actual numeric measurements:
    # if both rows have values, take the highest
    across(
      any_of(measurement_vars),
      ~ if (all(is.na(.x))) NA_real_ else max(.x, na.rm = TRUE)
    ),
    
    # For every other variable:
    # take the first available non-missing value
    across(
      -c(
        original_study_visit,
        any_of(measurement_vars)
      ),
      ~ {
        x <- .x[!is.na(.x)]
        
        if (length(x) == 0) {
          .x[NA_integer_][1]
        } else {
          x[1]
        }
      }
    ),
    
    .groups = "drop"
  )

#------Replace the two duplicated study visit row with the merged row
smart_analysis_final <- smart_analysis_final |>
  anti_join(
    tf_scheduled_pairs |>
      distinct(
        participant_id,
        study_day,
        study_visit
      ),
    by = c(
      "participant_id",
      "study_day",
      "study_visit"
    )
  ) |>
  bind_rows(
    merged_tf_scheduled
  ) |>
  select(- original_study_visit)

# Confirm there are no duplicates
smart_analysis_final |>
  count(
    participant_id,
    study_day,
    study_visit
  ) |>
  filter(n > 1) |>
  print(n = Inf)

View(smart_analysis_final)




#-------Remove rows without records after follow-up end-------------------------
#Checking if each row contains records 
check_followup_rows <- smart_analysis_final |>
  mutate(
    has_record = 
      !is.na(followup_date) |
      !is.na(temperature) |
      !is.na(haemoglobin) |
      !is.na(parasite_density) |
      !is.na(drug_administered) |
      !is.na(drug_dose) |
      !is.na(pk_concentration) |
      !is.na(pcr_conclusion_who) |
      !is.na(pcr_conclusion_2of3) |
      !is.na(new_complaint)
  )
# create a temporary order variable
check_followup_rows <- check_followup_rows |>
  mutate(
    visit_order = case_when(
      study_visit == "H0"  ~ 0,
      study_visit == "H8"  ~ 8,
      study_visit == "H24" ~ 24,
      study_visit == "H36" ~ 36,
      study_visit == "H48" ~ 48,
      study_visit == "H60" ~ 60,
      study_visit == "H72" ~ 72,
      
      str_detect(study_visit, "^D\\d+$") ~
        as.numeric(str_remove(study_visit, "D")) * 24,
      
      TRUE ~ NA_real_
    )
  )

#find thelast row that contains information for participants
last_real_followup <- check_followup_rows |>
  filter(
    has_record,
    !is.na(visit_order)
  ) |>
  group_by(participant_id) |>
  summarise(
    last_recorded_order = max(visit_order),
    .groups = "drop"
  )

#identify rows that would be removed
rows_after_followup <- check_followup_rows |>
  left_join(
    last_real_followup,
    by = "participant_id"
  ) |>
  filter(
    !has_record,
    !is.na(last_recorded_order),
    visit_order > last_recorded_order
  )

#check
rows_after_followup |>
  arrange(participant_id, visit_order) |>
  select(
    participant_id,
    study_day,
    study_visit,
    followup_date,
    temperature,
    haemoglobin,
    parasite_density,
    last_followup_visit,
    last_followup_date,
    study_end_reason
  ) |>
  print(n = Inf, width = Inf)

#-------Remove the identified rows without any records from the study visits----
smart_analysis_final <- smart_analysis_final |>
  anti_join(
    rows_after_followup |>
      distinct(
        participant_id,
        study_day,
        study_visit
      ),
    by = c(
      "participant_id",
      "study_day",
      "study_visit"
    )
  )

#==================MISSING PARASITE DENSITY RECORDS=============================

missing_parasite_density <- smart_analysis_final |>
  filter(
    is.na(parasite_density)
  ) |>
  select(
    participant_id,
    study_day,
    study_visit,
    parasite_density,
  )

missing_parasite_density |>
  print(n = Inf, width = Inf)

#write_csv(missing_parasite_density, "Discrepancy data/missing_parasite_density.csv")

missing_parasite_density |>
  count(
    study_day,
    study_visit,
    sort = TRUE
  ) |>
  print(n = Inf)

max(smart_analysis_final$haemoglobin, na.rm = TRUE)
max(smart_analysis_final$temperature, na.rm = TRUE) 


#====================ORDER FOLLOW-UP VISITS IN ANALYSIS DATASET=================
#                        
#===============================================================================
# 
# smart_analysis_final <- smart_analysis_final |>
#   mutate(
#     visit_order = case_when(
#       study_visit == "H0"  ~ 0,
#       study_visit == "H8"  ~ 8,
#       study_visit == "H24" ~ 24,
#       study_visit == "H36" ~ 36,
#       study_visit == "H48" ~ 48,
#       study_visit == "H60" ~ 60,
#       study_visit == "H72" ~ 72,
#       
#       # Any day-based visit, including derived days such as D5, D18, D31
#       str_detect(study_visit, "^D\\d+$") ~
#         as.numeric(str_remove(study_visit, "D")) * 24,
#       
#       TRUE ~ NA_real_
#     )
#   ) |>
#   arrange(
#     participant_id,
#     visit_order
#   )
# 
# #check for a participant
# smart_analysis_final |>
#   filter(participant_id == "S24-10061-AF") |>
#   select(
#     participant_id,
#     study_day,
#     study_visit,
#     visit_order,
#     temperature,
#     haemoglobin,
#     parasite_density
#   ) |>
#   print(n = Inf)


# #=============================================================================
# #                 CREATE DATA-QUERY WORKBOOK 
# #=============================================================================
# 
# # 1. Microscopy participants not found in SMART_analysis
# microscopy_not_in_smart <- unmatched_microscopy |>
#   anti_join(
#     SMART_analysis |>
#       distinct(participant_id),
#     by = "participant_id"
#   ) |>
#   select(
#     participant_id,
#     microscopy_sample_id,
#     study_visit,
#     unscheduled_visit,
#     parasite_density,
#     microscopy_positive
#   ) |>
#   arrange(participant_id, study_visit)
# 
# 
# # 2. Participant exists, but microscopy visit did not match
# microscopy_visit_query <- microscopy_visit_mismatch_check |>
#   select(
#     participant_id,
#     microscopy_sample_id,
#     study_visit,
#     unscheduled_visit,
#     last_followup_visit,
#     last_followup_date,
#     tf_visit,
#     tf_unscheduled_visit,
#     tf_visit_date
#   ) |>
#   arrange(participant_id)
# 
# # 3. Unusual microscopy sample IDs
# microscopy_unusual_id_query <- unusual_id |>
#   select(
#     sample_id,
#     reading_date_2,
#     reading_date_3,
#     reading_date_4,
#     detection_6,
#     detection_15,
#     detection_29,
#     detection_38,
#     x3rd_reading_done,
#     falciparum,
#     malariae,
#     ovale,
#     vivax,
#     final_parasitemia_over_all_readings_in_parasites_m_l
#   )
# 
# # 4. PCR visit discrepancies
# pcr_visit_query <- pcr_visit_mismatch_check |>
#   select(
#     participant_id,
#     pcr_failure_visit,
#     pcr_failure_sample_id,
#     pcr_conclusion_who,
#     pcr_conclusion_2of3,
#     last_followup_visit,
#     last_followup_date,
#     tf_visit,
#     tf_unscheduled_visit,
#     tf_visit_date,
#     study_end_date,
#     study_end_reason
#   ) |>
#   arrange(participant_id)
# 
# # Participants with no identifiable follow-up
# no_followup_query <- SMART_analysis |>
#   filter(
#     participant_id %in% c(
#       "S24-10062-AG",
#       "S24-10182-AF"
#     )
#   ) |>
#   select(
#     participant_id,
#     redcap_event_name,
#     study_visit,
#     dadmin,
#     dadmin1,
#     visitp_fup,
#     visitp_fupd,
#     a01_fuv,
#     a01_d_um,
#     a01_uv,
#     a02_uv,
#     d01_date_of_study_end,
#     d3_1_was_the_participant_fol,
#     d3_2_if_no_select_reason_f,
#     sremarks
#   ) |>
#   mutate(
#     across(
#       everything(),
#       ~ replace_na(as.character(.), "NA")
#     )
#   ) |>
#   arrange(participant_id)
# writexl::write_xlsx(
#   list(
#     Microscopy_Not_In_SMART = microscopy_not_in_smart,
#     Microscopy_Visit_Mismatch = microscopy_visit_query,
#     Microscopy_Unusual_IDs = microscopy_unusual_id_query,
#     PCR_Visit_Mismatch = pcr_visit_query,
#     No_followup = no_followup_query
#   ),
#   path = "Data/SMART_Data_Queries.xlsx"
# )
