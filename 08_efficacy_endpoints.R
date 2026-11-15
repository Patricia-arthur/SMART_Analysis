source("07_treatment_outcomes.R")
#==========================PRIMARY EFFICACY ANALYSIS============================
#----------------------------AL DAY 28 UNCORRECTED -----------------------------

al28_population <- eligible_participants |> 
  filter(treatment_arm == "AL") |> 
  select(
    participant_id,
    treatment_arm
  ) |> 
  distinct()

al28_population |> 
  ungroup() |> 
  summarise(
    n = n_distinct(participant_id)
  )

#----------------------AL day 28 LCF--------------------------------------------
al28_lcf <- bind_rows(
  
  #Criterion 1: danger sign/severe malaria + parasitaemia
  lcf_danger_severe |> 
    filter(
      treatment_arm == "AL",
      analysis_day <= 28
    ) |> 
    select(
      participant_id,
      treatment_arm
    ),
  
  #Criterion 2: fever >=37.5 + parasitaemia
  lcf_data |> 
    filter(
      treatment_arm =="AL",
      analysis_day <= 28,
      lcf_fever == TRUE
    ) |> 
    select(
      participant_id,
      treatment_arm
    )
) |> 
  distinct(
    participant_id,
    treatment_arm
  )

al28_lcf |> 
  ungroup() |> 
  summarise(
    n_lcf = n_distinct(participant_id)
  )

#----------------------AL day 28 LPF--------------------------------------------
#participants with no ETF and no LCF by day 28
al28_no_lcf <- lcf_data |> 
  filter(
    treatment_arm == "AL"
  ) |> 
  anti_join(
    al28_lcf,
    by = "participant_id"
  )

al28_lpf_data <- al28_no_lcf |>
  filter(
    analysis_day >= 7,
    analysis_day <= 28
  ) |>
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

#---------------Check lpf cases--------------
al28_lpf_data |>
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

#---------------- Day 28 AL LFP participants-----------------------------
al28_lpf <- al28_lpf_data |>
  filter(lpf == TRUE) |>
  distinct(
    participant_id,
    treatment_arm
  )

al28_lpf |>
  ungroup() |> 
  summarise(
    n_lpf = n_distinct(participant_id)
  )

#----------------------AL day 28 ACPR--------------------------------------------
al28_acpr_data <- al28_no_lcf |>
  anti_join(
    al28_lpf,
    by = "participant_id"
  ) |>
  filter(
    analysis_day == 28
  )

#Inspect day 28 records
al28_acpr_data |>
  select(
    participant_id,
    treatment_arm,
    analysis_day,
    study_visit,
    parasite_density,
    microscopy_positive,
    temperature
  ) |>
  print(n = Inf)

#-------acpr criterion: absence of parasitaemia
al28_acpr_data <- al28_acpr_data |>
  mutate(
    parasitaemia_absent = case_when(
      
      !is.na(parasite_density) &
        parasite_density == 0 ~ TRUE,
      
      microscopy_positive == 0 ~ TRUE,
      
      !is.na(parasite_density) &
        parasite_density > 0 ~ FALSE,
      
      microscopy_positive == 1 ~ FALSE,
      
      TRUE ~ NA
    )
  )

# check cases
al28_acpr_data |>
  count(
    parasitaemia_absent,
    .drop = FALSE
  )

#----------------Check indeterminate Day-28 AL cases----------------------------

al28_acpr_data |>
  filter(
    is.na(parasitaemia_absent)
  ) |>
  select(
    participant_id,
    treatment_arm,
    analysis_day,
    study_visit,
    parasite_density,
    microscopy_positive,
    temperature
  ) |>
  print(n = Inf)

# Confirm the sum is correct
# Number of AL Day-28 LCF participants
nrow(al28_lcf)

# Number of AL Day-28 LPF participants
nrow(al28_lpf)

# Number reaching the ACPR assessment stage
nrow(al28_acpr_data)

nrow(al28_lcf) +
  nrow(al28_lpf) +
  nrow(al28_acpr_data)

#---------------Check unaccounted AL day 28 participant------------------------

al28_accounted <- bind_rows(
  al28_lcf |>
    select(participant_id),
  
  al28_lpf |>
    select(participant_id),
  
  al28_acpr_data |>
    select(participant_id)
) |>
  distinct(participant_id)

al28_unaccounted <- al28_population |>
  anti_join(
    al28_accounted,
    by = "participant_id"
  )

al28_unaccounted

nrow(al28_unaccounted)

#----------------Inspect unaccounted participants around Day 28----------------

lcf_data |>
  semi_join(
    al28_unaccounted,
    by = "participant_id"
  ) |>
  filter(
    analysis_day >= 21,
    analysis_day <= 35
  ) |>
  select(
    participant_id, 
    analysis_day,
    study_visit,
    parasite_density,
    microscopy_positive,
    temperature
  ) |>
  arrange(
    participant_id,
    analysis_day
  ) |>
  print(n = Inf)

#----------------Last observed follow-up day-----------------------------------

al28_last_followup <- lcf_data |>
  semi_join(
    al28_unaccounted,
    by = "participant_id"
  ) |>
  group_by(participant_id) |>
  summarise(
    last_day = max(
      analysis_day,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  arrange(last_day)

al28_last_followup |>
  print(n = Inf)

#records after day 28
lcf_data |>
  semi_join(
    al28_unaccounted,
    by = "participant_id"
  ) |>
  select(
    participant_id,
    analysis_day,
    study_visit,
    parasite_density,
    microscopy_positive
  ) |>
  arrange(
    participant_id,
    analysis_day
  ) |>
  print(n = Inf)

#---------Participants among the 20 who do not appear in lcf_data
al28_no_lcf_records <- al28_unaccounted |>
  anti_join(
    lcf_data |>
      distinct(participant_id),
    by = "participant_id"
  )

nrow(al28_no_lcf_records)

#--------------Inspect all available records
smart_analysis_final |>
  semi_join(
    al28_no_lcf_records,
    by = "participant_id"
  ) |>
  select(
    participant_id,
    study_day,
    study_visit,
    parasite_density,
    microscopy_positive,
    temperature
  ) |>
  arrange(
    participant_id,
    study_day
  ) |>
  print(n = Inf)


# their last day of study
al28_no_lcf_last_followup <- smart_analysis_final |>
  semi_join(
    al28_no_lcf_records,
    by = "participant_id"
  ) |>
  mutate(
    day_number = readr::parse_number(study_day)
  ) |>
  group_by(participant_id) |>
  summarise(
    last_day = max(
      day_number,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) |>
  arrange(last_day)

al28_no_lcf_last_followup |>
  print(n = Inf)

eligible_participants |> filter(grepl("S24-10252-AG", participant_id))
etf_data |> filter(grepl("S24-10252-AG", participant_id))
lcf_data |> filter(grepl("S24-10247-AG", participant_id))

# Check ETF status of the 11 participants absent from lcf_data
etf_data |>
  semi_join(
    al28_no_lcf_records,
    by = "participant_id"
  ) |>
  select(
    participant_id,
    treatment_arm,
    starts_with("parasite_density"),
    starts_with("temperature"),
    etf_status
  ) |>
  print(n = Inf)

#count etf status for 20 unidentified participant without day 28 AL record
etf_data |>
  semi_join(
    al28_unaccounted,
    by = "participant_id"
  ) |>
  count(etf_status)

#
treatment_failure_check |>
  semi_join(
    eligible_participants,
    by = "participant_id"
  ) |>
  filter(
    tf_visit %in% c(
      "H24",
      "H48",
      "H72"
    ),
    tf_microscopy_result == 1,
    (
      tf_danger_sign == 1 |
        tf_severe_manifestation == 1
    )
  ) |>
  select(
    participant_id,
    tf_visit,
    tf_microscopy_result,
    tf_danger_sign,
    tf_severe_manifestation
  ) |>
  print(n = Inf)

#================AL DAY-28 UNCORRECTED FINAL OUTCOME============================

al28_outcome <- al28_population |>
  mutate(
    outcome = case_when(
      
      # ETF
      participant_id %in%
        etf_data$participant_id[
          etf_data$etf_status == "ETF"
        ] ~ "ETF",
      
      # LCF by Day 28
      participant_id %in%
        al28_lcf$participant_id ~ "LCF",
      
      # LPF by Day 28
      participant_id %in%
        al28_lpf$participant_id ~ "LPF",
      
      # ACPR
      participant_id %in%
        al28_acpr_data$participant_id[
          al28_acpr_data$parasitaemia_absent == TRUE
        ] ~ "ACPR",
      
      # Everything else cannot be evaluated
      TRUE ~ "Indeterminate"
    )
  )

al28_outcome |>
  count(outcome)

al28_outcome |>
  ungroup() |> 
  summarise(
    total_AL = n(),
    evaluable = sum(outcome != "Indeterminate"),
    not_evaluable = sum(outcome == "Indeterminate")
  )

#=====================PCR-uncorrected AL Day-28 outcomes===========================
al28_evaluable <- al28_outcome |>
  ungroup() |> 
  filter(outcome != "Indeterminate")

al28_evaluable |>
  ungroup() |> 
  count(outcome) |>
  mutate(
    total = sum(n),
    percentage = round(100 * n / total, 1)
  )

#-------------------AL day 28 uncorrected outcomes with 95% CI-----------------

al28_uncorrected_summary <- tibble(
  outcome = c("ETF", "LCF", "LPF", "ACPR"),
  n = c(0, 14, 9, 129),
  total = 152
)

# Function for exact binomial confidence interval
exact_ci <- function(x,n){
  ci <- binom.test(
    x = x,
    n = n,
    conf.level = 0.95
  )$ conf.int
  
  tibble(
    percentage = round(100* x/n, 1),
    lower_95_CI = round(100* ci[1], 1),
    upper_95_CI = round(100* ci[2], 1)
  )
}

al28_uncorrected_summary <- al28_uncorrected_summary |> 
  bind_cols(
    purrr::map2_dfr(
      al28_uncorrected_summary$n,
      al28_uncorrected_summary$total,
      exact_ci
    )
  )
al28_uncorrected_summary



# al28_uncorrected_summary <- al28_uncorrected_summary |>
#   mutate(
#     result = paste0(
#       n, "/", total,
#       " (", percentage, "%; 95% CI ",
#       lower_95_CI, "–", upper_95_CI, "%)"
#     )
#   )
# 
# al28_uncorrected_summary

#==========================AL day 28 PCR corrected outcome=======================
#--------AL lcf events

al28_lcf_events <- bind_rows(
  
  # Criterion 1: danger sign/severe malaria + parasitaemia
  lcf_danger_severe |>
    ungroup() |>
    filter(
      treatment_arm == "AL",
      analysis_day <= 28
    ) |>
    select(
      participant_id,
      treatment_arm,
      analysis_day,
      study_visit
    ),
  
  # Criterion 2: fever + parasitaemia
  lcf_data |>
    ungroup() |>
    filter(
      treatment_arm == "AL",
      analysis_day <= 28,
      lcf_fever == TRUE
    ) |>
    select(
      participant_id,
      treatment_arm,
      analysis_day,
      study_visit
    )
) |>
  group_by(participant_id) |>
  slice_min(
    analysis_day,
    n = 1,
    with_ties = FALSE
  ) |>
  ungroup() |>
  mutate(
    outcome = "LCF"
  )

al28_lcf_events

#-------AL lpf event

al28_lpf_events <- al28_lpf_data |>
  ungroup() |>
  filter(
    lpf == TRUE
  ) |>
  group_by(participant_id) |>
  slice_min(
    analysis_day,
    n = 1,
    with_ties = FALSE
  ) |>
  ungroup() |>
  select(
    participant_id,
    treatment_arm,
    analysis_day,
    study_visit
  ) |>
  mutate(
    outcome = "LPF"
  )
al28_lpf_events

#------------------------Total AL day 28 failures
al28_failures <- bind_rows(
  al28_lcf_events,
  al28_lpf_events
) |>
  arrange(
    analysis_day,
    participant_id
  )

al28_failures

#-----------------Identify failures that require PCR correction-----------------
al28_failures <- al28_failures |>
  mutate(
    pcr_required = case_when(
      
      outcome == "LCF" &
        analysis_day < 7 ~ FALSE,
      
      outcome == "LCF" &
        analysis_day >= 7 ~ TRUE,
      
      outcome == "LPF" ~ TRUE
    )
  )

al28_failures |>
  count(
    outcome,
    pcr_required
  )

#--------------Attach  PCR conclusion-------------------------------------------
al28_failures_pcr <- al28_failures |>
  left_join(
    SMART_pcr_analysis |>
      select(
        participant_id,
        pcr_failure_visit,
        pcr_failure_sample_id,
        pcr_conclusion_who,
        pcr_tf_visit_match
      ),
    by = "participant_id"
  )

#Inspect
al28_failures_pcr |>
  select(
    participant_id,
    outcome,
    analysis_day,
    study_visit,
    pcr_required,
    pcr_failure_visit,
    pcr_conclusion_who,
    pcr_tf_visit_match
  ) |>
  arrange(
    analysis_day,
    participant_id
  ) |>
  print(n = Inf)

# count pcr conclusion
al28_failures_pcr |>
  filter(
    pcr_required == TRUE
  ) |>
  count(
    outcome,
    pcr_conclusion_who,
    .drop = FALSE
  )

al28_failures_pcr |>
  filter(
    pcr_required == TRUE
  ) |>
  count(
    pcr_conclusion_who,
    .drop = FALSE
  )


al28_failures_pcr |>
  filter(
    pcr_required == TRUE,
    !is.na(pcr_conclusion_who),
    pcr_tf_visit_match == FALSE
  ) |>
  select(
    participant_id,
    outcome,
    analysis_day,
    study_visit,
    pcr_failure_visit,
    pcr_conclusion_who
  )

#================AL DAY-28 PCR-CORRECTED OUTCOME===============================

al28_corrected <- al28_outcome |>
  ungroup() |>
  left_join(
    al28_failures_pcr |>
      select(
        participant_id,
        outcome_uncorrected = outcome,
        pcr_required,
        pcr_conclusion_who
      ),
    by = "participant_id"
  ) |>
  mutate(
    corrected_outcome = case_when(
      
      # Original ACPR remains ACPR
      outcome == "ACPR" ~ "ACPR",
      
      # ETF would remain failure
      outcome == "ETF" ~ "ETF",
      
      # Recrudescence remains the same failure type
      outcome %in% c("LCF", "LPF") &
        pcr_conclusion_who == "Recrudescence OMS" ~ outcome,
      
      # New infection is excluded from corrected PP analysis
      outcome %in% c("LCF", "LPF") &
        pcr_conclusion_who == "New infection" ~
        "Reinfection",
      
      # PCR could not classify recurrence
      outcome %in% c("LCF", "LPF") &
        pcr_conclusion_who == "Indeterminate" ~
        "PCR indeterminate",
      
      # Already not evaluable in uncorrected analysis
      outcome == "Indeterminate" ~
        "Indeterminate",
      
      TRUE ~ "Indeterminate"
    )
  )

al28_corrected |>
  count(corrected_outcome)

al28_corrected |>
  summarise(
    total_AL = n(),
    total_accounted = sum(!is.na(corrected_outcome))
  )

#---------------Create AL D-28 PCR corrected population-------------------------
al28_corrected_population <- al28_corrected |>
  filter(
    corrected_outcome %in% c(
      "ETF",
      "LCF",
      "LPF",
      "ACPR"
    )
  )

al28_corrected_population |>
  count(corrected_outcome)

al28_corrected_population |>
  summarise(
    n_evaluable = n()
  )

#-----------------------AL D-28 PCR Corrected proportions-----------------------
al28_corrected_population|>
  count(corrected_outcome) |>
  mutate(
    total = sum(n),
    percentage = round(100 * n / total, 1)
  )

#-------------------AL day 28 corrected outcomes with 95% CI--------------------
al28_corrected_summary <- tibble(
  outcome = c("ETF", "LCF", "LPF", "ACPR"),
  n = c(0, 2, 1, 129),
  total = 132
)

al28_corrected_summary <- al28_corrected_summary |>
  bind_cols(
    purrr::map2_dfr(
      al28_corrected_summary$n,
      al28_corrected_summary$total,
      exact_ci
    )
  )

al28_corrected_summary

#====================AL DAY 28 OUTCOME TABLE====================================

al28_uncorrected_table <- al28_uncorrected_summary |>
  mutate(
    outcome = factor(
      outcome,
      levels = c("ETF", "LCF", "LPF", "ACPR")
    ),
    uncorrected_result = paste0(
      n, " (", percentage, "%; 95% CI ",
      lower_95_CI, "–", upper_95_CI, ")"
    )
  ) |>
  select(
    outcome,
    uncorrected_result
  )

al28_corrected_table <- al28_corrected_summary |>
  mutate(
    outcome = factor(
      outcome,
      levels = c("ETF", "LCF", "LPF", "ACPR")
    ),
    corrected_result = paste0(
      n," (", percentage, "%; 95% CI ",
      lower_95_CI, "–", upper_95_CI, ")"
    )
  ) |>
  select(
    outcome,
    corrected_result
  )

al28_outcome_table <- al28_uncorrected_table |>
  left_join(
    al28_corrected_table,
    by = "outcome"
  ) |>
  arrange(outcome)

al28_outcome_table

#------gt table

al28_outcome_gt <- al28_outcome_table |>
  mutate(
    outcome = recode(
      as.character(outcome),
      "ETF" = "Early treatment failure (ETF)",
      "LCF" = "Late clinical failure (LCF)",
      "LPF" = "Late parasitological failure (LPF)",
      "ACPR" = "Adequate clinical and parasitological response (ACPR)"
    )
  ) |>
  gt::gt() |>
  gt::cols_label(
    outcome = "Treatment outcome",
    uncorrected_result = gt::md(
      "**PCR-uncorrected**  \nN = 152"
    ),
    corrected_result = gt::md(
      "**PCR-corrected**  \nN = 132"
    )
  ) |>
  gt::tab_header(
    title = gt::md(
      "**Artemether-Lumefantrine (AL): day 28 treatment outcomes**"
    ),
    subtitle = "PCR-uncorrected and corrected efficacy outcomes"
  ) |>
  gt::tab_source_note(
    source_note = gt::md(
      "Values are presented as **n (%; 95% CI)**."
    )
  )

al28_outcome_gt

#save
gt::gtsave(
  al28_outcome_gt,
  "figures/AL28_efficacy_outcomes.png"
)

# #================CHECK AL DAY-28 ACPR AT DAY 42================================
# 
# al28_acpr_followup <- al28_outcome |>
#   ungroup() |>
#   filter(outcome == "ACPR") |>
#   select(
#     participant_id,
#     day28_outcome = outcome
#   ) |>
#   left_join(
#     treatment_outcome |>
#       ungroup() |>
#       filter(treatment_arm == "AL") |>
#       select(
#         participant_id,
#         day42_outcome = treatment_outcome
#       ),
#     by = "participant_id"
#   )
# 
# al28_acpr_followup |>
#   count(day42_outcome)

#================AP DAY-42 PCR-UNCORRECTED OUTCOME==============================

ap42_outcome <- treatment_outcome |>
  ungroup() |>
  filter(
    treatment_arm == "Artesunate-Pyranoridine"
  ) |>
  transmute(
    participant_id,
    treatment_arm,
    outcome = as.character(treatment_outcome)
  )

# Check
ap42_outcome |>
  count(outcome)

#---------------AP day 42 population--------------------------------------------

ap42_evaluable <- ap42_outcome |>
  filter(
    outcome != "Indeterminate"
  )

ap42_evaluable |>
  summarise(
    n_evaluable = n()
  )

#----------------AP day 42 failure participants---------------------------

ap42_failures <- ap42_outcome |>
  filter(
    outcome %in% c("LCF", "LPF")
  )

nrow(ap42_failures)

#-------------------AP day 42 uncorrected outcomes with 95% CI------------------
ap42_uncorrected_summary <- tibble(
  outcome = c("ETF", "LCF", "LPF", "ACPR"),
  n = c(0, 10, 10, 129),
  total = 149
) |>
  bind_cols(
    purrr::map2_dfr(
      c(0, 10, 10, 129),
      rep(149, 4),
      exact_ci
    )
  )

ap42_uncorrected_summary

#----------------------Attach PCR conclusion to failures-----------------------
ap42_failures_pcr <- ap42_failures |>
  left_join(
    SMART_pcr_analysis |>
      select(
        participant_id,
        pcr_failure_visit,
        pcr_conclusion_who,
        pcr_tf_visit_match
      ),
    by = "participant_id"
  )

# inspect
ap42_failures_pcr |>
  select(
    participant_id,
    outcome,
    pcr_failure_visit,
    pcr_conclusion_who,
    pcr_tf_visit_match
  ) |>
  print(n = Inf)

ap42_failures_pcr |>
  count(
    outcome,
    pcr_conclusion_who,
    .drop = FALSE
  )

ap42_failures_pcr |>
  count(
    pcr_conclusion_who,
    .drop = FALSE
  )

#================AP DAY-42 PCR-CORRECTED OUTCOME===============================

ap42_corrected <- ap42_outcome |>
  ungroup() |>
  left_join(
    ap42_failures_pcr |>
      select(
        participant_id,
        pcr_conclusion_who
      ),
    by = "participant_id"
  ) |>
  mutate(
    corrected_outcome = case_when(
      
      # ACPR remains ACPR
      outcome == "ACPR" ~ "ACPR",
      
      # ETF remains ETF
      outcome == "ETF" ~ "ETF",
      
      # Recrudescence remains treatment failure
      outcome == "LCF" &
        pcr_conclusion_who == "Recrudescence OMS" ~ "LCF",
      
      outcome == "LPF" &
        pcr_conclusion_who == "Recrudescence OMS" ~ "LPF",
      
      # Reinfection excluded
      outcome %in% c("LCF", "LPF") &
        pcr_conclusion_who == "New infection" ~ "Reinfection",
      
      # PCR could not classify recurrence
      outcome %in% c("LCF", "LPF") &
        pcr_conclusion_who == "Indeterminate" ~ "PCR indeterminate",
      
      # Missing PCR excluded
      outcome %in% c("LCF", "LPF") &
        is.na(pcr_conclusion_who) ~ "PCR NA",
      
      # Already indeterminate
      outcome == "Indeterminate" ~ "Indeterminate",
      
      TRUE ~ "Indeterminate"
    )
  )

ap42_corrected |>
  count(corrected_outcome)

#-----------------------AP D-42 PCR Corrected population-----------------------
ap42_corrected_population <- ap42_corrected |>
  filter(
    corrected_outcome %in% c(
      "ETF",
      "LCF",
      "LPF",
      "ACPR"  
    )
  )

ap42_corrected_population |>
  summarise(
    n_evaluable = n()
  )

#-----------------------AP D-42 PCR Corrected proportions-----------------------
ap42_corrected_population |>
  count(corrected_outcome) |>
  mutate(
    total = sum(n),
    percentage = round(100 * n / total, 1)
  )

#-------------------AP day 42 corrected outcomes with 95% CI--------------------

ap42_corrected_summary <- tibble(
  outcome = c("ETF", "LCF", "LPF", "ACPR"),
  n = c(0, 1, 2, 129),
  total = 132
)

ap42_corrected_summary <- ap42_corrected_summary |>
  bind_cols(
    purrr::map2_dfr(
      ap42_corrected_summary$n,
      ap42_corrected_summary$total,
      exact_ci
    )
  )

ap42_corrected_summary

#================AP DAY 42 OUTCOME TABLE=======================================

ap42_uncorrected_table <- ap42_uncorrected_summary |>
  mutate(
    outcome = factor(
      outcome,
      levels = c("ETF", "LCF", "LPF", "ACPR")
    ),
    uncorrected_result = paste0(
      n, " (", percentage, "%; 95% CI ",
      lower_95_CI, "–", upper_95_CI, ")"
    )
  ) |>
  select(
    outcome,
    uncorrected_result
  )

ap42_corrected_table <- ap42_corrected_summary |>
  mutate(
    outcome = factor(
      outcome,
      levels = c("ETF", "LCF", "LPF", "ACPR")
    ),
    corrected_result = paste0(
      n, " (", percentage, "%; 95% CI ",
      lower_95_CI, "–", upper_95_CI, ")"
    )
  ) |>
  select(
    outcome,
    corrected_result
  )

ap42_outcome_table <- ap42_uncorrected_table |>
  left_join(
    ap42_corrected_table,
    by = "outcome"
  ) |>
  arrange(outcome)

ap42_outcome_table

#------GT table

ap42_outcome_table <- ap42_outcome_table |>
  mutate(
    outcome = recode(
      as.character(outcome),
      "ETF" = "Early treatment failure (ETF)",
      "LCF" = "Late clinical failure (LCF)",
      "LPF" = "Late parasitological failure (LPF)",
      "ACPR" = "Adequate clinical and parasitological response (ACPR)"
    )
  )

ap42_outcome_gt <- ap42_outcome_table |>
  gt::gt() |>
  gt::cols_label(
    outcome = "Treatment outcome",
    uncorrected_result = gt::md(
      "**PCR-uncorrected**  \nN = 149"
    ),
    corrected_result = gt::md(
      "**PCR-corrected**  \nN = 132"
    )
  ) |>
  gt::tab_header(
    title = gt::md(
      "**Artesunate-Pyronaridine (AP): day 42 treatment outcomes**"
    ),
    subtitle = "PCR-uncorrected and corrected efficacy outcomes"
  ) |>
  gt::tab_source_note(
    source_note = gt::md(
      "Values are presented as **n (%; 95% CI)**."
    )
  )

ap42_outcome_gt

#save table
gt::gtsave(
  ap42_outcome_gt, "figures/AP42_efficacy_outcomes.png"
)
