source("08_efficacy_endpoints.R")


#============================PARASITE CLEARANCE ANALYSIS=========================

clearance_visits <- c(
  "H0", "H8", "H24", "H36",
  "H48", "H60", "H72"
)

clearance_data <- smart_analysis_final |>
  ungroup() |>
  semi_join(
    eligible_participants,
    by = "participant_id"
  ) |>
  filter(
    study_visit %in% clearance_visits
  ) |>
  mutate(
    clearance_hour = case_when(
      study_visit == "H0"  ~ 0,
      study_visit == "H8"  ~ 8,
      study_visit == "H24" ~ 24,
      study_visit == "H36" ~ 36,
      study_visit == "H48" ~ 48,
      study_visit == "H60" ~ 60,
      study_visit == "H72" ~ 72
    )
  ) |>
  select(
    participant_id,
    treatment_arm,
    site,
    study_visit,
    clearance_hour,
    parasite_density,
    microscopy_positive
  ) |>
  arrange(
    participant_id,
    clearance_hour
  )

clearance_data |>
  summarise(
    n_participants = n_distinct(participant_id)
  )

#-----------microscopy completeness by arm and time
clearance_data |>
  group_by(
    treatment_arm,
    study_visit,
    clearance_hour
  ) |>
  summarise(
    n_participants = n_distinct(participant_id),
    n_density_available = sum(!is.na(parasite_density)),
    n_positive = sum(parasite_density > 0, na.rm = TRUE),
    n_zero = sum(parasite_density == 0, na.rm = TRUE),
    n_missing = sum(is.na(parasite_density)),
    .groups = "drop"
  ) |>
  arrange(
    treatment_arm,
    clearance_hour
  ) |>
  print(n = Inf)

# Parasite clearance pattern
clearance_pattern_check <- clearance_data |>
  filter(
    !is.na(parasite_density)
  ) |>
  group_by(participant_id) |>
  arrange(
    clearance_hour,
    .by_group = TRUE
  ) |>
  mutate(
    cleared_once = cumany(parasite_density == 0),
    positive_after_zero = parasite_density > 0 & lag(cleared_once, default = FALSE)
  ) |>
  ungroup()

clearance_pattern_check |>
  filter(
    positive_after_zero == TRUE
  ) |>
  select(
    participant_id,
    treatment_arm,
    clearance_hour,
    study_visit,
    parasite_density
  ) |>
  print(n = Inf)


#-------------------------Check positive after zero parasites-------------------
#-------------participant with parasites after an observed clearance between visits
clearance_reversal_participants <- clearance_pattern_check |>
  filter(
    positive_after_zero == TRUE
  ) |>
  distinct(participant_id)

clearance_data |>
  semi_join(
    clearance_reversal_participants,
    by = "participant_id"
  ) |>
  select(
    participant_id,
    treatment_arm,
    clearance_hour,
    study_visit,
    parasite_density,
    microscopy_positive
  ) |>
  arrange(
    participant_id,
    clearance_hour
  ) |>
  print(n = Inf)

#-----------------Create a parasite status variable-------------------------
clearance_data <- clearance_data |>
  mutate(
    parasite_status = case_when(
      
      !is.na(parasite_density) &
        parasite_density > 0 ~ "Positive",
      
      !is.na(parasite_density) &
        parasite_density == 0 ~ "Negative",
      
      is.na(parasite_density) &
        microscopy_positive == 1 ~ "Positive",
      
      is.na(parasite_density) &
        microscopy_positive == 0 ~ "Negative",
      
      TRUE ~ NA_character_
    )
  )

clearance_data |>
  count(
    treatment_arm,
    clearance_hour,
    parasite_status,
    .drop = FALSE
  ) |>
  print(n = Inf)

#-----------------Derive participant level clearance time-----------------------
clearance_participant <- clearance_data |> 
  group_by(
    participant_id,
    treatment_arm,
  ) |> 
  arrange(
    clearance_hour,
    .by_group = TRUE
  ) |> 
  summarise(
    
    #Last observed positive microscopy result
    last_positive_hour = max(
      clearance_hour[!is.na(parasite_status) & parasite_status == "Positive"],
      na.rm = TRUE
    ),
    
    #first observed negative after the positive result
    clearance_time = {
      negative_hours <- clearance_hour[ !is.na(parasite_status) &
        parasite_status == "Negative" & clearance_hour > last_positive_hour
      ]
      
      if(length(negative_hours) > 0){
        min(negative_hours)
      } else{
        NA_real_
      }
    },
    
    #last hour at which microscopy status was actually known
    last_observed_hour = max(
      clearance_hour[!is.na(parasite_status)], na.rm = TRUE
    ),
    .groups = "drop"
  ) |> 
  mutate(
    # Event = observed parasite clearance
    clearance_event = !is.na(clearance_time),
    
    # KM time: event time if cleared;
    # otherwise censor at last observed microscopy
    time_to_clearance = if_else(
      clearance_event,
      clearance_time,
      last_observed_hour
    )
  )

clearance_participant |>
  print(n = Inf)

# count by treatment arm
clearance_participant |>
  count(
    treatment_arm,
    clearance_event
  )

#inspect unresolved parasite clearance time
clearance_participant |>
  filter(
    is.na(clearance_time)
  ) |>
  select(
    participant_id,
    treatment_arm,
    last_positive_hour,
    clearance_time,
    last_observed_hour
  ) |>
  print(n = Inf)

# summarise observed clearance time
clearance_participant |>
  count(
    treatment_arm,
    clearance_time
  ) |>
  arrange(
    treatment_arm,
    clearance_time
  )

#-------------------Median clearance time
clearance_participant |>
  filter(clearance_event == TRUE) |>
  group_by(treatment_arm) |>
  summarise(
    n = n(),
    median_clearance = median(clearance_time),
    q1 = quantile(clearance_time, 0.25),
    q3 = quantile(clearance_time, 0.75),
    min_clearance = min(clearance_time),
    max_clearance = max(clearance_time),
    .groups = "drop"
  )

#================PARASITE CLEARANCE CURVE======================================
#------------% of participant parasitaemic---------------------
parasite_clearance_summary <- clearance_data |>
  filter(
    !is.na(parasite_status)
  ) |>
  group_by(
    treatment_arm,
    clearance_hour
  ) |>
  summarise(
    n_observed = n(),
    n_positive = sum(parasite_status == "Positive"),
    percent_positive = 100 * n_positive / n_observed,
    .groups = "drop"
  )

parasite_clearance_summary

#----------plot the curve
parasite_clearance_plot <- parasite_clearance_summary |>
  ggplot(
    aes(
      x = clearance_hour,
      y = percent_positive,
      colour = treatment_arm,
      group = treatment_arm
    )
  ) +
  geom_line(
    linewidth = 1.1
  ) +
  geom_point(
    size = 3
  ) +
  scale_x_continuous(
    breaks = c(0, 8, 24, 36, 48, 60, 72),
    labels = c(
      "0", "8", "24",
      "36", "48", "60", "72"
    )
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    labels = function(x) paste0(x)
  ) +
  labs(
    title = "Parasite clearance by treatment arm",
    subtitle = "Percentage of participants with parasitaemia",
    x = "Follow-up time(hours)",
    y = "Participants with parasitaemia(%)",
    colour = "Treatment arm"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

parasite_clearance_plot

#----save the plot
ggsave(
  filename = "figures/percentage_parasite_clearance.png",
  plot = parasite_clearance_plot,
  height = 8,
  width = 10,
  dpi = 300
  )


#----------------Decline of parasite density------------------------------------
parasite_density_decline <- clearance_data |>
  filter(
    !is.na(parasite_density)
  ) |>
  group_by(
    treatment_arm,
    clearance_hour
  ) |>
  summarise(
    n = n(),
    
    # Geometric mean allowing zero parasite counts
    geometric_mean = exp(
      mean(log(parasite_density + 1))
    ) - 1,
    
    .groups = "drop"
  )

parasite_density_decline

#--------plot the curve

parasite_density_decline_plot <- parasite_density_decline |>
  ggplot(
    aes(
      x = clearance_hour,
      y = geometric_mean,
      colour = treatment_arm,
      group = treatment_arm
    )
  ) +
  geom_line(
    linewidth = 1.0
  ) +
  geom_point(
    size = 2
  ) +
  scale_x_continuous(
    breaks = c(0, 8, 24, 36, 48, 60, 72),
    labels = c(
      "0", "8", "24",
      "36", "48", "60", "72"
    )
  ) +
  scale_y_continuous(
    trans = scales::pseudo_log_trans(
      base = 10
    ),
    breaks = c(
      0, 10, 100, 1000,
      2000
    ),
    labels = scales::comma
  ) +
  labs(
    title = "Parasite clearance by treatment arm",
    subtitle = "Geometric mean decline in parasite density",
    x = "Follow-up time(hours)",
    y = "Parasite density (parasites/µL)",
    colour = "Treatment arm"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

parasite_density_decline_plot

#----------save the plot
ggsave(
  "figures/parasite_density_decline_log.png",
  parasite_density_decline_plot,
  width = 10,
  height = 8,
  dpi = 300
)

#====================== PARASITE CLEARANCE CURVE OVER 42 DAYS===================

parasite_42_visits <- c(
  "H0", "H24", "H48", "H72",
  "D7", "D14", "D21", "D28", "D35", "D42"
)

parasite_42_data <- smart_analysis_final |>
  ungroup() |>
  semi_join(
    eligible_participants,
    by = "participant_id"
  ) |>
  filter(
    study_visit %in% parasite_42_visits
  ) |>
  mutate(
    followup_day = case_when(
      study_visit == "H0"  ~ 0,
      study_visit == "H24" ~ 1,
      study_visit == "H48" ~ 2,
      study_visit == "H72" ~ 3,
      study_visit == "D7"  ~ 7,
      study_visit == "D14" ~ 14,
      study_visit == "D21" ~ 21,
      study_visit == "D28" ~ 28,
      study_visit == "D35" ~ 35,
      study_visit == "D42" ~ 42
    )
  )

parasite_42_summary <- parasite_42_data |>
  filter(
    !is.na(microscopy_positive)
  ) |>
  group_by(
    treatment_arm,
    followup_day
  ) |>
  summarise(
    n_observed = n(),
    n_positive = sum(microscopy_positive == 1),
    percent_positive = 100 * n_positive / n_observed,
    .groups = "drop"
  )

parasite_42_summary

#-------plot the clearance curve over 42 days
parasite_42_plot <- parasite_42_summary |>
  ggplot(
    aes(
      x = followup_day,
      y = percent_positive,
      colour = treatment_arm,
      group = treatment_arm
    )
  ) +
  geom_line(linewidth = 1) +
  geom_point(size = 2.0) +
  scale_x_continuous(
    breaks = c(0, 1, 2, 3, 7, 14, 21, 28, 35, 42)
  ) +
  labs(
    title = "Parasitaemia during the 42-day follow-up",
    subtitle = "Percentage of participants with parasitaemia by treatment arm",
    x = "Follow-up day",
    y = "Participants with parasitaemia (%)",
    colour = "Treatment arm"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

parasite_42_plot

#save plot
ggsave(
  filename = "figures/parasite_42_plot.png",
   plot = parasite_42_plot,
   width = 10,
   height = 8,
  dpi = 300
)


#=================KM: TIME TO PARASITE CLEARANCE================================
km_clearance_data <- clearance_participant |> 
  mutate(
    clearance_event_num = as.integer(clearance_event),
    treatment_arm = factor(
      treatment_arm,
      levels = c("AL", "AP")
    )
  )

km_clearance_data |> 
  count(
    treatment_arm,
    clearance_event_num
  )

#---------------------fit the km clearance curve--------------------------------
km_clearance_fit <- survfit2(
  Surv(time_to_clearance, clearance_event_num)~ treatment_arm, 
  data = km_clearance_data
  )

km_clearance_fit

#--------------------------------Log-rank test----------------------------------
pclearance_logrank <- survdiff(
  Surv(time_to_clearance, clearance_event_num) ~ treatment_arm, 
  data = km_clearance_data
)

pclearance_logrank

#---------pvalue
pclearance_pvalue <- pchisq(
  pclearance_logrank$chisq,
  df = length(pclearance_logrank$n) -1,
  lower.tail = FALSE
  )
pclearance_pvalue

plogrank_label <- paste0(
  "Log-rank test: \u03C7\u00B2 = ", round(pclearance_logrank$chisq, 2),
  ", p = ", format.pval(pclearance_pvalue, digits = 2, eps = 0.001)
)

#---------------------KM parasite clearance plot--------------------------------
pclearance_km_plot <- km_clearance_fit |> 
  ggsurvfit(linewidth = 0.8) +
  #add_confidence_interval() +
  add_censor_mark() +
  add_risktable(
    risktable_height = 0.25,
    risktable_group = "strata"
  ) +
  annotate(
    "text",
    x = 38,
    y = 0.25,
    label = plogrank_label,
    hjust = 0,
    size = 3.5
  )+
  labs(
    x = "Follow-up time(hours)",
    y = "Survival probability",
    title = "Kaplan-Meier: parasite clearance by treatment arm",
    caption = "Tick marks = censored observations"
  )+
  scale_x_continuous(
    breaks = c(0,8,24,36,48)
    # limits = c(0, 52),
    # expand = expansion(mult = c(0.01, 0.02))
  )+
  theme_ggsurvfit_default()

pclearance_km_plot

#-------save km clearance plot
ggsave(
  filename = "figures/pclearance_km_plot.png",
  plot = pclearance_km_plot,
  height = 8,
  width = 10,
  dpi = 300
)


#===================KM: TIME TO RECURRENT PARASITAEMIA==========================

#-----------------First recurrent parasitaemia----------------------------------
# Search all participants from Day 4 to Day 42. Event = parasite density > 0 OR microscopy positive.

recurrence_event <- smart_analysis_final |>
  ungroup() |>
  filter(
    readr::parse_number(study_day) >= 4,
    readr::parse_number(study_day) <= 42,
    (
      (!is.na(parasite_density) & parasite_density > 0) |
        microscopy_positive == 1 |
        tf_microscopy_result == 1 |
        (
          participant_id %in% treatment_failure_participant$participant_id &
            study_visit == tf_visit
        )
    )
  ) |>
  group_by(participant_id) |>
  slice_min(
    order_by = readr::parse_number(study_day),
    n = 1,
    with_ties = FALSE
  ) |>
  ungroup() |>
  transmute(
    participant_id,
    treatment_arm,
    recurrence_day = readr::parse_number(study_day)
  )

recurrence_event |>
  count(
    treatment_arm,
    recurrence_day
  ) |>
  print(n = Inf)

recurrence_event |>
  summarise(
    n_events = n(),
    n_participants = n_distinct(participant_id)
  )

#----------------------Attach PCR classification--------------------------------
recurrence_event <- recurrence_event |>
  left_join(
    SMART_pcr_analysis |>
      select(
        participant_id,
        pcr_conclusion_who
      ),
    by = "participant_id"
  )

recurrence_event |>
  count(
    treatment_arm,
    pcr_conclusion_who,
    .drop = FALSE
  )

View(recurrence_event)

#----------------Population for KM recurrence analysis--------------------------

recurrence_km_data <- smart_analysis_final |>
  ungroup() |>
  distinct(
    participant_id,
    treatment_arm
  )

recurrence_km_data <- recurrence_km_data |> 
  left_join(
    age_reference,
    by = "participant_id"
  )

#------------Attach recurrent events
recurrence_km_data <- recurrence_km_data |>
  left_join(
    recurrence_event |>
      select(
        participant_id,
        recurrence_day,
        pcr_conclusion_who
      ),
    by = "participant_id"
  )

#-----create event status
recurrence_km_data <- recurrence_km_data |> 
  mutate(
    recurrence_event = if_else(
      !is.na(recurrence_day),
      1L,
      0L
    )
  )

#------Get censoring time for those without recurrence

last_followup <- smart_analysis_final |>
  ungroup() |>
  mutate(
    followup_day = readr::parse_number(study_day)
  ) |>
  filter(
    followup_day <= 42,
    !is.na(parasite_density) |
      !is.na(microscopy_positive) |
      !is.na(tf_microscopy_result)
  ) |>
  group_by(participant_id) |>
  summarise(
    last_followup_day = max(
      followup_day,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# join to recurrence data
recurrence_km_data <- recurrence_km_data |>
  left_join(
    last_followup,
    by = "participant_id"
  ) |>
  mutate(
    time_to_recurrence = if_else(
      recurrence_event == 1,
      recurrence_day,
      last_followup_day
    )
  )

#----------check recurrence data

recurrence_km_data |>
  count(
    treatment_arm,
    recurrence_event
  )

recurrence_km_data |>
  ungroup() |> 
  summarise(
    n_participants = n(),
    n_unique = n_distinct(participant_id),
    n_recurrences = sum(recurrence_event == 1),
    n_censored = sum(recurrence_event == 0),
    missing_time = sum(is.na(time_to_recurrence))
  )

recurrence_km_data |>
  filter(
    recurrence_event == 0
  ) |>
  count(
    last_followup_day,
    sort = TRUE
  ) |>
  print(n = Inf)

#----check if participant has not event and time
recurrence_km_data |>
  filter(
    recurrence_event == 0,
    is.na(time_to_recurrence)
  ) |>
  select(
    participant_id,
    treatment_arm,
    age_group
  ) |>
  print(n = Inf)

#================ UNCORRECTED KAPLAN-MEIER ESTIMES OF RECURRENCE================

#===============KM: TIME TO RECURRENT PARASITAEMIA BY TREATMENT ARM=============

recurrence_km_fit_arm <- survfit2(
  Surv(time_to_recurrence, recurrence_event) ~ treatment_arm, 
  data = recurrence_km_data
)
recurrence_km_fit_arm

#----- estimated survival probability at Day 42
summary(
  recurrence_km_fit_arm,
  times = 42
)$surv


#-----------arm recurrence log-rank test-------------------
recurrence_logrank_arm <- survdiff(
  Surv(time_to_recurrence, recurrence_event) ~ treatment_arm, 
  data = recurrence_km_data
)
recurrence_logrank_arm

#-------recurrence p-value
recurrence_pval_arm <- pchisq(
  recurrence_logrank_arm$chisq,
  df = length(recurrence_logrank_arm$n) -1,
  lower.tail = FALSE
)
recurrence_pval_arm

recurrence_label_arm <- paste0(
  "Log-rank test: \u03C7\u00B2 = ", round(recurrence_logrank_arm$chisq, 2),
  ", p = ", format.pval(recurrence_pval_arm, digits = 3, eps = 0.001)
)

#---------------------KM recurrence by arm plot---------------------------------
recurrence_km_plot_arm <-recurrence_km_fit_arm |> 
  ggsurvfit(linewidth = 0.8) +
  #add_confidence_interval() +
  add_censor_mark()+
  add_risktable(
    risktable_group = "strata",
    risktable_height = 0.25
  ) +
  annotate(
    "text",
     x = 30,
     y = 0.25,
     label = recurrence_label_arm, 
     hjust = 0,
     size = 3.5,
  ) +
  labs(
    x = "Follow-up day",
    y = "Survival probability",
    title = "PCR-uncorrected Kaplan-Meier estimates of time to recurrent parasitaemia by treatment arm",
    caption = "Tick marks = censored observations" #Shaded bands = 95% CI ·
  )+
  scale_x_continuous(
    breaks = c(0, 1, 2, 3, 7, 14, 21, 28, 35, 42)
  )+
  theme_ggsurvfit_default()
recurrence_km_plot_arm

# save the arm km pcr-uncorrected plot
ggsave(
  filename = "figures/pcr-uncorrected_km_arm.png",
  plot = recurrence_km_plot_arm,
  height = 8,
  width = 10,
  dpi = 300
)


#===============KM: TIME TO RECURRENT PARASITAEMIA BY AGE GROUP=============
recurrence_km_fit_age <- survfit2(
  Surv(time_to_recurrence, recurrence_event) ~ age_group, 
  data = recurrence_km_data
)
recurrence_km_fit_age

#----- estimated survival probability at Day 42
summary(
  recurrence_km_fit_age,
  times = 42
)$surv


#-----------age_group recurrence log-rank test-------------------
recurrence_logrank_age <- survdiff(
  Surv(time_to_recurrence, recurrence_event) ~ age_group, 
  data = recurrence_km_data
)
recurrence_logrank_age

#-------recurrence p-value
recurrence_pval_age <- pchisq(
  recurrence_logrank_age$chisq,
  df = length(recurrence_logrank_age$n) -1,
  lower.tail = FALSE
)
recurrence_pval_age

recurrence_label_age <- paste0(
  "Log-rank test: \u03C7\u00B2 = ", round(recurrence_logrank_age$chisq, 2),
  ", p = ", format.pval(recurrence_pval_age, digits = 3, eps = 0.001)
)

#---------------------KM recurrence by age_group plot---------------------------
recurrence_km_plot_age <-recurrence_km_fit_age |> 
  ggsurvfit(linewidth = 0.8) +
  #add_confidence_interval() +
  add_censor_mark()+
  add_risktable(
    risktable_group = "strata",
    risktable_height = 0.25
  ) +
  annotate(
    "text",
    x = 30,
    y = 0.25,
    label = recurrence_label_age, 
    hjust = 0,
    size = 3.5,
  ) +
  labs(
    x = "Follow-up day",
    y = "Survival probability",
    title = "PCR-uncorrected Kaplan-Meier estimates of time to recurrent parasitaemia by age group",
    caption = "Tick marks = censored observations" #Shaded bands = 95% CI ·
  )+
  scale_x_continuous(
    breaks = c(0, 1, 2, 3, 7, 14, 21, 28, 35, 42)
  )+
  theme_ggsurvfit_default()
recurrence_km_plot_age

# save the age km pcr-uncorrected plot
ggsave(
  filename = "figures/pcr-uncorrected_km_age.png",
  plot = recurrence_km_plot_age,
  height = 8,
  width = 10,
  dpi = 300
)


#================ PCR CORRECTED KAPLAN-MEIER ESTIMES OF RECURRENCE==============

#--------------PCR corrected recurrence status----------------------
recurrence_km_data <- recurrence_km_data |> 
  mutate(
    
    pcr_recurrence_status = case_when(
      
      # No recurrence parasitaemia
      recurrence_event == 0 ~ "No recurrence",
      
      # PCR-confirmed recrudescence
      recurrence_event == 1 &
        pcr_conclusion_who == "Recrudescence OMS" ~ "Recrudescence",
      
      # PCR-confirmed new infection
      recurrence_event == 1 &
        pcr_conclusion_who == "New infection" ~ "Reinfection",
      
      # PCR couldn't classify the recurrence
      recurrence_event == 1 &
        (
          is.na(pcr_conclusion_who) |
          pcr_conclusion_who == "Indeterminate"
        ) ~ "PCR Indeterminate",
      
      TRUE ~NA_character_
    ),
    
    pcr_recurrence_event = case_when(
      
      # Recrudescence remains event
      pcr_recurrence_status == "Recrudescence" ~ 1L,
      
      # Reinfection is censored
      pcr_recurrence_status == "Reinfection" ~ 0L,
      
      # No recurrence is censored
      pcr_recurrence_status == "No recurrence" ~ 0L,
      
      # Missing/Indeterminate PCR is excluded
      pcr_recurrence_status == "PCR Indeterminate" ~ NA_integer_,
      
      TRUE ~ NA_integer_
    )
  )

# INSPECT IT 
recurrence_km_data |>
  count(
    treatment_arm,
    pcr_recurrence_status,
    pcr_recurrence_event,
    .drop = FALSE
  ) |>
  print(n = Inf)


recurrence_km_data |>
  filter(recurrence_event == 1) |>
  count(
    pcr_recurrence_status,
    .drop = FALSE
  )


#============ 1. KM Primary analysis (NA are dropped) ==========================

#===============KM: TIME TO RECURRENT PARASITAEMIA BY TREATMENT ARM=============
#cor means corrected 

cor_recurrence_km_fit_arm <- survfit2(
  Surv(time_to_recurrence, pcr_recurrence_event) ~ treatment_arm, 
  data = recurrence_km_data
)
cor_recurrence_km_fit_arm

#----- estimated survival probability at Day 42
summary(
  cor_recurrence_km_fit_arm,
  times = 42
)$surv


#-----------arm corrected recurrence log-rank test-------------------
cor_recurrence_logrank_arm <- survdiff(
  Surv(time_to_recurrence, pcr_recurrence_event) ~ treatment_arm, 
  data = recurrence_km_data
)
cor_recurrence_logrank_arm

#-------corrected recurrence p-value
cor_recurrence_pval_arm <- pchisq(
  cor_recurrence_logrank_arm$chisq,
  df = length(cor_recurrence_logrank_arm$n) -1,
  lower.tail = FALSE
)
cor_recurrence_pval_arm

cor_recurrence_label_arm <- paste0(
  "Log-rank test: \u03C7\u00B2 = ", round(cor_recurrence_logrank_arm$chisq, 2),
  ", p = ", format.pval(cor_recurrence_pval_arm, digits = 3, eps = 0.001)
)

#---------------------KM recurrence by arm plot---------------------------------
cor_recurrence_km_plot_arm <- cor_recurrence_km_fit_arm |> 
  ggsurvfit(linewidth = 0.8) +
  #add_confidence_interval() +
  add_censor_mark()+
  add_risktable(
    risktable_group = "strata",
    risktable_height = 0.25
  ) +
  annotate(
    "text",
    x = 30,
    y = 0.25,
    label = cor_recurrence_label_arm, 
    hjust = 0,
    size = 3.5
  ) +
  labs(
    x = "Follow-up day",
    y = "Survival probability",
    title = "PCR-corrected Kaplan-Meier estimates of time to recurrent parasitaemia by treatment arm",
    caption = "Tick marks = censored observations" #Shaded bands = 95% CI ·
  )+
  scale_x_continuous(
    breaks = c(0, 1, 2, 3, 7, 14, 21, 28, 35, 42)
  )+
  theme_ggsurvfit_default()
cor_recurrence_km_plot_arm

# save the arm km pcr-corrected plot
ggsave(
  filename = "figures/pcr-corrected_km_arm.png",
  plot = cor_recurrence_km_plot_arm,
  height = 8,
  width = 10,
  dpi = 300
)

#===============KM: TIME TO RECURRENT PARASITAEMIA BY AGE GROUP=============
cor_recurrence_km_fit_age <- survfit2(
  Surv(time_to_recurrence, pcr_recurrence_event) ~ age_group, 
  data = recurrence_km_data
  #na.action = na.omit
)
cor_recurrence_km_fit_age

#----- estimated survival probability at Day 42
summary(
  cor_recurrence_km_fit_age,
  times = 42
)$surv


#-----------corrected age_group recurrence log-rank test-------------------
cor_recurrence_logrank_age <- survdiff(
  Surv(time_to_recurrence, pcr_recurrence_event) ~ age_group, 
  data = recurrence_km_data
)
cor_recurrence_logrank_age

#------- corrected recurrence p-value
cor_recurrence_pval_age <- pchisq(
  cor_recurrence_logrank_age$chisq,
  df = length(cor_recurrence_logrank_age$n) -1,
  lower.tail = FALSE
)
cor_recurrence_pval_age

cor_recurrence_label_age <- paste0(
  "Log-rank test: \u03C7\u00B2 = ", round(cor_recurrence_logrank_age$chisq, 2),
  ", p = ", format.pval(cor_recurrence_pval_age, digits = 3, eps = 0.001)
)

#---------------------KM recurrence by age_group plot---------------------------
cor_recurrence_km_plot_age <-cor_recurrence_km_fit_age |> 
  ggsurvfit(linewidth = 0.8) +
  #add_confidence_interval() +
  add_censor_mark()+
  add_risktable(
    risktable_group = "strata",
    risktable_height = 0.25
  ) +
  annotate(
    "text",
    x = 30,
    y = 0.25,
    label = cor_recurrence_label_age, 
    hjust = 0,
    size = 3.5
  ) +
  labs(
    x = "Follow-up day",
    y = "Survival probability",
    title = "PCR-corrected Kaplan-Meier estimates of time to recurrent parasitaemia by age group",
    caption = "Tick marks = censored observations" #Shaded bands = 95% CI ·
  )+
  scale_x_continuous(
    breaks = c(0, 1, 2, 3, 7, 14, 21, 28, 35, 42)
  )+
  theme_ggsurvfit_default()
cor_recurrence_km_plot_age

# save the age km pcr-corrected plot
ggsave(
  filename = "figures/pcr-corrected_km_age.png",
  plot = cor_recurrence_km_plot_age,
  height = 8,
  width = 10,
  dpi = 300
)

#--------Check primary PCR-corrected analysis population---------------
recurrence_km_data |>
  ungroup() |> 
  summarise(
    total = n(),
    pcr_included = sum(!is.na(pcr_recurrence_event)),
    pcr_excluded = sum(is.na(pcr_recurrence_event)),
    recrudescence_events = sum(
      pcr_recurrence_event == 1,
      na.rm = TRUE
    ),
    missing_age_group = sum(is.na(age_group))
  )

cor_recurrence_km_fit_arm$na.action


#============= 42 DAYS SURVIVAL TABLE FOR PCR UNCORRECTED AND CORRECTED KM FIT====

#-------42 days survival table by arm-------------
#----pcr-uncorrected 

tbl1_arm_uncorrected <- survfit2(
  Surv(time_to_recurrence, recurrence_event) ~ treatment_arm, 
  data = recurrence_km_data
) |>
  tbl_survfit(
    times       = 42,
    label_header = "**PCR-uncorrected (95% CI)**"
  )

#----pcr-corrected 
tbl2_arm_corrected <- survfit2(
  Surv(time_to_recurrence, pcr_recurrence_event)~ treatment_arm,
  data = recurrence_km_data
) |> 
  tbl_survfit(
    times = 42,
    label_header = "**PCR-corrected (95% CI)**"
  )

tbl3_age_uncorrected <- survfit2(
  Surv(time_to_recurrence, recurrence_event) ~ age_group, 
  data = recurrence_km_data
) |>
  tbl_survfit(
    times       = 42,
    label_header = "**PCR-uncorrected (95% CI)**"
  )

#----pcr-corrected 
tbl4_age_corrected <- survfit2(
  Surv(time_to_recurrence, pcr_recurrence_event)~ age_group,
  data = recurrence_km_data
) |> 
  tbl_survfit(
    times = 42,
    label_header = "**PCR-corrected (95% CI)**"
  )

recurrence_by_arm <- tbl_merge(list(tbl1_arm_uncorrected, tbl2_arm_corrected ))

recurrence_by_age <- tbl_merge(list(tbl3_age_uncorrected, tbl4_age_corrected ))


recurrence_by_arm <- recurrence_by_arm |> 
  as_gt() |> 
  gt::tab_header(
    title = gt::md("**Day 42 survival by treatment arm**")
  ) 
recurrence_by_arm  

recurrence_by_age <- recurrence_by_age|> 
  as_gt() |>
  gt::tab_header(
    title = gt::md("**Day 42 survival by age group**")
  ) 
recurrence_by_age
 
#------------Save merged recurrence tables  
gt::gtsave(
  recurrence_by_arm,
  filename = "figures/km_survival_table_by_arm.png"
)

gt::gtsave(
  recurrence_by_age,
    filename = "figures/km_survival_table_by_age.png"
  )


#========2. Sensitivity 1: Treat indeterminate/missing PCR as censored (0) =====
recurrence_km_data$pcr_event_cens <- ifelse(
  is.na(recurrence_km_data$pcr_recurrence_event), 0L, 
  recurrence_km_data$pcr_recurrence_event
)

#============CENS KM: TIME TO RECURRENT PARASITAEMIA BY TREATMENT ARM=============
#cor means corrected, cens means censored 

cor_recurrence_km_fit_arm_cens <- survfit2(
  Surv(time_to_recurrence, pcr_event_cens) ~ treatment_arm, 
  data = recurrence_km_data
)
cor_recurrence_km_fit_arm_cens

#----- estimated survival probability at Day 42
summary(
  cor_recurrence_km_fit_arm_cens,
  times = 42
)$surv


#-----------cens arm corrected recurrence log-rank test-------------------
cor_recurrence_logrank_arm_cens <- survdiff(
  Surv(time_to_recurrence, pcr_event_cens) ~ treatment_arm, 
  data = recurrence_km_data
)
cor_recurrence_logrank_arm_cens

#-------cens corrected recurrence p-value
cor_recurrence_pval_arm_cens <- pchisq(
  cor_recurrence_logrank_arm_cens$chisq,
  df = length(cor_recurrence_logrank_arm_cens$n) -1,
  lower.tail = FALSE
)
cor_recurrence_pval_arm_cens

cor_recurrence_label_arm_cens <- paste0(
  "Log-rank test: \u03C7\u00B2 = ", round(cor_recurrence_logrank_arm_cens$chisq, 2),
  ", p = ", format.pval(cor_recurrence_pval_arm_cens, digits = 3, eps = 0.001)
)

#-----------------cens KM recurrence by arm plot---------------------------------
cor_recurrence_km_plot_arm_cens <- cor_recurrence_km_fit_arm_cens |> 
  ggsurvfit(linewidth = 0.8) +
  #add_confidence_interval() +
  add_censor_mark()+
  add_risktable(
    risktable_group = "strata",
    risktable_height = 0.25
  ) +
  annotate(
    "text",
    x = 30,
    y = 0.25,
    label = cor_recurrence_label_arm_cens, 
    hjust = 0,
    size = 3.5
  ) +
  labs(
    x = "Follow-up day",
    y = "Survival probability",
    title = "PCR-corrected Kaplan-Meier estimates of time to recurrent parasitaemia by treatment arm",
    subtitle = "Sensitivity analysis: missing/indeterminate PCR treated as censored",
    caption = "Tick marks = censored observations" #Shaded bands = 95% CI ·
  )+
  scale_x_continuous(
    breaks = c(0, 1, 2, 3, 7, 14, 21, 28, 35, 42)
  )+
  theme_ggsurvfit_default()
cor_recurrence_km_plot_arm_cens

# save the censored arm km pcr-corrected plot
ggsave(
  filename = "figures/pcr-corrected_km_arm_cens.png",
  plot = cor_recurrence_km_plot_arm_cens,
  height = 8,
  width = 10,
  dpi = 300
)

#===============CENS KM: TIME TO RECURRENT PARASITAEMIA BY AGE GROUP=============
cor_recurrence_km_fit_age_cens <- survfit2(
  Surv(time_to_recurrence, pcr_event_cens) ~ age_group, 
  data = recurrence_km_data
)
cor_recurrence_km_fit_age_cens

#----- estimated survival probability at Day 42
summary(
  cor_recurrence_km_fit_age_cens,
  times = 42
)$surv


#----------- cens corrected age_group recurrence log-rank test-------------------
cor_recurrence_logrank_age_cens <- survdiff(
  Surv(time_to_recurrence, pcr_event_cens) ~ age_group, 
  data = recurrence_km_data
)
cor_recurrence_logrank_age_cens

#-------cens corrected recurrence p-value
cor_recurrence_pval_age_cens <- pchisq(
  cor_recurrence_logrank_age_cens$chisq,
  df = length(cor_recurrence_logrank_age_cens$n) -1,
  lower.tail = FALSE
)
cor_recurrence_pval_age_cens

cor_recurrence_label_age_cens <- paste0(
  "Log-rank test: \u03C7\u00B2 = ", round(cor_recurrence_logrank_age_cens$chisq, 2),
  ", p = ", format.pval(cor_recurrence_pval_age_cens, digits = 3, eps = 0.001)
)

#-----------------cens KM recurrence by age_group plot---------------------------
cor_recurrence_km_plot_age_cens <-cor_recurrence_km_fit_age_cens |> 
  ggsurvfit(linewidth = 0.8) +
  #add_confidence_interval() +
  add_censor_mark()+
  add_risktable(
    risktable_group = "strata",
    risktable_height = 0.25
  ) +
  annotate(
    "text",
    x = 30,
    y = 0.25,
    label = cor_recurrence_label_age_cens, 
    hjust = 0,
    size = 3.5
  ) +
  labs(
    x = "Follow-up day",
    y = "Survival probability",
    title = "PCR-corrected Kaplan-Meier estimates of time to recurrent parasitaemia by age group",
    subtitle = "Sensitivity analysis: missing/indeterminate PCR treated as censored",
    caption = "Tick marks = censored observations" #Shaded bands = 95% CI ·
  )+
  scale_x_continuous(
    breaks = c(0, 1, 2, 3, 7, 14, 21, 28, 35, 42)
  )+
  theme_ggsurvfit_default()
cor_recurrence_km_plot_age_cens

# save the censored age km pcr-corrected plot
ggsave(
  filename = "figures/pcr-corrected_km_age_cens.png",
  plot = cor_recurrence_km_plot_age_cens,
  height = 8,
  width = 10,
  dpi = 300
)


#====3. Sensitivity 2: Treat indeterminate/missing PCR as recrudescence/event (1)====

recurrence_km_data$pcr_event_recrud <- ifelse(
  is.na(recurrence_km_data$pcr_recurrence_event), 1L, 
  recurrence_km_data$pcr_recurrence_event
)

#===========RECRUD KM: TIME TO RECURRENT PARASITAEMIA BY TREATMENT ARM=============
#cor means corrected, recrud means recrudescence

cor_recurrence_km_fit_arm_recrud <- survfit2(
  Surv(time_to_recurrence, pcr_event_recrud) ~ treatment_arm, 
  data = recurrence_km_data
)
cor_recurrence_km_fit_arm_recrud

#----- estimated survival probability at Day 42
summary(
  cor_recurrence_km_fit_arm_recrud,
  times = 42
)$surv


#-----------recrud arm corrected recurrence log-rank test-------------------
cor_recurrence_logrank_arm_recrud <- survdiff(
  Surv(time_to_recurrence, pcr_event_recrud) ~ treatment_arm, 
  data = recurrence_km_data
)
cor_recurrence_logrank_arm_recrud

#-------corrected recurrence p-value
cor_recurrence_pval_arm_recrud <- pchisq(
  cor_recurrence_logrank_arm_recrud$chisq,
  df = length(cor_recurrence_logrank_arm_recrud$n) -1,
  lower.tail = FALSE
)
cor_recurrence_pval_arm_recrud

cor_recurrence_label_arm_recrud <- paste0(
  "Log-rank test: \u03C7\u00B2 = ", round(cor_recurrence_logrank_arm_recrud$chisq, 2),
  ", p = ", format.pval(cor_recurrence_pval_arm_recrud, digits = 3, eps = 0.001)
)

#----------------recrud KM recurrence by arm plot---------------------------------
cor_recurrence_km_plot_arm_recrud <- cor_recurrence_km_fit_arm_recrud |> 
  ggsurvfit(linewidth = 0.8) +
  #add_confidence_interval() +
  add_censor_mark()+
  add_risktable(
    risktable_group = "strata",
    risktable_height = 0.25
  ) +
  annotate(
    "text",
    x = 30,
    y = 0.25,
    label = cor_recurrence_label_arm_recrud, 
    hjust = 0,
    size = 3.5
  ) +
  labs(
    x = "Follow-up day",
    y = "Survival probability",
    title = "PCR-corrected Kaplan-Meier estimates of time to recurrent parasitaemia by treatment arm",
    subtitle = "Sensitivity analysis: missing/indeterminate PCR treated as recrudescence",
    caption = "Tick marks = censored observations" #Shaded bands = 95% CI ·
  )+
  scale_x_continuous(
    breaks = c(0, 1, 2, 3, 7, 14, 21, 28, 35, 42)
  )+
  theme_ggsurvfit_default()
cor_recurrence_km_plot_arm_recrud

# save the recrudescence arm km pcr-corrected plot
ggsave(
  filename = "figures/pcr-corrected_km_arm_recrud.png",
  plot = cor_recurrence_km_plot_arm_recrud,
  height = 8,
  width = 10,
  dpi = 300
)


#=============RECRUD KM: TIME TO RECURRENT PARASITAEMIA BY AGE GROUP=============
cor_recurrence_km_fit_age_recrud <- survfit2(
  Surv(time_to_recurrence, pcr_event_recrud) ~ age_group, 
  data = recurrence_km_data
)
cor_recurrence_km_fit_age_recrud

#-----estimated survival probability at Day 42
summary(
  cor_recurrence_km_fit_age_recrud,
  times = 42
)$surv


#--------recrud corrected age_group recurrence log-rank test-------------------
cor_recurrence_logrank_age_recrud <- survdiff(
  Surv(time_to_recurrence, pcr_event_recrud) ~ age_group, 
  data = recurrence_km_data
)
cor_recurrence_logrank_age_recrud

#------- corrected recurrence p-value
cor_recurrence_pval_age_recrud <- pchisq(
  cor_recurrence_logrank_age_recrud$chisq,
  df = length(cor_recurrence_logrank_age_recrud$n) -1,
  lower.tail = FALSE
)
cor_recurrence_pval_age_recrud

cor_recurrence_label_age_recrud <- paste0(
  "Log-rank test: \u03C7\u00B2 = ", round(cor_recurrence_logrank_age_recrud$chisq, 2),
  ", p = ", format.pval(cor_recurrence_pval_age_recrud, digits = 3, eps = 0.001)
)

#-------------------recrud KM recurrence by age_group plot----------------------

cor_recurrence_km_plot_age_recrud <-cor_recurrence_km_fit_age_recrud |> 
  ggsurvfit(linewidth = 0.8) +
  #add_confidence_interval() +
  add_censor_mark()+
  add_risktable(
    risktable_group = "strata",
    risktable_height = 0.25
  ) +
  annotate(
    "text",
    x = 30,
    y = 0.25,
    label = cor_recurrence_label_age_recrud, 
    hjust = 0,
    size = 3.5
  ) +
  labs(
    x = "Follow-up day",
    y = "Survival probability",
    title = "PCR-corrected Kaplan-Meier estimates of time to recurrent parasitaemia by age group",
    subtitle = "Sensitivity analysis: missing/indeterminate PCR treated as recrudescence",
    caption = "Tick marks = censored observations" #Shaded bands = 95% CI ·
  )+
  scale_x_continuous(
    breaks = c(0, 1, 2, 3, 7, 14, 21, 28, 35, 42)
  )+
  theme_ggsurvfit_default()
cor_recurrence_km_plot_age_recrud

# save the recrudescence arm km pcr-corrected plot
ggsave(
  filename = "figures/pcr-corrected_km_age_recrud.png",
  plot = cor_recurrence_km_plot_age_recrud,
  height = 8,
  width = 10,
  dpi = 300
)


#============= 42 DAYS SURVIVAL TABLE FOR PCR CORRECTED, CENS AND RECRUD KM FIT====

#-------42 days survival table by arm-------------
#----pcr-uncorrected 

tbl5_arm_cens <- survfit2(
  Surv(time_to_recurrence, pcr_event_cens) ~ treatment_arm, 
  data = recurrence_km_data
) |>
  tbl_survfit(
    times       = 42,
    label_header = "**Indeterminate/Missing PCR as censored (95% CI)**"
  )

#----pcr-corrected 
tbl6_arm_recrud <- survfit2(
  Surv(time_to_recurrence, pcr_event_recrud)~ treatment_arm,
  data = recurrence_km_data
) |> 
  tbl_survfit(
    times = 42,
    label_header = "**Indeterminate/Missing PCR as recrudescence (95% CI)**"
  )

tbl7_age_cens <- survfit2(
  Surv(time_to_recurrence, pcr_event_cens) ~ age_group, 
  data = recurrence_km_data
) |>
  tbl_survfit(
    times       = 42,
    label_header = "**Indeterminate/Missing PCR as censored (95% CI)**"
  )

#----pcr-corrected 
tbl8_age_recrud <- survfit2(
  Surv(time_to_recurrence, pcr_event_recrud)~ age_group,
  data = recurrence_km_data
) |> 
  tbl_survfit(
    times = 42,
    label_header = "**Indeterminate/Missing PCR as recrudescence (95% CI)**"
  )

sensitivity_by_arm <- tbl_merge(list(tbl2_arm_corrected, tbl5_arm_cens,tbl6_arm_recrud))

sensitivity_by_age <- tbl_merge(list(tbl4_age_corrected, tbl7_age_cens, tbl8_age_recrud))


sensitivity_by_arm <- sensitivity_by_arm |> 
  as_gt() |> 
  gt::tab_header(
    title = gt::md("**PCR-corrected day 42 survival by treatment arm**")
  ) 
sensitivity_by_arm 

sensitivity_by_age  <- sensitivity_by_age |> 
  as_gt() |>
  gt::tab_header(
    title = gt::md("**PCR-corrected day 42 survival by age group**")
  ) 
sensitivity_by_age 

#------------Save merged recurrence tables  
gt::gtsave(
  sensitivity_by_arm,
  filename = "figures/sensitivity_survival_table_by_arm.png"
)

gt::gtsave(
  sensitivity_by_age ,
  filename = "figures/sensitivity_survival_table_by_age.png"
)


View(treatment_failure_participant)

treatment_failure_participant |> 
  anti_join(
    recurrence_event,
    by = "participant_id"
  ) 
  #print(n = Inf)

View(treatment_failure_participant)
treatment_failure_check|> 
  filter(participant_id %in% c("S24-10048-AG", "S24-10049-AG")) |> 
  select(
    participant_id,
    #study_visit,
    tf_visit,
    tf_microscopy_result
  )

#======================== PARASITE CLEARANCE ESTIMATOR==========================
# CALCULATING PARASITE CLEARANCE LIFE

pce_participants <- smart_analysis_final |>
  ungroup() |>
  semi_join(
    smart_analysis_final |>
      filter(
        study_visit == "H0",
        parasite_density > 0
      ) |>
      distinct(participant_id),
    by = "participant_id"
  )

pce_participants|> 
  select(
    participant_id,
    study_visit,
    treatment_arm,
    site,
    parasite_density
  ) |> 
  View()




# #================LOG-SCALE PARASITE DENSITY DECLINE====
# 
# parasite_density_log <- clearance_data |>
#   filter(
#     !is.na(parasite_density),
#     parasite_density > 0
#   ) |>
#   group_by(
#     treatment_arm,
#     clearance_hour
#   ) |>
#   summarise(
#     n_positive = n(),
#     geometric_mean = exp(mean(log(parasite_density))),
#     .groups = "drop"
#   )
# 
# parasite_density_log
# 
# parasite_density_log_plot <- parasite_density_log |>
#   ggplot(
#     aes(
#       x = clearance_hour,
#       y = geometric_mean,
#       colour = treatment_arm,
#       group = treatment_arm
#     )
#   ) +
#   geom_line(linewidth = 1.2) +
#   geom_point(size = 3) +
#   scale_x_continuous(
#     breaks = c(0, 8, 24, 36, 48, 60, 72),
#     labels = c("H0", "H8", "H24", "H36", "H48", "H60", "H72")
#   ) +
#   scale_y_log10(
#     labels = scales::comma
#   ) +
#   labs(
#     title = "Decline in parasite density during the first 72 hours after treatment",
#     subtitle = "Geometric mean parasite density among positive measurements",
#     x = "Time after treatment",
#     y = "Parasite density (parasites/µL, log scale)",
#     colour = "Treatment arm"
#   ) +
#   theme_minimal(base_size = 13) +
#   theme(
#     plot.title = element_text(face = "bold"),
#     legend.position = "bottom",
#     panel.grid.minor = element_blank()
#   )
# 
# parasite_density_log_plot


# #================RAW PARASITE DENSITY DECLINE==================================
# 
# raw_parasite_density_plot <- clearance_data |>
#   filter(
#     !is.na(parasite_density)
#   ) |>
#   ggplot(
#     aes(
#       x = clearance_hour,
#       y = parasite_density,
#       colour = treatment_arm
#     )
#   ) +
#   geom_point(
#     alpha = 0.25,
#     size = 1.5,
#     position = position_jitter(width = 1.2, height = 0)
#   ) +
#   scale_x_continuous(
#     breaks = c(0, 8, 24, 36, 48, 60, 72),
#     labels = c(
#       "H0", "H8", "H24",
#       "H36", "H48", "H60", "H72"
#     )
#   ) +
#   scale_y_continuous(
#     trans = "log1p",
#     labels = scales::comma
#   ) +
#   labs(
#     title = "Parasite density during the first 72 hours after treatment",
#     subtitle = "Observed parasite-density measurements by treatment arm",
#     x = "Time after treatment",
#     y = "Parasite density (parasites/µL, log scale)",
#     colour = "Treatment arm"
#   ) +
#   theme_minimal(base_size = 13) +
#   theme(
#     plot.title = element_text(face = "bold"),
#     legend.position = "bottom",
#     panel.grid.minor = element_blank()
#   )
# 
# raw_parasite_density_plot

# raw_parasite_density_plot +
#   facet_wrap(~ treatment_arm)
