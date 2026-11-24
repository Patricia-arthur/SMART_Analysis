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



#======================== PARASITE CLEARANCE ESTIMATOR =========================

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


#---------------- Prepare PCE upload file --------------------------------------
#----------AP participants
pce_AP_data <- pce_participants |>
  filter(
    treatment_arm == "AP",
    study_visit %in% c(
      "H0", "H8", "H24", "H36",
      "H48", "H60", "H72"
    )
  ) |>
  transmute(
    id = participant_id,
    time = parse_number(study_visit),
    count = as.numeric(parasite_density)
  ) |>
  arrange(
    id,
    time
  )


#check
pce_AP_data |>
  summarise(
    participants = n_distinct(id),
    records = n(),
    missing_count = sum(is.na(count))
  )

pce_AP_data |>
  group_by(time) |>
  summarise(
    participants = n_distinct(id),
    records = n(),
    missing_count = sum(is.na(count)),
    .groups = "drop"
  ) |>
  arrange(time)

pce_AP_data |>
  count(time)


pce_AP_data <-  pce_AP_data|>
  filter(
    !is.na(count) 
  ) 
write_csv(pce_AP_data, "Data/PCE_AP_data.csv")

#-----------AL Participants
pce_AL_data <- pce_participants |>
  filter(
    treatment_arm == "AL",
    study_visit %in% c(
      "H0", "H8", "H24", "H36",
      "H48", "H60", "H72"
    )
  ) |>
  transmute(
    id = participant_id,
    time = parse_number(study_visit),
    count = as.numeric(parasite_density)
  ) |>
  arrange(
    id,
    time
  )


#check
pce_AL_data |>
  summarise(
    participants = n_distinct(id),
    records = n(),
    missing_count = sum(is.na(count))
  )

pce_AL_data |>
  group_by(time) |>
  summarise(
    participants = n_distinct(id),
    records = n(),
    missing_count = sum(is.na(count)),
    .groups = "drop"
  ) |>
  arrange(time)

pce_AL_data |>
  count(time)

pce_AL_data <-  pce_AL_data|>
  filter(
    !is.na(count) 
  ) 
write_csv(pce_AL_data, "Data/PCE_AL_data.csv")

pce_AL_data |> 
  filter(id == "S24-10016-AG") 


#------------check for NA of counts at H8
H8_no_count <- SMARTUncomplicated |>
  filter(redcap_event_name == "hour_8_day_0_arm_1") |>
  select(
    id_um,
    redcap_event_name,
    da_ap_b01_1,
    da_ap_b01_2,
    da_al_c07_v2,
    da_al_b01_1,
    da_al_b01_2,
    da_al_c07
  ) |>
  mutate(
    across(id_um, as.character),
  ) |>
  rename(
    screening_id = id_um, 
    study_visit = redcap_event_name,
    study_drug_given_AP = da_ap_b01_1,
    if_no_specify_AP = da_ap_b01_2,
    remarks_AP = da_al_c07_v2,
    study_drug_given_AL = da_al_b01_1,
    if_no_specify_AL = da_al_b01_2,
    remarks_AL = da_al_c07
  ) |> 
  mutate(
    study_drug_given_AP = factor(
      study_drug_given_AP,
      levels = c(0,1),
      labels = c("NO", "Yes")
    ),
    
    study_drug_given_AL = factor(
      study_drug_given_AL,
      levels = c(0,1),
      labels = c("NO", "Yes")
    ),
    
    study_visit = factor(
      study_visit,
      levels = "hour_8_day_0_arm_1",
      labels = "H8"
    )
  )

H8_count <- pce_participants |> 
  filter(
    study_visit == "H8",
    is.na(parasite_density)
  ) |> 
  select(
    participant_id,
    screening_id,
    treatment_arm,
    parasite_density
  )

H8_no_count <- H8_no_count |>
  semi_join(
    H8_count, 
    by = "screening_id"
  )

H8_no_count <- H8_no_count |>
  left_join(
    H8_count, 
    by = "screening_id"
  )

pce_participants%>% has_name("screening_id")




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
