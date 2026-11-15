source("05_eda.R")
library(ggsurvfit)

#========================= FEVER RESOLUTION ANALYSIS============================

#check number of participants with fever at H0
smart_analysis_final |> 
  ungroup() |> 
  filter(study_visit == "H0") |> 
  summarise(
    n_participant = n_distinct(participant_id),
    n_febrile = sum(temperature >= 37.5, na.rm = TRUE),
    n_afebrile = sum(temperature < 37.5, na.rm = TRUE),
    n_missing_temperature = sum(is.na(temperature)),
    percentage = round(100 * n_febrile / n_participant, 2)
  )

#check within treatment arm
smart_analysis_final |> 
  ungroup() |> 
  filter(study_visit == "H0") |> 
  group_by(treatment_arm) |> 
  summarise(
    n = n(),
    n_febrile = sum(temperature >= 37.5, na.rm = TRUE),
    n_afebrile = sum(temperature < 37.5, na.rm = TRUE),
    .groups = "drop"
  )

#--------------------Fever clearance through H72--------------------------------
early_fever_visits <- tibble(
  study_visit = c(
    "H0", "H8", "H24", "H36",
    "H48", "H60", "H72"
  ),
  analysis_hour = c(
    0, 8, 24, 36,
    48, 60, 72
  )
)

#-----------------------IDs of participants febrile at H0------------------------------
febrile_participants <- smart_analysis_final |>
  ungroup() |>
  filter(
    study_visit == "H0",
    temperature >= 37.5
  ) |>
  distinct(
    participant_id,
    treatment_arm
  )

febrile_participants |> 
  count(treatment_arm)

#check how many participants still have fever at each visit
febrile_ids <- smart_analysis_final |> 
  ungroup() |> 
   filter(
     study_visit == "H0",
     temperature >= 37.5
   ) |> 
   pull(participant_id)

fever_group <- smart_analysis_final |>
  ungroup() |>
  filter(
    participant_id %in% febrile_ids
  )

fever_group |> 
  group_by(
    treatment_arm,
    study_visit
  ) |> 
  summarise(
    n_temperature = sum(!is.na(temperature)),
    n_febrile = sum(temperature >= 37.5, na.rm = TRUE),
    n_afebrile = sum(temperature < 37.5, na.rm = TRUE),
    .groups = "drop"
  ) |> 
  print(n =Inf)

#-----------------build  an early fever dataset---------------------------
fever_data <- tidyr::crossing(
  febrile_participants,
  early_fever_visits
) |>
  left_join(
    smart_analysis_final |>
      ungroup() |>
      select(
        participant_id,
        study_visit,
        temperature
      ),
    by = c(
      "participant_id",
      "study_visit"
    )
  ) |>
  arrange(
    participant_id,
    analysis_hour
  )

#check
fever_data |>
  count(participant_id) |>
  count(n)

fever_data |>
  filter(
    participant_id %in% c(
      "S24-10025-AG",
      "S24-10207-AG"
    )
  ) |>
  select(
    participant_id,
    treatment_arm,
    study_visit,
    analysis_hour,
    temperature
  ) |>
  print(n = Inf)

#---------------------Define fever clearance rule-------------------------------
find_fever_clearance <- function(
    analysis_hour, temperature
){
  #identify possible clearance visit(fever must remain resolve for at least 24hrs)
  candidate <- which(
    !is.na(temperature) &
      temperature < 37.5 &
      analysis_hour < max(analysis_hour) #<=48
  )
  
  if(length(candidate)== 0){
    return (NA_real_)
  }
  
  for (i in candidate) {
     target_hour <- analysis_hour[i] + 24
     
     # Available observations from candidate time onwards
     after_candidate <- which(
       analysis_hour >= analysis_hour[i] &
         !is.na(temperature)
     )
     
     # Need at least one observed temperature
     # 24 hours or more after the candidate
     confirmation <- after_candidate[
       analysis_hour[after_candidate] >= target_hour
     ]
     
     if (length(confirmation) == 0) {
       next
     }
     
     # No observed fever before confirmation
     observations_to_confirmation <- after_candidate[
       analysis_hour[after_candidate] <=
         analysis_hour[confirmation[1]]
     ]
     
     # Fever clearance occurs if all available temperatures
     # remain below 37.5°C
     if (
       all(
         temperature[observations_to_confirmation] < 37.5
       )
     ) {
       return(analysis_hour[i])
     }
  }
  
  NA_real_
}

#--------------Derive fever clearance time for each participant-----------------
fever_clearance <- fever_data |>
  group_by(
    participant_id,
    treatment_arm
  ) |>
  summarise(
    clearance_time = find_fever_clearance(
      analysis_hour,
      temperature
    ),
    .groups = "drop"
  )

# check clearance times
fever_clearance |>
  count(
    treatment_arm,
    clearance_time
  ) |>
  arrange(
    treatment_arm,
    clearance_time
  )

# Identify participants without clearance within H0 - D14
fever_clearance |>
  filter(is.na(clearance_time))

#---------------------Derive fever time and event-------------------------------
# 72 hours is used as censoring time for this current analysis,
fever_clearance <- fever_clearance |>
  mutate(
    fever_event = if_else(
      !is.na(clearance_time),
      1L,
      0L
    ),
    
    fever_time = if_else(
      fever_event == 1,
      clearance_time,
      72
    )
  )

#check
fever_clearance |>
  count(
    treatment_arm,
    fever_event,
    fever_time
  ) |>
  arrange(
    treatment_arm,
    fever_time
  )

#----------------KM fever resolution model--------------------------------------
#---------------------Kaplan-Meier model----------------------------------------
  
fever_km <- survfit2(
  Surv(fever_time, fever_event) ~ treatment_arm,
    data = fever_clearance
  )
  
#----------------------Log-rank test--------------------------------------------
fever_logrank <- survdiff(
  Surv(fever_time, fever_event) ~ treatment_arm,
  data = fever_clearance
)
fever_logrank

#---------------p-value
logrank_chisq <- fever_logrank$chisq

fever_pvalue <- pchisq(logrank_chisq, 
                       df = length(fever_logrank$n)-1,
                       lower.tail = FALSE
                      )

logrank_label <- paste0(
  "Log-rank test:  \u03C7\u00B2 = ",
  round(logrank_chisq, 2),
  ", p = ",
  format.pval(fever_pvalue, digits = 2, eps = 0.001)
)

#---------------------Kaplan-Meier plot-----------------------------------------

fever_km_plot <- fever_km |>
  ggsurvfit(linewidth = 1) +
  #add_confidence_interval() +
  add_censor_mark() +
  annotate(
    "text",
    x = 40,
    y = 0.25,
    label = logrank_label,
    hjust = 0,
    size = 3.5
  ) +
  add_risktable(
    times = c(0, 8, 24, 36, 48, 60, 72),
    #risktable_stats = "n.risk",
    risktable_group = "strata",
    # stats_label = list(
    #   n.risk = "Number at Risk"
    # ),
    risktable_height = 0.25
  ) +
  labs(
    x = "Follow-up time(hours)",
    y = "Survival probability",
    title = "Kaplan-Meier: fever resolution by treatment arm",
    caption = "Tick marks = censored observations" #Shaded bands = 95% CI ·
  ) +
  scale_x_continuous(
    breaks = c(0, 8, 24, 36, 48, 60, 72)
  ) +
  theme_ggsurvfit_default()

fever_km_plot

#------save the plot
ggsave(
  filename = "figures/fever_resolution_km.png",
  plot = fever_km_plot,
  width = 10,
  height = 8,
  dpi = 300
)
