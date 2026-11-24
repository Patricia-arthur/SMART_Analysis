source("10_recurrent_parasitaemia_KM.R")

#===================== COX REGRESSION ANALYSIS =================================

#--------------------------create Cox analysis dataset--------------------------

recurrence_cox_data <- recurrence_km_data |> 
  left_join(
    baseline_data |> 
      select(
        participant_id,
        site,
        sex,
        weight,
        baseline_temperature,
        baseline_haemoglobin,
        baseline_parasitaemia
      ),
    by = "participant_id"
  ) #|> 
  # mutate(
  #   site = factor(
  #     site,
  #     levels = c("Assin Fosu", "Agogo")
  #   )
  # )

#----check the cox dataset 
recurrence_cox_data |> 
  summarise(
    n_participants = n(),
    n_unique = n_distinct(participant_id),
    n_events = sum(recurrence_event == 1)
  )

# check missingness
recurrence_cox_data |>
  summarise(
    missing_treatment_arm = sum(is.na(treatment_arm)),
    missing_age = sum(is.na(age)),
    missing_age_group = sum(is.na(age_group)),
    missing_sex = sum(is.na(sex)),
    missing_weight = sum(is.na(weight)),
    missing_site = sum(is.na(site)),
    missing_temperature = sum(is.na(baseline_temperature)),
    missing_haemoglobin = sum(is.na(baseline_haemoglobin)),
    missing_parasitaemia = sum(is.na(baseline_parasitaemia)),
    missing_time = sum(is.na(time_to_recurrence)),
    missing_event = sum(is.na(recurrence_event))
  )

#-------check how variables are coded
recurrence_cox_data |>
  count(treatment_arm)

recurrence_cox_data |>
  count(age_group)

recurrence_cox_data |>
  count(sex)

recurrence_cox_data |>
  count(site)

View(recurrence_cox_data)


# check parasitaemia
recurrence_cox_data |>
  summarise(
    n_missing = sum(is.na(baseline_parasitaemia)),
    n_zero = sum(baseline_parasitaemia == 0, na.rm = TRUE),
    n_positive = sum(baseline_parasitaemia > 0, na.rm = TRUE)
  )

#================== Create a positive baseline parasitamia data=================
recurrence_cox_parasitaemia <- recurrence_cox_data |> 
  filter(
    !is.na(baseline_parasitaemia),
    baseline_parasitaemia > 0
  ) |> 
  mutate(
    log_parasitaemia = log10(baseline_parasitaemia)
  )

#---check
recurrence_cox_parasitaemia |>
  summarise(
    n_participants = n(),
    n_events = sum(recurrence_event == 1),
    min_parasitaemia = min(baseline_parasitaemia),
    median_parasitaemia = median(baseline_parasitaemia),
    max_parasitaemia = max(baseline_parasitaemia)
  )


#-------------visualise the raw versus log-transformed parasitaemia
ggplot(
  recurrence_cox_parasitaemia,
  aes(x = baseline_parasitaemia)
) +
  geom_histogram(bins = 30)+
  labs(
    x = "Baseline parasite density (parasite/µL)",
    y = "Frequency",
    title = "Distribution of baseline parasitaemia"
  ) +
  theme_minimal()

ggplot(
  recurrence_cox_parasitaemia,
  aes(x = log_parasitaemia)
) +
  geom_histogram(bins = 30)+
  labs(
    x = "Baseline parasite density (parasite/µL)",
    y = "Frequency",
    title = "Distribution of baseline parasitaemia"
  ) +
  theme_minimal()

#============================ THE COX MODEL=====================================

#============================= UNIVARIABLE COX =================================

#----------treatment arm
cox_arm <- coxph(
  Surv(time_to_recurrence, recurrence_event) ~ treatment_arm,
  data = recurrence_cox_data
)
summary(cox_arm)

# fitted to valid parasitaemia
cox_arm_valid <- coxph(
  Surv(time_to_recurrence, recurrence_event) ~ treatment_arm,
  data = recurrence_cox_parasitaemia
)
summary(cox_arm_valid)

#----------age
cox_age <- coxph(
  Surv(time_to_recurrence, recurrence_event) ~ age,
  data = recurrence_cox_data
)
summary(cox_age)

#----------age_group
cox_age_group <- coxph(
  Surv(time_to_recurrence, recurrence_event) ~ age_group,
  data = recurrence_cox_data
)
summary(cox_age_group)

#--------------parasitaemia
cox_parasitaemia <- coxph(
  Surv(time_to_recurrence, recurrence_event) ~ log_parasitaemia,
  data = recurrence_cox_parasitaemia
)
summary(cox_parasitaemia)

#----------haemoglobin
cox_haemoglobin <- coxph(
  Surv(time_to_recurrence, recurrence_event) ~ baseline_haemoglobin,
  data = recurrence_cox_data
)
summary(cox_haemoglobin)

#----------site
cox_site <- coxph(
  Surv(time_to_recurrence, recurrence_event) ~ site,
  data = recurrence_cox_data
)
summary(cox_site)

cox_site_valid <- coxph(
  Surv(time_to_recurrence, recurrence_event) ~ site,
  data = recurrence_cox_parasitaemia
)
summary(cox_site_valid)

#----------sex
cox_sex <- coxph(
  Surv(time_to_recurrence, recurrence_event) ~ sex,
  data = recurrence_cox_data
)
summary(cox_sex)

#-------------temperature
cox_temperature <- coxph(
  Surv(time_to_recurrence, recurrence_event) ~ baseline_temperature,
  data = recurrence_cox_data
)
summary(cox_temperature)

cox_weight <- coxph(
  Surv(time_to_recurrence, recurrence_event) ~ weight,
  data = recurrence_cox_data
)
summary(cox_weight)


#===================== MULTIVARIABLE COX MODEL =================================

#---------------- Model 1: Full cohort -----------------------------------------
# Uses all participants  but baseline parasitaemia excluded 

cox_full<- coxph(
  Surv(time_to_recurrence, recurrence_event) ~
    treatment_arm +
    age +
    baseline_haemoglobin +
    site,
  data = recurrence_cox_data
)

summary(cox_full)

#----------------------- Model 2: restricted cohort ----------------------------
# restricted to participant with parasitaemia >0 but parasitaemia still excluded

cox_restricted <- coxph(
  Surv(time_to_recurrence, recurrence_event) ~
    treatment_arm +
    age +
    baseline_haemoglobin +
    site,
  data = recurrence_cox_parasitaemia
)

summary(cox_restricted)

#---------------- Model 3: Valid baseline parasitaemia -------------------------
# Uses participants with baseline_parasitaemia > 0

cox_multi <- coxph(
  Surv(time_to_recurrence, recurrence_event) ~
    treatment_arm +
    age +
    log_parasitaemia +
    baseline_haemoglobin +
    site,
  data = recurrence_cox_parasitaemia
)

summary(cox_multi)

#----------check total participants and events in each model
cox_multi$n
cox_multi$nevent

cox_full$n
cox_full$nevent

cox_restricted$n
cox_restricted$nevent


#========================== THE PH ASSUMPTION ==================================

#------------------------------PH for model 1-----------------------------------
# all participants but parasitaemia excluded
ph_cox_full <- cox.zph(cox_full)
ph_cox_full 

#--------------------------------PH for model 2---------------------------------
# restricited participants butc parasitaemia excluded
ph_cox_restricted <- cox.zph(cox_restricted)
ph_cox_restricted

#--------------------------------PH for model 3---------------------------------
#only valid parasitaemia included
ph_cox_multi <-cox.zph(cox_multi)
ph_cox_multi

#--------------------Visualise the PH assumptions ------------------------------
#----model 1
ggcoxzph(
  ph_cox_full,
  point.alpha = 0.3,
  point.col = "red",
  caption =  "Shaded band = 95% CI on Schoenfeld residuals · Flat line supports PH assumption"
)

#-------model 2
ggcoxzph(
  ph_cox_restricted,
  point.alpha = 0.3,
  point.col = "blue",
  caption =  "Shaded band = 95% CI on Schoenfeld residuals · Flat line supports PH assumption"
)

#-------model 3
ggcoxzph(
  ph_cox_multi,
  point.alpha = 0.3,
  point.col = "violet",
  caption =  "Shaded band = 95% CI on Schoenfeld residuals · Flat line supports PH assumption"
)


#--------------Test if a variable improves the model
anova(
  cox_restricted,
  cox_multi,
  test = "LRT"
)


#===================== TABLES for HAZARD RATIOS ================================

#------------HR for model 3
hr_cox_multi <- cox_multi |> 
  tbl_regression(
    exp = TRUE,
    label = list(
      treatment_arm ~ "Treatment arm",
      age ~ "Age(months)",
      log_parasitaemia ~ "Parasite density (/µL)",
      baseline_haemoglobin ~ "Haemoglobin (g/dL)"
    )
  ) |> 
  bold_p(t = 0.05)

hr_cox_multi <- hr_cox_multi |> 
  as_gt() |> 
  tab_header(
    title = md("**Adjusted hazard ratio estimates**")
  )
hr_cox_multi 

#----save the table
gtsave(
  hr_cox_multi,
  filename = "figures/hr_cox_multi.png"
)

#------------HR for model 2
hr_cox_restricted <- cox_restricted |> 
  tbl_regression(
    exp = TRUE,
    label = list(
      treatment_arm ~ "Treatment arm",
      age ~ "Age(months)",
      baseline_haemoglobin ~ "Haemoglobin (g/dL)"
    )
  ) |> 
  bold_p(t = 0.05)
hr_cox_restricted

hr_cox_restricted <- hr_cox_restricted |> 
  as_gt() |> 
  tab_header(
    title = md("**Adjusted hazard ratio estimates**")
  )
hr_cox_restricted 

#----save the table
gtsave(
  hr_cox_restricted,
  filename = "figures/hr_cox_restricted.png"
)

#------------HR for model 1
hr_cox_full <- cox_full |> 
  tbl_regression(
    exp = TRUE,
    label = list(
      treatment_arm ~ "Treatment arm",
      age ~ "Age(months)",
      baseline_haemoglobin ~ "Haemoglobin (g/dL)"
    )
  ) |> 
  bold_p(t = 0.05)
hr_cox_full

hr_cox_full <- hr_cox_full|> 
  as_gt() |> 
  tab_header(
    title = md("**Adjusted hazard ratio estimates**")
  )
hr_cox_full

#----save the table
gtsave(
  hr_cox_full,
  filename = "figures/hr_cox_full.png"
)


#------------univarable models
hr_cox_arm_valid <- cox_arm_valid |> 
  tbl_regression(
    exp = TRUE,
    include = treatment_arm,
    label = list(
      treatment_arm ~ "Treatment arm"
    )
  ) |> 
  bold_p(t = 0.05)
hr_cox_arm_valid

hr_cox_parasitaemia <- cox_parasitaemia |> 
  tbl_regression(
    exp = TRUE,
    include = log_parasitaemia,
    label = list(
      log_parasitaemia ~ "Parasite density (/µL)"
    )
  ) |> 
  bold_p(t = 0.05)
hr_cox_parasitaemia

hr_cox_site <- cox_site_valid |> 
  tbl_regression(
    exp = TRUE,
    include = site,
  ) |> 
  bold_p(t = 0.05)
hr_cox_site

#-----------------Combine Adjusted and unadjusted tables-------------------------
#---------------------with treatment arm
tbl_adj_arm <- cox_multi |>
  tbl_regression(
    exp = TRUE,
    include = treatment_arm,
    label = list(
      treatment_arm ~ "Treatment arm"
    )
  ) |> 
  bold_p(t = 0.05) |>
  modify_header(estimate = "**Adjusted HR**")

tbl_unadj_arm <- hr_cox_arm_valid |> 
  modify_header(estimate = "**Unadjusted HR**")

# Merge the tables
tbl_arm_comb <- tbl_merge(
  tbls = list(tbl_unadj_arm, tbl_adj_arm),
  tab_spanner = c("**Unadjusted**", "**Adjusted**")
)

tbl_arm_comb <- tbl_arm_comb |> 
  as_gt() |> 
  tab_header(
    title = md("**Unadjusted and adjusted hazard ratio estimates by treatment arm**")
  )
tbl_arm_comb 

gtsave(
  tbl_arm_comb,
  filename = "figures/combined_adj_unadj_hr_arm.png"
)

#--------------------with parasitaemia
tbl_adj_para <- cox_multi |>
  tbl_regression(
    exp = TRUE,
    include = log_parasitaemia,
    label = list(
      log_parasitaemia ~ "Parasite density (/µL)"
    )
  ) |> 
  bold_p(t = 0.05) |>
  modify_header(estimate = "**Adjusted HR**")

tbl_unadj_para <- hr_cox_parasitaemia |> 
  modify_header(estimate = "**Unadjusted HR**")

# Merge the tables
tbl_para_comb <- tbl_merge(
  tbls = list(tbl_unadj_para, tbl_adj_para),
  tab_spanner = c("**Unadjusted**", "**Adjusted**")
)

tbl_para_comb <- tbl_para_comb |> 
  as_gt() |> 
  tab_header(
    title = md("**Unadjusted and adjusted hazard ratio estimates by parasite density**")
  )
tbl_para_comb 

#save the table
gtsave(
  tbl_para_comb,
  filename = "figures/combined_adj_unadj_hr_pasitaemia.png"
)

#------------with site
tbl_adj_site <- cox_multi |> 
  tbl_regression(
    exp = TRUE,
    include = site
  ) |> 
  bold_p(t = 0.05) |> 
  modify_header(estimate = "**Unadjusted HR**")

tbl_unadj_site <- hr_cox_site |> 
  modify_header(estimate = "**Unadjusted HR**")

# Merge the tables
tbl_site_comb <- tbl_merge(
  tbls = list(tbl_unadj_site, tbl_adj_site ),
  tab_spanner = c("**Unadjusted**", "**Adjusted**")
)

tbl_site_comb <- tbl_site_comb |> 
  as_gt() |> 
  tab_header(
    title = md("**Unadjusted and adjusted hazard ratio estimates by site**")
  )
tbl_site_comb

#save the table
gtsave(
  tbl_site_comb,
  filename = "figures/combined_adj_unadj_hr_site.png"
)

#================ TIME-VARYING COX MODEL =======================================

cox_multi_tv <- coxph(
  Surv(time_to_recurrence, recurrence_event) ~
    arm_AP +
    age +
    log_parasitaemia +
    baseline_haemoglobin +
    site +
    tt(arm_AP),
  
  data = recurrence_cox_parasitaemia |>
    mutate(
      arm_AP = if_else(treatment_arm == "AP", 1, 0)
    ),
  
  tt = function(x, t, ...) {
    x * log(t + 1) 
  }
)

summary(cox_multi_tv)


#---------------- Treatment HR over follow-up ----------------------------------

b <- coef(cox_multi_tv)
V <- vcov(cox_multi_tv)

hr_days <- tibble(
  day = c(4, 7, 14, 21, 28, 35, 42)
) |>
  mutate(
    time_term = log(day + 1),
    
    log_hr =
      b["arm_AP"] +
      b["tt(arm_AP)"] * time_term,
    
    se_log_hr = sqrt(
      V["arm_AP", "arm_AP"] +
        time_term^2 * V["tt(arm_AP)", "tt(arm_AP)"] +
        2 * time_term * V["arm_AP", "tt(arm_AP)"]
    ),
    
    HR = exp(log_hr),
    lower_95_CI = exp(log_hr - 1.96 * se_log_hr),
    upper_95_CI = exp(log_hr + 1.96 * se_log_hr)
  ) |>
  select(
    day,
    HR,
    lower_95_CI,
    upper_95_CI
  )

hr_days

#-----------l
#just a check with the tv model
cox_multi_tv_linear <- coxph(
  Surv(time_to_recurrence, recurrence_event) ~
    arm_AP +
    age +
    log_parasitaemia +
    baseline_haemoglobin +
    site +
    tt(arm_AP),
  
  data = recurrence_cox_parasitaemia |>
    mutate(
      arm_AP = if_else(treatment_arm == "AP", 1, 0)
    ),
  
  tt = function(x, t, ...) {
    x * t
  }
)
cox_multi_tv_linear


AIC(
  cox_multi_tv,
  cox_multi_tv_linear
)

#================ MODEL VALIDATION =============================================

#---------------- 1. Functional form of continuous variables ------------------

# Age
cox_age_spline <- update(
  cox_multi,
  . ~ . - age + pspline(age, df = 4)
)

anova(
  cox_multi,
  cox_age_spline,
  test = "LRT"
)

# Log10 baseline parasitaemia
cox_parasitaemia_spline <- update(
  cox_multi,
  . ~ . - log_parasitaemia + pspline(log_parasitaemia, df = 4)
)

anova(
  cox_multi,
  cox_parasitaemia_spline,
  test = "LRT"
)

# Baseline haemoglobin
cox_hb_spline <- update(
  cox_multi,
  . ~ . - baseline_haemoglobin +
    pspline(baseline_haemoglobin, df = 4)
)

anova(
  cox_multi,
  cox_hb_spline,
  test = "LRT"
)

#----------------Suggest Cox Model--------------------------------
# Age doesn't appear to be linear from the model validation

cox_multi_tv_final <- coxph(
  Surv(time_to_recurrence, recurrence_event) ~
    arm_AP +
    pspline(age, df = 4) +
    log_parasitaemia +
    baseline_haemoglobin +
    site +
    tt(arm_AP),
  
  data = recurrence_cox_parasitaemia |>
    mutate(
      arm_AP = if_else(treatment_arm == "AP", 1, 0)
    ),
  
  tt = function(x, t, ...) {
    x * log(t + 1)
  }
)
cox_multi_tv_final

# check AIC
AIC(
  cox_multi_tv,
  cox_multi_tv_final
)

#---------------- Treatment HR over follow-up from new tv model

b2 <- coef(cox_multi_tv)
V2 <- vcov(cox_multi_tv)

hr_days2 <- tibble(
  day = c(4, 7, 14, 21, 28, 35, 42)
) |>
  mutate(
    time_term = log(day + 1),
    
    log_hr =
      b2["arm_AP"] +
      b2["tt(arm_AP)"] * time_term,
    
    se_log_hr = sqrt(
      V2["arm_AP", "arm_AP"] +
        time_term^2 * V2["tt(arm_AP)", "tt(arm_AP)"] +
        2 * time_term * V2["arm_AP", "tt(arm_AP)"]
    ),
    
    HR = exp(log_hr),
    lower_95_CI = exp(log_hr - 1.96 * se_log_hr),
    upper_95_CI = exp(log_hr + 1.96 * se_log_hr)
  ) |>
  select(
    day,
    HR,
    lower_95_CI,
    upper_95_CI
  )

hr_days2

#=============== FOREST PLOT: SIMPLE MULTIVARIABLE COX MODEL ===================
cox_multi_plot_data <- broom::tidy(
  cox_multi,
  exponentiate = TRUE,
  conf.int = TRUE
) |> 
  mutate(
    variable = case_when(
      term == "treatment_armAP" ~ "Treatment: AP vs AL",
      term == "age" ~ "Age",
      term == "log_parasitaemia" ~ "Baseline parasitaemia(10-fold increase)",
      term == "siteAgogo" ~ "Site: Agogo vs Assin Fosu",
      TRUE ~ term
    )
  ) 

#plot the curve
cox_multi_plot <- ggplot(
  cox_multi_plot_data,
  aes(
    x = estimate,
    y = reorder(variable, estimates)
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed"
  ) +
  geom_errorbar(
    aes(
      xmin = conf.low,
      xmax = conf.high
    ),
    width = 0.15,
    oreintation = "y"
  )  +
  geom_point(
    size = 3
  ) +
  scale_x_log10() +
  labs(
    x = "Adjusted hazard ratio (95% CI)",
    y = NULL,
    title   = "Adjusted Hazard Ratios — Multivariable Cox Model",
    caption = "Horizontal bars = 95% CI · Dashed line = HR of 1 (no effect)"
  ) +
  theme_minimal(base_size = 12)

cox_multi_plot

# #---------------- 2. Influential observations ---------------------------------
# 
# dfbetas_cox <- residuals(
#   cox_multi,
#   type = "dfbetas"
# )
# 
# influence_check <- bind_cols(
#   recurrence_cox_parasitaemia |>
#     select(participant_id),
#   as.data.frame(dfbetas_cox)
# )
# 
# # Screening threshold
# dfbeta_threshold <- 2 / sqrt(cox_multi$n)
# 
# dfbeta_threshold
# 
# # check
# influential_participants <- influence_check |>
#   pivot_longer(
#     cols = -participant_id,
#     names_to = "coefficient",
#     values_to = "dfbeta"
#   ) |>
#   filter(
#     abs(dfbeta) > dfbeta_threshold
#   ) |>
#   arrange(
#     desc(abs(dfbeta))
#   )
# 
# influential_participants |>
#   print(n = Inf)
# 
# # check
# deviance_check <- recurrence_cox_parasitaemia |>
#   select(participant_id) |>
#   mutate(
#     deviance_residual = residuals(
#       cox_multi,
#       type = "deviance"
#     )
#   ) |>
#   arrange(
#     desc(abs(deviance_residual))
#   )
# 
# deviance_check |>
#   filter(
#     abs(deviance_residual) > 3
#   ) |>
#   print(n = Inf)

# names(baseline_data)

# baseline_data |>
#   filter(
#     #site == "Agogo",
#     baseline_parasitaemia > 0
#   ) |>
#   count(treatment_arm)
# 
# View(baseline_data)
# 
#  baseline_data |> 
#    count(site)
#  