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
    baseline_parasitaemia = log10(baseline_parasitaemia)
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
  recurrence_cox_data,
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
  aes(x = baseline_parasitaemia)
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

cox_age_valid <- coxph(
  Surv(time_to_recurrence, recurrence_event) ~ age,
  data = recurrence_cox_parasitaemia
)
summary(cox_age_valid)

#----------age_group
cox_age_group <- coxph(
  Surv(time_to_recurrence, recurrence_event) ~ age_group,
  data = recurrence_cox_data
)
summary(cox_age_group)

#--------------parasitaemia
cox_parasitaemia <- coxph(
  Surv(time_to_recurrence, recurrence_event) ~ baseline_parasitaemia,
  data = recurrence_cox_parasitaemia
)
summary(cox_parasitaemia)

#----------haemoglobin
cox_haemoglobin <- coxph(
  Surv(time_to_recurrence, recurrence_event) ~ baseline_haemoglobin,
  data = recurrence_cox_data
)
summary(cox_haemoglobin)

cox_haemoglobin_valid <- coxph(
  Surv(time_to_recurrence, recurrence_event) ~ baseline_haemoglobin,
  data = recurrence_cox_parasitaemia
)
summary(cox_haemoglobin_valid)

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

#---------------- Model 3: Valid baseline parasitaemia--Primary model--------------------------
# Uses participants with baseline_parasitaemia > 0

cox_multi <- coxph(
  Surv(time_to_recurrence, recurrence_event) ~
    treatment_arm +
    age +
    baseline_parasitaemia +
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

#--------------------------------PH for model 3- Primary model------------------
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

#---------------------Visualise PH assumption for model 3:primary model --------
ph_plot <- ggcoxzph(
  ph_cox_multi,
  point.alpha = 0.3,
  point.col = "violet",
  caption =  "Shaded band = 95% CI on Schoenfeld residuals · Flat line supports PH assumption"
)
ph_plot

#------save the ph plot
global_p <- ph_cox_multi$table["GLOBAL", "p"]

ph_plot_combined <- ggpubr::ggarrange(
  plotlist = ph_plot,
  ncol = 2,
  nrow = 3
)

ph_plot_final <- ggpubr::annotate_figure(
  ph_plot_combined,
  top = ggpubr::text_grob(
    paste0(
      "Global Schoenfeld Test p: ",
      round(global_p, 4)
    ),
    size = 12
  ),
  bottom = ggpubr::text_grob(
    "Dashed bands = approximate 95% CI · Flat line supports PH assumption",
    size = 10
  )
)

ph_plot_final

ggsave(
  filename = "figures/cox_ph_assumption_plot.png",
  plot = ph_plot_final,
  width = 10,
  height = 7,
  dpi = 300,
  bg = "white"
)

#--------------Test if a variable improves the model
anova(
  cox_restricted,
  cox_multi,
  test = "LRT"
)

#================ PH ASSUMPTION TABLE: Primary model ===========================
ph_table <- ph_cox_multi$table |>
  as.data.frame() |>
  tibble::rownames_to_column("Characteristic") |>
  rename(
    `Chi-square` = chisq,
    `df` = df,
    `p-value` = p
  ) |>
  mutate(
    Characteristic = recode(
      Characteristic,
      treatment_arm = "Treatment arm",
      age = "Age",
      log_parasitaemia = "Baseline parasitaemia",
      baseline_haemoglobin = "Baseline haemoglobin",
      site = "Study site",
      GLOBAL = "Global test"
    ),
    `PH assumption` = if_else(
      `p-value` < 0.05,
      "Violated",
      "Satisfied"
    )
  )

ph_table_gt <- ph_table |>
  gt() |>
  fmt_number(
    columns = `Chi-square`,
    decimals = 2
  ) |>
  fmt_number(
    columns = `p-value`,
    decimals = 3
  ) |>
  tab_header(
    title = md("**Proportional hazards assumption**")
  ) |> 
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  )

ph_table_gt

#save table
gt::gtsave(
  ph_table_gt,
  filename = "figures/cox_ph_assumption_table.png"
)

#===================== HAZARD RATIOS TABLES  ================================

#------------HR for model 3
hr_cox_multi <- cox_multi |> 
  tbl_regression(
    exp = TRUE,
    label = list(
      treatment_arm ~ "Treatment arm",
      age ~ "Age(months)",
      baseline_parasitaemia ~ "Baseline parasitaemia (10-fold increase)",
      baseline_haemoglobin ~ "Baseline haemoglobin (g/dL)"
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
      baseline_haemoglobin ~ "Baseline haemoglobin (g/dL)"
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
      baseline_haemoglobin ~ "Baseline haemoglobin (g/dL)"
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


#------------ Hazard ratios for univarable models-------------------------------
#-------treatment arm
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

hr_cox_arm_valid <- hr_cox_arm_valid |> 
  as_gt() |> 
  tab_header(
    title = md("**Undjusted hazard ratio estimate**")
  )
hr_cox_arm_valid

#----save the table
gtsave(
  hr_cox_arm_valid,
  filename = "figures/hr_cox_arm_valid.png"
)

#---------parasitaemia
hr_cox_parasitaemia <- cox_parasitaemia |> 
  tbl_regression(
    exp = TRUE,
    include = baseline_parasitaemia,
    label = list(
      baseline_parasitaemia ~ "Baseline parasitaemia (10-fold increase)"
    )
  ) |> 
  bold_p(t = 0.05)
hr_cox_parasitaemia

hr_cox_parasitaemia <- hr_cox_parasitaemia |> 
  as_gt() |> 
  tab_header(
    title = md("**Undjusted hazard ratio estimate**")
  )
hr_cox_arm_valid

#----save the table
gtsave(
  hr_cox_parasitaemia,
  filename = "figures/hr_cox_parasitaemia.png"
)

#----------site
hr_cox_site <- cox_site_valid |> 
  tbl_regression(
    exp = TRUE,
    include = site,
  ) |> 
  bold_p(t = 0.05)
hr_cox_site

hr_cox_site <- hr_cox_site |> 
  as_gt() |> 
  tab_header(
    title = md("**Undjusted hazard ratio estimate**")
  )
hr_cox_site

#----save the table
gtsave(
  hr_cox_site,
  filename = "figures/hr_cox_site.png"
)


#-----------age
hr_cox_age <- cox_age_valid|> 
  tbl_regression(
    exp = TRUE,
    include = age,
  ) |> 
  bold_p(t = 0.05)
hr_cox_age

hr_cox_age <- hr_cox_age |> 
  as_gt() |> 
  tab_header(
    title = md("**Undjusted hazard ratio estimate**")
  )
hr_cox_age

#----save the table
gtsave(
  hr_cox_age,
  filename = "figures/hr_cox_age.png"
)

#-----------haemoglobin
hr_cox_haemoglobin <- cox_haemoglobin_valid|> 
  tbl_regression(
    exp = TRUE,
    include = baseline_haemoglobin,
    label = list(
      baseline_haemoglobin ~ "Baseline haemoglobin (g/dL)"
    )
  ) |> 
  bold_p(t = 0.05)
hr_cox_haemoglobin

hr_cox_haemoglobin <- hr_cox_haemoglobin |> 
  as_gt() |> 
  tab_header(
    title = md("**Undjusted hazard ratio estimate**")
  )
hr_cox_haemoglobin

#----save the table
gtsave(
  hr_cox_haemoglobin,
  filename = "figures/hr_cox_haemoglobin.png"
)

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
    include = baseline_parasitaemia,
    label = list(
      baseline_parasitaemia ~ "Baseline parasitaemia (10-fold increase)"
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
# treatment arm violates the ph assumption and that we include time-varying model
cox_multi_tv <- coxph(
  Surv(time_to_recurrence, recurrence_event) ~
    arm_AP +
    age +
    baseline_parasitaemia +
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

#================ HAZARD RATIO TABLE: TIME-VARYING COX MODEL ==================

hr_cox_multi_tv <- cox_multi_tv |>
  tbl_regression(
    exp = TRUE,
    label = list(
      arm_AP ~ "Treatment: AP vs AL",
      age ~ "Age (months)",
      baseline_parasitaemia ~ "Baseline parasitaemia (10-fold increase)",
      baseline_haemoglobin ~ "Baseline haemoglobin (g/dL)",
      site ~ "Study site"
    )
  ) |>
  bold_p(t = 0.05)

hr_cox_multi_tv

hr_cox_multi_tv <- hr_cox_multi_tv |> 
  as_gt() |> 
  tab_header(
    title = md("**Adjusted hazard ratio estimates**")
  )
hr_cox_multi_tv 

#----save the table
gtsave(
  hr_cox_multi_tv,
  filename = "figures/hr_cox_multi_tv.png"
)

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
  . ~ . - baseline_parasitaemia + pspline(baseline_parasitaemia, df = 4)
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
# Age doesn't appear to be linear from the model validation and has been included 
# in the time-varying model

cox_multi_tv_final <- coxph(
  Surv(time_to_recurrence, recurrence_event) ~
    arm_AP +
    pspline(age, df = 4) +
    baseline_parasitaemia +
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
summary(cox_multi_tv_final)

# check AIC
AIC(
  cox_multi_tv,
  cox_multi_tv_final
)

#================ HAZARD RATIO TABLE: FINAL TIME-VARYING COX MODEL =============

hr_cox_multi_tv_final <- cox_multi_tv_final |>
  tbl_regression(
    exp = TRUE,
    label = list(
      arm_AP ~ "Treatment: AP vs AL",
      baseline_parasitaemia ~ "Baseline parasitaemia (10-fold increase)",
      baseline_haemoglobin ~ "Baseline haemoglobin (g/dL)",
      site ~ "Study site"
    )
  ) |>
  modify_table_body(
    ~ .x |>
      add_row(
        label = "Age spline (linear)",
        estimate = NA_real_,
        conf.low = NA_real_,
        conf.high = NA_real_,
        p.value = 0.95943,
        .after = 1
      ) |>
      add_row(
        label = "Age spline (non-linear)",
        estimate = NA_real_,
        conf.low = NA_real_,
        conf.high = NA_real_,
        p.value = 0.00344,
        .after = 2
      ) |>
      mutate(
        label = if_else(
          label == "tt(arm_AP)",
          "Treatment × log(time + 1)",
          label
        )
      )
  ) |>
  bold_p(t = 0.05)

hr_cox_multi_tv_final

hr_cox_multi_tv_final <- hr_cox_multi_tv_final |> 
  as_gt() |> 
  tab_header(
    title = md("**Adjusted hazard ratio estimates**")
  )
hr_cox_multi_tv_final 

#----save the table
gtsave(
  hr_cox_multi_tv_final,
  filename = "figures/hr_cox_multi_tv_final.png"
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
      term == "age" ~ "Age (months)",
      term == "baseline_parasitaemia" ~ "Baseline parasitaemia (10-fold increase)",
      term == "baseline_haemoglobin" ~ "Haemoglobin(g/dL)",
      term == "siteAgogo" ~ "Site: Agogo vs Assin Fosu",
      TRUE ~ term
    )
  ) 

#plot the curve
cox_multi_plot <- ggplot(
  cox_multi_plot_data,
  aes(
    x = estimate,
    y = reorder(variable, estimate)
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed"
  ) +
  geom_errorbar(
    aes(
      xmin = conf.low,
      xmax = conf.high,
      colour = p.value < 0.05
    ),
    width = 0.15,
    orientation = "y"
  ) +
  geom_point(
    aes(colour = p.value < 0.05),
    size = 3
  ) +
  scale_color_manual(
    name = "Statistical significance",
    values = c("TRUE" = "#008080", "FALSE" = "gray60"),
    labels = c("TRUE" = "p < 0.05", "FALSE" = "p ≥ 0.05"),
    breaks = c("TRUE", "FALSE")
  ) +
  scale_x_log10(
    breaks = c(0.5, 0.75, 1, 1.25, 1.5, 2, 2.5, 3),
    labels = c("0.50", "0.75", "1.00", "1.25", "1.50", "2.00", "2.50", "3.00")
  ) +
  labs(
    x = "Adjusted hazard ratio (95% CI)",
    y = NULL,
    title   = "Multivariable Cox proportional hazards model",
    caption = "Horizontal bars = 95% CI · Dashed line = HR of 1 (no effect)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top") 
  

cox_multi_plot

ggsave(
  filename = "figures/cox_multi_forest_ggplot.png",
  plot = cox_multi_plot,
  width = 9,
  height = 5,
  dpi = 300
)

#-------------- Using ggforest for the plot------------------------------------

cox_plot <- ggforest(
  model = cox_multi,
  data = model.frame(cox_multi),
  main = "Adjusted hazard ratios",
  fontsize = 0.85
)

cox_plot

# save 
ggsave(
  filename = "figures/cox_multi_ggforest.png",
  plot = cox_plot,
  width = 9,
  height = 5,
  dpi = 300
)


#================ FOREST PLOT: TIME-VARYING COX MODEL ========================

# extract model results
cox_tv_plot_data <- broom::tidy(
  cox_multi_tv,
  exponentiate = TRUE,
  conf.int = TRUE
) |>
  mutate(
    variable = case_when(
      term == "arm_AP" ~ "Treatment: AP vs AL",
      term == "age" ~ "Age (months)",
      term == "baseline_parasitaemia" ~ "Baseline parasitaemia (10-fold increase)",
      term == "baseline_haemoglobin" ~ "Baseline haemoglobin (g/dL)",
      term == "siteAgogo" ~ "Site: Agogo vs Assin Fosu",
      term == "tt(arm_AP)" ~ "Time-varying treatment effect",
      TRUE ~ term
    )
  )

# plot
cox_multi_tv_plot <- ggplot(
  cox_tv_plot_data,
  aes(
    x = estimate,
    y = reorder(variable, estimate)
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed"
  ) +
  geom_errorbar(
    aes(
      xmin = conf.low,
      xmax = conf.high,
      colour = p.value < 0.05
    ),
    width = 0.15,
    orientation = "y"
  ) +
  geom_point(
    aes(colour = p.value < 0.05),
    size = 3
  ) +
  scale_color_manual(
    name = "Statistical significance",
    values = c("TRUE" = "#008080", "FALSE" = "gray60"),
    labels = c("TRUE" = "p < 0.05", "FALSE" = "p ≥ 0.05"),
    breaks = c("TRUE", "FALSE")
  ) +
  scale_x_log10() +
  labs(
    x = "Adjusted hazard ratio (95% CI)",
    y = NULL,
    title = "Time-varying Cox proportional hazards model",
    caption = "Horizontal bars = 95% CI · Dashed line = HR of 1 (no effect)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

cox_multi_tv_plot 

# save 
ggsave(
  filename = "figures/cox_multi_tv_forest_ggplot.png",
  plot = cox_multi_tv_plot,
  width = 9,
  height = 5,
  dpi = 300
)

#-------------- Using ggforest for the plot------------------------------------

cox_plot_tv <- ggforest(
  model = cox_multi_tv,
  data = model.frame(cox_multi_tv),
  main = "Adjusted hazard ratios",
  fontsize = 0.85
)
cox_plot_tv

# save 
ggsave(
  filename = "figures/cox_multi_tv_ggforest.png",
  plot = cox_plot_tv,
  width = 9,
  height = 5,
  dpi = 300
)

#================ FOREST PLOT: FINAL TIME-VARYING COX MODEL ==================

tv_final_sum <- summary(cox_multi_tv_final)

cox_tv_final_plot_data <- tibble(
  term = rownames(tv_final_sum$coefficients),
  estimate = exp(tv_final_sum$coefficients[, "coef"]),
  conf.low = exp(tv_final_sum$coefficients[, "coef"] - 1.96 * tv_final_sum$coefficients[, "se(coef)"]),
  conf.high = exp(tv_final_sum$coefficients[, "coef"] + 1.96 * tv_final_sum$coefficients[, "se(coef)"]),
  p.value = tv_final_sum$coefficients[, "p"]
) |>
  mutate(
    variable = case_when(
      term == "arm_AP" ~ "Treatment: AP vs AL",
      term == "pspline(age, df = 4), linear" ~ "Age spline (linear)",
      term == "pspline(age, df = 4), nonlin" ~ "Age spline (non-linear)",
      term == "baseline_parasitaemia" ~ "Baseline parasitaemia (10-fold increase)",
      term == "baseline_haemoglobin" ~ "Baseline haemoglobin (g/dL)",
      term == "siteAgogo" ~ "Site: Agogo vs Assin Fosu",
      term == "tt(arm_AP)" ~ "Time-varying treatment effect",
      TRUE ~ term
    )
  )

cox_multi_tv_final_plot <- ggplot(
  cox_tv_final_plot_data,
  aes(
    x = estimate,
    y = reorder(variable, estimate)
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed"
  ) +
  geom_errorbar(
    aes(
      xmin = conf.low,
      xmax = conf.high,
      colour = p.value < 0.05
    ),
    width = 0.15,
    orientation = "y"
  ) +
  geom_point(
    aes(colour = p.value < 0.05),
    size = 3
  ) +
  scale_color_manual(
    name = "Statistical significance",
    values = c("TRUE" = "#008080", "FALSE" = "gray60"),
    labels = c("TRUE" = "p < 0.05", "FALSE" = "p ≥ 0.05"),
    breaks = c("TRUE", "FALSE")
  ) +
  scale_x_log10() +
  labs(
    x = "Adjusted hazard ratio (95% CI)",
    y = NULL,
    title = "Final time-varying Cox proportional hazards model",
    caption = "Horizontal bars = 95% CI · Dashed line = HR of 1 (no effect)"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

cox_multi_tv_final_plot

# save 
ggsave(
  filename = "figures/cox_multi_tv_final_forest_ggplot.png",
  plot = cox_multi_tv_final_plot,
  width = 9,
  height = 5,
  dpi = 300
)

#-------------- Using ggforest for the plot------------------------------------

cox_plot_tv_final <- ggforest(
  model = cox_multi_tv_final,
  data = model.frame(cox_multi_tv_final),
  main = "Adjusted hazard ratios",
  fontsize = 0.85
)
cox_plot_tv_final 

#save
ggsave(
  filename = "figures/cox_multi_tv_final_ggforest.png",
  plot = cox_plot_tv_final,
  width = 9,
  height = 5,
  dpi = 300
)

#=================== FIT COX CURVE =============================================

cox_adjusted_curves <- ggadjustedcurves(
  fit = cox_multi,
  data = as.data.frame(recurrence_cox_parasitaemia),
  variable = "treatment_arm",
  method = "conditional"
) +
  scale_x_continuous(
    breaks = c(0, 1, 2, 3, 4, 7, 14, 21, 28, 35, 42),
    limits = c(0, 42)
  )+
  labs(
    x = "Follow-up day",
    y = "Adjusted probability of no recurrence",
    title = "Adjusted survival curves by treatment arm",
    #subtitle = "Adjusted for age, baseline parasitaemia, haemoglobin and site"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

cox_adjusted_curves

#save
ggsave(
  filename = "figures/cox_adjusted_curves.png",
  plot = cox_adjusted_curves,
  width = 10,
  height = 7,
  dpi = 300
)
