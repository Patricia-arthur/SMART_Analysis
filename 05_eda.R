#============================EXPLORATORY DATA ANALYSIS==========================
source("04_analysis_data.R")
library(gtsummary)
library(gt)
library(webshot2)

#==============================BASELINE DATASET=================================

baseline_data <- smart_analysis_final |>
  ungroup()|> 
  filter(study_visit == "H0") |>
  select(
    participant_id,
    screening_id,
    treatment_arm,
    site,
    age,
    sex,
    weight,
    baseline_temperature = temperature,
    baseline_haemoglobin = haemoglobin,
    baseline_parasitaemia = parasite_density
  ) |>
  mutate(
    age_group = case_when(
      age < 60  ~ "<5 years",
      age >= 60 ~ ">=5 years",
      TRUE ~ NA_character_
    )
  )

#check there is one row per participant
baseline_data |>
  count(participant_id) |>
  filter(n > 1)

#check number of participants at baseline
baseline_data |>
  summarise(
    n_participants = n(),
    n_AL = sum(treatment_arm == "AL", na.rm = TRUE),
    n_AP = sum(treatment_arm == "AP", na.rm = TRUE)
  )

#check the range of age,weight, temp, hb
baseline_data |>
  summarise(
    min_age_months = min(age, na.rm = TRUE),
    max_age_months = max(age, na.rm = TRUE),
    median_age_months = median(age, na.rm = TRUE),
    min_weight = min(weight, na.rm = TRUE),
    max_weight = max(weight, na.rm = TRUE),
    median_weight = median(weight, na.rm = TRUE)
  )

baseline_data |>
  summarise(
    min_temperature = min(baseline_temperature, na.rm = TRUE),
    max_temperature = max(baseline_temperature, na.rm = TRUE),
    median_temperature = median(baseline_temperature, na.rm = TRUE),
    min_haemoglobin = min(baseline_haemoglobin, na.rm = TRUE),
    max_haemoglobin = max(baseline_haemoglobin, na.rm = TRUE),
    median_haemoglobin = median(baseline_haemoglobin, na.rm = TRUE)
  )

#check age group
baseline_data |>
  count(
    treatment_arm,
    age_group
  )

# Check missingness at baseline
baseline_data |>
  summarise(
    missing_age = sum(is.na(age)),
    missing_sex = sum(is.na(sex)),
    missing_weight = sum(is.na(weight)),
    missing_baseline_temperature = sum(is.na(baseline_temperature)),
    missing_baseline_haemoglobin = sum(is.na(baseline_haemoglobin)),
    missing_baseline_parasitaemia = sum(is.na(baseline_parasitaemia))
  )

#parasite density at baseline
SMART_analysis |>
  ungroup() |>
  filter(study_visit == "H0") |>
  summarise(
    n_participants = n_distinct(participant_id),
    n_with_parasite_density = sum(!is.na(parasite_density)),
    n_without_parasite_density = sum(is.na(parasite_density))
  )

baseline_data |>
  filter(age < 6) |>
  select(
    participant_id,
    treatment_arm,
    site,
    age
  )

#--------------------------------Check variable coding--------------------------
baseline_data |>
  count(sex)

baseline_data |>
  count(site)

baseline_data |>
  count(treatment_arm)

#--------------------------------------Label variables--------------------------
baseline_data <- baseline_data |> 
  mutate(
    sex = factor(
      sex,
      levels = c(0,1),
      labels = c("Female", "Male")
    ),
    
    site = factor(
      site,
      levels = c(1,2),
      labels = c("Assin Fosu", "Agogo")
    ),
    
    treatment_arm = factor(
      treatment_arm,
      levels = c("AL", "AP"),
      labels = c("Artemether-Lumefantrine", "Artesunate-Pyranoridine")
    ),
    
    age_group = factor(
      age_group,
      levels = c ("<5 years", ">=5 years")
    )
  )

# Count participants according to site and sex
baseline_data |>
  count(sex)

baseline_data |>
  count(site)

#===================BASELINE CHARACTERISTICS TABLES=============================
#=================1. BASELINE CHARACTERISTICS BY ARM=============================
table1_arm <- baseline_data |> 
  select(
    treatment_arm,
    age,
    age_group,
    sex,
    weight,
    baseline_temperature,
    baseline_haemoglobin,
    #baseline_parasitaemia,
  ) |> 
  tbl_summary(
    by = treatment_arm,
    
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    
    missing = "ifany",
    
    label = list(
      age ~ "Age (months)",
      age_group ~ "Age group",
      sex ~ "Sex",
      weight ~ "Weight (kg)",
      baseline_temperature ~ "Temperature (°C)",
      baseline_haemoglobin ~ "Haemoglobin (g/dL)"
      #baseline_parasitaemia ~ "Parasite density (/µL)"
    )
  ) |> 
  #add_n() |>
  add_overall() |> 
  bold_labels()

table_arm_gt <- table1_arm |> 
  as_gt() |> 
  gt:: tab_header(
    title = gt::md("**Baseline characteristics by treatment arm**")
  )
table_arm_gt

#--- Save the table---
gt:: gtsave(
  table_arm_gt,
  "figures/baseline_charac_arm.png"
)


#================BASELINE PARASITAEMIA ISSUE====================================
#---check participants with baseline parasitaemia = o
baseline_data |>
  filter(!is.na(baseline_parasitaemia)) |>
  count(
    baseline_parasitaemia == 0
  )

baseline_data |>
  summarise(
    n_available = sum(!is.na(baseline_parasitaemia)),
    n_zero = sum(baseline_parasitaemia == 0, na.rm = TRUE),
    n_positive = sum(baseline_parasitaemia > 0, na.rm = TRUE)
  )

baseline_zero_parasitaemia <- baseline_data |>
  filter(
    baseline_parasitaemia == 0
  ) |>
  select(
    participant_id,
    baseline_parasitaemia
  )
write_csv(baseline_zero_parasitaemia, "Discrepancy data/H0_parasitaemia.csv")

baseline_missing_parasitaemia <- baseline_data |>
  filter(
    is.na(baseline_parasitaemia)
  ) |>
  select(
    participant_id,
    baseline_parasitaemia
  )
write_csv(baseline_missing_parasitaemia, "Discrepancy data/missing_baseline_parasitaemia.csv")

#=================2. BASELINE CHARACTERISTICS BY SITE===========================
table1_site <- baseline_data |> 
  select(
    site,
    age,
    age_group,
    sex,
    weight,
    baseline_temperature,
    baseline_haemoglobin,
    #baseline_parasitaemia
  ) |> 
  tbl_summary(
    by = site,
    
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    
    missing = "ifany",
    
    label = list(
      age ~ "Age (months)",
      age_group ~ "Age_group",
      sex ~ "Sex",
      weight ~ "Weight (kg)",
      baseline_temperature ~ "Temperature (°C)",
      baseline_haemoglobin ~ "Haemoglobin (g/dL)"
      #baseline_parasitaemia ~ "Parasite density (/µL)"
    )
  ) |> 
  #add_n() |> 
  add_overall() |> 
  bold_labels()

table_site_gt <- table1_site |> 
  as_gt() |> 
  gt:: tab_header(
    title = gt::md("**Baseline characteristics by study site**")
  )
table_site_gt

#--- Save the table---
gt:: gtsave(
  table_site_gt,
  "figures/baseline_charac_site.png"
)

# Treatment allocation by site
baseline_data |>
  count(
    site,
    treatment_arm
  )

#================3. BASELINE CHARACTERISTICS BY ARM UNDER SITE==================
table_site_arm <- baseline_data |> 
  select(
    site,
    treatment_arm,
    age,
    age_group,
    sex,
    weight,
    baseline_temperature,
    baseline_haemoglobin
    #baseline_parasitaemia
  ) |> 
  tbl_strata(
    strata = site,
    .tbl_fun = ~ .x |> 
      tbl_summary(
        by = treatment_arm,
        
        # type = list(
        #   baseline_parasitaemia ~ "continuous"
        # ),
        
        statistic = list(
          age ~ "{median} ({p25}, {p75})",
          weight ~ "{mean} ({sd})",
          baseline_temperature ~ "{mean} ({sd})",
          baseline_haemoglobin ~ "{mean} ({sd})",
          #baseline_parasitaemia ~ "{median} ({p25}, {p75})",
          all_categorical() ~ "{n} ({p}%)"
        ),
        
        digits = list(
          age ~ 0,
          weight ~ 1,
          baseline_temperature ~ 1,
          baseline_haemoglobin ~ 1,
          #baseline_parasitaemia ~ 0,
          all_categorical() ~ c(0,1)
        ),
        
        missing = "ifany",
        
        label = list(
          age ~ "Age (months)",
          age_group ~ "Age_group",
          sex ~ "Sex",
          weight ~ "Weight (kg)",
          baseline_temperature ~ "Temperature (°C)",
          baseline_haemoglobin ~ "Haemoglobin (g/dL)"
          #baseline_parasitaemia ~ "Parasite density (/µL)"
        )
      ) |> 
      bold_labels(),
    .combine_with = "tbl_merge"
  )

table_site_arm_gt <- table_site_arm |> 
  as_gt() |> 
  gt:: tab_header(
    title = gt::md("**Baseline characteristics by study site and treatment arm**")
  )
table_site_arm_gt

#--- Save the table---
gt:: gtsave(
  table_site_arm_gt,
  "figures/baseline_charac_site_arm.png"
)

#===================BASELINE DESCRIPTIVE CURVES=================================
#--------------descriptive table for temperature------------
temperature_summary <- smart_analysis_final |> 
  ungroup() |> 
  group_by(
    treatment_arm,
    study_day,
    study_visit
  ) |> 
  summarise(
    n_rows = n(),
    
    n_temperature = sum(!is.na(temperature)),
    
    n_missing_temperature = sum(is.na(temperature)),
    
    mean_temperature = if(
      all(is.na(temperature))
    ) NA_real_ else mean(temperature, na.rm = TRUE),
    
    sd_temperature = if(
      all(is.na(temperature))
    ) NA_real_ else sd(temperature, na.rm = TRUE),
    
    median_temperature = if(
      all(is.na(temperature))
    ) NA_real_ else median(temperature, na.rm = TRUE),
    
    q1_temperature = if(
      all(is.na(temperature))
    ) NA_real_ else quantile(temperature, 0.25, na.rm = TRUE),
    
    q3_temperature = if(
      all(is.na(temperature))
    ) NA_real_ else quantile(temperature, 0.75, na.rm = TRUE),
    
    .groups = "drop"
  )

#check the table
temperature_summary |> 
  print(n = Inf)

#===================TEMPERATURE CURVE AT SCHEDULED VISITS=======================
scheduled_visits <- c(
  "H0", "H8", "H24", "H36", "H48", "H60","H72","D7","D14","D21","D28","D35","D42"
)

temperature_curve <- smart_analysis_final |> 
  ungroup() |> 
  filter(
    study_visit %in% scheduled_visits
  ) |> 
  mutate(
    study_visit = factor(
      study_visit,
      levels = scheduled_visits
    )
  ) |> 
  group_by(
    treatment_arm,
    study_visit
  ) |> 
  summarise(
    n = sum(!is.na(temperature)),
    mean_temperature = mean(
      temperature, na.rm = TRUE
    ),
    sd_temperature = sd(temperature, na.rm = TRUE),
    .groups = "drop"
  )

temperature_curve |> 
  print(n = Inf)

#-------------Draw the descriptive curve----------------------------------------
temperature_plot <- ggplot(
  temperature_curve,
  aes(
    x = study_visit,
    y = mean_temperature,
    group = treatment_arm,
    color = treatment_arm
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  labs(
    x = "Study visit",
    y = "Mean temperature (°C)",
    color = "Treatment arm",
    title = "Temperature curves by treatment arm"
  ) +
  theme_minimal()

temperature_plot

#---Save the plot
ggsave(
  "figures/temperature_curve.png",
  plot = temperature_plot,
  width = 9,
  height = 6,
  dpi = 300
)

#--------Interpretation of the temperature curve--------------------------------
# Mean temperature declined rapidly in both treatment arms during the early 
# follow-up period. At H0, participants in both the artemether-lumefantrine (AL) 
# and artesunate-pyronaridine (AP) arms had the highest mean temperatures, 
# approximately 37.3°C. Mean temperature decreased substantially by H8 and reached 
# approximately 36.5°C by H24. Thereafter, temperatures remained relatively stable 
# within the normal range throughout follow-up to Day 42. The trajectories for AL 
# and AP were very similar during the first 72 hours, with only small differences
# observed at later follow-up visits. Overall, the descriptive pattern suggests 
# rapid resolution of elevated temperature following treatment in both groups, 
# with no marked difference in temperature trajectory between the two treatment arms.


#--------Descriptive table for haemoglobin--------------------------------------
haemoglobin_summary <- smart_analysis_final |>
  ungroup() |>
  filter(study_visit %in% scheduled_visits) |>
  group_by(
    treatment_arm,
    study_day,
    study_visit
  ) |>
  summarise(
    n_rows = n(),
    
    n_haemoglobin = sum(!is.na(haemoglobin)),
    
    n_missing_haemoglobin = sum(is.na(haemoglobin)),
    
    mean_haemoglobin = if (
      all(is.na(haemoglobin))
    ) NA_real_ else mean(haemoglobin, na.rm = TRUE),
    
    sd_haemoglobin = if (
      sum(!is.na(haemoglobin)) <= 1
    ) NA_real_ else sd(haemoglobin, na.rm = TRUE),
    
    median_haemoglobin = if (
      all(is.na(haemoglobin))
    ) NA_real_ else median(haemoglobin, na.rm = TRUE),
    
    q1_haemoglobin = if (
      all(is.na(haemoglobin))
    ) NA_real_ else quantile(haemoglobin, 0.25, na.rm = TRUE),
    
    q3_haemoglobin = if (
      all(is.na(haemoglobin))
    ) NA_real_ else quantile(haemoglobin, 0.75, na.rm = TRUE),
    
    .groups = "drop"
  )

haemoglobin_summary |> 
  print(n = Inf)

#check missing haemoglobin
smart_analysis_final |> 
  filter( is.na(haemoglobin)) |> 
  count(
    n_row = n(),
    study_visit
  )

#===================HAEMOGLOBIN CURVE AT SPECIFIC VISITS========================
hb_visits  <- c(
  "H0","H72","D7","D14","D28","D42"
)

haemoglobin_curve <- smart_analysis_final |> 
  ungroup() |> 
  filter(
    study_visit %in% hb_visits
  ) |> 
  mutate(
    study_visit = factor(
      study_visit,
      levels = hb_visits
    )
  ) |> 
  group_by(
    treatment_arm,
    study_visit
  ) |> 
    summarise(
      n = sum(!is.na(haemoglobin)),
      
      mean_haemogloblin = mean(haemoglobin, na.rm = TRUE),
      
      sd_haemoglobin = sd(haemoglobin, na.rm = TRUE),
      
      .groups = "drop"
    )

haemoglobin_curve |> 
  print(n = Inf)

#-------------------------------plot the curve----------------------------------
haemoglobin_plot <- ggplot(
  haemoglobin_curve,
  aes(
    x = study_visit,
    y = mean_haemogloblin,
    group = treatment_arm,
    color = treatment_arm
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.5) +
  labs(
    x = "Study visits",
    y = " Mean haemoglobin (g/dL)",
    color = " Treatment arm",
    title =  "Haemoglobin curves by treatment arm"
  ) +
  theme_minimal()
haemoglobin_plot

ggsave(
  "figures/haemoglobin_curve.png",
  plot = haemoglobin_plot,
  width = 9,
  height = 6,
  dpi = 300
)

#-----------Interpretation of the HAEMOGLOBIN curve-----------------------------
# Mean haemoglobin concentration followed a similar pattern in both treatment arms
# over the follow-up period. In both the artemether-lumefantrine (AL) and 
# artesunate-pyronaridine (AP) groups, mean haemoglobin decreased slightly from
# baseline (H0) to H72. This was followed by a gradual increase from Day 7 through
# Day 28, with mean haemoglobin concentrations rising above baseline levels in both
# treatment arms. By Day 42, haemoglobin remained higher than at baseline, although
# a small decline from the Day 28 level was observed in both groups. The AL and AP
# trajectories were broadly similar throughout follow-up, with AP showing slightly 
# higher mean haemoglobin concentrations at H0, Day 28, and Day 42. The 
# descriptive pattern indicates an initial decline in haemoglobin during the early 
# post-treatment period followed by progressive recovery over subsequent follow-up.

#--------Descriptive table for biochemistry-------------------------------------
# check visits where data exists

biochemistry_availability <- smart_analysis_final |>
  ungroup() |>
  filter(study_visit %in% scheduled_visits) |>
  select(
    treatment_arm,
    study_day,
    study_visit,
    total_bilirubin,
    direct_bilirubin,
    albumin,
    alt,
    ast,
    alp,
    ggt,
    total_protein
  ) |>
  pivot_longer(
    cols = c(
      total_bilirubin,
      direct_bilirubin,
      albumin,
      alt,
      ast,
      alp,
      ggt,
      total_protein
    ),
    names_to = "biochemistry_marker",
    values_to = "value"
  ) |>
  group_by(
    treatment_arm,
    study_day,
    study_visit,
    biochemistry_marker
  ) |>
  summarise(
    n_rows = n(),
    n_available = sum(!is.na(value)),
    n_missing = sum(is.na(value)),
    .groups = "drop"
  )
  
biochemistry_availability |>
  filter(n_available > 0) |>
  print(n = Inf) |> 
  distinct(study_visit)

#===================HAEMOGLOBIN CURVE AT SPECIFIC VISITS========================
biochem_visits <- c(
  "H0", "H72", "D7", "D14", "D28", "D42"
)

# Put all biochemical markers into one long dataset
biochemistry_long <- smart_analysis_final |>
  ungroup() |>
  filter(
    study_visit %in% biochem_visits
  ) |>
  select(
    treatment_arm,
    study_visit,
    total_bilirubin,
    direct_bilirubin,
    albumin,
    alt,
    ast,
    alp,
    ggt,
    total_protein
  ) |>
  pivot_longer(
    cols = c(
      total_bilirubin,
      direct_bilirubin,
      albumin,
      alt,
      ast,
      alp,
      ggt,
      total_protein
    ),
    names_to = "biochemistry_marker",
    values_to = "value"
  ) |>
  mutate(
    study_visit = factor(
      study_visit,
      levels = biochem_visits
    )
  )

#Summarise each marker by treatment arm and visit
biochemistry_summary <- biochemistry_long |>
  group_by(
    treatment_arm,
    study_visit,
    biochemistry_marker
  ) |>
  summarise(
    n = sum(!is.na(value)),
    
    median_value = if (
      all(is.na(value))
    ) NA_real_ else median(value, na.rm = TRUE),
    
    q1 = if (
      all(is.na(value))
    ) NA_real_ else quantile(value, 0.25, na.rm = TRUE),
    
    q3 = if (
      all(is.na(value))
    ) NA_real_ else quantile(value, 0.75, na.rm = TRUE),
    
    .groups = "drop"
  )

biochemistry_summary |>
  print(n = Inf)

#Give the markers readable names
biochemistry_summary <- biochemistry_summary |>
  mutate(
    biochemistry_marker = recode(
      biochemistry_marker,
      total_bilirubin = "Total bilirubin",
      direct_bilirubin = "Direct bilirubin",
      albumin = "Albumin",
      alt = "ALT",
      ast = "AST",
      alp = "ALP",
      ggt = "GGT",
      total_protein = "Total protein"
    )
  )

#-------------------plot the curve----------------------------------------------
biochemistry_plot <- ggplot(
  biochemistry_summary,
  aes(
    x = study_visit,
    y = median_value,
    group = treatment_arm,
    color = treatment_arm
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(
    ~ biochemistry_marker,
    scales = "free_y"
  ) +
  labs(
    x = "Study visit",
    y = "Median value",
    color = "Treatment arm",
    title = "Biochemistry trajectories by treatment arm"
  ) +
  theme_minimal()
biochemistry_plot

#---save the plot
ggsave(
  "figures/biochemistry_curve.png",
  plot = biochemistry_plot,
  width = 12,
  height = 9,
  dpi = 300
)

#---------------------Interpretation of the BIOCHEMISTRY curves-----------------
# Biochemical parameters showed broadly similar temporal patterns between the 
# artemether-lumefantrine (AL) and artesunate-pyronaridine (AP) treatment arms, 
# although some differences were observed for individual markers. Median albumin
# initially declined slightly by H72 before increasing during follow-up in both 
# treatment groups. Alkaline phosphatase (ALP) generally increased over time, with
# the largest values observed at Day 42, particularly in the AP arm.
# 
# Alanine aminotransferase (ALT) and aspartate aminotransferase (AST) showed some 
# fluctuations over follow-up rather than a consistent increasing or decreasing pattern.
# ALT values were generally higher in the AP arm at several visits, while AST 
# trajectories were broadly similar between the two groups.
# 
# Both total and direct bilirubin declined markedly from baseline during the early
# follow-up period, reaching their lowest levels around Days 7–14 before showing a 
# modest increase towards Day 42. This declining pattern was observed in both treatment
# arms. Gamma-glutamyl transferase (GGT) remained relatively stable over follow-up, 
# although median values were consistently somewhat higher in the AP arm than in the 
# AL arm.
# 
# Total protein showed relatively small fluctuations from H0 through Day 28, 
# followed by a noticeable decline at Day 42 in both treatment groups. The 
# biochemical trajectories were generally comparable between AL and AP, with no 
# obvious sustained divergence between the treatment arms based on the descriptive
# curves alone. These findings should be interpreted cautiously because biochemical
# measurements were not available for all participants at every scheduled visit.