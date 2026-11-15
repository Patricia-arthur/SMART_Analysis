source("09_parasite_clearance.R")

#========STUDY TOWNS WITH BASELINE AVAILABLE PARASITE DENSITY===================
study_towns <- SMARTUncomplicated |> 
  filter(!is.na(h01_3um)) |> 
  select(
    h01_3um,
    b06um
  ) |> 
  distinct() |> 
  rename(
    participant_id = h01_3um,
    study_villages = b06um
  ) |> 
  mutate(
    participant_id = as.character(participant_id)
  )

available_parasite_count <- baseline_data |> 
  filter(baseline_parasitaemia > 0 ) |> 
  select(
    participant_id,
    baseline_parasitaemia
  ) |> 
  rename(
    parasite_density_H0 = baseline_parasitaemia
  )

study_towns <- study_towns |> 
  left_join(
    available_parasite_count,
    by = "participant_id"
  ) 

study_towns <- study_towns |> 
  filter(!is.na(parasite_density_H0))

View(study_towns)
write_csv(study_towns, "Discrepancy data/study_towns.csv" )
