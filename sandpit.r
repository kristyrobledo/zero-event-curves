## check lacc data
source("utils.R")
library(readxl)
library(tidyverse)

shape <- readxl::read_excel("C:/Users/kman2804/OneDrive - The University of Sydney (Staff)/1. RESEARCH PROJECTS/3. GYNO ONC/CRE ALIVEWELL/LACC_shape/Shape/SHAPE_balance_table_KM_coordinates.xlsx",
  sheet = "SHAPEM_by_surgery group"
)

# remove row at time 0
shape <- shape[shape$Time_months > 0, ]

# run for shape by trt
# Curve for each stratum (adds a `group` column automatically)
shape_by_stratum <- event_free_curve(shape$Time_years, group = shape$Stratum)

# Bootstrap point estimate + 95% CI at t = 3 years, one row per stratum
shape_est_3yr <- bootstrap_event_free_band(
  shape$Time_years,
  group      = shape$Stratum,
  times_grid = 3,
  B          = 2000,
  seed       = 1
)

shape_est_3yr

## run for LACC
lacc <- haven::read_sas("C:/Users/kman2804/OneDrive - The University of Sydney (Staff)/1. RESEARCH PROJECTS/3. GYNO ONC/CRE ALIVEWELL/LACC_shape/LACC/patients.sas7bdat")
head(lacc)

lacc$dfsdays
lacc$AlloTreat
lacc_by_stratum <- event_free_curve(lacc$dfsdays, group = lacc$AlloTreat)

# Bootstrap point estimate + 95% CI at t = 3 years, one row per stratum
lacc_est_3yr <- bootstrap_event_free_band(
  lacc$dfsdays,
  group      = lacc$AlloTreat,
  times_grid = 3 * 365.25,
  B          = 2000,
  seed       = 1
)

lacc_est_3yr

## run for combined
# combine datasets
shape |>
  mutate(
    dataset = "shape",
    trt = Stratum,
    dfsyrs = Time_years
  ) |>
  select(dataset, trt, dfsyrs) -> shape2
lacc |>
  mutate(
    dataset = "lacc",
    trt = AlloTreat,
    dfsyrs = dfsdays / 365.25
  ) |>
  select(dataset, trt, dfsyrs) -> lacc2

# combine
all <- rbind(shape2, lacc2)
head(all)

all_by_trt <- event_free_curve(all$dfsyrs, group = all$trt)

# Bootstrap point estimate + 95% CI at t = 3 years, one row per stratum
all_est_3yr <- bootstrap_event_free_band(
  all$dfsyrs,
  group      = all$trt,
  times_grid = 3,
  B          = 2000,
  seed       = 1
)

all_est_3yr
table(all$trt)
## oops align trt coding, one has MIS/open the other 1/2
all |>
  mutate(trt = case_when(
    trt == "MIS" ~ "MIS",
    trt == "open" ~ "open",
    trt == 1 ~ "open",
    trt == 2 ~ "MIS"
  )) -> allfix

allfix_est_3yr <- bootstrap_event_free_band(
  allfix$dfsyrs,
  group      = allfix$trt,
  times_grid = 3,
  B          = 2000,
  seed       = 1
)

allfix_est_3yr

# save two tables for manuscript
write.csv(allfix_est_3yr, "allfix_est_3yr.csv", row.names = FALSE)
write.csv(shape_est_3yr, "shape_est_3yr.csv", row.names = FALSE)
write.csv(lacc_est_3yr, "lacc_est_3yr.csv", row.names = FALSE)

## check 10 000
allfix_est_3yr <- bootstrap_event_free_band(
  allfix$dfsyrs,
  group      = allfix$trt,
  times_grid = 3,
  B          = 100000,
  seed       = 1
)

allfix_est_3yr

## look at difference in 3 year estimates
bootstrap_event_free_diff(
  allfix$dfsyrs,
  group = allfix$trt,
  group1 = "open", group2 = "MIS", # use your actual stratum labels
  times_grid = 3,
  B = 2000, seed = 1
)

bootstrap_event_free_diff(
  allfix$dfsyrs,
  group = allfix$trt,
  group1 = "open", group2 = "MIS", # use your actual stratum labels
  conf_level = 0.975,
  times_grid = 3,
  B = 2000, seed = 1
)
