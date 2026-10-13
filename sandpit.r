## check lacc data
source("utils.R")
library(readxl)
library(tidyverse)
library(rlang)

shape <- readxl::read_excel("C:/Users/kman2804/OneDrive - The University of Sydney (Staff)/1. RESEARCH PROJECTS/3. GYNO ONC/CRE ALIVEWELL/LACC_shape/Shape/SHAPE_balance_table_KM_coordinates.xlsx",
  sheet = "SHAPEM_by_surgery group"
)
shape$Stratum <- recode(as.character(shape$Stratum),
  "MIS" = "MIS",
  "open" = "open"
)

# remove row at time 0
shape <- shape[shape$Time_months > 0, ]
table(shape$Stratum)
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
lacc$AlloTreat <- recode(as.character(lacc$AlloTreat),
  "1" = "open",
  "2" = "MIS"
)
head(lacc)

table(lacc$AlloTreat)
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


# save two tables for manuscript
write.csv(all_est_3yr, "all_est_3yr.csv", row.names = FALSE)
write.csv(shape_est_3yr, "shape_est_3yr.csv", row.names = FALSE)
write.csv(lacc_est_3yr, "lacc_est_3yr.csv", row.names = FALSE)

## get median FU time from inverse KM by trial
library(survival)

Surv_lacc <- Surv(lacc$dfsdays / 365.25, rep(1, nrow(lacc))) # Assuming all events are observed
survfit(Surv_lacc ~ 1)

# get median survival

Surv_shape <- Surv(shape$Time_years, rep(1, nrow(shape))) # Assuming all events are observed
survfit(Surv_shape ~ 1)


## check 10 000
all_est_3yr10000 <- bootstrap_event_free_band(
  all$dfsyrs,
  group      = all$trt,
  times_grid = 3,
  B          = 100000,
  seed       = 1
)

all_est_3yr10000

## maturity of data
followup_maturity(
  data = lacc,
  followup = dfsdays / 365.25,
  horizon = 3,
  group = AlloTreat
)

followup_maturity(
  data = shape,
  followup = Time_years,
  horizon = 3,
  group = Stratum
)

## look at difference in 3 year estimates
diff.all <- bootstrap_event_free_diff(
  all$dfsyrs,
  group = all$trt,
  group1 = "MIS", group2 = "open", # use your actual stratum labels
  times_grid = 3,
  conf_level = 0.95,
  B = 2000, seed = 1
)


## look at difference in 3 year estimates for shape
shape.diff <- bootstrap_event_free_diff(
  shape$Time_years,
  group = shape$Stratum,
  group1 = "MIS", group2 = "open", # use your actual stratum labels
  times_grid = 3,
  conf_level = 0.95,
  B = 2000, seed = 1
)

# look at difference in 3 year estimates for lacc
lacc.diff <- bootstrap_event_free_diff(
  lacc$dfsdays,
  group = lacc$AlloTreat,
  group1 = "MIS", group2 = "open",
  times_grid = 3 * 365.25,
  conf_level = 0.95,
  B = 2000, seed = 1
)

## forestplot for shape and lacc
library(forestplot)
## stick together estimates for forest plot
forest_data <- data.frame(
  cohort = c("SHAPE", "LACC", "Combined"),
  MIS = c(shape.diff$estimate_group1, lacc.diff$estimate_group1, diff.all$estimate_group1),
  open = c(shape.diff$estimate_group2, lacc.diff$estimate_group2, diff.all$estimate_group2),
  estimate = c(shape.diff$diff, lacc.diff$diff, diff.all$diff),
  se = c(shape.diff$se, lacc.diff$se, diff.all$se),
  ci_lower = c(shape.diff$ci_lower, lacc.diff$ci_lower, diff.all$ci_lower),
  ci_upper = c(shape.diff$ci_upper, lacc.diff$ci_upper, diff.all$ci_upper)
)
write.csv(forest_data, "forest_data-3yrs.csv", row.names = FALSE)

# rescale to percentages and round to 1 decimal place
forest_data$MIS <- round(forest_data$MIS * 100, 1)
forest_data$open <- round(forest_data$open * 100, 1)
forest_data$estimate <- round(forest_data$estimate * 100, 5)
forest_data$ci_lower <- round(forest_data$ci_lower * 100, 5)
forest_data$ci_upper <- round(forest_data$ci_upper * 100, 5)


forest_data$estimate_ci <- sprintf(
  "%.1f (%.1f to %.1f)",
  forest_data$estimate,
  forest_data$ci_lower,
  forest_data$ci_upper
)

forest_labels <- rbind(
  c("Cohort", "Open", "MIS", "Difference (95% CI)"),
  as.matrix(forest_data[, c("cohort", "open", "MIS", "estimate_ci")])
)

forestplot(
  labeltext = forest_labels,
  mean = c(NA, forest_data$estimate),
  lower = c(NA, forest_data$ci_lower),
  upper = c(NA, forest_data$ci_upper),
  is.summary = c(TRUE, rep(FALSE, nrow(forest_data))),
  zero = 0,
  xticks = seq(-20, 30, by = 5),
  xlab = "Difference in 3-year Event-free probability (MIS - open)",
  boxsize = 0.2,
  lwd.ci = 2,
  ci.vertices = TRUE,
  ci.vertices.height = 0.1,
  col = fpColors(box = "#2C7FB8", line = "#2C7FB8", summary = "#253494"),
  txt_gp = fpTxtGp(
    label = gpar(cex = 0.9),
    ticks = gpar(cex = 0.8),
    xlab = gpar(cex = 0.9),
    title = gpar(cex = 1)
  )
)
