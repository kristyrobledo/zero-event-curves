#' Event free curve
#'
#' Reverse Kaplan-Meier censoring distribution and expected event-free curve
#'
#' Implements the two quantities described for studies with zero observed
#' events:
#'   G(t_i) = G(t_{i-1}) * (1 - c_i / n_i)
#'                       [reverse-KM censoring distribution]
#'   E(t_i) = G(t_{i-1}) * (1 - 1 / (2.5 * n_i))
#'                       [expected event-free curve, Quigley adjustment]
#'
#' where n_i is the number of subjects still at risk (under follow-up)
#' IMMEDIATELY BEFORE t_i, and c_i is the number of subjects censored
#' (lost to follow-up / reaching end of observation) exactly at t_i.
#'
#' The function does not assume the input is sorted, does not assume the
#' data starts at time 0, and aggregates ties automatically (multiple
#' subjects censored at the same time are treated as a single c_i > 1).
#'
#' @param time A numeric vector of each subject's last-observed (censoring)
#'   time. One value per subject -- do not pre-aggregate. Every subject in
#'   the dataset is assumed to be censored (this method is for the
#'   zero-events case), so `length(time)` is taken as the total sample size N
#'   (or the size of each group, when `group` is supplied).
#' @param group Optional grouping vector, same length as `time`. When
#'   supplied, the curve is computed separately within each group and the
#'   results are row-bound together with an added `group` column.
#'
#' @return A data.frame with one row per distinct observed time (plus a
#'   leading t = 0 anchor row when the data does not already start at 0),
#'   containing:
#'     time              - the distinct censoring time
#'     n_at_risk_before  - n_i, subjects at risk immediately before this time
#'     n_censored        - c_i, subjects censored at this time
#'     censoring_dist    - G(t_i), the reverse-KM censoring distribution
#'     event_free_curve  - E(t_i), the expected event-free probability
#'     group             - only present when `group` is supplied
#'
#' @examples
#' # Paper's own illustrative example: 15 subjects, 3 lost at 4mo,
#' # 7 at 12mo, 2 at 19mo, 3 still under follow-up (administratively
#' # censored) at 24mo. Deliberately unsorted / not starting at 0.
#' toy_times <- c(12, 12, 12, 12, 12, 12, 12, 4, 4, 4, 24, 24, 24, 19, 19)
#' event_free_curve(toy_times)
event_free_curve <- function(time, group = NULL) {
  if (!is.null(group)) {
    if (length(group) != length(time)) { # nolint: indentation_linter.
      stop("`time` and `group` must be the same length.")
    }
    splits <- split(time, group)
    out <- lapply(names(splits), function(g) {
      d <- event_free_curve(splits[[g]])
      d$group <- g
      d
    })
    return(do.call(rbind, out))
  }

  if (!is.numeric(time)) stop("`time` must be a numeric vector.")
  if (anyNA(time)) stop("`time` contains missing values.")
  if (any(time < 0)) stop("`time` values must be non-negative.")

  utime_orig <- sort(unique(time))
  core <- .curve_core(time, utime_orig)

  out <- data.frame(
    time             = core$time,
    n_at_risk_before = core$n_before,
    n_censored       = core$cens,
    censoring_dist   = core$G,
    event_free_curve = core$E
  )

  # Prepend a t = 0 reference row (G = E = 1) only if the data doesn't
  # already start there -- if a subject is genuinely censored at time 0,
  # that time point is itself t_1 and already carries n_i = N.
  if (min(core$time) > 0) {
    anchor <- data.frame( # nolint: indentation_linter.
      time = 0, n_at_risk_before = length(time), n_censored = 0,
      censoring_dist = 1, event_free_curve = 1
    )
    out <- rbind(anchor, out)
  }

  rownames(out) <- NULL
  out
}

#' Internal, fully-vectorized core computation (no data.frame overhead)
#'
#' Given raw subject times and a fixed sorted grid of the distinct times
#' that can occur (`utime_grid`), returns the observed distinct times in
#' `time` (a subset of `utime_grid`) with their n_before/cens/G/E vectors.
#' Used by both `event_free_curve()` and the bootstrap function, so the
#' bootstrap never has to repeat `sort(unique(...))` or use `table()`
#' (both of which are the main bottlenecks in a naive R implementation).
.curve_core <- function(time, utime_grid) {
  N <- length(time)
  idx <- match(time, utime_grid)
  counts <- tabulate(idx, nbins = length(utime_grid))
  pos <- counts > 0
  rt <- utime_grid[pos]
  rc <- counts[pos]
  k <- length(rt)

  if (sum(rc) != N) {
    stop("Censoring counts do not sum to the total sample size (check for duplicate IDs / NAs).") # nolint: line_length_linter.
  }

  n_before <- N - c(0, cumsum(rc)[-k])
  Gstep <- 1 - rc / n_before
  Gcum <- cumprod(Gstep)
  Gprev <- c(1, Gcum[-k])
  E <- Gprev * (1 - 1 / (2.5 * n_before))

  list(time = rt, n_before = n_before, cens = rc, G = Gcum, E = E)
}

#' Evaluate a fitted curve (censoring_dist or event_free_curve) at a
#' specific time point
#'
#' Both curves are step functions, constant on [t_i, t_{i+1}). To get the
#' value at an arbitrary time (e.g. 3 years), find the last observed step
#' at or before that time and carry it forward.
#'
#' @param curve_df Output of `event_free_curve()` (single group).
#' @param time_point The time at which to read off the curve.
#' @param column Which column to evaluate: "event_free_curve" (default) or
#'   "censoring_dist".
#'
#' @return A single numeric value, or NA (with a warning) if `time_point`
#'   is beyond the last observed time in `curve_df` -- the method gives no
#'   information past the last censoring time, so this is extrapolation
#'   and is deliberately not supported.
evaluate_curve_at <- function(curve_df, time_point, column = "event_free_curve") {
  if (time_point < min(curve_df$time)) {
    return(1) # before anything has happened, everyone is still event-free
  }
  if (time_point > max(curve_df$time)) {
    warning("`time_point` is beyond the last observed time; returning NA (no extrapolation).") # nolint: line_length_linter.
    return(NA_real_)
  }
  idx <- max(which(curve_df$time <= time_point))
  curve_df[[column]][idx]
}

#' Bootstrap confidence band (or single-point CI) for the event-free curve
#'
#' Fast nonparametric bootstrap: resamples subjects with replacement,
#' recomputes the curve on each resample using the same vectorized
#' `.curve_core()` used by `event_free_curve()`, and reads off ALL
#' requested grid times from each resample in one vectorized lookup
#' (`findInterval`) rather than looping over times. This means a full
#' confidence band -- or a single time point -- costs the same B
#' resamples; the cost does not scale with the number of grid points.
#'
#' Efficiency notes:
#'   - `table()` is replaced with `match()` + `tabulate()` against a grid
#'     computed once, up front, outside the bootstrap loop.
#'   - No data.frame / rbind is built inside the loop.
#'   - The cumulative product is vectorized (`cumprod`) instead of an
#'     explicit per-time-step for() loop.
#'   - All requested times (`times_grid`) are evaluated in one
#'     `findInterval()` call per resample.
#'
#' @param time Numeric vector of each subject's censoring time.
#' @param group Optional grouping vector, same length as `time`. When
#'   supplied, resampling and the curve/CI are computed separately within
#'   each group (each group is bootstrapped from its own subjects only),
#'   and the results are row-bound together with an added `group` column.
#'   `times_grid = NULL` then defaults to each group's own observed times.
#' @param times_grid Vector of times at which to estimate the curve /
#'   band. Defaults to every distinct observed time in `time` (i.e. the
#'   natural step points) if left NULL. Pass a single value for a
#'   single-time CI (e.g. `times_grid = 3`).
#' @param column "event_free_curve" (default) or "censoring_dist".
#' @param B Number of bootstrap resamples (default 2000).
#' @param conf_level Confidence level for the percentile interval (default 0.95).
#' @param seed Optional seed for reproducibility.
#'
#' @return A data.frame with one row per element of `times_grid` (per
#'   group, if `group` is supplied):
#'   time, estimate (from the original data), se, ci_lower, ci_upper,
#'   n_boot_used (resamples whose own follow-up reached that time --
#'   this can vary across grid points, and is worth checking near the
#'   end of follow-up where it may drop off), group (only if supplied).
bootstrap_event_free_band <- function(time, group = NULL, times_grid = NULL,
                                      column = "event_free_curve",
                                      B = 2000, conf_level = 0.95, seed = NULL) {
  if (!is.null(group)) {
    if (length(group) != length(time)) {
      stop("`time` and `group` must be the same length.")
    }
    splits <- split(time, group)
    out <- lapply(names(splits), function(g) {
      d <- bootstrap_event_free_band(splits[[g]],
        times_grid = times_grid, column = column,
        B = B, conf_level = conf_level, seed = seed
      )
      d$group <- g
      d
    })
    return(do.call(rbind, out))
  }

  if (!is.null(seed)) set.seed(seed)

  N <- length(time)
  utime_orig <- sort(unique(time))
  if (is.null(times_grid)) times_grid <- utime_orig
  times_grid <- sort(times_grid)

  point_curve <- event_free_curve(time)
  point_est <- vapply(times_grid, function(tp) {
    evaluate_curve_at(point_curve, tp, column)
  }, numeric(1))

  ng <- length(times_grid)
  mat <- matrix(NA_real_, nrow = B, ncol = ng)

  for (b in seq_len(B)) {
    samp <- sample(time, size = N, replace = TRUE)
    core <- .curve_core(samp, utime_orig)

    vals <- if (column == "event_free_curve") core$E else core$G
    times_ext <- c(0, core$time)
    vals_ext <- c(1, vals)

    idx <- findInterval(times_grid, times_ext)
    idx[idx == 0] <- 1
    row_vals <- vals_ext[idx]
    row_vals[times_grid > max(core$time)] <- NA_real_ # no extrapolation
    mat[b, ] <- row_vals
  }

  n_dropped <- colSums(is.na(mat))
  if (any(n_dropped > 0)) {
    warning(sprintf(
      "Some bootstrap resamples did not reach some grid times (max dropped: %d of %d, at t = %s). Their CI is based on fewer resamples -- see n_boot_used.", # nolint: line_length_linter.
      max(n_dropped), B, times_grid[which.max(n_dropped)]
    ))
  }

  alpha <- 1 - conf_level
  se <- apply(mat, 2, stats::sd, na.rm = TRUE)
  ci <- apply(mat, 2, stats::quantile, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
  n_used <- apply(mat, 2, function(x) sum(!is.na(x)))

  data.frame(
    time        = times_grid,
    estimate    = point_est,
    se          = se,
    ci_lower    = ci[1, ],
    ci_upper    = ci[2, ],
    n_boot_used = n_used
  )
}

#' Convenience wrapper: bootstrap CI at a single time point
#'
#' Thin wrapper around `bootstrap_event_free_band()` for the common case
#' of wanting one number (e.g. "what's the estimate and CI at 3 years?").
#' Returns a list rather than a one-row data.frame for easy `$` access.
#' Not vectorized over `group` -- with `group` supplied this returns a
#' data.frame (one row per group) instead, since there's more than one
#' number to return.
bootstrap_event_free <- function(time, time_point, group = NULL,
                                 column = "event_free_curve",
                                 B = 2000, conf_level = 0.95, seed = NULL) {
  row <- bootstrap_event_free_band(time,
    group = group, times_grid = time_point, column = column,
    B = B, conf_level = conf_level, seed = seed
  )
  if (!is.null(group)) {
    return(row)
  }
  list(
    estimate = row$estimate, se = row$se,
    ci_lower = row$ci_lower, ci_upper = row$ci_upper,
    n_boot_used = row$n_boot_used
  )
}

# ---- quick self-test against the paper's own worked example ----------------
if (sys.nframe() == 0) {
  toy_times <- c(12, 12, 12, 12, 12, 12, 12, 4, 4, 4, 24, 24, 24, 19, 19)
  print(event_free_curve(toy_times))
  # Expected censoring_dist at t=4,12,19,24: 0.8, 0.333333, 0.2, 0
  # Expected event_free_curve at same times: 0.973333, 0.773333, 0.306667, 0.173333

  # Single time point
  boot <- bootstrap_event_free(toy_times, time_point = 20, B = 2000, seed = 1)
  cat(sprintf(
    "\nAt t = 20: estimate = %.4f, SE = %.4f, 95%% CI = (%.4f, %.4f), n_boot_used = %d\n",
    boot$estimate, boot$se, boot$ci_lower, boot$ci_upper, boot$n_boot_used
  ))

  # Vector of time points
  band <- bootstrap_event_free_band(toy_times, times_grid = c(5, 10, 15, 20), B = 2000, seed = 1)
  print(band)

  # Grouped example: curve and bootstrap band computed within each group
  set.seed(2)
  grp_time <- c(toy_times, toy_times * 1.3)
  grp <- rep(c("Control", "Treatment"), each = length(toy_times))

  print(event_free_curve(grp_time, group = grp))
  print(bootstrap_event_free_band(grp_time, group = grp, B = 1000, seed = 2))
  print(bootstrap_event_free(grp_time, time_point = 15, group = grp, B = 1000, seed = 2))
}

#' Internal: run B bootstrap resamples of `time` and read off `times_grid`
#'
#' Shared by `bootstrap_event_free_band()` and `bootstrap_event_free_diff()`
#' so the resampling loop only exists once. Returns the raw B x
#' length(times_grid) matrix of resampled curve values (NA past a given
#' resample's own last observed time), left unsummarized so callers can
#' either take per-column CIs (`bootstrap_event_free_band()`) or difference
#' two such matrices column-by-column (`bootstrap_event_free_diff()`).
.bootstrap_matrix <- function(time, times_grid, column, B) {
  N <- length(time)
  utime_orig <- sort(unique(time))
  ng <- length(times_grid)
  mat <- matrix(NA_real_, nrow = B, ncol = ng)

  for (b in seq_len(B)) {
    samp <- sample(time, size = N, replace = TRUE)
    core <- .curve_core(samp, utime_orig)

    vals <- if (column == "event_free_curve") core$E else core$G
    times_ext <- c(0, core$time)
    vals_ext <- c(1, vals)

    idx <- findInterval(times_grid, times_ext)
    idx[idx == 0] <- 1
    row_vals <- vals_ext[idx]
    row_vals[times_grid > max(core$time)] <- NA_real_ # no extrapolation
    mat[b, ] <- row_vals
  }
  mat
}

#' Bootstrap the difference between two groups' curves, with a 95% CI
#'
#' Independently bootstraps `group1` and `group2` at the requested time(s)
#' (same B for both, drawn from each group's own subjects) and takes the
#' difference resample-by-resample -- the standard percentile-bootstrap
#' approach for a difference between two independent estimators. Pairing
#' resample b of group1 with resample b of group2 is arbitrary (they're
#' independent draws), but that's fine: what matters is the resulting
#' distribution of differences, not any particular pairing.
#'
#' @param time Numeric vector of each subject's censoring time (both groups).
#' @param group Grouping vector, same length as `time`.
#' @param group1,group2 The two levels of `group` to compare; the
#'   reported difference is `group1` minus `group2`.
#' @param times_grid Vector of times to evaluate the difference at.
#'   Defaults to the distinct observed times pooled across both groups.
#'   Pass a single value for a difference at one time point (e.g. 3 years).
#' @param column "event_free_curve" (default) or "censoring_dist".
#' @param B Number of bootstrap resamples per group (default 2000).
#' @param conf_level Confidence level for the percentile interval (default 0.95).
#' @param seed Optional seed for reproducibility (set once, shared across
#'   both groups' resampling).
#'
#' @return A data.frame with one row per element of `times_grid`:
#'   time, estimate_group1, estimate_group2, diff, se, ci_lower, ci_upper,
#'   n_boot_used (resamples where BOTH groups' curves reached that time --
#'   worth checking near the end of follow-up, same as `n_boot_used` in
#'   `bootstrap_event_free_band()`).
bootstrap_event_free_diff <- function(time, group, group1, group2,
                                      times_grid = NULL,
                                      column = "event_free_curve",
                                      B = 2000, conf_level = 0.95, seed = NULL) {
  if (length(group) != length(time)) {
    stop("`time` and `group` must be the same length.")
  }
  if (!(group1 %in% group) || !(group2 %in% group)) {
    stop("`group1` and `group2` must both be levels present in `group`.")
  }
  if (identical(group1, group2)) stop("`group1` and `group2` must differ.")

  t1 <- time[group == group1]
  t2 <- time[group == group2]

  if (is.null(times_grid)) times_grid <- sort(unique(time))
  times_grid <- sort(times_grid)

  curve1 <- event_free_curve(t1)
  curve2 <- event_free_curve(t2)
  est1 <- vapply(times_grid, function(tp) evaluate_curve_at(curve1, tp, column), numeric(1))
  est2 <- vapply(times_grid, function(tp) evaluate_curve_at(curve2, tp, column), numeric(1))

  if (!is.null(seed)) set.seed(seed)
  mat1 <- .bootstrap_matrix(t1, times_grid, column, B)
  mat2 <- .bootstrap_matrix(t2, times_grid, column, B)
  diff_mat <- mat1 - mat2

  n_dropped <- colSums(is.na(diff_mat))
  if (any(n_dropped > 0)) {
    warning(sprintf(
      "Some bootstrap resamples did not reach some grid times in one or both groups (max dropped: %d of %d, at t = %s). The CI there is based on fewer resamples -- see n_boot_used.",
      max(n_dropped), B, times_grid[which.max(n_dropped)]
    ))
  }

  alpha <- 1 - conf_level
  se <- apply(diff_mat, 2, stats::sd, na.rm = TRUE)
  ci <- apply(diff_mat, 2, stats::quantile, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
  n_used <- apply(diff_mat, 2, function(x) sum(!is.na(x)))

  data.frame(
    time            = times_grid,
    estimate_group1 = est1,
    estimate_group2 = est2,
    diff            = est1 - est2,
    se              = se,
    ci_lower        = ci[1, ],
    ci_upper        = ci[2, ],
    n_boot_used     = n_used
  )
}

#' Convenience wrapper: bootstrap CI at a single time point
#'
#' Thin wrapper around `bootstrap_event_free_band()` for the common case
#' of wanting one number (e.g. "what's the estimate and CI at 3 years?").
#' Returns a list rather than a one-row data.frame for easy `$` access.
#' Not vectorized over `group` -- with `group` supplied this returns a
#' data.frame (one row per group) instead, since there's more than one
#' number to return.
bootstrap_event_free <- function(time, time_point, group = NULL,
                                 column = "event_free_curve",
                                 B = 2000, conf_level = 0.95, seed = NULL) {
  row <- bootstrap_event_free_band(time,
    group = group, times_grid = time_point, column = column,
    B = B, conf_level = conf_level, seed = seed
  )
  if (!is.null(group)) {
    return(row)
  }
  list(
    estimate = row$estimate, se = row$se,
    ci_lower = row$ci_lower, ci_upper = row$ci_upper,
    n_boot_used = row$n_boot_used
  )
}
