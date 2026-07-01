# ============================================================
# MASLD,S1-3 Prevalence, 95% CI and chi-squared tests for forest plots
# Hospital = 1; Community = 2
# ============================================================
library(dplyr)
data <- read.csv(
  file.choose(),
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

subgroup_vars <- c(
  "Age_group", "Sex", "Ethnicity", "BMI_group",
  "Glu_group", "ALT_group", "AST_group"
)

outcome_vars <- c("S1", "S2", "S3", "MASLD")

data <- data[, !is.na(names(data)) & names(data) != "", drop = FALSE]

hospital_data <- data[
  as.character(data$Source.of.population) == "1",
  ,
  drop = FALSE
]

community_data <- data[
  as.character(data$Source.of.population) == "2",
  ,
  drop = FALSE
]

subgroup_vars <- c(
  "Age_group", "Sex", "Ethnicity", "BMI_group",
  "Glu_group", "ALT_group", "AST_group"
)

outcome_vars <- c("S1", "S2", "S3", "MASLD")

# Calculate prevalence and exact 95% binomial CI
calculate_prevalence <- function(dat, subgroup_var, outcome_var, source_name) {
  
  dat %>%
    transmute(
      subgroup_level = .data[[subgroup_var]],
      outcome = .data[[outcome_var]]
    ) %>%
    filter(!is.na(subgroup_level), outcome %in% c(0, 1)) %>%
    group_by(subgroup_level) %>%
    summarise(
      cases = sum(outcome == 1),
      total = n(),
      .groups = "drop"
    ) %>%
    rowwise() %>%
    mutate(
      prevalence = cases / total,
      CI_lower = binom.test(cases, total)$conf.int[1],
      CI_upper = binom.test(cases, total)$conf.int[2]
    ) %>%
    ungroup() %>%
    transmute(
      source = source_name,
      subgroup = subgroup_var,
      subgroup_level = as.character(subgroup_level),
      outcome = outcome_var,
      cases = cases,
      total = total,
      prevalence_percent = prevalence * 100,
      CI_lower_percent = CI_lower * 100,
      CI_upper_percent = CI_upper * 100
    )
}

# Format very small chi-squared P values
format_chisq_p <- function(statistic, df, p_value) {
  
  if (is.na(p_value)) return(NA_character_)
  
  if (p_value > 0) {
    return(format(p_value, scientific = TRUE, digits = 12, trim = TRUE))
  }
  
  log_p <- pchisq(
    statistic,
    df = df,
    lower.tail = FALSE,
    log.p = TRUE
  )
  
  exponent <- floor(log_p / log(10))
  mantissa <- exp(log_p - exponent * log(10))
  
  paste0(formatC(mantissa, format = "f", digits = 8), "e", exponent)
}

# Pearson chi-squared test for prevalence differences across subgroups
run_chisq <- function(dat, subgroup_var, outcome_var, source_name) {
  
  d <- data.frame(
    subgroup = dat[[subgroup_var]],
    outcome = dat[[outcome_var]]
  )
  
  d <- d[
    !is.na(d$subgroup) & d$outcome %in% c(0, 1),
    ,
    drop = FALSE
  ]
  
  tab <- table(d$subgroup, d$outcome)
  
  if (nrow(tab) < 2 || ncol(tab) < 2) {
    return(data.frame(
      source = source_name,
      subgroup = subgroup_var,
      outcome = outcome_var,
      test = "Pearson chi-squared test",
      statistic = NA,
      df = NA,
      p_value = NA,
      p_value_display = NA,
      n_total = nrow(d)
    ))
  }
  
  result <- suppressWarnings(
    chisq.test(tab, correct = FALSE)
  )
  
  data.frame(
    source = source_name,
    subgroup = subgroup_var,
    outcome = outcome_var,
    test = "Pearson chi-squared test",
    statistic = unname(result$statistic),
    df = unname(result$parameter),
    p_value = unname(result$p.value),
    p_value_display = format_chisq_p(
      unname(result$statistic),
      unname(result$parameter),
      unname(result$p.value)
    ),
    n_total = nrow(d)
  )
}

# Calculate prevalence data for one population source
make_prevalence_data <- function(dat, source_name) {
  
  bind_rows(lapply(subgroup_vars, function(g) {
    bind_rows(lapply(outcome_vars, function(o) {
      calculate_prevalence(dat, g, o, source_name)
    }))
  }))
}

# Calculate chi-squared results for one population source
make_chisq_data <- function(dat, source_name) {
  
  bind_rows(lapply(subgroup_vars, function(g) {
    bind_rows(lapply(outcome_vars, function(o) {
      run_chisq(dat, g, o, source_name)
    }))
  }))
}

# Community and Hospital prevalence data
forest_data <- bind_rows(
  make_prevalence_data(community_data, "Community"),
  make_prevalence_data(hospital_data, "Hospital")
)

# Community and Hospital chi-squared results
chisq_results <- bind_rows(
  make_chisq_data(community_data, "Community"),
  make_chisq_data(hospital_data, "Hospital")
)

# Add chi-squared statistics and P values to plotting data
forest_data <- forest_data %>%
  left_join(
    chisq_results %>%
      select(
        source, subgroup, outcome,
        statistic, df, p_value_display
      ),
    by = c("source", "subgroup", "outcome")
  )

write.csv(
  forest_data,
  file = "R_forest_plot_data_S.csv",
  row.names = FALSE
)

write.csv(
  chisq_results,
  file = "R_subgroup_chisq_results_S.csv",
  row.names = FALSE
)

message("Completed successfully.")
message("Created files:")
message("  - R_forest_plot_data_S.csv")
message("  - R_subgroup_chisq_results_S.csv")


# ============================================================
# Fibrosis prevalence data and chi-squared tests
# Hospital = 1; Community = 2
# Outcomes: F1, F2 and F3
# ============================================================

subgroup_vars <- c(
  "Age_group", "Sex", "Ethnicity", "BMI_group",
  "Glu_group", "ALT_group", "AST_group"
)

outcome_vars <- c("F1", "F2", "F3")

hospital_data <- data[
  as.character(data$Source.of.population) == "1",
  ,
  drop = FALSE
]

community_data <- data[
  as.character(data$Source.of.population) == "2",
  ,
  drop = FALSE
]

# Calculate prevalence and exact 95% binomial confidence intervals
calculate_prevalence <- function(dat, subgroup_var, outcome_var, source_name) {
  
  d <- dat[, c(subgroup_var, outcome_var), drop = FALSE]
  names(d) <- c("subgroup_level", "outcome")
  
  d$subgroup_level <- as.character(d$subgroup_level)
  d$outcome <- suppressWarnings(as.numeric(as.character(d$outcome)))
  
  d <- d[
    !is.na(d$subgroup_level) &
      d$subgroup_level != "" &
      d$outcome %in% c(0, 1),
    ,
    drop = FALSE
  ]
  
  level_list <- unique(d$subgroup_level)
  
  result <- lapply(level_list, function(level_name) {
    
    x <- d[d$subgroup_level == level_name, , drop = FALSE]
    
    cases <- sum(x$outcome == 1)
    N <- nrow(x)
    ci <- binom.test(cases, N)$conf.int * 100
    
    data.frame(
      source = source_name,
      subgroup = subgroup_var,
      subgroup_level = level_name,
      outcome = outcome_var,
      cases = cases,
      N = N,
      prevalence_percent = cases / N * 100,
      CI_lower_percent = ci[1],
      CI_upper_percent = ci[2],
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, result)
}

# Full-precision P value in scientific notation
format_chisq_p <- function(statistic, df) {
  
  log_p <- pchisq(
    statistic,
    df = df,
    lower.tail = FALSE,
    log.p = TRUE
  )
  
  if (is.na(log_p) || !is.finite(log_p)) {
    return(NA_character_)
  }
  
  exponent <- floor(log_p / log(10))
  mantissa <- exp(log_p - exponent * log(10))
  
  sprintf("%.12fe%+d", mantissa, as.integer(exponent))
}

# Figure-friendly P value
format_figure_p <- function(statistic, df) {
  
  log_p <- pchisq(
    statistic,
    df = df,
    lower.tail = FALSE,
    log.p = TRUE
  )
  
  if (is.na(log_p) || !is.finite(log_p)) {
    return(NA_character_)
  }
  
  if (log_p < log(0.001)) {
    return("<0.001")
  }
  
  sprintf("%.3f", exp(log_p))
}

# Pearson chi-squared test for prevalence differences across subgroups
run_chisq <- function(dat, subgroup_var, outcome_var, source_name) {
  
  d <- dat[, c(subgroup_var, outcome_var), drop = FALSE]
  names(d) <- c("subgroup_level", "outcome")
  
  d$subgroup_level <- as.character(d$subgroup_level)
  d$outcome <- suppressWarnings(as.numeric(as.character(d$outcome)))
  
  d <- d[
    !is.na(d$subgroup_level) &
      d$subgroup_level != "" &
      d$outcome %in% c(0, 1),
    ,
    drop = FALSE
  ]
  
  tab <- table(
    d$subgroup_level,
    factor(d$outcome, levels = c(0, 1))
  )
  
  if (nrow(tab) < 2 || any(colSums(tab) == 0)) {
    return(data.frame(
      source = source_name,
      subgroup = subgroup_var,
      outcome = outcome_var,
      test = "Pearson chi-squared test",
      statistic = NA,
      df = NA,
      p_value = NA,
      p_value_full = NA,
      p_value_figure = NA,
      n_total = nrow(d),
      min_expected_count = NA,
      any_expected_count_lt_5 = NA,
      stringsAsFactors = FALSE
    ))
  }
  
  test <- suppressWarnings(chisq.test(tab, correct = FALSE))
  
  data.frame(
    source = source_name,
    subgroup = subgroup_var,
    outcome = outcome_var,
    test = "Pearson chi-squared test",
    statistic = unname(test$statistic),
    df = unname(test$parameter),
    p_value = unname(test$p.value),
    p_value_full = format_chisq_p(
      unname(test$statistic),
      unname(test$parameter)
    ),
    p_value_figure = format_figure_p(
      unname(test$statistic),
      unname(test$parameter)
    ),
    n_total = nrow(d),
    min_expected_count = min(test$expected),
    any_expected_count_lt_5 = any(test$expected < 5),
    stringsAsFactors = FALSE
  )
}

# Run all outcomes and subgroups for one population source
make_prevalence_data <- function(dat, source_name) {
  
  do.call(rbind, lapply(subgroup_vars, function(g) {
    do.call(rbind, lapply(outcome_vars, function(o) {
      calculate_prevalence(dat, g, o, source_name)
    }))
  }))
}

make_chisq_data <- function(dat, source_name) {
  
  do.call(rbind, lapply(subgroup_vars, function(g) {
    do.call(rbind, lapply(outcome_vars, function(o) {
      run_chisq(dat, g, o, source_name)
    }))
  }))
}

# Community and Hospital prevalence data
fibrosis_prevalence_data <- rbind(
  make_prevalence_data(community_data, "Community"),
  make_prevalence_data(hospital_data, "Hospital")
)

# Community and Hospital chi-squared results
fibrosis_chisq_results <- rbind(
  make_chisq_data(community_data, "Community"),
  make_chisq_data(hospital_data, "Hospital")
)

write.csv(
  fibrosis_prevalence_data,
  "R_forest_plot_data_F.csv",
  row.names = FALSE
)

write.csv(
  fibrosis_chisq_results,
  "R_subgroup_chisq_results_F.csv",
  row.names = FALSE
)

message("Completed successfully.")
message("Created files:")
message("1. R_forest_plot_data_F.csv")
message("2. R_subgroup_chisq_results_F.csv")

# ============================================================
# Total population MASLD prevalence
# Hospital + Community combined
# ============================================================

subgroup_vars <- c(
  "Age_group", "Sex", "Ethnicity", "BMI_group",
  "Glu_group", "ALT_group", "AST_group"
)

# Keep Hospital and Community participants with valid MASLD status
total_data <- data[
  as.character(data$Source.of.population) %in% c("1", "2") &
    as.character(data$MASLD) %in% c("0", "1"),
  ,
  drop = FALSE
]

total_data$MASLD <- suppressWarnings(
  as.numeric(as.character(total_data$MASLD))
)

# Calculate overall MASLD prevalence
total_cases <- sum(total_data$MASLD == 1, na.rm = TRUE)
total_N <- sum(!is.na(total_data$MASLD))
total_CI <- binom.test(total_cases, total_N)$conf.int * 100

overall_total_MASLD <- data.frame(
  subgroup = "Overall",
  subgroup_level = "Overall",
  cases = total_cases,
  N = total_N,
  prevalence_percent = total_cases / total_N * 100,
  CI_lower_percent = total_CI[1],
  CI_upper_percent = total_CI[2],
  stringsAsFactors = FALSE
)

# Calculate MASLD prevalence for each subgroup
calculate_total_MASLD <- function(dat, subgroup_var) {
  
  d <- dat[, c(subgroup_var, "MASLD"), drop = FALSE]
  names(d) <- c("subgroup_level", "MASLD")
  
  d$subgroup_level <- as.character(d$subgroup_level)
  d$MASLD <- suppressWarnings(as.numeric(as.character(d$MASLD)))
  
  d <- d[
    !is.na(d$subgroup_level) &
      d$subgroup_level != "" &
      d$MASLD %in% c(0, 1),
    ,
    drop = FALSE
  ]
  
  result <- lapply(unique(d$subgroup_level), function(level_name) {
    
    x <- d[d$subgroup_level == level_name, , drop = FALSE]
    
    cases <- sum(x$MASLD == 1)
    N <- nrow(x)
    CI <- binom.test(cases, N)$conf.int * 100
    
    data.frame(
      subgroup = subgroup_var,
      subgroup_level = level_name,
      cases = cases,
      N = N,
      prevalence_percent = cases / N * 100,
      CI_lower_percent = CI[1],
      CI_upper_percent = CI[2],
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, result)
}

# Pearson chi-squared test for subgroup differences
run_total_chisq <- function(dat, subgroup_var) {
  
  d <- dat[, c(subgroup_var, "MASLD"), drop = FALSE]
  names(d) <- c("subgroup_level", "MASLD")
  
  d$subgroup_level <- as.character(d$subgroup_level)
  d$MASLD <- suppressWarnings(as.numeric(as.character(d$MASLD)))
  
  d <- d[
    !is.na(d$subgroup_level) &
      d$subgroup_level != "" &
      d$MASLD %in% c(0, 1),
    ,
    drop = FALSE
  ]
  
  tab <- table(
    d$subgroup_level,
    factor(d$MASLD, levels = c(0, 1))
  )
  
  test <- chisq.test(tab, correct = FALSE)
  
  data.frame(
    subgroup = subgroup_var,
    test = "Pearson chi-squared test",
    statistic = unname(test$statistic),
    df = unname(test$parameter),
    p_value = unname(test$p.value),
    p_value_display = format(
      test$p.value,
      scientific = TRUE,
      digits = 12,
      trim = TRUE
    ),
    n_total = nrow(d),
    stringsAsFactors = FALSE
  )
}

# Generate prevalence data
total_MASLD_data <- do.call(
  rbind,
  lapply(subgroup_vars, function(v) {
    calculate_total_MASLD(total_data, v)
  })
)

total_MASLD_data <- rbind(
  overall_total_MASLD,
  total_MASLD_data
)

# Generate chi-squared results
total_MASLD_chisq <- do.call(
  rbind,
  lapply(subgroup_vars, function(v) {
    run_total_chisq(total_data, v)
  })
)

# Export files
write.csv(
  total_MASLD_data,
  "R_forest_plot_data_total_MASLD.csv",
  row.names = FALSE
)

write.csv(
  total_MASLD_chisq,
  "R_subgroup_chisq_results_total_MASLD.csv",
  row.names = FALSE
)


# ============================================================
# Hospital ethnicity comparison excluding ethnicity code 5
# Outcomes: MASLD, S1, S2, S3, F1, F2, F3
# Output: chi-squared statistic, df, exact P value and N
# ============================================================

outcome_vars <- c("MASLD", "S1", "S2", "S3", "F1", "F2", "F3")

# Keep Hospital participants only and exclude Ethnicity code 5
hospital_data <- data[
  as.character(data$Source.of.population) == "1" &
    as.character(data$Ethnicity) != "5" &
    !is.na(data$Ethnicity),
  ,
  drop = FALSE
]

# Convert outcomes to numeric
for (v in outcome_vars) {
  hospital_data[[v]] <- suppressWarnings(
    as.numeric(as.character(hospital_data[[v]]))
  )
}

# Format extremely small P values without showing 0
format_exact_p <- function(statistic, df) {
  
  log_p <- pchisq(
    statistic,
    df = df,
    lower.tail = FALSE,
    log.p = TRUE
  )
  
  if (is.na(log_p) || !is.finite(log_p)) {
    return(NA_character_)
  }
  
  exponent <- floor(log_p / log(10))
  mantissa <- exp(log_p - exponent * log(10))
  
  sprintf("%.12fe%d", mantissa, as.integer(exponent))
}

# Pearson chi-squared test across the remaining ethnicity groups
run_ethnicity_chisq <- function(outcome_var) {
  
  d <- hospital_data[, c("Ethnicity", outcome_var), drop = FALSE]
  names(d) <- c("Ethnicity", "Outcome")
  
  d$Ethnicity <- as.character(d$Ethnicity)
  d$Outcome <- suppressWarnings(as.numeric(as.character(d$Outcome)))
  
  d <- d[
    !is.na(d$Ethnicity) &
      d$Ethnicity != "" &
      d$Outcome %in% c(0, 1),
    ,
    drop = FALSE
  ]
  
  tab <- table(
    d$Ethnicity,
    factor(d$Outcome, levels = c(0, 1))
  )
  
  test <- chisq.test(tab, correct = FALSE)
  
  data.frame(
    source = "Hospital",
    subgroup = "Ethnicity excluding code 5",
    outcome = outcome_var,
    test = "Pearson chi-squared test",
    statistic = unname(test$statistic),
    df = unname(test$parameter),
    p_value = unname(test$p.value),
    p_value_exact = format_exact_p(
      unname(test$statistic),
      unname(test$parameter)
    ),
    n_total = nrow(d),
    min_expected_count = min(test$expected),
    any_expected_count_lt_5 = any(test$expected < 5),
    stringsAsFactors = FALSE
  )
}

hospital_ethnicity_statistics <- do.call(
  rbind,
  lapply(outcome_vars, run_ethnicity_chisq)
)

write.csv(
  hospital_ethnicity_statistics,
  "Hospital_ethnicity_excluding5_statistics.csv",
  row.names = FALSE
)
