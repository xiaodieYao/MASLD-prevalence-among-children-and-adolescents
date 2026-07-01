# ============================================================
# Figure 4 data:
## Community vs Hospital MASLD prevalence by subgroup
# Source.of.population: 1 = Hospital; 2 = Community
# MASLD: 0 = No; 1 = Yes
# ============================================================

rm(list = ls())

# Import data
data <- read.csv(
  file.choose(),
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Variables used in Figure 4
subgroup_vars <- c(
  "Age_group",
  "Sex",
  "Ethnicity",
  "BMI_group",
  "Glu_group",
  "ALT_group",
  "AST_group",
  "Region1"
)

panel_names <- c(
  Age_group = "Age (year)",
  Sex = "Sex",
  Ethnicity = "Ethnicity",
  BMI_group = "Weight status",
  Glu_group = "Glucose",
  ALT_group = "ALT",
  AST_group = "AST",
  Region1 = "Region"
)

# Check required variables
required_vars <- c(
  "Source.of.population",
  "MASLD",
  subgroup_vars
)

missing_vars <- setdiff(required_vars, names(data))

if (length(missing_vars) > 0) {
  stop(
    "Missing variable(s): ",
    paste(missing_vars, collapse = ", ")
  )
}

# Prepare variables
data$Source.of.population <- as.character(data$Source.of.population)
data$MASLD <- suppressWarnings(as.numeric(as.character(data$MASLD)))

for (v in subgroup_vars) {
  data[[v]] <- as.character(data[[v]])
}

# Retain valid Community and Hospital participants
data <- data[
  data$Source.of.population %in% c("1", "2") &
    data$MASLD %in% c(0, 1),
  ,
  drop = FALSE
]

# Calculate exact binomial 95% confidence interval
get_prevalence <- function(outcome_vector) {
  
  N <- length(outcome_vector)
  cases <- sum(outcome_vector == 1)
  
  if (N == 0) {
    return(list(
      cases = NA,
      N = NA,
      prevalence = NA,
      lower = NA,
      upper = NA
    ))
  }
  
  ci <- binom.test(cases, N)$conf.int * 100
  
  list(
    cases = cases,
    N = N,
    prevalence = cases / N * 100,
    lower = ci[1],
    upper = ci[2]
  )
}

# Full precision Pearson chi-squared P value
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
format_figure_p <- function(p_value) {
  
  if (is.na(p_value)) {
    return("")
  }
  
  if (p_value < 0.001) {
    return("<0.001")
  }
  
  sprintf("%.3f", p_value)
}

# Compare Community versus Hospital within one subgroup level
compare_prevalence <- function(community_outcome, hospital_outcome) {
  
  community_N <- length(community_outcome)
  hospital_N <- length(hospital_outcome)
  
  if (community_N == 0 || hospital_N == 0) {
    return(list(
      test = NA,
      statistic = NA,
      df = NA,
      p_value = NA,
      p_value_exact = NA,
      p_value_figure = "",
      min_expected_count = NA,
      expected_count_lt_5 = NA
    ))
  }
  
  community_cases <- sum(community_outcome == 1)
  hospital_cases <- sum(hospital_outcome == 1)
  
  test_table <- matrix(
    c(
      community_cases,
      community_N - community_cases,
      hospital_cases,
      hospital_N - hospital_cases
    ),
    nrow = 2,
    byrow = TRUE
  )
  
  rownames(test_table) <- c("Community", "Hospital")
  colnames(test_table) <- c("MASLD", "Non_MASLD")
  
  # No comparison is possible if all participants have the same outcome
  if (any(colSums(test_table) == 0)) {
    return(list(
      test = "Not estimable",
      statistic = NA,
      df = NA,
      p_value = NA,
      p_value_exact = NA,
      p_value_figure = "",
      min_expected_count = NA,
      expected_count_lt_5 = NA
    ))
  }
  
  chisq_result <- suppressWarnings(
    chisq.test(test_table, correct = FALSE)
  )
  
  min_expected <- min(chisq_result$expected)
  low_expected <- any(chisq_result$expected < 5)
  
  # Use Fisher's exact test when expected counts are small
  if (low_expected) {
    
    fisher_result <- fisher.test(
      test_table,
      alternative = "two.sided"
    )
    
    return(list(
      test = "Fisher's exact test",
      statistic = NA,
      df = NA,
      p_value = fisher_result$p.value,
      p_value_exact = format(
        fisher_result$p.value,
        scientific = TRUE,
        digits = 12,
        trim = TRUE
      ),
      p_value_figure = format_figure_p(fisher_result$p.value),
      min_expected_count = min_expected,
      expected_count_lt_5 = TRUE
    ))
  }
  
  list(
    test = "Pearson chi-squared test",
    statistic = unname(chisq_result$statistic),
    df = unname(chisq_result$parameter),
    p_value = unname(chisq_result$p.value),
    p_value_exact = format_chisq_p(
      unname(chisq_result$statistic),
      unname(chisq_result$parameter)
    ),
    p_value_figure = format_figure_p(chisq_result$p.value),
    min_expected_count = min_expected,
    expected_count_lt_5 = FALSE
  )
}

# Generate results for all panels
figure4_data <- list()
row_id <- 1

for (v in subgroup_vars) {
  
  valid_levels <- unique(
    data[[v]][
      !is.na(data[[v]]) &
        data[[v]] != ""
    ]
  )
  
  for (level in valid_levels) {
    
    community_outcome <- data$MASLD[
      data$Source.of.population == "2" &
        data[[v]] == level
    ]
    
    hospital_outcome <- data$MASLD[
      data$Source.of.population == "1" &
        data[[v]] == level
    ]
    
    community_result <- get_prevalence(community_outcome)
    hospital_result <- get_prevalence(hospital_outcome)
    
    comparison_result <- compare_prevalence(
      community_outcome,
      hospital_outcome
    )
    
    figure4_data[[row_id]] <- data.frame(
      panel = panel_names[v],
      subgroup_variable = v,
      subgroup_level = level,
      
      community_cases = community_result$cases,
      community_N = community_result$N,
      community_prevalence = community_result$prevalence,
      community_lower_95CI = community_result$lower,
      community_upper_95CI = community_result$upper,
      
      hospital_cases = hospital_result$cases,
      hospital_N = hospital_result$N,
      hospital_prevalence = hospital_result$prevalence,
      hospital_lower_95CI = hospital_result$lower,
      hospital_upper_95CI = hospital_result$upper,
      
      comparison = "Community vs Hospital",
      test = comparison_result$test,
      statistic = comparison_result$statistic,
      df = comparison_result$df,
      p_value = comparison_result$p_value,
      p_value_exact = comparison_result$p_value_exact,
      p_value_figure = comparison_result$p_value_figure,
      min_expected_count = comparison_result$min_expected_count,
      expected_count_lt_5 = comparison_result$expected_count_lt_5,
      
      stringsAsFactors = FALSE
    )
    
    row_id <- row_id + 1
  }
}

figure4_data <- do.call(rbind, figure4_data)
# ============================================================
# Add overall Community vs Hospital MASLD prevalence
# ============================================================

community_overall <- data$MASLD[
  data$Source.of.population == "2"
]

hospital_overall <- data$MASLD[
  data$Source.of.population == "1"
]

community_overall_result <- get_prevalence(community_overall)
hospital_overall_result <- get_prevalence(hospital_overall)

overall_comparison_result <- compare_prevalence(
  community_overall,
  hospital_overall
)

overall_row <- data.frame(
  panel = "Overall",
  subgroup_variable = "Overall",
  subgroup_level = "Overall",
  
  community_cases = community_overall_result$cases,
  community_N = community_overall_result$N,
  community_prevalence = community_overall_result$prevalence,
  community_lower_95CI = community_overall_result$lower,
  community_upper_95CI = community_overall_result$upper,
  
  hospital_cases = hospital_overall_result$cases,
  hospital_N = hospital_overall_result$N,
  hospital_prevalence = hospital_overall_result$prevalence,
  hospital_lower_95CI = hospital_overall_result$lower,
  hospital_upper_95CI = hospital_overall_result$upper,
  
  comparison = "Community vs Hospital",
  test = overall_comparison_result$test,
  statistic = overall_comparison_result$statistic,
  df = overall_comparison_result$df,
  p_value = overall_comparison_result$p_value,
  p_value_exact = overall_comparison_result$p_value_exact,
  p_value_figure = overall_comparison_result$p_value_figure,
  min_expected_count = overall_comparison_result$min_expected_count,
  expected_count_lt_5 = overall_comparison_result$expected_count_lt_5,
  
  stringsAsFactors = FALSE
)

# Put Overall row at the top
figure4_data <- rbind(
  overall_row,
  figure4_data
)

write.csv(
  figure4_data,
  "Figure4_Community_vs_Hospital_prevalence.csv",
  row.names = FALSE
)


# ============================================================
# Figure 5: Annual MASLD prevalence among adolescents aged 12-18 years
# Age_group = 3
# Year variable: Year1 (binary variable)
#
# Output:
# Figure5_annual_MASLD_age12_18.csv
# ============================================================

rm(list = ls())

# ------------------------------------------------------------
# Import data
# ------------------------------------------------------------

data <- read.csv(
  file.choose(),
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# ------------------------------------------------------------
# Variables
# ------------------------------------------------------------

subgroup_vars <- c(
  "Source.of.population",
  "Sex",
  "Race",
  "BMI_group",
  "ALT_group",
  "AST_group"
)

panel_names <- c(
  "Population source",
  "Sex",
  "Race",
  "Weight status",
  "ALT",
  "AST"
)

required_vars <- c(
  "Age_group",
  "Year1",
  "MASLD",
  subgroup_vars
)

missing_vars <- setdiff(required_vars, names(data))

if (length(missing_vars) > 0) {
  stop(
    "Missing variable(s): ",
    paste(missing_vars, collapse = ", ")
  )
}

# ------------------------------------------------------------
# Prepare data
# ------------------------------------------------------------

data$Age_group <- trimws(as.character(data$Age_group))
data$Year1 <- trimws(as.character(data$Year1))
data$MASLD <- suppressWarnings(as.numeric(as.character(data$MASLD)))

for (v in subgroup_vars) {
  data[[v]] <- trimws(as.character(data[[v]]))
}

# Keep adolescents aged 12-18 years only.
# Remove missing Year1 values and invalid/missing MASLD values.
adolescent_data <- data[
  data$Age_group == "3" &
    !is.na(data$Year1) &
    data$Year1 != "" &
    data$MASLD %in% c(0, 1),
  ,
  drop = FALSE
]

if (nrow(adolescent_data) == 0) {
  stop("No valid observations remained after filtering Age_group = 3.")
}

# Identify the two Year1 categories automatically
year_levels <- unique(adolescent_data$Year1)

# Sort Year1 values if they are numeric
year_numbers <- suppressWarnings(as.numeric(year_levels))

if (all(!is.na(year_numbers))) {
  year_levels <- year_levels[order(year_numbers)]
}

if (length(year_levels) != 2) {
  stop(
    "Year1 must contain exactly two non-missing categories after filtering. Found: ",
    paste(year_levels, collapse = ", ")
  )
}

# ------------------------------------------------------------
# Helper functions
# ------------------------------------------------------------

format_subgroup_label <- function(variable, value) {
  
  if (variable == "Source.of.population") {
    if (value == "1") return("Hospital")
    if (value == "2") return("Community")
  }
  
  return(value)
}

calculate_prevalence <- function(outcome_vector) {
  
  N <- length(outcome_vector)
  cases <- sum(outcome_vector == 1)
  
  if (N == 0) {
    return(list(
      cases = NA,
      N = NA,
      prevalence = NA,
      CI_lower = NA,
      CI_upper = NA
    ))
  }
  
  CI <- binom.test(cases, N)$conf.int * 100
  
  list(
    cases = cases,
    N = N,
    prevalence = cases / N * 100,
    CI_lower = CI[1],
    CI_upper = CI[2]
  )
}

# Full scientific-notation P value from Pearson chi-squared test
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
  
  sprintf("%.12fe%+d", mantissa, as.integer(exponent))
}

# P value for plotting
format_figure_p <- function(statistic, df) {
  
  log_p <- pchisq(
    statistic,
    df = df,
    lower.tail = FALSE,
    log.p = TRUE
  )
  
  if (is.na(log_p) || !is.finite(log_p)) {
    return("")
  }
  
  if (log_p < log(0.001)) {
    return("<0.001")
  }
  
  sprintf("%.3f", exp(log_p))
}

# Pearson chi-squared test between the two Year1 categories
run_year_chisq <- function(dat) {
  
  d <- dat[
    dat$Year1 %in% year_levels &
      dat$MASLD %in% c(0, 1),
    c("Year1", "MASLD"),
    drop = FALSE
  ]
  
  tab <- table(
    factor(d$Year1, levels = year_levels),
    factor(d$MASLD, levels = c(0, 1))
  )
  
  if (
    any(rowSums(tab) == 0) ||
    any(colSums(tab) == 0)
  ) {
    return(list(
      test = "Not estimable",
      chi2_statistic = NA,
      df = NA,
      p_value_exact = NA,
      p_value_display = "",
      min_expected_count = NA,
      any_expected_count_lt_5 = NA
    ))
  }
  
  test_result <- suppressWarnings(
    chisq.test(tab, correct = FALSE)
  )
  
  statistic <- unname(test_result$statistic)
  df <- unname(test_result$parameter)
  
  list(
    test = "Pearson chi-squared test",
    chi2_statistic = statistic,
    df = df,
    p_value_exact = format_exact_p(statistic, df),
    p_value_display = format_figure_p(statistic, df),
    min_expected_count = min(test_result$expected),
    any_expected_count_lt_5 = any(test_result$expected < 5)
  )
}

# Create two Year1 rows for one line in the figure
analyse_line <- function(
    dat,
    panel,
    subgroup,
    variable = NULL,
    subgroup_value = NULL
) {
  
  d <- dat
  
  if (!is.null(variable)) {
    d <- d[
      d[[variable]] == subgroup_value,
      ,
      drop = FALSE
    ]
  }
  
  stats <- run_year_chisq(d)
  
  output <- lapply(seq_along(year_levels), function(i) {
    
    current_year <- year_levels[i]
    outcome <- d$MASLD[d$Year1 == current_year]
    result <- calculate_prevalence(outcome)
    
    # Put statistics only on the first Year1 row
    if (i == 1) {
      test <- stats$test
      chi2_statistic <- stats$chi2_statistic
      df <- stats$df
      p_value_exact <- stats$p_value_exact
      P_value <- stats$p_value_display
      min_expected_count <- stats$min_expected_count
      any_expected_count_lt_5 <- stats$any_expected_count_lt_5
    } else {
      test <- ""
      chi2_statistic <- NA
      df <- NA
      p_value_exact <- ""
      P_value <- ""
      min_expected_count <- NA
      any_expected_count_lt_5 <- NA
    }
    
    data.frame(
      Panel = panel,
      Subgroup = subgroup,
      Year1 = current_year,
      Cases = result$cases,
      N = result$N,
      Prevalence = result$prevalence,
      CI_lower = result$CI_lower,
      CI_upper = result$CI_upper,
      test = test,
      Chi2_statistic = chi2_statistic,
      df = df,
      p_value_exact = p_value_exact,
      P_value = P_value,
      min_expected_count = min_expected_count,
      any_expected_count_lt_5 = any_expected_count_lt_5,
      stringsAsFactors = FALSE
    )
  })
  
  do.call(rbind, output)
}

# ------------------------------------------------------------
# Generate one combined CSV
# ------------------------------------------------------------

result_list <- list()
result_id <- 1

for (i in seq_along(subgroup_vars)) {
  
  variable <- subgroup_vars[i]
  panel <- panel_names[i]
  
  # Add overall prevalence to the first panel
  if (variable == "Source.of.population") {
    
    result_list[[result_id]] <- analyse_line(
      dat = adolescent_data,
      panel = panel,
      subgroup = "Overall"
    )
    
    result_id <- result_id + 1
  }
  
  subgroup_levels <- unique(
    adolescent_data[[variable]][
      !is.na(adolescent_data[[variable]]) &
        adolescent_data[[variable]] != ""
    ]
  )
  
  for (level in subgroup_levels) {
    
    result_list[[result_id]] <- analyse_line(
      dat = adolescent_data,
      panel = panel,
      subgroup = format_subgroup_label(variable, level),
      variable = variable,
      subgroup_value = level
    )
    
    result_id <- result_id + 1
  }
}

figure5_data <- do.call(rbind, result_list)

# ------------------------------------------------------------
# Export one CSV file
# ------------------------------------------------------------

write.csv(
  figure5_data,
  "Figure5_annual_MASLD_age12_18.csv",
  row.names = FALSE,
  na = ""
)

message("Completed successfully.")
message("Created file: Figure5_annual_MASLD_age12_18.csv")

