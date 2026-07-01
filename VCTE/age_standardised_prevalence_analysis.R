# Age-standardised prevalence analysis for MASLD, steatosis and fibrosis outcomes
#
# Purpose:
#   1. Retrieve 2020 World Bank population weights for ages 6-18 years.
#   2. Calculate age-standardised MASLD prevalence overall and by population source.
#   3. Calculate age-standardised prevalence of S1-S3 and F1-F3 overall and by source.
#   4. Estimate 95% confidence intervals using bootstrap resampling.
#
# Required input variables in the selected CSV file:
#   Age, MASLD, Source.of.population, S1, S2, S3, F1, F2, F3
#
# Source.of.population coding:
#   1 = Hospital
#   2 = Community
#
# Output files:
#   age_standardized_prevalence_results.csv
#   age_standardized_prevalence_SF_results.csv
#
# Note:
#   This script does not include or upload participant-level data.

# -----------------------------
# 0. Packages
# -----------------------------
required_packages <- c("WDI", "dplyr", "tidyr", "stringr")

missing_packages <- setdiff(required_packages, rownames(installed.packages()))
if (length(missing_packages) > 0) {
  stop(
    "Please install the following R packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(WDI)
  library(dplyr)
  library(tidyr)
  library(stringr)
})

# -----------------------------
# 1. Settings
# -----------------------------
ref_year <- 2020
ages <- 6:18
bootstrap_replicates <- 1000
set.seed(123)

outcomes <- c("S1", "S2", "S3", "F1", "F2", "F3")

# -----------------------------
# 2. Retrieve World Bank standard population weights
# -----------------------------
female_indicators <- sprintf("SP.POP.AG%02d.FE.IN", ages)
male_indicators   <- sprintf("SP.POP.AG%02d.MA.IN", ages)
population_indicators <- c(female_indicators, male_indicators)

wb_raw <- WDI(
  country = "all",
  indicator = population_indicators,
  start = ref_year,
  end = ref_year,
  extra = TRUE
)

std_pop <- wb_raw %>%
  filter(region != "Aggregates") %>%
  select(all_of(population_indicators)) %>%
  pivot_longer(
    cols = everything(),
    names_to = "indicator",
    values_to = "population"
  ) %>%
  mutate(
    age = as.integer(str_extract(indicator, "(?<=AG)\\d{2}"))
  ) %>%
  group_by(age) %>%
  summarise(
    std_population = sum(population, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(age) %>%
  mutate(
    weight = std_population / sum(std_population)
  )

if (!isTRUE(all.equal(sum(std_pop$weight), 1))) {
  stop("Standard population weights do not sum to 1.")
}

# -----------------------------
# 3. Import and validate study data
# -----------------------------
message("Select the de-identified CSV analysis dataset.")
input_file <- file.choose()

df <- read.csv(
  input_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_columns <- c(
  "Age", "MASLD", "Source.of.population",
  "S1", "S2", "S3", "F1", "F2", "F3"
)

missing_columns <- setdiff(required_columns, names(df))
if (length(missing_columns) > 0) {
  stop(
    "The input dataset is missing the following required columns: ",
    paste(missing_columns, collapse = ", ")
  )
}

df <- df %>%
  mutate(
    Age = suppressWarnings(as.numeric(Age)),
    Source.of.population = suppressWarnings(as.numeric(Source.of.population)),
    across(all_of(c("MASLD", outcomes)), ~ suppressWarnings(as.numeric(.x)))
  )

# -----------------------------
# 4. Helper functions
# -----------------------------

# Calculate age-standardised prevalence for MASLD.
calc_age_std_masld <- function(data, std_pop) {
  age_prev <- data %>%
    group_by(age) %>%
    summarise(
      n = n(),
      cases = sum(MASLD, na.rm = TRUE),
      prev = cases / n,
      .groups = "drop"
    )

  std_pop %>%
    left_join(age_prev, by = "age") %>%
    summarise(
      age_std_prev = sum(weight * prev, na.rm = TRUE)
    ) %>%
    pull(age_std_prev)
}

# Bootstrap 95% CI for age-standardised MASLD prevalence.
boot_age_std_masld <- function(data, std_pop, B = 1000) {
  point_estimate <- calc_age_std_masld(data, std_pop)

  bootstrap_estimates <- replicate(B, {
    bootstrap_data <- data[
      sample(seq_len(nrow(data)), size = nrow(data), replace = TRUE),
    ]
    calc_age_std_masld(bootstrap_data, std_pop)
  })

  ci <- quantile(
    bootstrap_estimates,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  )

  data.frame(
    age_standardized_prevalence = point_estimate,
    lower_95CI = as.numeric(ci[1]),
    upper_95CI = as.numeric(ci[2])
  )
}

# Calculate age-standardised prevalence for a binary outcome.
calc_age_std_outcome <- function(data, outcome, std_pop) {
  data_use <- data %>%
    filter(!is.na(.data[[outcome]]))

  age_prev <- data_use %>%
    group_by(age) %>%
    summarise(
      n = n(),
      cases = sum(.data[[outcome]] == 1, na.rm = TRUE),
      prev = cases / n,
      .groups = "drop"
    )

  std_pop %>%
    left_join(age_prev, by = "age") %>%
    summarise(
      age_std_prev = sum(weight * prev, na.rm = FALSE)
    ) %>%
    pull(age_std_prev)
}

# Stratified bootstrap 95% CI for an age-standardised binary outcome.
boot_age_std_outcome <- function(data, outcome, std_pop, B = 1000) {
  data_use <- data %>%
    filter(!is.na(.data[[outcome]]))

  point_estimate <- calc_age_std_outcome(data_use, outcome, std_pop)

  bootstrap_estimates <- replicate(B, {
    bootstrap_data <- data_use %>%
      group_by(age) %>%
      slice_sample(prop = 1, replace = TRUE) %>%
      ungroup()

    calc_age_std_outcome(bootstrap_data, outcome, std_pop)
  })

  ci <- quantile(
    bootstrap_estimates,
    probs = c(0.025, 0.975),
    na.rm = TRUE
  )

  data.frame(
    outcome = outcome,
    age_standardized_prevalence = point_estimate,
    lower_95CI = as.numeric(ci[1]),
    upper_95CI = as.numeric(ci[2])
  )
}

# -----------------------------
# 5. MASLD age-standardised prevalence
# -----------------------------
df_masld <- df %>%
  filter(
    !is.na(Age),
    !is.na(MASLD),
    Age >= 6,
    Age < 19
  ) %>%
  mutate(
    age = floor(Age),
    source = case_when(
      Source.of.population == 1 ~ "Hospital",
      Source.of.population == 2 ~ "Community",
      TRUE ~ NA_character_
    )
  )

# Overall estimate includes all eligible participants with non-missing MASLD.
overall_masld_result <- boot_age_std_masld(
  data = df_masld,
  std_pop = std_pop,
  B = bootstrap_replicates
) %>%
  mutate(group = "Overall")

# Source-specific estimates include participants with a defined population source.
source_masld_result <- df_masld %>%
  filter(!is.na(source)) %>%
  group_by(source) %>%
  group_modify(~ boot_age_std_masld(
    data = .x,
    std_pop = std_pop,
    B = bootstrap_replicates
  )) %>%
  ungroup() %>%
  transmute(
    group = source,
    age_standardized_prevalence,
    lower_95CI,
    upper_95CI
  )

final_age_std_result <- bind_rows(
  overall_masld_result,
  source_masld_result
) %>%
  mutate(
    age_standardized_prevalence_percent = age_standardized_prevalence * 100,
    lower_95CI_percent = lower_95CI * 100,
    upper_95CI_percent = upper_95CI * 100
  )

# -----------------------------
# 6. S1-S3 and F1-F3 age-standardised prevalence
# -----------------------------
df_outcomes <- df %>%
  filter(
    !is.na(Age),
    Age >= 6,
    Age < 19
  ) %>%
  mutate(
    age = floor(Age),
    source = case_when(
      Source.of.population == 1 ~ "Hospital",
      Source.of.population == 2 ~ "Community",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(source))

final_outcome_age_std <- bind_rows(
  lapply(outcomes, function(outcome_name) {
    overall_result <- boot_age_std_outcome(
      data = df_outcomes,
      outcome = outcome_name,
      std_pop = std_pop,
      B = bootstrap_replicates
    ) %>%
      mutate(group = "Overall")

    source_result <- df_outcomes %>%
      group_by(source) %>%
      group_modify(~ boot_age_std_outcome(
        data = .x,
        outcome = outcome_name,
        std_pop = std_pop,
        B = bootstrap_replicates
      )) %>%
      ungroup() %>%
      transmute(
        outcome,
        group = source,
        age_standardized_prevalence,
        lower_95CI,
        upper_95CI
      )

    bind_rows(overall_result, source_result)
  })
) %>%
  mutate(
    age_standardized_prevalence_percent = age_standardized_prevalence * 100,
    lower_95CI_percent = lower_95CI * 100,
    upper_95CI_percent = upper_95CI * 100
  )

# -----------------------------
# 7. Save results
# -----------------------------
output_dir <- getwd()

write.csv(
  final_age_std_result,
  file = file.path(output_dir, "age_standardized_prevalence_results.csv"),
  row.names = FALSE
)

write.csv(
  final_outcome_age_std,
  file = file.path(output_dir, "age_standardized_prevalence_SF_results.csv"),
  row.names = FALSE
)

message("Analysis completed successfully.")
message("Output files were saved to: ", output_dir)
