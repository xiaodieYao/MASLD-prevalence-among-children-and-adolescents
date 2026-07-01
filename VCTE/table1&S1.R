rm(list = ls())

library(tableone)

# Import data
input_file <- file.choose()
data <- read.csv(
  input_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# Display p values without showing zero
format_p <- function(p) {
  if (is.na(p)) return(NA_character_)
  if (p == 0) return("<2.2e-16")
  format(p, scientific = TRUE, digits = 12, trim = TRUE)
}

# ============================================================
# Table 1: Community versus Hospital
# ============================================================

m1 <- CreateTableOne(
  vars = c(
    "Age", "Age_group", "Sex", "Ethnicity", "LSM", "LSM_group",
    "CAP", "CAP_group", "BMI_zscore", "BMI_group",
    "Blood.pressure", "Glucose", "Glu_group", "ALT", "ALT_group",
    "AST", "AST_group", "GGT", "PLT", "Total.cholesterol",
    "HDL.cholesterol", "Triglycerides", "LDL.cholesterol", "MASLD"
  ),
  strata = "Source.of.population",
  data = data,
  factorVars = c(
    "Age_group", "Sex", "Ethnicity", "LSM_group", "CAP_group",
    "BMI_group", "Blood.pressure", "Glu_group", "ALT_group",
    "AST_group", "MASLD"
  ),
  addOverall = TRUE,
  testApprox = chisq.test,
  argsApprox = list(correct = FALSE),
  testNonNormal = wilcox.test,
  argsNonNormal = list(
    alternative = "two.sided",
    exact = FALSE,
    correct = FALSE
  )
)

m1_csv <- print(
  m1,
  smd = TRUE,
  showAllLevels = TRUE,
  quote = FALSE,
  nonnormal = TRUE,
  noSpaces = TRUE,
  printToggle = FALSE
)

write.csv(m1_csv, file = "m1.csv")

# Full statistics for Table 1
# Hospital = 1; Community = 2

continuous_vars_m1 <- c(
  "Age", "LSM", "CAP", "BMI_zscore", "Glucose", "ALT", "AST",
  "GGT", "PLT", "Total.cholesterol", "HDL.cholesterol",
  "Triglycerides", "LDL.cholesterol"
)

categorical_vars_m1 <- c(
  "Age_group", "Sex", "Ethnicity", "LSM_group", "CAP_group",
  "BMI_group", "Blood.pressure", "Glu_group", "ALT_group",
  "AST_group", "MASLD"
)

stats_m1 <- data[
  as.character(data$Source.of.population) %in% c("1", "2"),
  ,
  drop = FALSE
]

cat_results_m1 <- lapply(categorical_vars_m1, function(v) {
  
  d <- stats_m1[, c("Source.of.population", v), drop = FALSE]
  d <- d[complete.cases(d), , drop = FALSE]
  
  test <- chisq.test(
    table(d[[v]], d$Source.of.population),
    correct = FALSE
  )
  
  data.frame(
    variable = v,
    test = "Pearson chi-squared test",
    statistic_name = "X-squared",
    statistic = unname(test$statistic),
    df = unname(test$parameter),
    p_value = unname(test$p.value),
    p_value_display = format_p(test$p.value),
    n_community = sum(as.character(d$Source.of.population) == "2"),
    n_hospital = sum(as.character(d$Source.of.population) == "1"),
    n_total = nrow(d),
    stringsAsFactors = FALSE
  )
})

cont_results_m1 <- lapply(continuous_vars_m1, function(v) {
  
  d <- stats_m1[, c("Source.of.population", v), drop = FALSE]
  d[[v]] <- suppressWarnings(as.numeric(as.character(d[[v]])))
  d <- d[complete.cases(d), , drop = FALSE]
  
  community_values <- d[[v]][as.character(d$Source.of.population) == "2"]
  hospital_values <- d[[v]][as.character(d$Source.of.population) == "1"]
  
  test <- wilcox.test(
    community_values,
    hospital_values,
    alternative = "two.sided",
    exact = FALSE,
    correct = FALSE
  )
  
  data.frame(
    variable = v,
    test = "Wilcoxon rank-sum test",
    statistic_name = "W",
    statistic = unname(test$statistic),
    df = NA,
    p_value = unname(test$p.value),
    p_value_display = format_p(test$p.value),
    n_community = length(community_values),
    n_hospital = length(hospital_values),
    n_total = nrow(d),
    stringsAsFactors = FALSE
  )
})

full_statistics_m1 <- do.call(
  rbind,
  c(cat_results_m1, cont_results_m1)
)

write.csv(
  full_statistics_m1,
  file = "m1_full_statistics.csv",
  row.names = FALSE
)

# ============================================================
# Table S1: Non-MASLD versus MASLD
# ============================================================

m2 <- CreateTableOne(
  vars = c(
    "Age", "Age_group", "Sex", "Ethnicity", "Source.of.population",
    "LSM", "LSM_group", "CAP", "CAP_group", "BMI_zscore",
    "BMI_group", "Blood.pressure", "Glucose", "Glu_group",
    "ALT", "ALT_group", "AST", "AST_group", "GGT", "PLT",
    "Total.cholesterol", "HDL.cholesterol", "Triglycerides",
    "LDL.cholesterol"
  ),
  strata = "MASLD",
  data = data,
  factorVars = c(
    "Age_group", "Sex", "Ethnicity", "Source.of.population",
    "LSM_group", "CAP_group", "BMI_group", "Blood.pressure",
    "Glu_group", "ALT_group", "AST_group"
  ),
  addOverall = TRUE,
  testApprox = chisq.test,
  argsApprox = list(correct = FALSE),
  testNonNormal = wilcox.test,
  argsNonNormal = list(
    alternative = "two.sided",
    exact = FALSE,
    correct = FALSE
  )
)

m2_csv <- print(
  m2,
  smd = TRUE,
  showAllLevels = TRUE,
  quote = FALSE,
  nonnormal = TRUE,
  noSpaces = TRUE,
  printToggle = FALSE
)

write.csv(m2_csv, file = "m2.csv")

# Full statistics for Table S1
# MASLD = 0: Non-MASLD; MASLD = 1: MASLD

continuous_vars_m2 <- c(
  "Age", "LSM", "CAP", "BMI_zscore", "Glucose", "ALT", "AST",
  "GGT", "PLT", "Total.cholesterol", "HDL.cholesterol",
  "Triglycerides", "LDL.cholesterol"
)

categorical_vars_m2 <- c(
  "Age_group", "Sex", "Ethnicity", "Source.of.population",
  "LSM_group", "CAP_group", "BMI_group", "Blood.pressure",
  "Glu_group", "ALT_group", "AST_group"
)

stats_m2 <- data[
  as.character(data$MASLD) %in% c("0", "1"),
  ,
  drop = FALSE
]

cat_results_m2 <- lapply(categorical_vars_m2, function(v) {
  
  d <- stats_m2[, c("MASLD", v), drop = FALSE]
  d <- d[complete.cases(d), , drop = FALSE]
  
  test <- chisq.test(
    table(d[[v]], d$MASLD),
    correct = FALSE
  )
  
  data.frame(
    variable = v,
    test = "Pearson chi-squared test",
    statistic_name = "X-squared",
    statistic = unname(test$statistic),
    df = unname(test$parameter),
    p_value = unname(test$p.value),
    p_value_display = format_p(test$p.value),
    n_non_MASLD = sum(as.character(d$MASLD) == "0"),
    n_MASLD = sum(as.character(d$MASLD) == "1"),
    n_total = nrow(d),
    stringsAsFactors = FALSE
  )
})

cont_results_m2 <- lapply(continuous_vars_m2, function(v) {
  
  d <- stats_m2[, c("MASLD", v), drop = FALSE]
  d[[v]] <- suppressWarnings(as.numeric(as.character(d[[v]])))
  d <- d[complete.cases(d), , drop = FALSE]
  
  non_masld_values <- d[[v]][as.character(d$MASLD) == "0"]
  masld_values <- d[[v]][as.character(d$MASLD) == "1"]
  
  test <- wilcox.test(
    non_masld_values,
    masld_values,
    alternative = "two.sided",
    exact = FALSE,
    correct = FALSE
  )
  
  data.frame(
    variable = v,
    test = "Wilcoxon rank-sum test",
    statistic_name = "W",
    statistic = unname(test$statistic),
    df = NA,
    p_value = unname(test$p.value),
    p_value_display = format_p(test$p.value),
    n_non_MASLD = length(non_masld_values),
    n_MASLD = length(masld_values),
    n_total = nrow(d),
    stringsAsFactors = FALSE
  )
})

full_statistics_m2 <- do.call(
  rbind,
  c(cat_results_m2, cont_results_m2)
)

write.csv(
  full_statistics_m2,
  file = "m2_full_statistics.csv",
  row.names = FALSE
)

message("Completed successfully.")
message("Created files: m1.csv, m1_full_statistics.csv, m2.csv, m2_full_statistics.csv")

