# ============================================================
# Age-specific prevalence data for Figure S6
# Figure A: Community population (Source.of.population = 2)
# Figure B: Hospital population (Source.of.population = 1)
#
# Age variable: Age1
# Outcomes: MASLD, S1, S2, S3, F1, F2, F3
#
# Each outcome uses its own non-missing denominator.
#
# Output:
# Figure_age_specific_prevalence_outcome_specific_N.csv
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
# Required variables
# ------------------------------------------------------------

required_vars <- c(
  "Age1",
  "Source.of.population",
  "MASLD",
  "S1",
  "S2",
  "S3",
  "F1",
  "F2",
  "F3"
)

missing_vars <- setdiff(required_vars, names(data))

if (length(missing_vars) > 0) {
  stop(
    "Missing variable(s): ",
    paste(missing_vars, collapse = ", ")
  )
}

# ------------------------------------------------------------
# Prepare variables
# ------------------------------------------------------------

data$Age1 <- trimws(as.character(data$Age1))

data$Source.of.population <- trimws(
  as.character(data$Source.of.population)
)

outcomes_left <- c("S1", "S2", "S3", "MASLD")
outcomes_right <- c("F1", "F2", "F3")

all_outcomes <- c(outcomes_left, outcomes_right)

for (v in all_outcomes) {
  data[[v]] <- suppressWarnings(
    as.numeric(trimws(as.character(data[[v]])))
  )
}

# Remove missing age-group values
data <- data[
  !is.na(data$Age1) &
    data$Age1 != "",
  ,
  drop = FALSE
]

if (nrow(data) == 0) {
  stop("No valid observations remained after removing missing Age1 values.")
}

# ------------------------------------------------------------
# Define age-group order and labels
# ------------------------------------------------------------

age_levels <- unique(data$Age1)

age_numbers <- suppressWarnings(as.numeric(age_levels))

if (all(!is.na(age_numbers))) {
  age_levels <- age_levels[order(age_numbers)]
}

# Supports  Age1 = 6/10/13/16
age_label_map <- c(
  "6" = "6–9",
  "10" = "10–12",
  "13" = "13–15",
  "16" = "16–18"
)

get_age_label <- function(x) {
  
  if (x %in% names(age_label_map)) {
    return(age_label_map[x])
  }
  
  return(x)
}

# ------------------------------------------------------------
# Calculate prevalence and exact binomial 95% CI
# ------------------------------------------------------------

calculate_prevalence <- function(outcome_vector) {
  
  valid_values <- outcome_vector[
    !is.na(outcome_vector)
  ]
  
  N <- length(valid_values)
  
  if (N == 0) {
    return(data.frame(
      Cases = NA,
      N = NA,
      Prevalence = NA,
      CI_lower = NA,
      CI_upper = NA,
      stringsAsFactors = FALSE
    ))
  }
  
  cases <- sum(valid_values == 1)
  
  CI <- binom.test(
    cases,
    N
  )$conf.int * 100
  
  data.frame(
    Cases = cases,
    N = N,
    Prevalence = cases / N * 100,
    CI_lower = CI[1],
    CI_upper = CI[2],
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# Create data for one population
# ------------------------------------------------------------

create_population_data <- function(
    population_data,
    population_name,
    figure_label
) {
  
  result_list <- list()
  result_id <- 1
  
  panel_info <- list(
    "Steatosis and MASLD" = outcomes_left,
    "Fibrosis" = outcomes_right
  )
  
  for (panel_name in names(panel_info)) {
    
    outcome_list <- panel_info[[panel_name]]
    
    for (outcome in outcome_list) {
      
      for (current_age in age_levels) {
        
        outcome_vector <- population_data[[outcome]][
          population_data$Age1 == current_age
        ]
        
        prevalence_result <- calculate_prevalence(
          outcome_vector
        )
        
        result_list[[result_id]] <- data.frame(
          Figure = figure_label,
          Population = population_name,
          Panel = panel_name,
          Outcome = outcome,
          Age1 = current_age,
          Age_group = get_age_label(current_age),
          N_definition = paste0(
            "Non-missing ",
            outcome,
            " values in this age group"
          ),
          prevalence_result,
          stringsAsFactors = FALSE
        )
        
        result_id <- result_id + 1
      }
    }
  }
  
  do.call(rbind, result_list)
}

# ------------------------------------------------------------
# Figure A: Community
# ------------------------------------------------------------

community_data <- data[
  data$Source.of.population == "2",
  ,
  drop = FALSE
]

if (nrow(community_data) == 0) {
  stop(
    "No Community observations found. Please check Source.of.population coding."
  )
}

figure_a_data <- create_population_data(
  population_data = community_data,
  population_name = "Community",
  figure_label = "A"
)

# ------------------------------------------------------------
# Figure B: Hospital
# ------------------------------------------------------------

hospital_data <- data[
  data$Source.of.population == "1",
  ,
  drop = FALSE
]

if (nrow(hospital_data) == 0) {
  stop(
    "No Hospital observations found. Please check Source.of.population coding."
  )
}

figure_b_data <- create_population_data(
  population_data = hospital_data,
  population_name = "Hospital",
  figure_label = "B"
)

# ------------------------------------------------------------
# Combine and export
# ------------------------------------------------------------

figure_data <- rbind(
  figure_a_data,
  figure_b_data
)

figure_data$Figure <- factor(
  figure_data$Figure,
  levels = c("A", "B")
)

figure_data$Panel <- factor(
  figure_data$Panel,
  levels = c("Steatosis and MASLD", "Fibrosis")
)

figure_data$Outcome <- factor(
  figure_data$Outcome,
  levels = c(
    "S1", "S2", "S3", "MASLD",
    "F1", "F2", "F3"
  )
)

write.csv(
  figure_data,
  "Figure_age_specific_prevalence_outcome_specific_N.csv",
  row.names = FALSE,
  na = ""
)

message("Completed successfully.")
message(
  "Created file: Figure_age_specific_prevalence_outcome_specific_N.csv"
)

# ============================================================
#FigS7
#============================================================
# Figure: Liver fibrosis prevalence across liver steatosis grades
#
# Figure A: Community population (Source.of.population = 2)
# Figure B: Hospital population (Source.of.population = 1)
#
# CAP_group coding:
# 1 = S0
# 2  = S1
# 3  = S2
# 4  = S3
#
# Outcomes:
# F1, F2, F3
#
# Output:
# Figure_steatosis_fibrosis_data.csv
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
# Check required variables
# ------------------------------------------------------------

required_vars <- c(
  "Source.of.population",
  "CAP_group",
  "F1",
  "F2",
  "F3"
)

missing_vars <- setdiff(required_vars, names(data))

if (length(missing_vars) > 0) {
  stop(
    "Missing variable(s): ",
    paste(missing_vars, collapse = ", ")
  )
}

# ------------------------------------------------------------
# Prepare variables
# ------------------------------------------------------------

data$Source.of.population <- trimws(
  as.character(data$Source.of.population)
)

data$CAP_group <- trimws(
  as.character(data$CAP_group)
)

for (v in c("F1", "F2", "F3")) {
  data[[v]] <- suppressWarnings(
    as.numeric(as.character(data[[v]]))
  )
}

# ------------------------------------------------------------
# Convert CAP_group into S0-S3 labels
# ------------------------------------------------------------

data$CAP_group_plot <- NA_character_

data$CAP_group_plot[
  data$CAP_group == "1"
] <- "S0"

data$CAP_group_plot[
  data$CAP_group == "2"
] <- "S1"

data$CAP_group_plot[
  data$CAP_group == "3"
] <- "S2"

data$CAP_group_plot[
  data$CAP_group == "4"
] <- "S3"

cap_levels <- c("S0", "S1", "S2", "S3")

# ------------------------------------------------------------
# Function: calculate prevalence and exact 95% CI
# ------------------------------------------------------------

calculate_prevalence <- function(outcome_vector) {
  
  # Only retain valid binary outcome values
  valid_values <- outcome_vector[
    outcome_vector %in% c(0, 1)
  ]
  
  N <- length(valid_values)
  
  if (N == 0) {
    return(
      data.frame(
        Cases = NA,
        N = NA,
        Prevalence = NA,
        CI_lower = NA,
        CI_upper = NA,
        stringsAsFactors = FALSE
      )
    )
  }
  
  cases <- sum(valid_values == 1)
  
  CI <- binom.test(
    cases,
    N
  )$conf.int * 100
  
  data.frame(
    Cases = cases,
    N = N,
    Prevalence = cases / N * 100,
    CI_lower = CI[1],
    CI_upper = CI[2],
    stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------
# Function: calculate data for one population
# ------------------------------------------------------------

create_population_data <- function(
    input_data,
    population_code,
    population_name,
    figure_label
) {
  
  population_data <- input_data[
    input_data$Source.of.population == population_code &
      !is.na(input_data$CAP_group_plot),
    ,
    drop = FALSE
  ]
  
  if (nrow(population_data) == 0) {
    stop(
      "No data found for ",
      population_name,
      ". Please check Source.of.population coding."
    )
  }
  
  # Number of participants in each CAP group
  steatosis_distribution <- table(
    factor(
      population_data$CAP_group_plot,
      levels = cap_levels
    )
  )
  
  outcome_map <- c(
    "F1" = "F1",
    "F2" = "F2",
    "F3" = "≥F3"
  )
  
  result_list <- list()
  result_id <- 1
  
  for (outcome_var in names(outcome_map)) {
    
    outcome_label <- outcome_map[[outcome_var]]
    
    for (cap_level in cap_levels) {
      
      outcome_vector <- population_data[[outcome_var]][
        population_data$CAP_group_plot == cap_level
      ]
      
      prevalence_result <- calculate_prevalence(
        outcome_vector
      )
      
      result_list[[result_id]] <- data.frame(
        Figure = figure_label,
        Population = population_name,
        CAP_group = cap_level,
        Steatosis_distribution_N = as.numeric(
          steatosis_distribution[cap_level]
        ),
        Outcome = outcome_label,
        Cases = prevalence_result$Cases,
        N = prevalence_result$N,
        Prevalence = prevalence_result$Prevalence,
        CI_lower = prevalence_result$CI_lower,
        CI_upper = prevalence_result$CI_upper,
        N_definition = paste0(
          "Participants with valid binary ",
          outcome_var,
          " data within this CAP group"
        ),
        stringsAsFactors = FALSE
      )
      
      result_id <- result_id + 1
    }
  }
  
  do.call(rbind, result_list)
}

# ------------------------------------------------------------
# Figure A: Community population
# ------------------------------------------------------------

figure_A <- create_population_data(
  input_data = data,
  population_code = "2",
  population_name = "Community",
  figure_label = "A"
)

# ------------------------------------------------------------
# Figure B: Hospital population
# ------------------------------------------------------------

figure_B <- create_population_data(
  input_data = data,
  population_code = "1",
  population_name = "Hospital",
  figure_label = "B"
)

# ------------------------------------------------------------
# Combine and sort data
# ------------------------------------------------------------

figure_data <- rbind(
  figure_A,
  figure_B
)

figure_data$Figure <- factor(
  figure_data$Figure,
  levels = c("A", "B")
)

figure_data$CAP_group <- factor(
  figure_data$CAP_group,
  levels = cap_levels
)

figure_data$Outcome <- factor(
  figure_data$Outcome,
  levels = c("F1", "F2", "≥F3")
)

figure_data <- figure_data[
  order(
    figure_data$Figure,
    figure_data$Outcome,
    figure_data$CAP_group
  ),
]

# ------------------------------------------------------------
# Export one CSV file
# ------------------------------------------------------------

write.csv(
  figure_data,
  "Figure_steatosis_fibrosis_data.csv",
  row.names = FALSE,
  na = ""
)

message("Completed successfully.")
message("Created file: Figure_steatosis_fibrosis_data.csv")
