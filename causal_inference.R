# Load necessary libraries
library(dplyr)
library(ggplot2)
library(randomForest)
library(grf)
library(causalweight)

# Convert variables to appropriate types
data_sample$A <- as.numeric(as.character(data_sample$A))
data_out_of_sample$A <- as.numeric(as.character(data_out_of_sample$A))
data_sample$group <- as.factor(data_sample$group)
data_out_of_sample$group <- as.factor(data_out_of_sample$group)


# Function using logistic regression and linear regression
estimate_nuisance_params_logistic_linear <- function(data_sample, data_out_of_sample) {
  set.seed(sample.int(1e9, 1))
  # Combine labeled and unlabeled datasets
  combined_data <- rbind(
    data_sample[, setdiff(names(data_sample), "y")],
    data_out_of_sample
  )
  
  
  # Propensity score estimation using logistic regression
  covariate_names <- paste0("X", 1:11)
  formula_ps <- as.formula(paste("A ~", paste(covariate_names, collapse = " + ")))
  
  # Fit logistic regression on combined data
  model_ps <- glm(formula_ps, data = combined_data, family = binomial(link = "logit"))
  
  # Get propensity score predictions for both datasets
  propensity_scores <- predict(model_ps, newdata = combined_data, type = "response")
  
  # Outcome model estimation using linear regression
  formula_om <- as.formula(paste("y ~ A +", paste(covariate_names, collapse = " + ")))
  
  # Fit linear regression on data_sample
  model_om <- lm(formula_om, data = data_sample)
  
  # Get outcome predictions for combined data
  combined_data_no_y <- combined_data
  combined_data_no_y$y <- NULL  # Remove 'y' if present
  predicted_outcomes_all <- predict(model_om, newdata = combined_data_no_y)
  
  # Get outcome predictions for unlabeled data only
  n_sample <- nrow(data_sample)
  n_combined <- nrow(combined_data)
  predicted_outcomes_unlabeled <- predicted_outcomes_all[(n_sample + 1):n_combined]
  
  # Return propensity scores and predicted outcomes
  return(list(
    propensity_scores = propensity_scores,
    predicted_outcomes_all = predicted_outcomes_all,
    predicted_outcomes_unlabeled = predicted_outcomes_unlabeled
  ))
}


# Function using random forests
estimate_nuisance_params_random_forest <- function(data_sample, data_out_of_sample) {
  # Ensure required packages are installed and loaded
  if (!require("randomForest", quietly = TRUE)) install.packages("randomForest")
  library(randomForest)
  
  # Combine labeled and unlabeled datasets
  combined_data <- rbind(
    data_sample[, setdiff(names(data_sample), "y")],
    data_out_of_sample
  )
  
  # Propensity score estimation using random forest
  covariate_names <- paste0("X", 1:11)
  formula_ps <- as.formula(paste("factor(A) ~", paste(covariate_names, collapse = " + ")))
  
  # Convert A to factor for classification
  combined_data$A <- as.factor(combined_data$A)
  
  # Fit random forest for propensity score estimation
  model_ps <- randomForest(formula_ps, data = combined_data)
  
  # Get propensity score predictions (probabilities) for both datasets
  propensity_scores <- predict(model_ps, newdata = combined_data, type = "prob")[, 2]  # Probability of A=1
  
  # Outcome model estimation using random forest
  formula_om <- as.formula(paste("y ~ A +", paste(covariate_names, collapse = " + ")))
  
  # Ensure that A is numeric in data_sample
  data_sample$A <- as.numeric(as.character(data_sample$A))
  
  # Fit random forest regression model on data_sample
  model_om <- randomForest(formula_om, data = data_sample)
  
  # Get outcome predictions for combined data
  combined_data$A <- as.numeric(as.character(combined_data$A))
  predicted_outcomes_all <- predict(model_om, newdata = combined_data)
  
  # Get outcome predictions for unlabeled data only
  n_sample <- nrow(data_sample)
  n_combined <- nrow(combined_data)
  predicted_outcomes_unlabeled <- predicted_outcomes_all[(n_sample + 1):n_combined]
  
  # Return propensity scores and predicted outcomes
  return(list(
    propensity_scores = propensity_scores,
    predicted_outcomes_all = predicted_outcomes_all,
    predicted_outcomes_unlabeled = predicted_outcomes_unlabeled
  ))
}

# MAIPW estimator per group
calculate_MAIPW_per_group <- function(data_sample, data_out_of_sample, 
                                      propensity_scores, predicted_outcomes_unlabelled) {
  
  # Split propensity scores for labeled and unlabeled data
  n_sample <- nrow(data_sample)
  n_out <- nrow(data_out_of_sample)
  ps_sample <- propensity_scores[1:n_sample]
  ps_out <- propensity_scores[(n_sample + 1):(n_sample + n_out)]
  
  # Initialize a vector to store ATE per group
  groups <- unique(data_sample$group)
  ate_maipw_per_group <- numeric(length(groups))
  names(ate_maipw_per_group) <- groups
  
  for (group in groups) {
    # Filter data for each group
    sample_group <- data_sample[data_sample$group == group, ]
    out_group <- data_out_of_sample[data_out_of_sample$group == group, ]
    
    # Get propensity scores for this group
    ps_sample_group <- ps_sample[data_sample$group == group]
    ps_out_group <- ps_out[data_out_of_sample$group == group]
    
    # First term (treated group) for labeled and unlabeled data
    num1_labeled <- sum(sample_group$A * sample_group$y / ps_sample_group)
    num1_unlabeled <- sum(out_group$A * predicted_outcomes_unlabelled[data_out_of_sample$group == group] / ps_out_group)
    denom1 <- sum(c(sample_group$A / ps_sample_group, out_group$A / ps_out_group))
    
    # Second term (control group) for labeled and unlabeled data
    num2_labeled <- sum((1 - sample_group$A) * sample_group$y / (1 - ps_sample_group))
    num2_unlabeled <- sum((1 - out_group$A) * predicted_outcomes_unlabelled[data_out_of_sample$group == group] / (1 - ps_out_group))
    denom2 <- sum(c((1 - sample_group$A) / (1 - ps_sample_group), (1 - out_group$A) / (1 - ps_out_group)))
    
    # Calculate MAIPW estimator for this group
    ate_maipw_per_group[group] <- (num1_labeled + num1_unlabeled) / denom1 - 
      (num2_labeled + num2_unlabeled) / denom2
  }
  
  return(ate_maipw_per_group)
}


# Doubly Robust Estimator
calculate_DR_per_group <- function(data_sample, data_out_of_sample, 
                                             propensity_scores, predicted_outcomes_all) {
  
  # Split propensity scores and predicted outcomes
  n_sample <- nrow(data_sample)
  ps_sample <- propensity_scores[1:n_sample]
  
  # Initialize a vector to store ATE per group
  groups <- unique(data_sample$group)
  ate_dr_per_group <- numeric(length(groups))
  names(ate_dr_per_group) <- groups
  
  for (group in groups) {
    # Filter data for each group
    sample_group <- data_sample[data_sample$group == group, ]
    out_group <- data_out_of_sample[data_out_of_sample$group == group, ]
    
    # Get propensity scores and predicted outcomes for this group
    ps_sample_group <- ps_sample[data_sample$group == group]
    predicted_outcomes_sample_group <- predicted_outcomes_all[1:n_sample][data_sample$group == group]
    predicted_outcomes_all_group <- c(predicted_outcomes_sample_group, 
                                      predicted_outcomes_all[(n_sample + 1):length(predicted_outcomes_all)][data_out_of_sample$group == group])
    
    # First term: average of outcome model predictions for both labeled and unlabeled data
    first_term <- mean(predicted_outcomes_all_group)
    
    # Second term: adjustment using labeled data only
    second_term <- mean(
      (sample_group$A / ps_sample_group) * 
        (sample_group$y - predicted_outcomes_sample_group) -
        ((1 - sample_group$A) / (1 - ps_sample_group)) * 
        (sample_group$y - predicted_outcomes_sample_group)
    )
    
    # Calculate DR estimator for this group
    ate_dr_per_group[group] <- first_term + second_term
  }
  
  return(ate_dr_per_group)
}



# Function to calculate the true ATE per group using plasmode simulation
calculate_true_ATE_per_group <- function(data_sample, data_out_of_sample) {
  
  # Set a random seed for reproducibility
  set.seed(sample.int(1e9, 1))
  
  # Ensure 'group' is a factor with consistent levels across datasets
  all_groups <- union(levels(data_sample$group), levels(data_out_of_sample$group))
  data_sample$group <- factor(data_sample$group, levels = all_groups)
  data_out_of_sample$group <- factor(data_out_of_sample$group, levels = all_groups)
  
  # Combine labeled and unlabeled datasets
  combined_data <- rbind(
    data_sample[, setdiff(names(data_sample), "y")],
    data_out_of_sample
  )
  
  # Define covariate names
  covariate_names <- paste0("X", 1:11)
  
  # Outcome model estimation using linear regression on labeled data
  # Including 'group' as a fixed effect
  formula_om <- as.formula(paste("y ~ A +", paste(covariate_names, collapse = " + "), "+ group"))
  
  # Fit linear regression model for outcomes
  model_om <- tryCatch(
    lm(formula_om, data = data_sample),
    error = function(e) { 
      stop("Outcome model failed to converge: ", e$message) 
    }
  )
  
  # Extract outcome model coefficients
  alpha_hat <- coef(model_om)
  
  # Extract residual standard error (sigma) from the outcome model
  sigma_squared <- summary(model_om)$sigma^2
  
  # Ensure that the coefficients include the intercept and treatment effect
  if (!all(c("(Intercept)", "A") %in% names(alpha_hat))) {
    stop("Outcome model must include intercept and treatment indicator 'A'.")
  }
  
  # Extract intercept and treatment effect
  alpha0 <- alpha_hat["(Intercept)"]
  alpha1 <- alpha_hat["A"]
  
  # Extract covariate coefficients (excluding intercept, treatment, and group effects)
  alpha_covariates <- alpha_hat[covariate_names]
  
  # Extract group coefficients (excluding the reference group absorbed in the intercept)
  group_levels <- levels(combined_data$group)
  group_coefficients <- alpha_hat[grep("^group", names(alpha_hat))]
  
  # Create a mapping from group levels to coefficients
  group_coefs <- numeric(length(group_levels))
  names(group_coefs) <- group_levels
  
  # Assign coefficients to groups (reference group coefficient is zero)
  group_coefs[names(group_coefficients)] <- group_coefficients
  group_coefs[is.na(group_coefs)] <- 0  # Reference group
  
  # Extract covariate matrix
  Xi <- as.matrix(combined_data[, covariate_names])
  
  # Simulate error terms for potential outcomes
  epsilon_0 <- rnorm(n = nrow(combined_data), mean = 0, sd = sqrt(sigma_squared))
  epsilon_1 <- rnorm(n = nrow(combined_data), mean = 0, sd = sqrt(sigma_squared))
  
  # Get group effect for each individual
  group_effects <- group_coefs[as.character(combined_data$group)]
  
  # Simulate potential outcomes under control (Y0) and treatment (Y1)
  Y0 <- alpha0 + Xi %*% alpha_covariates + group_effects + epsilon_0
  Y1 <- alpha0 + alpha1 + Xi %*% alpha_covariates + group_effects + epsilon_1
  
  # Create a data frame to hold results
  sim_data <- data.frame(
    group = combined_data$group,
    Y1 = as.numeric(Y1),
    Y0 = as.numeric(Y0)
  )
  
  # Compute the true ATE per group
  ATE_per_group <- sim_data %>%
    group_by(group) %>%
    summarise(ATETrue = mean(Y1 - Y0))
  
  # Return the true ATE per group
  return(ATE_per_group)
}
# Function to perform multiple true ATE simulations per district
simulate_true_ATE_per_group <- function(data_sample, data_out_of_sample, nSim = 1000) {
  
  # Run the first simulation to get district (group) identifiers
  initial_result <- calculate_true_ATE_per_group(data_sample, data_out_of_sample)
  group_ids <- initial_result$group
  
  # Preallocate a matrix to store ATE values for each district and simulation
  ATE_values <- matrix(NA, nrow = nSim, ncol = length(group_ids))
  colnames(ATE_values) <- group_ids
  
  # Run simulations
  for (i in 1:nSim) {
    # Calculate true ATE per district for this simulation
    sim_result <- calculate_true_ATE_per_group(data_sample, data_out_of_sample)
    
    # Store the ATE values in the matrix for each group
    ATE_values[i, ] <- sim_result$ATETrue
    
    # Optional: Print progress every 100 simulations
    if (i %% 100 == 0) {
      cat("Simulation", i, "completed.\n")
    }
  }
  
  # Calculate cumulative mean and standard deviation of ATE per district
  mean_ATE_per_group <- colMeans(ATE_values, na.rm = TRUE)
  sd_ATE_per_group <- apply(ATE_values, 2, sd, na.rm = TRUE)
  
  # Return results as a list with ATE values, mean ATE, and standard deviation per group
  return(list(
    mean_ATE_per_group = mean_ATE_per_group,
    sd_ATE_per_group = sd_ATE_per_group
  ))
}

# Modified function to save the plot
plot_simulation_results <- function(simulation_results, filename = "simulation_results.png") {
  library(ggplot2)
  
  # Prepare data frame from simulation_results
  df <- data.frame(
    group = names(simulation_results$mean_ATE_per_group),
    mean_ATE = simulation_results$mean_ATE_per_group,
    sd_ATE = simulation_results$sd_ATE_per_group
  )
  
  # Ensure 'group' is treated as a factor and levels are in order
  df$group <- factor(df$group, levels = df$group)
  
  # Create the plot
  p <- ggplot(df, aes(x = group, y = mean_ATE)) +
    geom_errorbar(aes(ymin = mean_ATE - sd_ATE, ymax = mean_ATE + sd_ATE), 
                  width = 0.2, color = "red") +
    geom_point(color = "black") + 
    xlab("District") +
    ylab("Average Treatment Effect (ATE)") +
    ylim(1, 1.5) + # Set y-axis limits from 0 to 1.5 
    ggtitle("Mean ATE per District with Standard Deviation") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  # Save the plot as a .png file
  ggsave(filename, plot = p, width = 10, height = 6, dpi = 300)
  
  # Optionally, display the plot
  print(p)
}


simulation_results <- simulate_true_ATE_per_group(data_sample, data_out_of_sample, nSim = 1000)
plot_simulation_results(simulation_results)

# Step 1: Estimate nuisance parameters (propensity scores and outcome model predictions)
nuisance_params <- estimate_nuisance_params_logistic_linear(data_sample, data_out_of_sample)

# Extract propensity scores and predicted outcomes from nuisance parameters
propensity_scores <- nuisance_params$propensity_scores
predicted_outcomes_unlabeled <- nuisance_params$predicted_outcomes_unlabeled
predicted_outcomes_all <- nuisance_params$predicted_outcomes_all

# Step 2: Calculate ATE per group using the MAIPW estimator
ate_maipw_per_group <- calculate_MAIPW_per_group(
  data_sample = data_sample,
  data_out_of_sample = data_out_of_sample,
  propensity_scores = propensity_scores,
  predicted_outcomes_unlabelled = predicted_outcomes_unlabeled
)

# Step 3: Calculate ATE per group using the DR estimator
ate_dr_per_group <- calculate_DR_per_group(
  data_sample = data_sample,
  data_out_of_sample = data_out_of_sample,
  propensity_scores = propensity_scores,
  predicted_outcomes_all = predicted_outcomes_all
)

# Step 1: Estimate nuisance parameters using Random Forest
nuisance_params_rf <- estimate_nuisance_params_random_forest(data_sample, data_out_of_sample)

# Extract propensity scores and predicted outcomes from random forest nuisance parameters
propensity_scores_rf <- nuisance_params_rf$propensity_scores
predicted_outcomes_unlabeled_rf <- nuisance_params_rf$predicted_outcomes_unlabeled
predicted_outcomes_all_rf <- nuisance_params_rf$predicted_outcomes_all

# Step 2: Calculate ATE per group using the MAIPW estimator with Random Forest nuisance parameters
ate_maipw_per_group_rf <- calculate_MAIPW_per_group(
  data_sample = data_sample,
  data_out_of_sample = data_out_of_sample,
  propensity_scores = propensity_scores_rf,
  predicted_outcomes_unlabelled = predicted_outcomes_unlabeled_rf
)

# Step 3: Calculate ATE per group using the DR estimator with Random Forest nuisance parameters
ate_dr_per_group_rf <- calculate_DR_per_group(
  data_sample = data_sample,
  data_out_of_sample = data_out_of_sample,
  propensity_scores = propensity_scores_rf,
  predicted_outcomes_all = predicted_outcomes_all_rf
)

# Function to calculate Mean and Median MSE for each estimator
calculate_mse_summary <- function(ate_dr_per_group, ate_maipw_per_group, simulation_results) {

  # Extract the true ATE values from simulation results
  true_ate_values <- simulation_results$mean_ATE_per_group

  # Calculate MSE for each district for the DR estimator
  mse_dr <- sapply(names(ate_dr_per_group), function(group) {
    true_ate <- true_ate_values[group]
    estimated_ate <- ate_dr_per_group[group]
    mean((estimated_ate - true_ate)^2)
  })

  # Calculate MSE for each district for the MAIPW estimator
  mse_maipw <- sapply(names(ate_maipw_per_group), function(group) {
    true_ate <- true_ate_values[group]
    estimated_ate <- ate_maipw_per_group[group]
    mean((estimated_ate - true_ate)^2)
  })

  # Compute mean and median MSE for each estimator
  mean_mse_dr <- mean(mse_dr)
  median_mse_dr <- median(mse_dr)
  mean_mse_maipw <- mean(mse_maipw)
  median_mse_maipw <- median(mse_maipw)

  # Return results as two lists: one for DR and one for MAIPW
  return(list(
    DR = list(mean_mse = mean_mse_dr, median_mse = median_mse_dr),
    MAIPW = list(mean_mse = mean_mse_maipw, median_mse = median_mse_maipw)
  ))
}

# Usage with logistic/linear regression nuisance parameters
mse_summary_logistic <- calculate_mse_summary(
  ate_dr_per_group = ate_dr_per_group,
  ate_maipw_per_group = ate_maipw_per_group,
  simulation_results = simulation_results
)

# Usage with random forest nuisance parameters
mse_summary_rf <- calculate_mse_summary(
  ate_dr_per_group = ate_dr_per_group_rf,
  ate_maipw_per_group = ate_maipw_per_group_rf,
  simulation_results = simulation_results
)

# Print results
cat("Logistic/Linear Regression Nuisance Parameters - DR Estimator MSE (Mean, Median): ",
    mse_summary_logistic$DR$mean_mse, mse_summary_logistic$DR$median_mse, "\n")
cat("Logistic/Linear Regression Nuisance Parameters - MAIPW Estimator MSE (Mean, Median): ",
    mse_summary_logistic$MAIPW$mean_mse, mse_summary_logistic$MAIPW$median_mse, "\n")

cat("Random Forest Nuisance Parameters - DR Estimator MSE (Mean, Median): ",
    mse_summary_rf$DR$mean_mse, mse_summary_rf$DR$median_mse, "\n")
cat("Random Forest Nuisance Parameters - MAIPW Estimator MSE (Mean, Median): ",
    mse_summary_rf$MAIPW$mean_mse, mse_summary_rf$MAIPW$median_mse, "\n")
