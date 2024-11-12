# Install and load necessary packages
library(ggplot2)
library(cowplot)

n_values <- rep(c(100, 200, 500, 1000), each = 4)
Semi_Parametric_Method <- rep(c('Doubly Robust', 'Doubly Robust', 'MAIPW', 'MAIPW'), times = 4)
Nuisance_Parameter_Estimator <- rep(c('Logistic/Linear Regression', 'Random Forest', 'Logistic/Linear Regression', 'Random Forest'), times = 4)

Mean_MSE <- c(
  4.075102, 3.508268, 0.000765701, 0.02001351,   # n = 100
  4.078445, 3.50001, 0.000843975, 0.02016137,     # n = 200
  4.076975, 3.491045, 0.000797057, 0.02200245,    # n = 500
  4.07806, 3.514266, 0.000797057, 0.02212761      # n = 1000
)

Median_MSE <- c(
  3.119898, 3.350692, 0.000330206, 0.006280697,   # n = 100
  3.093296, 3.362853, 0.000420843, 0.006256257,   # n = 200
  3.093079, 3.387544, 0.000369306, 0.00819386,    # n = 500
  3.111192, 3.43054, 0.000369306, 0.007530356     # n = 1000
)

data <- data.frame(
  n = n_values,
  Semi_Parametric_Method = Semi_Parametric_Method,
  Nuisance_Parameter_Estimator = Nuisance_Parameter_Estimator,
  Mean_MSE = Mean_MSE,
  Median_MSE = Median_MSE
)

# Define colors for the bars
color_values <- c('Logistic/Linear Regression' = 'purple', 'Random Forest' = 'orange')

# Function to create Mean MSE plot for a given Semi-Parametric Method
plot_mean_mse <- function(data, method) {
  data_subset <- subset(data, Semi_Parametric_Method == method)
  
  p <- ggplot(data_subset, aes(x = factor(n), y = Mean_MSE, fill = Nuisance_Parameter_Estimator)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +  # Adjust bar width
    scale_fill_manual(values = color_values) +
    labs(title = paste("Mean MSE for", method, "Estimator"),
         x = "Number of Simulation Iterations (m)",
         y = "Mean MSE",
         fill = "Nuisance Parameter Estimator") +
    theme_minimal()
  
  return(p)
}

# Function to create Median MSE plot for a given Semi-Parametric Method
plot_median_mse <- function(data, method) {
  data_subset <- subset(data, Semi_Parametric_Method == method)
  
  p <- ggplot(data_subset, aes(x = factor(n), y = Median_MSE, fill = Nuisance_Parameter_Estimator)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +  # Adjust bar width
    scale_fill_manual(values = color_values) +
    labs(title = paste("Median MSE for", method, "Estimator"),
         x = "Number of Simulation Iterations (m)",
         y = "Median MSE",
         fill = "Nuisance Parameter Estimator") +
    theme_minimal()
  
  return(p)
}

# Generate plots for Doubly Robust estimator
p1 <- plot_mean_mse(data, 'Doubly Robust')
p2 <- plot_median_mse(data, 'Doubly Robust')

# Generate plots for MAIPW estimator
p3 <- plot_mean_mse(data, 'MAIPW')
p4 <- plot_median_mse(data, 'MAIPW')

# Arrange the plots into two separate 1x2 grids
grid_doubly_robust <- plot_grid(
  p1, p2,
  labels = c('A', 'B'),
  ncol = 2
)

grid_maipw <- plot_grid(
  p3, p4,
  labels = c('C', 'D'),
  ncol = 2
)

# Save the separate grids
ggsave("Doubly_Robust_Plots.png", grid_doubly_robust, width = 12, height = 5)
ggsave("MAIPW_Plots.png", grid_maipw, width = 12, height = 5)

# Display the grids
print(grid_doubly_robust)
print(grid_maipw)
