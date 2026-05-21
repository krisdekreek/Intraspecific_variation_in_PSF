    ### This is an R script with functions that can be sourced into your R script.
    ### This can reduce the lines of code in your script
    ### Use the function source() to access this script. Run this line before using the pre-made functions
    ### For example: source("C:/Users/kreek001/OneDrive - Wageningen University & Research/R_Scripts")
    ### Autor: Kris de Kreek


### Function to check model assumptions with the DHARMa package
library(DHARMa)
DHARMa.sum <- function(Model) {
  simulationOutput <- simulateResiduals(Model, n = 1000, seed = 256)
  par(mfrow = c(2, 2))
  cat("1. Kolmogorov-Smirnov test indicates deviations of the residuals from uniform distribution. H0: The model fits the data well; Ha: The model does not fit the data well")
  print(testUniformity(simulationOutput, plot = TRUE))
  plotResiduals(simulationOutput)               # Plot residuals
  cat("2. Test for over-dispersion. H0: No over-dispersion in the model; Ha: Over-dispersion in the model")
  print(testDispersion(simulationOutput, plot = TRUE)) # Tests over-dispersion: is observed data more/less dispersed than expected under fitted model
  cat("3. Test for outliers. H0: No ouliers in the model; Ha: Outliers in the model")
  print(testOutliers(simulationOutput, plot = TRUE))   # Tests for outliers
  par(mfrow = c(1, 1))
}

### Made when last three plots do give errors.
DHARMa.sum2 <- function(Model) {
  simulationOutput <- simulateResiduals(Model, n = 1000, seed = 256)
  #par(mfrow = c(2, 1))
  cat("1. Kolmogorov-Smirnov test indicates deviations of the residuals from uniform distribution. H0: The model fits the data well; Ha: The model does not fit the data well")
  print(testUniformity(simulationOutput, plot = TRUE))
  #plotResiduals(simulationOutput)               # Plot residuals
  #cat("2. Test for over-dispersion. H0: No over-dispersion in the model; Ha: Over-dispersion in the model")
  #print(testDispersion(simulationOutput, plot = TRUE)) # Tests over-dispersion: is observed data more/less dispersed than expected under fitted model
  #cat("3. Test for outliers. H0: No ouliers in the model; Ha: Outliers in the model")
  #print(testOutliers(simulationOutput, plot = TRUE))   # Tests for outliers
  #par(mfrow = c(1, 1))
}

#### This is specifically when 5 plots are in the output. the last plot is not very informative an is left out. This happens in the case of binomial and beta-binomial distributions
DHARMa.sum.3 <- function(Model) {
  simulationOutput <- simulateResiduals(Model, n = 1000, seed = 256)
  par(mfrow = c(1, 3))
  cat("1. Kolmogorov-Smirnov test indicates deviations of the residuals from uniform distribution. H0: The model fits the data well; Ha: The model does not fit the data well")
  print(testUniformity(simulationOutput, plot = TRUE))
  plotResiduals(simulationOutput)               # Plot residuals
  cat("2. Test for over-dispersion. H0: No over-dispersion in the model; Ha: Over-dispersion in the model")
  print(testDispersion(simulationOutput, plot = TRUE)) # Tests over-dispersion: is observed data more/less dispersed than expected under fitted model
  cat("3. Test for outliers. H0: No ouliers in the model; Ha: Outliers in the model")
  print(testOutliers(simulationOutput, plot = FALSE))   # Tests for outliers
  par(mfrow = c(1, 1))
}

#### Levene test with DHARMA
# # DHARMa object
# sim_res <- simulateResiduals(glmg1)
# 
# # Levene-type test on DHARMa residuals
# leveneTest(sim_res$scaledResiduals ~ data_RI$Soil_conditioning)