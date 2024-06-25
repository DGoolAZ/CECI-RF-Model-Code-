Species Distribution Modeling and Analysis
This repository contains R code for performing Species Distribution Modeling (SDM) using a Random Forest model. The script covers data preprocessing, model training, prediction, and visualization, including the creation of Partial Dependence Plots (PDPs) and ROC curves.

Prerequisites
Make sure you have the following R packages installed:

r
install.packages(c("terra", "sf", "data.table", "randomForest", "viridis", "ggplot2", "vip", "pROC", "yardstick", "pdp"))
Description of the Code


1. Load Required Libraries
The script begins by loading the necessary libraries for geospatial data handling, data manipulation, modeling, and visualization.

2. Define Helper Functions
A helper function impute_na_with_mean is defined to impute missing values in raster layers with the mean value.

3. Load and Preprocess Data
Presence Data: Presence data is loaded from a CSV file and converted into an sf object.
Raster Data: Environmental raster data is loaded, aligned, and missing values are imputed. The rasters are then combined into a stack.

4. Generate Background Data
Background data points are generated to provide a contrast to the presence data.

5. Combine Data
Presence and background data are combined into a single dataset for modeling. The response variable is set to 1 for presence data and 0 for background data.

6. Train the Random Forest Model
A Random Forest model is trained using the combined dataset. Predictor variables are extracted, and the model is fitted.

7. Predict on Raster Data
The trained model is used to predict species presence probabilities across the landscape. The predictions are saved as a raster file.

8. Variable Importance and Partial Dependence Plots
Permutation Importance: The script calculates permutation importance for the predictor variables and plots the results.
Partial Dependence Plots: PDPs are generated for each predictor variable and displayed one at a time, prompting the user to press Enter to proceed to the next plot.

9. ROC Curve and AUC
An ROC curve is plotted for the model, and the AUC value is displayed on the plot.

