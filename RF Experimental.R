# Load the necessary libraries
library(terra)
library(sf)
library(data.table)
library(randomForest)
library(viridis)
library(ggplot2)
library(vip)
library(pROC)
library(yardstick)
library(pdp)
library(caret)
library(spThin)
library(sp)
library(terra)
library(sf)
library(data.table)
library(sp)


# Define a function to manually thin the points based on a minimum distance
thin_points_sf <- function(points, min_dist) {
  points_sf <- st_as_sf(points, coords = c("longitude", "latitude"), crs = 4326)
  thinned_sf <- points_sf[1, ]  # Start with the first point
  
  for (i in 2:nrow(points_sf)) {
    point <- points_sf[i, ]
    distances <- st_distance(thinned_sf, point)
    if (all(as.numeric(distances) >= min_dist * 1000)) {  # Convert km to meters
      thinned_sf <- rbind(thinned_sf, point)
    }
  }
  
  return(thinned_sf)
}

# Load presence data
presence_data <- fread("F:\\Modeling\\CECI Update\\Presance\\Pima_Presaince.csv")

# Adjust column names based on actual data structure
presence_data_clean <- presence_data[, .(latitude = get("latitude"), longitude = get("longitude"))]
presence_data_clean <- na.omit(presence_data_clean)

# Debugging: Check if presence_data_clean is empty
if (nrow(presence_data_clean) == 0) {
  stop("Presence data is empty after cleaning.")
}

# Define minimum distance in kilometers
min_dist_km <- 4  # Example: 4 km

# Apply the custom thinning function
presence_data_thinned_sf <- thin_points_sf(presence_data_clean, min_dist_km)

# Convert the thinned sf object back to a data table
presence_data_thinned <- as.data.table(st_coordinates(presence_data_thinned_sf))
setnames(presence_data_thinned, c("X", "Y"), c("longitude", "latitude"))

# Visualize presence data points after thinning
ggplot(presence_data_thinned, aes(x = longitude, y = latitude)) +
  geom_point(color = "red") +
  ggtitle("Presence Data Points After Thinning") +
  xlab("Longitude") +
  ylab("Latitude")


# Convert to sf object
presence_data_sf <- st_as_sf(presence_data_thinned, coords = c("longitude", "latitude"), crs = 4326)

# Read in raster data
raster_directory <- "F:\\Modeling\\CECI Update\\ENV Rast Pima"
dem_file <- file.path(raster_directory, "DEM.tif")

# Load the DEM raster as a template for resampling
template_raster <- rast(dem_file)

# Load other rasters, align them, and impute NAs
raster_files <- list.files(path = raster_directory, pattern = "\\.tif$", full.names = TRUE)
raster_files <- raster_files[!grepl("DEM.tif|VegType.tif", raster_files)]
aligned_rasters <- list()
for (file in raster_files) {
  r <- rast(file)
  r <- if (!compareGeom(r, template_raster, stopOnError = FALSE)) project(r, template_raster) else r
  r <- resample(r, template_raster, method = "bilinear")
  r <- impute_na_with_mean(r)
  raster_name <- tools::file_path_sans_ext(basename(file))
  aligned_rasters[[raster_name]] <- r
  
  # Debugging: Print summary of each raster
  cat("Summary of raster", raster_name, ":\n")
  print(summary(r))
  cat("\n")
}

# Add DEM to the list, also impute NAs
aligned_rasters[["DEM"]] <- impute_na_with_mean(template_raster)

# Combine aligned and imputed rasters into a stack
raster_stack <- rast(aligned_rasters)
names(raster_stack) <- make.names(names(raster_stack))

# Ensure all raster layers have no NA values
for (i in seq_along(raster_stack)) {
  if (sum(is.na(values(raster_stack[[i]]))) > 0) {
    raster_stack[[i]] <- impute_na_with_mean(raster_stack[[i]])
  }
  
  # Debugging: Check for NAs after imputation
  cat("Number of NAs in raster", names(raster_stack)[i], "after imputation:", sum(is.na(values(raster_stack[[i]]))), "\n")
}

# Extract raster values for presence data
presence_values <- terra::extract(raster_stack, presence_data_sf)
presence_values <- presence_values[, -1]  # Remove ID column
presence_df <- as.data.table(presence_data_sf)
presence_df <- cbind(presence_df, presence_values)
presence_df$response <- 1

# Debugging: Print summary of presence data
cat("Summary of presence data:\n")
print(summary(presence_df))
cat("\n")

# Generate background data
set.seed(17172)
background_points <- spatSample(raster_stack, size = nrow(presence_data_thinned) * 10, xy = TRUE, as.df = TRUE)  # Generate more background points initially
background_sf <- st_as_sf(background_points, coords = c("x", "y"), crs = 4326)
background_values <- terra::extract(raster_stack, background_sf)
background_values <- background_values[, -1]
background_df <- as.data.table(background_sf)
background_df <- cbind(background_df, background_values)
background_df$response <- 0

# Debugging: Print summary of background data
cat("Summary of background data:\n")
print(summary(background_df))
cat("\n")

# Ensure presence_df and background_df have the same columns
common_cols <- intersect(names(presence_df), names(background_df))
presence_df <- presence_df[, ..common_cols]
background_df <- background_df[, ..common_cols]

# Combine presence and background data
combined_data <- rbind(presence_df, background_df, fill = TRUE)

# Remove geometry column if it exists
combined_data <- combined_data[, !"geometry", with = FALSE]

# Remove rows with any remaining NA values
combined_data_clean <- na.omit(combined_data)
combined_data_clean$response <- as.factor(combined_data_clean$response)

# Rename levels of the response variable to valid R variable names
levels(combined_data_clean$response) <- make.names(levels(combined_data_clean$response))

# Debugging: Print summary of combined data
cat("Summary of combined data:\n")
print(summary(combined_data_clean))
cat("\n")

# Define response and predictors
response_var <- "response"
predictors <- setdiff(names(combined_data_clean), response_var)

# Train the Random Forest model
rf_formula <- as.formula(paste(response_var, "~", paste(predictors, collapse = " + ")))

# Set up cross-validation
control <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary)

# Train the Random Forest model with cross-validation
rf_model_cv <- train(
  rf_formula,
  data = combined_data_clean,
  method = "rf",
  trControl = control,
  metric = "ROC"
)

# Print cross-validated model results
print(rf_model_cv)

# Extract and plot cross-validated AUC
cv_results <- rf_model_cv$results
cv_auc <- cv_results[which.max(cv_results$ROC), "ROC"]

# Predict raster using the trained model
predictions <- terra::predict(raster_stack, rf_model_cv$finalModel, type = "prob", index = 2, filename="F:\\Modeling\\CECI Update\\RF_ModelR1.tif", na.rm=TRUE, progress="text", overwrite=TRUE)
predicted_raster <- rast("F:\\Modeling\\CECI Update\\RF_ModelR1.tif")
plot(predicted_raster, main="Predicted Probability of Presence", col=viridis::viridis(100))

# Define a custom prediction function for permutation importance
predict_function <- function(object, newdata) {
  predict(object, newdata, type = "prob")[, 2]
}

# Calculate permutation importance with AUC as the metric
perm_importance <- vip::vi(
  object = rf_model_cv$finalModel,
  method = "permute",
  target = "response",
  train = combined_data_clean,
  metric = yardstick::roc_auc_vec,
  pred_wrapper = predict_function,
  nsim = 50,  # Number of permutations
  sample_frac = 1,  # Use the whole dataset for each permutation
  smaller_is_better = FALSE
)

# Plot permutation importance using ggplot2
vip::vip(perm_importance, num_features = length(predictors), geom = "col") +
  theme_minimal() +
  ggtitle("Permutation Variable Importance Plot") +
  xlab("Variables") +
  ylab("Importance")

# Create PDPs for all predictor variables
for (var in predictors) {
  pd <- partial(rf_model_cv$finalModel, pred.var = var, train = combined_data_clean, prob = TRUE)
  p <- autoplot(pd) +
    ggtitle(paste("Partial Dependence Plot for", var)) +
    xlab(var) +
    ylab("Predicted Probability of Presence") +
    theme_minimal()
  print(p)
  cat("Press Enter to continue to the next plot...")
  readline()
}

# Extract the predictions and plot ROC curves for each fold
predictions <- rf_model_cv$pred
folds <- unique(predictions$Resample)

# Initialize plot
plot(0, 0, type = "n", xlab = "1 - Specificity", ylab = "Sensitivity", xlim = c(0, 1), ylim = c(0, 1), main = "ROC Curves for Each Fold")

# Plot ROC curve for each fold
for (fold in folds) {
  fold_predictions <- predictions[predictions$Resample == fold, ]
  roc_obj <- roc(fold_predictions$obs, fold_predictions$yes)
  plot(roc_obj, col = "blue", add = TRUE)
}

# Add a legend with the fold names
legend("bottomright", legend = paste("Fold", folds), col = "blue", lty = 1)

# Extract and plot cross-validated AUC
cv_results <- rf_model_cv$results
cv_auc <- cv_results[which.max(cv_results$ROC), "ROC"]

# Generate the overall ROC curve
train_predictions <- predict(rf_model_cv, combined_data_clean, type = "prob")[, 2]
roc_obj <- roc(combined_data_clean$response, train_predictions)
plot(roc_obj, main = "Cross-Validated ROC Curve for Random Forest Model", col = "blue")
mtext(paste("AUC:", round(cv_auc, 3)), side = 1, line = 2.5, at = 0.5, col = "blue")