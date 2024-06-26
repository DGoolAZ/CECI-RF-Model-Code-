# Load the necessary libraries
library(terra)
library(sf)
library(data.table)
library(randomForest)
library(viridis)
library(ggplot2)
library(vip)
library(pROC)
library(pdp)

# STEP 1: Load in Data ----------------------------------------------------------

# Define a function to impute NA values with the mean for a raster layer
impute_na_with_mean <- function(r) {
  vals <- values(r)
  na_mean <- mean(vals, na.rm = TRUE)
  vals[is.na(vals)] <- na_mean
  values(r) <- vals
  return(r)
}

# Load presence data
presence_data <- fread("F:\\Modeling\\CECI Update\\Presance\\Pima_Presaince.csv")

# Adjust column names based on actual data structure
presence_data_clean <- presence_data[, .(latitude = get("latitude"), longitude = get("longitude"))]
presence_data_clean <- na.omit(presence_data_clean)

# Convert to sf object
presence_data_sf <- st_as_sf(presence_data_clean, coords = c("longitude", "latitude"), crs = 4326)

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

###################################################
#               STEP 1 Explanation                #
###################################################
# In Step 1, we prepare our data for modeling. We start by loading presence data from a CSV file and adjusting column names to match the structure of the data. 
# Next, we convert this data into a spatial points dataframe (sf object) using longitude and latitude coordinates. 
# We then read in raster data stored in a specified directory, ensuring all rasters align with a template raster (DEM.tif) by resampling and imputing missing values using a bilinear method. 
# Finally, we extract values from these rasters corresponding to the presence data locations, preparing them for subsequent modeling steps.

# STEP 2: GENERATE BACKGROUND DATA ----------------------------------------------

# Generate background data
set.seed(17172)
background_points <- spatSample(raster_stack, size = nrow(presence_data_clean), xy = TRUE, as.df = TRUE)
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

# Debugging: Print summary of combined data
cat("Summary of combined data:\n")
print(summary(combined_data_clean))
cat("\n")

#########################################
#        STEP 2 Explanation             #
#########################################
# In Step 2, we create background data for our modeling process. Using spatial sampling (`spatSample`), we generate random points across our raster stack, ensuring they align with our presence data's spatial extent. 
# These points are converted into a spatial points dataframe (`sf object`) and values are extracted from our raster stack at these locations. 
# After combining these background values with our presence data, we ensure they share the same columns and clean the combined dataset by removing any remaining missing values, preparing it for model training.

# STEP 3: TRAIN RF MODEL --------------------------------------------------------

# Define response and predictors
response_var <- "response"
predictors <- setdiff(names(combined_data_clean), response_var)

# Train the Random Forest model
rf_formula <- as.formula(paste(response_var, "~", paste(predictors, collapse = " + ")))

rf_model <- randomForest(rf_formula,
                         data = combined_data_clean,
                         ntree = 1000,
                         mtry = 2, 
                         importance = TRUE,
                         type = "prob")

# Debugging: Print the random forest model
print(rf_model)

#######################################
#          STEP 3 Explanation         #
#######################################
# In Step 3, we train a Random Forest (RF) model using our prepared data. We define the response variable (`response`) and predictors by removing the response variable from the dataset columns. 
# The RF model is then built with 1000 trees (`ntree = 1000`) and 2 variables (`mtry = 2`) considered at each split, utilizing the entire dataset (`type = "prob"` for probability predictions). 
# This model will be used to predict the probability of presence across our raster stack based on the environmental predictors identified.

# STEP4: PREDICT OVER RASTER ----------------------------------------------------
# Predict raster using the trained model
predictions <- terra::predict(raster_stack, rf_model, type = "prob", index = 2, filename="F:\\Modeling\\CECI Update\\RF_ModelR1.tif", na.rm=TRUE, progress="text", overwrite=TRUE)

# Load the predictions as a raster
predicted_raster <- rast("F:\\Modeling\\CECI Update\\RF_ModelR1.tif")

# Debugging: Print summary of predicted raster
cat("Summary of predicted raster:\n")
print(summary(predicted_raster))
cat("\n")

# Plot the predicted raster
plot(predicted_raster, main="Predicted Probability of Presence", col=viridis::viridis(100))

#################################
#         STEP 4 Explanation    #
#################################
# In Step 4, we apply our trained RF model to predict the probability of species presence across our entire raster stack. 
# The model generates predictions and saves them as a raster file (`RF_ModelR1.tif`). 
# After loading these predictions into R as a raster object, we summarize its properties and visualize the predicted probabilities using a color gradient (`viridis::viridis(100)`).

# STEP 5 ANALYSIS ---------------------------------------------------------------

predict_function <- function(object, newdata) {
  predict(object, newdata, type = "prob")[, 2]
}

# Predict probabilities on the training data
train_predictions <- predict(rf_model, combined_data_clean, type = "prob")[, 2]

# Create a variable importance plot
importance_df <- data.frame(Variable = rownames(rf_model$importance),
                            Importance = rf_model$importance[, "MeanDecreaseGini"])

# Plot using ggplot2
ggplot(importance_df, aes(x = reorder(Variable, Importance), y = Importance)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  xlab("Variables") +
  ylab("Importance") +
  ggtitle("Variable Importance Plot") +
  theme_minimal()

# Define a custom function to calculate AUC for the RF model
calculate_auc <- function(actual, predicted) {
  roc_obj <- roc(actual, predicted)
  auc(roc_obj)
}

# Calculate AUC for the RF model
auc_value <- calculate_auc(combined_data_clean$response, train_predictions)
cat("AUC for the RF model:", auc_value, "\n")

#################################
#         STEP 5 Explanation    #
#################################
# In Step 5, we analyze the performance and insights from our trained Random Forest (RF) model. 
# We first define a prediction function to compute probabilities using our model. 
# Then, we predict probabilities on the training data (`combined_data_clean`) and visualize the variable importance using a bar plot with ggplot2. 
# Additionally, we calculate the Area Under the Curve (AUC) to assess the model's predictive accuracy, providing a quantitative measure (`auc_value`) of its performance.