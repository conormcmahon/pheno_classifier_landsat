
library(tidyverse)
library(raster)
#library(sp)
library(rgdal)
library(here)
library(caret)
library(snow)
library(randomForest)
library(rlist)
library(rgeos)
# Reprojects, resamples, and crops a raster image to match another raster
source(here::here("aggregate_custom.R"))

directory <- "D:/SERDP/Vandenberg/good_data/indices/"

set.seed(7)

# ***** Load Data for Model ****

# Load Terrain Information
slope <- raster("D:/SERDP/Vandenberg/hydrology/slope.tif")
slope <- aggregate_custom(slope, ndvi_dry, reproj_method="bilinear")
dem <- raster("D:/SERDP/Vandenberg/hydrology/dem_resamp.tif")
dem <- aggregate_custom(dem, ndvi_dry, reproj_method="bilinear")
channel_dist <- raster("D:/SERDP/Vandenberg/hydrology/phenology_classifier/distance_to_channel_ln.tif")
channel_dist <- aggregate_custom(channel_dist, ndvi_dry, reproj_method="bilinear")
# Combine Phenology and Terrain information
#   Something seems to be wrong with some of the December phenology data... not sure why yet, just excluding them for now
# Load phenology rasters
load_phenoseries <- function(year)
{
  ndvi <- brick(paste(directory,"pheno_NDVI_",year,".tif",sep=""))
  wetness <- brick(paste(directory,"pheno_tc_wet_",year,".tif",sep=""))
  brightness <- brick(paste(directory,"pheno_tc_bright_",year,".tif",sep=""))
  phenoseries <- stack(ndvi[[1:11]],wetness[[1:11]],brightness[[1:11]])
  return(phenoseries)
}
phenoseries_dry <- stack(load_phenoseries(2015),slope,dem,channel_dist)
phenoseries_wet <- stack(load_phenoseries(2019),slope,dem,channel_dist)

# **** Load training polygons ****
segments <- readOGR(here::here("polygons","training_classes.shp"))
segments <- spTransform(segments, crs(phenoseries))
segments$area <- area(segments)/900
as.data.frame(segments) %>% group_by(class) %>% summarize(total_area = sum(area)) %>% arrange(-total_area)


generate_training_data <- function(phenoseries)
{
  
  # Generate Ground Truth Segment Groups
  # --- Riparian ---
  riparian_segments <- segments[segments$class=="Riparian",]
  riparian_df <- raster::extract(phenoseries, riparian_segments, cellnumbers=TRUE)
  riparian_df <- as.data.frame(list.rbind(riparian_df))
  riparian_df <- cbind(riparian_df, coordinates(phenoseries)[riparian_df[,1],])
  riparian_df$class <- "Riparian"
  # --- Water ---
  oak_segments <- segments[segments$class=="Oak",]
  oak_df <- raster::extract(phenoseries, oak_segments, cellnumbers=TRUE)
  oak_df <- as.data.frame(list.rbind(oak_df))
  oak_df <- cbind(oak_df, coordinates(phenoseries)[oak_df[,1],])
  oak_df$class <- "Oak"
  # --- Chaparral ---
  chaparral_segments <- segments[segments$class=="Chaparral",]
  chaparral_df <- raster::extract(phenoseries, chaparral_segments, cellnumbers=TRUE)
  chaparral_df <- as.data.frame(list.rbind(chaparral_df))
  chaparral_df <- cbind(chaparral_df, coordinates(phenoseries)[chaparral_df[,1],])
  chaparral_df$class <- "Chaparral"
  # --- Annual ---
  annual_segments <- segments[segments$class=="Annual",]
  annual_df <- raster::extract(phenoseries, annual_segments, cellnumbers=TRUE)
  annual_df <- as.data.frame(list.rbind(annual_df))
  annual_df <- cbind(annual_df, coordinates(phenoseries)[annual_df[,1],])
  annual_df$class <- "Annual"
  # --- EucCup ---
  euccup_segments <- segments[segments$class=="EucCup",]
  euccup_df <- raster::extract(phenoseries, euccup_segments, cellnumbers=TRUE)
  euccup_df <- as.data.frame(list.rbind(euccup_df))
  euccup_df <- cbind(euccup_df, coordinates(phenoseries)[euccup_df[,1],])
  euccup_df$class <- "Evergreen Nonnative"
  # --- Agriculture ---
  agriculture_segments <- segments[segments$class=="Agriculture",]
  agriculture_df <- raster::extract(phenoseries, agriculture_segments, cellnumbers=TRUE)
  agriculture_df <- as.data.frame(list.rbind(agriculture_df))
  agriculture_df <- cbind(agriculture_df, coordinates(phenoseries)[agriculture_df[,1],])
  agriculture_df$class <- "Agriculture"
  # --- Soil ---
  soil_segments <- segments[segments$class=="Soil",]
  soil_df <- raster::extract(phenoseries, soil_segments, cellnumbers=TRUE)
  soil_df <- as.data.frame(list.rbind(soil_df))
  soil_df <- cbind(soil_df, coordinates(phenoseries)[soil_df[,1],])
  soil_df$class <- "Soil"
  # --- Water ---
  water_segments <- segments[segments$class=="Water",]
  water_df <- raster::extract(phenoseries, water_segments, cellnumbers=TRUE)
  water_df <- as.data.frame(list.rbind(water_df))
  water_df <- cbind(water_df, coordinates(phenoseries)[water_df[,1],])
  water_df$class <- "Water"
  # --- Pavement ---
  urban_segments <- segments[segments$class=="Urban",]
  urban_df <- raster::extract(phenoseries, urban_segments, cellnumbers=TRUE)
  urban_df <- as.data.frame(list.rbind(urban_df))
  urban_df <- cbind(urban_df, coordinates(phenoseries)[urban_df[,1],])
  urban_df$class <- "Urban"
  # --- All Label Data ---
  label_df <- rbind(riparian_df, oak_df, chaparral_df, annual_df, euccup_df, agriculture_df, soil_df, water_df, urban_df)
  label_df <- label_df %>% 
    dplyr::select(2:ncol(label_df)) %>%
    drop_na()
  label_df$class <- factor(label_df$class, 
                           levels = c("Riparian", "Oak", "Chaparral", "Annual", "Evergreen Nonnative", "Agriculture", "Soil", "Water", "Urban"))
  names(label_df) <- c(as.character(1:33),"slope","dem","channel_dist","x","y","class")
  
  
  
  # **** Visualize Greenness Phenology **** 
  pheno_df_long <- label_df %>% pivot_longer((1:11), names_to="month", values_to="ndvi")
  pheno_df_long$month <- as.numeric(pheno_df_long$month)
  pheno_stats <- pheno_df_long %>% 
    group_by(class, month) %>%
    summarise(p_05=quantile(ndvi, probs=c(0.05)), 
              p_25=quantile(ndvi, probs=c(0.25)),
              p_50=quantile(ndvi, probs=c(0.50)),
              p_75=quantile(ndvi, probs=c(0.75)),
              p_95=quantile(ndvi, probs=c(0.95)))
  # Natural Vegetation Phenology
  veg_pheno_plot <- ggplot(data=pheno_stats[pheno_stats$class %in% factor(c("Riparian","Oak","Chaparral","Annual","Evergreen Nonnative","Agriculture")),], aes(x=month)) + 
    geom_line(aes(y=p_05), linetype="dashed", color="red") + 
    geom_line(aes(y=p_25), color = "cyan3") + 
    geom_line(aes(y=p_50)) + 
    geom_line(aes(y=p_75), color = "cyan3") + 
    geom_line(aes(y=p_95), linetype="dashed", color="red") + 
    geom_hline(yintercept=0, color="gray") +
    facet_wrap(~class, ncol=1) + 
    labs(x = "Month of Year",
         y = "NDVI") + 
    ggtitle("Greenness Phenology by Vegetation Type")
  
  # **** Visualize Wetness Phenology **** 
  pheno_wetness_df_long <- label_df %>% pivot_longer((12:22), names_to="month", values_to="ndvi")
  pheno_wetness_df_long$month <- as.numeric(pheno_wetness_df_long$month)
  pheno_wetness_stats <- pheno_wetness_df_long %>% 
    group_by(class, month) %>%
    summarise(p_05=quantile(ndvi, probs=c(0.05)), 
              p_25=quantile(ndvi, probs=c(0.25)),
              p_50=quantile(ndvi, probs=c(0.50)),
              p_75=quantile(ndvi, probs=c(0.75)),
              p_95=quantile(ndvi, probs=c(0.95)))
  # Natural Vegetation Phenology
  veg_wetness_plot <- ggplot(data=pheno_wetness_stats[pheno_wetness_stats$class %in% factor(c("Riparian","Oak","Chaparral","Annual","Evergreen Nonnative","Agriculture")),], 
                             aes(x=month-11)) + 
    geom_line(aes(y=p_05), linetype="dashed", color="red") + 
    geom_line(aes(y=p_25), color = "cyan3") + 
    geom_line(aes(y=p_50)) + 
    geom_line(aes(y=p_75), color = "cyan3") + 
    geom_line(aes(y=p_95), linetype="dashed", color="red") + 
    geom_hline(yintercept=0, color="gray") +
    facet_wrap(~class, ncol=1) + 
    labs(x = "Month of Year",
         y = "Tasselled Cap Wetness") + 
    ggtitle("Wetness Phenology by Vegetation Type")
  
  
  # **** Visualize Brightness Phenology **** 
  pheno_brightness_df_long <- label_df %>% pivot_longer((23:33), names_to="month", values_to="ndvi")
  pheno_brightness_df_long$month <- as.numeric(pheno_brightness_df_long$month)
  pheno_brightness_stats <- pheno_brightness_df_long %>% 
    group_by(class, month) %>%
    summarise(p_05=quantile(ndvi, probs=c(0.05)), 
              p_25=quantile(ndvi, probs=c(0.25)),
              p_50=quantile(ndvi, probs=c(0.50)),
              p_75=quantile(ndvi, probs=c(0.75)),
              p_95=quantile(ndvi, probs=c(0.95)))
  # Natural Vegetation Phenology
  veg_brightness_plot <- ggplot(data=pheno_brightness_stats[pheno_brightness_stats$class %in% factor(c("Riparian","Oak","Chaparral","Annual","Evergreen Nonnative","Agriculture")),], 
                                aes(x=month-22)) + 
    geom_line(aes(y=p_05), linetype="dashed", color="red") + 
    geom_line(aes(y=p_25), color = "cyan3") + 
    geom_line(aes(y=p_50)) + 
    geom_line(aes(y=p_75), color = "cyan3") + 
    geom_line(aes(y=p_95), linetype="dashed", color="red") + 
    geom_hline(yintercept=0, color="gray") +
    facet_wrap(~class, ncol=1) + 
    labs(x = "Month of Year",
         y = "Tasselled Cap Brightness") + 
    ggtitle("Brightness Phenology by Vegetation Type")
  
  
  
  
  
  # Rename month columns to have text names (not numbers) 
  names(label_df) <- c(paste("NDVI_",(month.name)[1:11],sep=""),
                       paste("Wetness_",(month.name)[1:11],sep=""),
                       paste("Brightness_",(month.name)[1:11],sep=""), 
                       "slope","dem","channel_dist","x","y","class")
  
  # ------ Generate Training Data ------
  num_samples <- 50
  # --- Riparian ---
  riparian_df <- label_df[label_df$class=="Riparian",]
  riparian_training_indices <- sample(1:nrow(riparian_df),num_samples)
  train <- riparian_df[riparian_training_indices,]
  leftover <- riparian_df[-riparian_training_indices,]
  # --- Oak ---
  oak_df <- label_df[label_df$class=="Oak",]
  oak_training_indices <- sample(1:nrow(oak_df),num_samples)
  train <- rbind(train, oak_df[oak_training_indices,])
  leftover <- rbind(leftover, oak_df[-oak_training_indices,])
  # --- Chaparral ---
  chaparral_df <- label_df[label_df$class=="Chaparral",]
  chaparral_training_indices <- sample(1:nrow(chaparral_df),num_samples)
  train <- rbind(train, chaparral_df[chaparral_training_indices,])
  leftover <- rbind(leftover, chaparral_df[-chaparral_training_indices,])
  # --- Annual ---
  annual_df <- label_df[label_df$class=="Annual",]
  annual_training_indices <- sample(1:nrow(annual_df),num_samples)
  train <- rbind(train, annual_df[annual_training_indices,])
  leftover <- rbind(leftover, annual_df[-annual_training_indices,])
  # --- Evergreen Nonnative ---
  euccup_df <- label_df[label_df$class=="Evergreen Nonnative",]
  euccup_training_indices <- sample(1:nrow(euccup_df),num_samples)
  train <- rbind(train, euccup_df[euccup_training_indices,])
  leftover <- rbind(leftover, euccup_df[-euccup_training_indices,])
  # --- Agriculture ---
  agriculture_df <- label_df[label_df$class=="Agriculture",]
  agriculture_training_indices <- sample(1:nrow(agriculture_df),num_samples)
  train <- rbind(train, agriculture_df[agriculture_training_indices,])
  leftover <- rbind(leftover, agriculture_df[-agriculture_training_indices,])
  # --- Soil ---
  soil_df <- label_df[label_df$class=="Soil",]
  soil_training_indices <- sample(1:nrow(soil_df),num_samples)
  train <- rbind(train, soil_df[soil_training_indices,])
  leftover <- rbind(leftover, soil_df[-soil_training_indices,])
  # --- Water ---
  water_df <- label_df[label_df$class=="Water",]
  water_training_indices <- sample(1:nrow(water_df),num_samples)
  train <- rbind(train, water_df[water_training_indices,])
  leftover <- rbind(leftover, water_df[-water_training_indices,])
  # --- Urban ---
  urban_df <- label_df[label_df$class=="Urban",]
  urban_training_indices <- sample(1:nrow(urban_df),num_samples)
  train <- rbind(train, urban_df[urban_training_indices,])
  leftover <- rbind(leftover, urban_df[-urban_training_indices,])
  
  
  
  # Generate some validation data
  #   to keep it balanced, we'll figure out the minimum number of remaining samples across all classes:
  remaining_per_class <- leftover %>% group_by(class) %>% summarise(num=n())
  min_remaining <- min(remaining_per_class$num)
  # --- Riparian ---
  riparian_df <- leftover[leftover$class=="Riparian",]
  riparian_training_indices <- sample(1:nrow(riparian_df),min_remaining)
  valid <- riparian_df[riparian_training_indices,]
  leftover_final <- riparian_df[-riparian_training_indices,]
  # --- Oak ---
  oak_df <- leftover[leftover$class=="Oak",]
  oak_training_indices <- sample(1:nrow(oak_df),min_remaining)
  valid <- rbind(valid, oak_df[oak_training_indices,])
  leftover_final <- rbind(leftover_final, oak_df[-oak_training_indices,])
  # --- Chaparral ---
  chaparral_df <- leftover[leftover$class=="Chaparral",]
  chaparral_training_indices <- sample(1:nrow(chaparral_df),min_remaining)
  valid <- rbind(valid, chaparral_df[chaparral_training_indices,])
  leftover_final <- rbind(leftover_final, chaparral_df[-chaparral_training_indices,])
  # --- Annual ---
  annual_df <- leftover[leftover$class=="Annual",]
  annual_training_indices <- sample(1:nrow(annual_df),min_remaining)
  valid <- rbind(valid, annual_df[annual_training_indices,])
  leftover_final <- rbind(leftover_final, annual_df[-annual_training_indices,])
  # --- Evergreen Nonnative ---
  euccup_df <- leftover[leftover$class=="Evergreen Nonnative",]
  euccup_training_indices <- sample(1:nrow(euccup_df),min_remaining)
  valid <- rbind(valid, euccup_df[euccup_training_indices,])
  leftover_final <- rbind(leftover_final, euccup_df[-euccup_training_indices,])
  # --- Agriculture ---
  agriculture_df <- leftover[leftover$class=="Agriculture",]
  agriculture_training_indices <- sample(1:nrow(agriculture_df),min_remaining)
  valid <- rbind(valid, agriculture_df[agriculture_training_indices,])
  leftover_final <- rbind(leftover_final, agriculture_df[-agriculture_training_indices,])
  # --- Soil ---
  soil_df <- leftover[leftover$class=="Soil",]
  soil_training_indices <- sample(1:nrow(soil_df),min_remaining)
  valid <- rbind(valid, soil_df[soil_training_indices,])
  leftover_final <- rbind(leftover_final, soil_df[-soil_training_indices,])
  # --- Urban ---
  urban_df <- leftover[leftover$class=="Urban",]
  urban_training_indices <- sample(1:nrow(urban_df),min_remaining)
  valid <- rbind(valid, urban_df[urban_training_indices,])
  leftover_final <- rbind(leftover_final, urban_df[-urban_training_indices,])
  # --- Water ---
  water_df <- leftover[leftover$class=="Water",]
  water_training_indices <- sample(1:nrow(water_df),min_remaining)
  valid <- rbind(valid, water_df[water_training_indices,])
  leftover_final <- rbind(leftover_final, water_df[-water_training_indices,])

  return(list(train,valid,label_df,leftover))
}

data_wet <- generate_training_data(phenoseries_wet)
data_dry <- generate_training_data(phenoseries_dry)
train_all <- rbind(data_wet[[1]],data_dry[[1]])
valid_all <- rbind(data_wet[[2]],data_dry[[2]])
leftover_all <- rbind(data_wet[[4]],data_dry[[4]])

# --- Train RF Model ---
modFit_rf <- train(as.factor(class) ~ NDVI_January+NDVI_February+NDVI_March+NDVI_April+NDVI_May+NDVI_June+NDVI_July+NDVI_August+NDVI_September+NDVI_October+NDVI_November+
                                      Brightness_January+Brightness_February+Brightness_March+Brightness_April+Brightness_May+Brightness_June+Brightness_July+Brightness_August+Brightness_September+Brightness_October+Brightness_November+
                                      Wetness_January+Wetness_February+Wetness_March+Wetness_April+Wetness_May+Wetness_June+Wetness_July+Wetness_August+Wetness_September+Wetness_October+Wetness_November+
                                      slope+dem+channel_dist, 
                   method="rf", data = train_all)
# --- Validation  Model (Balanced) ---
validation_results <- raster::predict(modFit_rf, newdata=valid_all)
# Confusion Matrix with Validation Data
confusionMatrix(validation_results, as.factor(valid_all$class))
# --- Validation Model (All) ---
validation_results_all <- raster::predict(modFit_rf, newdata=leftover)
confusionMatrix(validation_results_all, as.factor(leftover$class))
# --- Variable Importance ---
varImp(modFit_rf)

target_years <- c(2013,2017)#seq(2013,2019,2)

predict_raster <- function(year)
{
  phenoseries <- stack(load_phenoseries(year),slope,dem,channel_dist)
  # --- Generate Output Raster ---
  names(phenoseries) <- c(paste("NDVI_",(month.name)[1:11],sep=""),
                              paste("Wetness_",(month.name)[1:11],sep=""),
                              paste("Brightness_",(month.name)[1:11],sep=""), 
                              "slope","dem","channel_dist")
  pheno_df <- as.data.frame(phenoseries, xy=TRUE)
  predictions_all <- raster::predict(modFit_rf, newdata=pheno_df)
  predictions_all_num <- 1:length(predictions_all)
  predictions_all_num[predictions_all=="Riparian"] <- 1 
  predictions_all_num[predictions_all=="Oak"] <- 2
  predictions_all_num[predictions_all=="Chaparral"] <- 3
  predictions_all_num[predictions_all=="Annual"] <- 4
  predictions_all_num[predictions_all=="Evergreen Nonnative"] <- 5
  predictions_all_num[predictions_all=="Agriculture"] <- 6
  predictions_all_num[predictions_all=="Soil"] <- 7
  predictions_all_num[predictions_all=="Water"] <- 8
  predictions_all_num[predictions_all=="Urban"] <- 9
  predictions_df <- pheno_df %>% drop_na() %>% dplyr::select(1,2)
  predictions_df$class <- predictions_all_num
  predictions_raster <- rasterFromXYZ(predictions_df)
  crs(predictions_raster) <- crs(phenoseries)
  writeRaster(predictions_raster,
              paste("D:/SERDP/Vandenberg/Landsat/classified/classes_",year,".tif",sep="")
              ,overwrite=TRUE)
}

