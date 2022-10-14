
library(tidyverse)
library(raster)

# Directory in which to run phenology program
directory <- "D:/landsat_vafb/data/stacked_scenes"

# Get all files in directory
files <- list.files(path=directory, pattern="")
files_full <- list.files(path=directory, full.names=TRUE)
# Filter for only large data files (no metadata, headers, etc.)
files <- files[sapply(files_full,file.size)>1000]
# Filter out input .tif stacked rasters and MESMA .csv files
files <- files[!substr(files,nchar(files)-3,nchar(files))==".tif"]
files <- files[!substr(files,nchar(files)-3,nchar(files))==".csv"]

# Range of years to analyze
# Can skip years if desired (e.g every two years) by changing last number below
target_years <- seq(1984,2011,1)

# Information on dates covered by image scenes
dates <- data.frame(year = as.numeric(as.character(substr(files,18,21))),
                    month = as.numeric(as.character(substr(files,22,23))),
                    day = as.numeric(as.character(substr(files,24,25))),
                    index = 1:length(files))
dates <- dates %>%
  mutate(doy = lubridate::yday(paste(year,month,day,sep="-")))

# Loads a .tif raster file
# Transforms it into a dataframe format
# Adds information about overpass date and image quality (cloud cover, etc.)
load_file <- function(index, band)
{
  # Get Dataset
  collection <- substr(files[[index]],6,9)
  dataset <- substr(files[[index]],1,4)
  filename_base_end <- str_locate(files[[1]], "_stacked.tif")[[1]]
  filename_base <- substr((files[[index]]),1,filename_base_end[[1]])
  # Load scene 
  scene <- raster(paste(directory,"/",files[[index]],sep=""),band=band)
  print(paste(directory,"/",filename_base,"quality.tif",sep=""))
  quality <- raster(paste(directory,"/",filename_base,"quality.tif",sep=""))
  #  if(collection=="L2SP")
  #  {
  #    good <- (quality==0)
  #  }
  if(dataset %in% c("LT05","LT04"))
  {
    good <- (quality %in% c(0, 672,676,680,684, 5440))
  }
  else if(dataset=="LC08")
  {
    good <- (quality %in% c(21824))
  }
  else
  {
    print(paste("ERROR: dataset ",dataset," not in expected values."),sep="")
  }
  
  # Return raster as a data frame 
  data <- as.data.frame(scene, xy=TRUE)
  names(data) <- c("easting","northing","value")
  data$doy <- rep(dates[index,]$doy, nrow(data))
  
  # Pixelwise image quality information (e.g. if clouds cover areas of the scene)
  goodness <- as.data.frame(good, xy=TRUE)
  names(goodness) <- c("easting","northing","good")
  bad_points <- sum(goodness$good==0)
  
  # Some debugging (report how much of scene is cloudy)
  print(paste("  For scene ",paste(directory,"/",files[[index]],sep="")), sep="")
  print(paste("  And quality raster ", paste(directory,"/",filename_base,"quality.tif",sep="")), sep="")
  print(paste("  Dropping ",
              bad_points," 
              points out of ",
              nrow(goodness),
              " total points, or ",
              bad_points/nrow(goodness)*100,"%.",
              sep=""))
  print("")
  
  # Remove all 'bad' pixels
  temp <- cbind(data,goodness)
  temp <- temp[temp$good > 0,]
  data <- temp[,c(1,2,3,4,7)] %>%
    drop_na()
  
  return(data)
}

# Loads individual raster files
#   converts to dataframe, 
#   rejects points with cloud clover, 
#   adds information about date
aggregate_scenes <- function(indices, band)
{
  scenes <- lapply(indices, load_file, band=band)
  year_data <- bind_rows(scenes)
}

# Get raster size, list of all Northing/Easting values, and map between pixel indices and coordinates
scene_for_crs <- raster(paste(directory,"/",files[[1]],sep=""))
scene_for_crs_df <- as.data.frame(scene_for_crs, xy=TRUE)
northings <- sort(unique(scene_for_crs_df$y))
num_rows <- length(northings)
eastings <- sort(unique(scene_for_crs_df$x))
num_cols <- length(eastings)
north_map <- 1:length(northings)
names(north_map) <- northings
east_map <- 1:length(eastings)
names(east_map) <- eastings
# Initialize Data Structures
output <- array(rep(0,num_rows*num_cols*12), c(num_rows,num_cols,12))
scene_count_raster <- scene_for_crs*0
scene_count <- array(rep(0,num_rows*num_cols), c(num_rows,num_cols))

# Actually process data!
generate_pheno_raster <- function(band, name)
{
  for(current_year in target_years)
  {
    # Get dates for all scenes within this or the next year
    #year_dates <- dates %>% filter(year>=current_year & year<=current_year+1)
    # Get dates for all scenes within this year
    year_dates <- dates %>% filter(year==current_year)
    print(current_year)
    
    # Reset Output Matrix
    output <- output*0
    
    # Get all scenes within the target year
    year_data <- aggregate_scenes(year_dates$index, band)
    
    # Setup output filename
    output_filename <- paste(directory,"/pheno_",name,"_",current_year,".tif",sep="")
    output_count_filename <- paste(directory,"/pheno_",name,"_",current_year,"_scene_count_annual.tif",sep="")
    
    print(paste("Starting work on ", output_filename, ", which is ",num_cols,"x",num_rows,"in size.", sep=""))
    
    year_data$northing <- as.character(year_data$northing)
    year_data$easting <- as.character(year_data$easting)
    ind_i <- 0
    for(northing_i in unique(northings))
    {
      ind_j <- 0
      data_i <- filter(year_data, northing==northing_i)
      data_i_list <- data_i %>% group_split(easting) #split(data_i, data_i[1])
      for(pixel in data_i_list)
      {
        if(nrow(pixel) < 2)
        {
          ind_j <- ind_j+1
          next
        }
        # Bookend phenoseries with year-shifted first and last values, to prevent erroneous forecasting/hindcasting
        early_point <- pixel[nrow(pixel),]
        early_point$doy <- early_point$doy-365
        late_point <- pixel[1,]
        late_point$doy <- late_point$doy+365
        pixel <- rbind(early_point, pixel, late_point)
        # Estimate phenology curve - linear interpolation between datapoints
        new_series <- approx(pixel$doy,pixel$value,xout=((1:12)*30))
        # Populate phenology data matrix
        output[north_map[[pixel[1,]$northing]],east_map[[pixel[1,]$easting]],1:12] <- new_series$y
        # Update raster which keeps track of how many good scenes there are for each pixel 
        scene_count[ind_i+1,ind_j+1] <- nrow(pixel) - 2
        ind_j <- ind_j+1
      }
      print(paste(ind_i/num_rows*100, "% done with file ", output_filename, sep=""))
      ind_i <- ind_i+1
    }
    # Populate raster using data matrix
    pheno_raster <- flip(brick(output),direction=2)
    scene_count_raster <- flip(raster(scene_count),direction=2)
    # Georeference rasters
    crs(pheno_raster) <- crs(scene_for_crs)
    origin(pheno_raster) <- origin(scene_for_crs)
    extent(pheno_raster) <- extent(scene_for_crs)
    crs(scene_count_raster) <- crs(scene_for_crs)
    origin(scene_count_raster) <- origin(scene_for_crs)
    extent(scene_count_raster) <- extent(scene_for_crs)
    # Output raster to disk 
    writeRaster(pheno_raster, filename=output_filename,overwrite=TRUE)
    writeRaster(scene_count_raster, filename=output_count_filename,overwrite=TRUE)
  }
}

generate_pheno_raster(1,"GV")
generate_pheno_raster(2,"NPV")
generate_pheno_raster(3,"Shade")
generate_pheno_raster(4,"Soil")