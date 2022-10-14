
# Get NDVI phenology by stream in the VSFB area


# Masks INPUT_IMAGE to riparian vegetation around one target stream
#   STREAM_NAME - string name for target stream, assumes in NHD+ Flowlines variable "GNIS_Name"
#   RIPARIAN_MASK - 1 at riparian cells, NA at non-riparian cells
#   FLOWLINES - SpatialLineDataFrame object containing vector info on streams, assuming NHD+ format
maskToStream <- function(stream_name, riparian_mask, flowlines, input_image, buffer_width)
{
  target_stream <- flowlines[flowlines$GNIS_Name == stream_name,]
  target_stream_buffer <- buffer(target_stream, width=buffer_width)
  output_image <- mask(input_image, target_stream_buffer) * riparian_mask
  output_image <- crop(output_image, target_stream_buffer)
  return(output_image)
}

# Load Phenology Data
phenoseries_2014 <- stack("D:/SERDP/Vandenberg/good_data/indices/pheno_NDVI_2014_annual.tif")
phenoseries_2015 <- stack("D:/SERDP/Vandenberg/good_data/indices/pheno_NDVI_2015_annual.tif")
phenoseries_2016 <- stack("D:/SERDP/Vandenberg/good_data/indices/pheno_NDVI_2016_annual.tif")
phenoseries_2017 <- stack("D:/SERDP/Vandenberg/good_data/indices/pheno_NDVI_2017_annual.tif")
phenoseries_2018 <- stack("D:/SERDP/Vandenberg/good_data/indices/pheno_NDVI_2018_annual.tif")
phenoseries <- stack(phenoseries_2014, phenoseries_2015, phenoseries_2016, phenoseries_2017, phenoseries_2018)


# Load Hydrology data
flowlines <- readOGR("D:/SERDP/Vandenberg/hydrology/derived/flowlines_clipped_narrow.shp")
flowlines <- spTransform(flowlines, crs(phenoseries_2016))
# Subset to a few major streams we are interested in
target_streams <- c("Santa Ynez River", "Santa Maria River","San Antonio Creek")
flowlines <- flowlines[flowlines$GNIS_Name %in% target_streams,]

# Choose a Riparian Mask 
#   Classifier was run in timeseries
#   Keep something as 'riparian' if it was classified that way in at least 90% of models
classified_years <- c(2013,2015,2017,2019)
riparian_mask <- (raster(paste("D:/SERDP/Vandenberg/Landsat/classified/classes_",
                               classified_years[1],
                               ".tif",
                               sep=""))==1)

for(ind in classified_years)
{
  riparian_mask <- riparian_mask +
    (raster(paste("D:/SERDP/Vandenberg/Landsat/classified/classes_",
                  ind,
                  ".tif",
                  sep=""))==1)
}
riparian_mask <- riparian_mask >= 3
riparian_mask[riparian_mask == 0] <- NA

# Generate Image for Santa Ynez River
phenoseries_syr <- maskToStream("Santa Ynez River", riparian_mask, flowlines, phenoseries, buffer_width=1000)
names(phenoseries_syr) <- c(paste("ndvi_2014_",1:12,sep=""), 
                            paste("ndvi_2015_",1:12,sep=""), 
                            paste("ndvi_2016_",1:12,sep=""), 
                            paste("ndvi_2017_",1:12,sep=""), 
                            paste("ndvi_2018_",1:12,sep=""))
phenoseries_syr_df <- as.data.frame(phenoseries_syr) %>%
  drop_na %>%
  pivot_longer(1:60, names_to="raster_name", values_to="ndvi") %>%
  mutate(year = as.numeric(substr(raster_name, 6,9)),
         month = as.numeric(substr(raster_name, 11,12)))
phenoseries_syr_summary <- phenoseries_syr_df %>%
  group_by(year, month) %>%
  summarize(ndvi = median(ndvi))
phenoseries_syr_summary$stream <- rep("Santa Ynez River", nrow(phenoseries_syr_summary))

# Generate Image for Santa Maria River
phenoseries_smr <- maskToStream("Santa Maria River", riparian_mask, flowlines, phenoseries, buffer_width=1000)
names(phenoseries_smr) <- c(paste("ndvi_2014_",1:12,sep=""), 
                            paste("ndvi_2015_",1:12,sep=""), 
                            paste("ndvi_2016_",1:12,sep=""), 
                            paste("ndvi_2017_",1:12,sep=""), 
                            paste("ndvi_2018_",1:12,sep=""))
phenoseries_smr_df <- as.data.frame(phenoseries_smr) %>%
  drop_na %>%
  pivot_longer(1:60, names_to="raster_name", values_to="ndvi") %>%
  mutate(year = as.numeric(substr(raster_name, 6,9)),
         month = as.numeric(substr(raster_name, 11,12)))
phenoseries_smr_summary <- phenoseries_smr_df %>%
  group_by(year, month) %>%
  summarize(ndvi = median(ndvi))
phenoseries_smr_summary$stream <- rep("Santa Maria River", nrow(phenoseries_smy_summary))

# Generate Image for San Antonio Creek
phenoseries_sac <- maskToStream("San Antonio Creek", riparian_mask, flowlines, phenoseries, buffer_width=1000)
names(phenoseries_sac) <- c(paste("ndvi_2014_",1:12,sep=""), 
                            paste("ndvi_2015_",1:12,sep=""), 
                            paste("ndvi_2016_",1:12,sep=""), 
                            paste("ndvi_2017_",1:12,sep=""), 
                            paste("ndvi_2018_",1:12,sep=""))
phenoseries_sac_df <- as.data.frame(phenoseries_sac) %>%
  drop_na %>%
  pivot_longer(1:60, names_to="raster_name", values_to="ndvi") %>%
  mutate(year = as.numeric(substr(raster_name, 6,9)),
         month = as.numeric(substr(raster_name, 11,12)))
phenoseries_sac_summary <- phenoseries_sac_df %>%
  group_by(year, month) %>%
  summarize(ndvi = median(ndvi))
phenoseries_sac_summary$stream <- rep("San Antonio Creek", nrow(phenoseries_sac_summary))

phenoseries_summary <- rbind(phenoseries_syr_summary, phenoseries_smr_summary, phenoseries_sac_summary)

ggplot(data=phenoseries_summary) + 
  geom_line(aes(x=month, y=ndvi, group=year, col=year)) + 
  facet_wrap(~stream) + 
  xlab("Month") + 
  scale_x_continuous(breaks=seq(1,6)*2, expand=c(0,0)) + 
  ylab("Median NDVI") + 
  theme_bw() + 
  ggtitle("Seasonal Phenology of Riparian Plants")

