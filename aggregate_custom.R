
aggregate_custom <- function(input_raster, target_raster, reproj_method="bilinear")
{
  output_raster <- projectRaster(input_raster, target_raster, method=reproj_method)
  output_raster <- crop(output_raster, projectExtent(input_raster, crs(output_raster)))
  
  return(output_raster)
}