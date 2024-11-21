
## data from isimip https://data.isimip.org/
# Install the ncdf4 and raster packages if not already installed
install.packages("ncdf4")
install.packages("raster")
install.packages("RColorBrewer")

# Load the libraries
library(ncdf4)
library(raster)
library(RColorBrewer)

# calculate gridcell area
    # Define the grid dimensions: 360 rows (latitudes) and 720 columns (longitudes)
    nrows <- 360
    ncols <- 720

    # Latitude and longitude resolution (degrees)
    res_deg <- 0.5

    # Create the grid for latitude and longitude
    lats <- seq(-89.75, 89.75, by=res_deg)  # Latitude grid (-90 to 90 degrees)
    lons <- seq(-179.75, 179.75, by=res_deg)  # Longitude grid (-180 to 180 degrees)

    # Calculate the latitudinal distance (constant)
    lat_dist <- 111000 * res_deg  # 111 km per degree of latitude, multiplied by resolution

    # Calculate the longitudinal distance that varies with latitude
    lon_dist <- 111000 * cos(lats * pi / 180) * res_deg  # Varies with latitude

    # Create a matrix to hold the area values for each grid cell
    area_matrix <- matrix(0, nrow=nrows, ncol=ncols)

    # Calculate the area for each grid cell
      for (i in 1:nrows) {
        for (j in 1:ncols) {
          # Multiply the latitudinal distance by the longitudinal distance for each grid cell
            area_matrix[i, j] <- lat_dist * lon_dist[i]
          }
        }

    # Convert the matrix to a raster
    cell_area_raster <- raster(area_matrix)

    # Set the extent and CRS for the raster (global grid)
    extent(cell_area_raster) <- extent(-180, 180, -90, 90)
    crs(cell_area_raster) <- CRS("+proj=longlat +datum=WGS84")

    # Plot the cell area raster
    plot(cell_area_raster, main="Grid Cell Areas (m²)", col=terrain.colors(100))


# Open NetCDF files
nc_factual <- nc_open("/Users/stijn.hantson/Documents/Documents/clases/Bosque_cambio_climatio/2024/semana13/orchidee-mict_gswp3-w5e5_obsclim_histsoc_default_gpp-total_global_monthly_1901_2019.nc")
nc_counterfactual <- nc_open("/Users/stijn.hantson/Documents/Documents/clases/Bosque_cambio_climatio/2024/semana13/orchidee-mict_gswp3-w5e5_counterclim_histsoc_default_gpp-total_global_monthly_1901_2019.nc")

# Display metadata
print(nc_factual)
print(nc_counterfactual)

# Close NetCDF files after extracting data
nc_close(nc_factual)
nc_close(nc_counterfactual)

# Load each NetCDF file as a raster stack
factual_stack <- stack("/Users/stijn.hantson/Documents/Documents/clases/Bosque_cambio_climatio/2024/semana13/orchidee-mict_gswp3-w5e5_obsclim_histsoc_default_gpp-total_global_monthly_1901_2019.nc", varname = "gpp-total")
counterfactual_stack <- stack("/Users/stijn.hantson/Documents/Documents/clases/Bosque_cambio_climatio/2024/semana13/orchidee-mict_gswp3-w5e5_counterclim_histsoc_default_gpp-total_global_monthly_1901_2019.nc", varname = "gpp-total")

# Define the number of years
months_in_year <- 12
num_years <- nlayers(factual_stack) / months_in_year
seconds_per_month <- 30 * 24 * 60 * 60  # For a 365-day calendar

# Initialize stacks for storing annual GPP for each scenario
annual_gpp_factual <- stack()
annual_gpp_counterfactual <- stack()

# Loop through each year
for (year in 1:num_years) {
  # Define the start and end layers for the 12-month period
  start_idx <- (year - 1) * months_in_year + 1
  end_idx <- year * months_in_year
  
  # Sum the layers over the 12-month period for both factual and counterfactual stacks
  annual_sum_factual <- sum(factual_stack[[start_idx:end_idx]], na.rm = TRUE) * seconds_per_month
  annual_sum_counterfactual <- sum(counterfactual_stack[[start_idx:end_idx]], na.rm = TRUE) * seconds_per_month
  
  # Add the annual sum as a new layer in the annual GPP stacks
  annual_gpp_factual <- addLayer(annual_gpp_factual, annual_sum_factual)
  annual_gpp_counterfactual <- addLayer(annual_gpp_counterfactual, annual_sum_counterfactual)
}

# Calculate mean annual GPP for each grid cell across the last 10 years
mean_annual_gpp_factual <- calc(annual_gpp_factual[[109:119]], mean, na.rm = TRUE)
mean_annual_gpp_counterfactual <- calc(annual_gpp_counterfactual[[109:119]], mean, na.rm = TRUE)

plot(mean_annual_gpp_factual, main="Mean Annual GPP (Factual) [kg m -2 year-1]")
plot(mean_annual_gpp_counterfactual, main="Mean Annual GPP (Factual) [kg m -2 year-1]")


# Calculate the difference in mean annual GPP
gpp_difference <- mean_annual_gpp_factual - mean_annual_gpp_counterfactual
gpp_difference[gpp_difference == 0] <- NA

palette <- brewer.pal(11, "RdYlGn")

plot(gpp_difference, col=palette, zlim=c(-1.5, 1.5))

# calculate global mean GPP

per_gridcell_gpp_factual = mean_annual_gpp_factual*cell_area_raster
per_gridcell_gpp_counterfactual = mean_annual_gpp_counterfactual*cell_area_raster

total_gpp_factual <- cellStats(per_gridcell_gpp_factual, stat = "sum", na.rm = TRUE)
total_gpp_counterfactual <- cellStats(per_gridcell_gpp_counterfactual, stat = "sum", na.rm = TRUE)

total_gpp_factual_pg = total_gpp_factual* 1e-12
total_gpp_counterfactual_pg = total_gpp_counterfactual* 1e-12

print(total_gpp_factual_pg)
print(total_gpp_counterfactual_pg)
