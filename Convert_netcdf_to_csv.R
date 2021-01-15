# This file takes netcdf files with min and max daily temperatures on a UK grid and combines them to give individual
# csv files for each grid square with a time series of temperatures.

# load packages
library(ncdf4)
library(fields)
library(sp)
library(rgdal)
library(abind)

setwd("N:/WNV_Paper/Data/UKCP18/")

# Load data files.  Can be found at address below.  In this case I have loaded 2070-2080 for climate run 13.
# Separate netcdf files for minimum and maximum daily temperature which I want to combine into a series of
# .csv files for each grid square because that's what the WNV model code is set up to work with.
#Citable as:  Met Office Hadley Centre (2018): UKCP18 Regional Projections on a 12km grid over the UK for 1980-2080.
#Centre for Environmental Data Analysis, 15/01/2021. 
#https://catalogue.ceda.ac.uk/uuid/589211abeb844070a95d061c8cc7f604
min_file <- nc_open("tasmin_rcp85_land-rcm_uk_12km_13_day_20701201-20801130.nc")
max_file <- nc_open("tasmax_rcp85_land-rcm_uk_12km_13_day_20701201-20801130.nc")

# Extract temperature variables
temp_min <- ncvar_get(min_file,"tasmin")
temp_max <- ncvar_get(max_file,"tasmax")

# Extract longitude and latitude.  There is some inconsistency in the way different versions of the climate files
# are presented so the commented out code may be required.
# #lon <- min_file$var$grid_longitude$dim[[1]]$vals
# #lat <- min_file$var$grid_longitude$dim[[2]]$vals
lon <- min_file$var$longitude$dim[[1]]$vals
lat <- min_file$var$longitude$dim[[2]]$vals
# 

# Load the outline of the UK
setwd("N:/WNV_Paper/Data/UKCP18")
GBR.outline <- readOGR(dsn="GBR_adm", layer="GBR_adm0")

# Subset so that only grid squares over land within the UK are run
lon.lat.grid <- expand.grid(lon,lat)
temp.grid <- SpatialPoints(lon.lat.grid, proj4string=CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000"))
temp.grid.conv <- spTransform(temp.grid, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
points.in.gbr <- over(temp.grid.conv, GBR.outline)
coords.in.gbr <- which(!is.na(points.in.gbr$SOVEREIGN), arr.ind=TRUE)

# Subset to the time period desired
temp_min <- temp_min[,,2190:3270]
temp_max <- temp_max[,,2190:3270]

# Bind together minimum and maximum temperatures to create a single time series
temp_minmax <- abind(temp_min[,,1], temp_max[,,1],along=3)
temp_minmax[,,1] <- 10
for(i in 2:1080){
  temp_minmax <- abind(temp_minmax,temp_min[,,i], temp_max[,,i])
}

# Write a .csv file for each grid square
for(i in 1:length(lat)){
  for(j in 1:length(lon)){
    if(((i-1)*length(lon)+j) %in% coords.in.gbr){
      current_coord <- SpatialPoints(matrix(c(lon[j],lat[i]),ncol=2), proj4string=CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000"))
      trans.coord <- spTransform(current_coord, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
      write.csv(c(trans.coord@coords[2],temp_minmax[j,i,]),file=paste("C:/Users/davewi/Documents/Temperature_data/2077-2079/Run 13/tasminmax_rcp85_land_rcm_uk_12km_13_day_2077_2079_lon",j,"lat",i,".csv",sep=""),row.names=FALSE)
    }
  }
}