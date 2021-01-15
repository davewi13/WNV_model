# This file takes the WNV model (run in fortran) output as an input and creates dataframes summarising peak MIR
# and vector-host ratios on the UK grid for plotting.  Note that the VHRs are not required for any of the plots
# in the paper, however they were used to help tune the vector-host ratio to a sensible level.

# load packages
library(ncdf4)
library(sp)
library(rgeos)
library(rgdal)

# Will need the latitude and longitude vectors again (originally used in "Convert_netcdf_to_csv.R")
setwd("N:/WNV_Paper/Data/UKCP18/")
min_file <- nc_open("tasmin_rcp85_land-rcm_uk_12km_13_day_20701201-20801130.nc")
max_file <- nc_open("tasmax_rcp85_land-rcm_uk_12km_13_day_20701201-20801130.nc")
lon <- min_file$var$longitude$dim[[1]]$vals
lat <- min_file$var$longitude$dim[[2]]$vals
 
# Create vectors listing the possible run numbers and the time periods considered.  The years actually refer to
# the last year in their decade e.g. 2070 refers to climate data from 2077-2079 with WNV introduction in 2079
runsvec <- c("01","04","05","06","07","08","09","10","11","12","13","15")
yearsvec <- c("2010", "2020", "2030", "2040", "2050", "2060", "2070")
counter1 <- 1
counter2 <- 1
# Loop through all years and climate runs for which the WNV model was run
for(year in yearsvec){
  for(run in runsvec){
    # This assumes that the set of output files for each combination of year and climate model run is stored
    # in its own folder (and then within a subfolder called output files, though this could easily be omitted)
    setwd(paste0("C:/Users/davewi/Documents/Results/WNV_UKCP18_",run,"_",year,"s_120intro/Output_files/"))
    
    # List all the output files
    file.names <- list.files()
    
    # Create space for results
    lat.lon.mat <- matrix(NA,nrow=length(file.names),ncol=2)
    peak.mir <- numeric(length(file.names))
    peak.vhr <- numeric(length(file.names))
    mean.vhr <- numeric(length(file.names))
    colnames(lat.lon.mat) <- c("lon","lat")
    
    # Loop through all output files
    for(i in 1:length(file.names)){
      # Extract the relevant latitude and longitude from the filename
      if(nchar(file.names[i]) == 64){
        lat.lon.vec <- as.numeric(c(substr(file.names[i], 54, 55), substr(file.names[i], 59, 60)))
      } else {
        lat.lon.vec <- as.numeric(c(substr(file.names[i], 54, 55), substr(file.names[i], 59, 61)))
      }
      # Load the file
      temp.file <- read.table(file.names[i])
      # Calculate the MIR at each time point.  Column 7 corresponds to infectious mosquitoes, 6 is exposed and 5 is
      # susceptible.  The multiplier of 0.5 is because MIR is only based on females but the model predicts total
      # adult abundance.  Drop row 1 to avoid ever dividing by zero.  Peak MIR should never be on day 1 anyway.
      mir.vec <- (temp.file[-1,7]/((temp.file[-1,5]+temp.file[-1,6])*0.5))*1000
      # Take the max value
      peak.mir[i] <- max(mir.vec)
      # Calculate vector-host ratio over the biting season in the year of WNV introduction
      vhr.vec <- ((temp.file[-c(1:880,971:1080),7]+temp.file[-c(1:880,971:1080),6]+temp.file[-c(1:880,971:1080),5])*0.5)/(temp.file[-c(1:880,971:1080),20]+temp.file[-c(1:880,971:1080),21]+temp.file[-c(1:880,971:1080),22])
      # Take the peak and the mean
      peak.vhr[i] <- max(vhr.vec)
      mean.vhr[i] <- mean(vhr.vec)
      # Store latitude and longitude
      lat.lon.mat[i,] <- c(lon[lat.lon.vec[1]], lat[lat.lon.vec[2]])
    }
    # Combine latitude and longitude with MIR and VHR predictions
    sp.lat.lon.mat <- SpatialPoints(lat.lon.mat, proj4string=CRS("+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +no_defs 
                                                             +ellps=airy +datum=OSGB36 +units=m"))
    lat.lon.trans.coord <- spTransform(sp.lat.lon.mat, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    lat.lon.trans.coord$mir <- peak.mir
    lat.lon.trans.coord$vhr <- peak.vhr
    lat.lon.trans.coord$meanvhr <- mean.vhr
    
    # Save the output file
    save(lat.lon.trans.coord, file=paste0("C:/Users/davewi/Documents/Results/120intro_dfs/results_dataframe_",run,"_",year,"s_120intro.Rdata"))
  }
}
