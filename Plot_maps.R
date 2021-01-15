# This code is used to generate the maps of the spatial distribution of WNV risk (as seen in Figures 1, 2, 4, SF1 and SF2).
# Note that part 2 of this code (the legend functions) must be run first before attempting to create the maps.
# The code here takes the dataframes created by "Convert_WNVmodel_output_to_dataframes.R" and the files in the 
# Datafiles folder as input and produces maps.

# Load packages
library(rgdal)
library(raster)
library(rgeos)
library(sp)
library(RColorBrewer)

# Create vectors of all possible climate runs and years
runsvec <- c("01","04","05","06","07","08","09","10","11","12","13","15")
yearsvec <- c("2010", "2020", "2030", "2040", "2050", "2060", "2070")
# This yearstext vector is used to annotate the maps
yearstext <- c("2019", "2029", "2039", "2049", "2059", "2069", "2079")
counter1 <- 1
counter2 <- 1
# Loop through all combinations of years and runs
for(year in yearsvec){
  for(run in runsvec){
    # Load the dataframe created previously
    filename <- paste0("C:/Users/davewi/Documents/Results/Tmin_dfs/results_dataframe_",run,"_",year,"s_Tmins.Rdata")
    load(file=filename)
    
    options(scipen=22)
    
    # should contain 'temp.grid'; load grid cell centre coordinates
    setwd("N:/WNV_Paper/Files from Stephen")
    load("Grid_in_trans_merc.Rdata")
    # number of grid cells
    n=dim(temp.grid@coords)[1]
    # attributes for OS National Grid
    bng = '+proj=tmerc +lat_0=49 +lon_0=-2 +k=0.9996012717 +x_0=400000 +y_0=-100000 +ellps=airy +datum=OSGB36 +units=m +no_defs'
    # should contain 'map' and 'ix'
    load("map.RData")
    map@proj4string<-CRS(bng)
    # transform to OS grid
    MIRos<-spTransform(lat.lon.trans.coord,CRS(bng))
    
    # initially risk is set to zero for every grid cell NB this includes cells in the sea and Republic of Ireland i.e. all those in temp.grid
    risk.mir <- rep(0,n) 
    risk.vhr <- rep(0,n)
    risk.meanvhr <- rep(0,n)
    MIRdata<-as.data.frame(MIRos)
    x=temp.grid@coords[,1]
    y=temp.grid@coords[,2]

    # There seems to be an issue with the masking on some of the coastal grid squares due to grid squares with a large proportion of sea giving unreliable results.  Remove these.
    rem.vec <- c(70, 97, 111, 112, 128, 143, 144, 158, 159, 165, 197, 206, 209, 236, 237, 241, 278, 279, 304, 360, 423, 485, 520, 545, 571, 691, 771, 814, 929, 1248, 1292, 1332, 1368, 1438, 1465, 1472, 1557, 1581, 1619, 1644, 1645)
    MIRdata <- MIRdata[-rem.vec,]
    MIRdata <- na.omit(MIRdata)
    
    for(i in 1:dim(MIRdata)[1]) # assign i'th risk value to appropriate grid cell
    {
      index=(1:n)[abs(MIRdata$lon[i]-x)<0.1 & abs(MIRdata$lat[i]-y)<0.1] # allow some tolerance in case of discrepancies due to the projection transformations
      if(length(index)!=1) {cat("\nError!\n"); browser()} # something's gone wrong
      risk.mir[index]=MIRdata$mir[i]
      risk.vhr[index]=MIRdata$vhr[i]
      risk.meanvhr[index]=MIRdata$meanvhr[i]
    }
    
    # Set up colours for map
    colors <- brewer.pal(4, "Reds") # for risk maps, data maps
    pal <- colorRampPalette(colors)
    
    n.breaks <- 21
    
    # Give each grid cell a colour
    breakvector.mir <- c(-0.001, seq(0.25,5,by=0.25))
    if(max(MIRdata$mir) > breakvector.mir[length(breakvector.mir)]){breakvector.mir[n.breaks]=(max(MIRdata$mir)+0.01)}
    
    breakvector.vhr <- ((-0.001-(1/(n.breaks))+(1/(n.breaks))*(1:n.breaks))*10)
    breakvector.vhr[n.breaks] = (max(MIRdata$vhr)+0.01)
    
    breakvector.meanvhr<- ((-0.001-(1/(n.breaks-1))+(1/(n.breaks-1))*(1:n.breaks))*5)
    breakvector.meanvhr[n.breaks]=(max(MIRdata$meanvhr)+0.01)
    
    lablist<-c("white",pal(n.breaks-2)); # thresholds for categorising risks
    mc.mir <- as.character(cut(risk.mir, breaks=breakvector.mir, labels=lablist))
    mc.vhr <- as.character(cut(risk.vhr, breaks=breakvector.vhr, labels=lablist))
    mc.meanvhr <- as.character(cut(risk.meanvhr, breaks=breakvector.meanvhr, labels=lablist))
    
    # Give each polygon a colour e.g. a white cell on the England/Scotland border will be split into 2 polygons and both of these polygons need to be given a white colour
    fullmc.mir <- rep(0,sum(ix))
    fullmc.vhr<-rep(0,sum(ix))
    fullmc.meanvhr <- rep(0,sum(ix))
    ptr=1
    for(j in 1:n)
    {
      num=sum(ix[,j])
      if(num>0) {fullmc.mir[ptr:(ptr-1+num)]<-mc.mir[j];ptr<-ptr+num;}
      if(num>0) {fullmc.vhr[ptr:(ptr-1+num)]<-mc.vhr[j];ptr<-ptr+num;}
      if(num>0) {fullmc.meanvhr[ptr:(ptr-1+num)]<-mc.meanvhr[j];ptr<-ptr+num;}
    }
    
    # Plot map: remember to run legend functions below before doing this
    mapname=paste0("C:/Users/davewi/Documents/WNV_Paper_Plots/",year,"s_",run,"_Tmin.png")
    png(mapname,width=1230,height=1405)
    plot(map,col=fullmc.mir,axes=F,border=NA,bg="lightgrey",main="",cex.main=4,xaxs="i",yaxs="i")
    xx=c(0,40000,40000,0)
    yy=c(0,400000,400000,0)
    text(x=0, y=1100000, yearstext[counter1], cex= 7)
    # Legend is commented out here to create multiple maps with one legend for all of them
    # Uncomment for standalone maps
    #legend.gradient(cbind(x = xx-200000, y = yy+50000), cols = c("white",pal(21)), title = "", limits = c("0","7.5","2.5","5"))
    dev.off()
    
    # Repeat code from mapname declaration to dev.off() with minor tweaks for VHR maps.
    counter2 <- counter2 + 1
  }
  counter1 <- counter1 + 1
  
}


##############################################################################################################################################################
### PART 2: LEGEND FUNCIONS                        ###########################################################################################################
### RUN THESE BEFORE PLOTTING MAPS                 ###########################################################################################################
###                                                ###########################################################################################################
##############################################################################################################################################################
legend.gradient<-function (pnts, cols = heat.colors(100), limits = c(0, 1), title = "Legend", 
                           ...) 
{
  pnts = try(as.matrix(pnts), silent = T)
  
  
  yvals = seq(min(pnts[, 2]), max(pnts[, 2]), length = length(cols) + 
                1)
  for (i in 1:length(cols)) {
    polygon(x = pnts[, 1], y = c(yvals[i], yvals[i], yvals[i + 
                                                             1], yvals[i + 1]), col = cols[i], border = F)
  }
  text(max(pnts[, 1]), min(pnts[, 2]), labels = limits[1], 
       pos = 4, ...)
  text(max(pnts[, 1]), max(pnts[, 2]), labels = limits[2], 
       pos = 4, ...)
  text(max(pnts[, 1]), min(pnts[, 2])+(1/3)*(max(pnts[, 2])-min(pnts[, 2])), labels = limits[3], 
       pos = 4, ...)
  text(max(pnts[, 1]), min(pnts[, 2])+(2/3)*(max(pnts[, 2])-min(pnts[, 2])), labels = limits[4], 
       pos = 4, ...)
  text(min(pnts[, 1]), max(pnts[, 2]), labels = title, adj = c(0, 
                                                               -1), ...)
}

legend.groups<-function (pnts, cols = heat.colors(100), limits = c(0, 1), title = "Legend", 
                         ...) 
{
  pnts = try(as.matrix(pnts), silent = T)
  
  
  yvals = seq(min(pnts[, 2]), max(pnts[, 2]), length = length(cols) + 
                1)
  for (i in 1:length(cols)) {
    polygon(x = pnts[, 1], y = c(yvals[i], yvals[i], yvals[i + 
                                                             1], yvals[i + 1]), col = cols[i], border = F)
  }
  #text(max(pnts[, 1]), min(pnts[, 2]), labels = limits[1], 
  #pos = 4, ...)
  #text(max(pnts[, 1]), max(pnts[, 2]), labels = limits[2], 
  #pos = 4, ...)
  text(max(pnts[, 1]), min(pnts[, 2])+(1/6)*(max(pnts[, 2])-min(pnts[, 2])), labels = limits[1], 
       pos = 4, ...)
  text(max(pnts[, 1]), min(pnts[, 2])+(3/6)*(max(pnts[, 2])-min(pnts[, 2])), labels = limits[2], 
       pos = 4, ...)
  text(max(pnts[, 1]), min(pnts[, 2])+(5/6)*(max(pnts[, 2])-min(pnts[, 2])), labels = limits[3], 
       pos = 4, ...)
  text(min(pnts[, 1]), max(pnts[, 2]), labels = title, adj = c(0, 
                                                               -1), ...)
}