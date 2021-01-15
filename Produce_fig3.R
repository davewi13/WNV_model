# Function for use later
whichpart <- function(x, n=1664) {
  nx <- length(x)
  p <- nx-n
  print(nx)
  print(p)
  xp <- sort(x, partial=p)[p]
  which(x > xp)
}

# Biting rate relationship with temperature function
biting_rate <- function(temp, q=7.83*10^-5, Tmin=11.4, Tmax=45.2, calc_med=TRUE){
  temp[temp > 45.2] <- 45.2
  biting_rate <- q*temp*(temp-Tmin)*sqrt(Tmax-temp)
  biting_rate[biting_rate < 0.005] <- 0.005
  if((length(biting_rate) > 1) && (calc_med==TRUE)){biting_rate <- median(biting_rate)}
  return(biting_rate)
}

# Create vectors of time periods and runs both for referencing and for plotting
yearvec <- c("2010s","2020s","2030s","2040s","2050s","2060s","2070s")
runsvec <- c("01","04","05","06","07","08","09","10","11","12","13","15")
yearvectemp <- c("2017-2019","2027-2029","2037-2039","2047-2049","2057-2059","2067-2069","2077-2079")
runsvectemp <- c("Run 1","Run 4","Run 5","Run 6","Run 7","Run 8","Run 9","Run 10","Run 11","Run 12","Run 13","Run 15")
# There seems to be an issue with the masking on some of the coastal grid squares due to grid squares with a large proportion of sea giving unreliable results.  Remove these.
droplon <- c(30,33,35,38,39,40,42,42,45,46,47,49,55,56,57,58,60,60,61,64,65,66,68,68)
droplat <- c(12,52,13,14,25,37,13,41,17,42,25,16,17,16,17,18,18,45,18,18,18,38,24,25)
droplon <- as.character(droplon)
droplat <- as.character(droplat)
droplonlat <- as.character(length(droplon))
for(a in 1:length(droplon)){
  droplonlat[a] <- paste0(droplon[a],droplat[a])
}

# Create space for results
results.mat <- matrix(c(rep(yearvec,each=20424), rep(rep(runsvec, each=1702),7), rep(1:1702, 84), rep(NA, 142968*12)), byrow=F, nrow=142968) 
colnames(results.mat) <- c("Year","Run.No","Replicate","Lon","Lat","EIP.peak","EIP.avg","EIP.sum","BR.peak","BR.avg","BR.sum","VHR.peak","VHR.avg","VHR.sum","MIR")
counter <- 1

# Loop through all combinations of years and climate runs
for(i in 1:length(yearvec)){
  for(j in 1:length(runsvec)){
    # Set directory to location of file storage and drop the files that have masking issues
    filepath <- paste0("C:/Users/davewi/Documents/Results/WNV_UKCP18_",runsvec[j],"_",yearvec[i],"_150intro/Output_files")
    setwd(filepath)
    file.names <- list.files()
    file.inds <- which(paste0(substr(file.names,54,55),substr(file.names,59,60)) %in% droplonlat)
    file.names <- file.names[-file.inds]
    # Create space for outputs
    peak.mir <- numeric(length(file.names))
    temp.names <- numeric(length(file.names))
    
    # Loop through all files and determine peak mir for each one
    for(k in 1:length(file.names)){
      temp.file <- read.table(file.names[k], skip=869)
      mir.vec <- (temp.file[,7]/((temp.file[,5]+temp.file[,6])*0.5))*1000
      peak.mir[k] <- max(mir.vec)
      if(nchar(file.names[k]) == 65){
        temp.names[k] <- paste0("C:/Users/davewi/Documents/Temperature_data/",yearvectemp[i],"/",runsvectemp[j],"/tasminmax_rcp85_land_rcm_uk_12km_",runsvec[j],"_day_",substr(yearvectemp[i],1,4),"_",substr(yearvectemp[i],6,9),substr(file.names[k],50,61),".csv")
      } else {
        temp.names[k] <- paste0("C:/Users/davewi/Documents/Temperature_data/",yearvectemp[i],"/",runsvectemp[j],"/tasminmax_rcp85_land_rcm_uk_12km_",runsvec[j],"_day_",substr(yearvectemp[i],1,4),"_",substr(yearvectemp[i],6,9),substr(file.names[k],50,60),".csv")
      }
    }
    # Create vector where we'll store the time of the peak MIR
    peak.mir.inds <- peak.mir
    
    # Loop through length of all files
    for(l in 1:length(peak.mir.inds)){
      # Load results
      temp.file <- read.table(file.names[l], skip=869)
      # Load temperature data for this file
      temperature.file <- read.table(temp.names[l], skip=1739)
      # Populate the results matrix with EIP, VHR and BR values
      temp.eip <- 1:length(temp.file[,19])-temp.file[,19]
      results.mat[counter,1] <- substr(file.names[l], 46, 49)
      results.mat[counter,2] <- substr(file.names[l], 34, 35)
      results.mat[counter,4] <- lon[as.numeric(substr(file.names[l], 54, 55))]
      if(nchar(file.names[l]) == 65){
        results.mat[counter,5] <- lat[as.numeric(substr(file.names[l], 59, 61))]
      } else {
        results.mat[counter,5] <- lat[as.numeric(substr(file.names[l], 59, 60))]
      }
      results.mat[counter,6] <- temp.file[sum(temp.eip < 0)+1,19]
      results.mat[counter,7] <- median(temp.file[1:90,19])
      results.mat[counter,8] <- sum(median(temp.file[1:90,19]))
      results.mat[counter,9] <- biting_rate(temperature.file[1:4,1])
      results.mat[counter,10] <- biting_rate(temperature.file[1:180,1])
      results.mat[counter,11] <- sum(biting_rate(temperature.file[1:180,1], calc_med=FALSE))
      vhr.vec <- (rowSums(temp.file[1:90,5:7])*0.5)/rowSums(temp.file[1:90,20:22])
      results.mat[counter,12] <- vhr.vec[1]
      results.mat[counter,13] <- median(vhr.vec)
      results.mat[counter,14] <- sum(vhr.vec)
      results.mat[counter,15] <- peak.mir[l]
      counter <- counter+1
    }
    print(paste0("Year ",yearvec[i]," run number ", runsvec[j], " completed"))
  }
}

# Convert columns to numeric
results.mat <- data.frame(results.mat)
results.mat$Lat <- as.numeric(as.character(results.mat$Lat))
results.mat$Lon <- as.numeric(as.character(results.mat$Lon))
results.mat$EIP.peak <- as.numeric(as.character(results.mat$EIP.peak))
results.mat$BR.peak <- as.numeric(as.character(results.mat$BR.peak))
results.mat$VHR.peak <- as.numeric(as.character(results.mat$VHR.peak))
results.mat$EIP.avg <- as.numeric(as.character(results.mat$EIP.avg))
results.mat$BR.avg <- as.numeric(as.character(results.mat$BR.avg))
results.mat$VHR.avg <- as.numeric(as.character(results.mat$VHR.avg))
results.mat$EIP.sum <- as.numeric(as.character(results.mat$EIP.sum))
results.mat$BR.sum <- as.numeric(as.character(results.mat$BR.sum))
results.mat$VHR.sum <- as.numeric(as.character(results.mat$VHR.sum))
results.mat$MIR <- as.numeric(as.character(results.mat$MIR))

# The above is slow so there's an intermediate save option here
#save(results.mat, file="C:/Users/davewi/Documents/Results/WNV_allresults_withlatlon.Rdata")

# Load packages required for plot
library(ggplot2)
library(colorspace)
library(gridExtra)
library(cowplot)
library(RColorBrewer)

# Set up colours to match colour palette from maps
colors <- brewer.pal(9, "Reds") # for risk maps, data maps
pal <- colorRampPalette(colors)

# Subset to only plot certain years and order so that higher MIRs are plotted last to avoid being masked
# Not ideal but I couldn't find a better solution to the overplotting issues
results.for.plotting <- results.mat
results.for.plotting <- subset(results.for.plotting, Year %in% c("2049","2059","2069","2079"))
results.for.plotting <- results.for.plotting[order(results.for.plotting$MIR),]

# Make Figure 3
p1 <- ggplot(results.for.plotting, aes(x=VHR.avg, y=BR.avg)) + 
  geom_point(aes(color=MIR), size=0.9, alpha=0.8) +
  scale_color_gradient(low=pal(20)[1], high=pal(20)[20]) + #"#FFF5F0", high="#67000D") +
  theme_bw() +
  theme(panel.border = element_blank(), 
        axis.line = element_line(),
        legend.title = element_text(size=12),
        axis.title = element_text(size=12)) +
  facet_wrap(~Year, nrow=1) +
  ylab("Biting rate") +
  xlab("Vector-host ratio")

MIR.leg <- get_legend(p1)
p1 <- p1 + theme(legend.position = "none")

p2 <- ggplot(results.for.plotting, aes(x=VHR.avg, y=1/EIP.avg)) + 
  geom_point(aes(color=MIR), size=0.9, alpha=0.8) +
  scale_color_gradient(low=pal(20)[1], high=pal(20)[20]) + #"#FFF5F0", high="#67000D") +
  theme_bw() +
  theme(panel.border = element_blank(), 
        axis.line = element_line(),
        axis.title = element_text(size=12),
        legend.position = "none",
        strip.background.x = element_blank(),
        strip.text.x = element_blank()) +
  facet_wrap(~Year, nrow=1) +
  ylab("Viral replication rate") +
  xlab("Vector-host ratio")

p3 <- ggplot(results.for.plotting, aes(x=BR.avg, y=1/EIP.avg)) + 
  geom_point(aes(color=MIR), size=0.9, alpha=0.8) +
  scale_color_gradient(low=pal(20)[1], high=pal(20)[20]) + #"#FFF5F0", high="#67000D") +
  theme_bw() +
  theme(panel.border = element_blank(), 
        axis.line = element_line(),
        axis.title = element_text(size=12),
        legend.position = "none",
        strip.background.x = element_blank(),
        strip.text.x = element_blank()) +
  facet_wrap(~Year, nrow=1) +
  ylab("Viral replciation rate") +
  xlab("Biting rate") +
  scale_x_continuous(breaks=c(0.00, 0.03,0.06,0.09))

grid.arrange(p1,p2,p3,MIR.leg, nrow=3, ncol=2, widths=c(10,1), layout_matrix=rbind(c(1,4),c(2,4),c(3,4)))
