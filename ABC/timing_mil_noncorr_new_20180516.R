###################################################################################################################################
#
# Copyright 2018 IRD & Grenoble-Alpes University
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/> or
# write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston,
# MA 02110-1301, USA.
#
# You should have received a copy of the CeCILL-C license with this program.
#If not see <http://www.cecill.info/licences/Licence_CeCILL-C_V1-en.txt>
#
# Intellectual property belongs to IRD 
#
# Written by Philippe Cubry, Yves Vigouroux, Olivier François
#
###################################################################################################################################

#### Functions ####
library(ggplot2)
library(raster)
library(rgdal)
library(cowplot)
library(MASS)
library(coda)

nat.earth <- stack('/media/philippe/OSDisk/Users/cubry/Documents/GIS/naturalearth/HYP_50M_SR_W.tif')
setwd("/media/philippe/OSDisk/Users/cubry/Documents/")
lakes <- readOGR("/media/philippe/OSDisk/Users/cubry/Documents/GIS/naturalearth",layer = "ne_10m_lakes")
rivers <- readOGR("/media/philippe/OSDisk/Users/cubry/Documents/GIS/naturalearth",layer = "ne_50m_rivers_lake_centerlines")
playas <- readOGR("/media/philippe/OSDisk/Users/cubry/Documents/GIS/naturalearth",layer = "ne_50m_playas")
ocean <- readOGR("/media/philippe/OSDisk/Users/cubry/Documents/GIS/naturalearth",layer = "ne_10m_ocean")
shapefile <- readOGR("/media/philippe/OSDisk/Users/cubry/Documents/GIS/naturalearth",layer = "ne_50m_admin_0_countries")
data <- fortify(shapefile)
data.lakes <- fortify(lakes)
data.playas <- fortify(playas)
data.ocean <- fortify(ocean)
data.rivers <- fortify(rivers)

# Define extent and load raster for the map
extent <- c(-25,57,-5,40)
nat.crop <- crop(nat.earth, y=extent(extent))
rast.table <- data.frame(xyFromCell(nat.crop, 1:ncell(nat.crop)),
                         getValues(nat.crop/255))
rast.table$rgb <- with(rast.table, rgb(HYP_50M_SR_W.1,
                                       HYP_50M_SR_W.2,
                                       HYP_50M_SR_W.3,
                                      1))

rm(nat.earth,nat.crop)

# Rationale
# Posterior to the Splatche inference of the onset of the cultivated pearl millet expansion, we performed additional 5000 simulations were we recovered time of arrival in locations were we had real time indication for the presence of cultivated millet. The idea beyond this experiment was to try infer the timing for the onset of expansion.

### Calibration dates available UPDATED 16052018
#Following the indications found in bibliography by Christian Dupuy, we identified the following candidate times: 

#Country/Region|Approximative calibrated date (BC)  
#--------------|-----------------------  
#Mali|-2579 ? -2409  
# Mauritania|-2134 ? -1562  
# Congo|-392 ? -115 
# South Africa|300  
# India|-1800  
# Nigeria|-1616 ? -1506 
# Ghana|-2401 ? -1300 
# Cameroon|-756 ? -197  
# Burkina|-1191 ? -851
# Kenya|775-980
# Zanzibar|605-760
# Tanzania|640-770
# 
# Corresponding mean and sd
# Pop|mean|sd
# Mali|4494|42.5
# Mauritania|3848|143
# Congo|2253.5|69.25
# Nigeria|3561|27.5
# Ghana|3850.5|275.25
# Cameroon|2476.5|139.75
# Burkina|3021|85
# Kenya|1122.5|51.25
# Zanzibar|1317.5|38.75
# Tanzania|1295|32.5

### Simulations
map_archeo_retained <-  ggplot()+  
    geom_polygon(aes(x = long, y = lat, group = group),
                 fill="lightblue", data = data.ocean, size = .3) +
    geom_polygon(aes(x = long, y = lat, group = group), data = data,
                 colour ='antiquewhite4', fill = "white", size = .3) +
    geom_polygon(aes(x = long, y = lat, group = group), data = data,
                 colour ='antiquewhite4', fill = NA, size = .3) +
    geom_polygon(aes(x = long, y = lat, group = group), data = data.lakes,
                 colour ='lightblue',fill="lightblue", size = .3) +
    geom_path(aes(x = long, y = lat, group = group), data = data.rivers,
              colour ='lightblue') +
    geom_point(
      aes(x = 0.23,y = 16.8))+
      geom_point(
      aes(x = 14.18,y = 1.13))+
      geom_point(
      aes(x = 7.7,y = 9.62))+
      geom_point(
      aes(x = -0.78,y = 10.53))+
      geom_point(
      aes(x = 9.92,y = 2.9))+
      geom_point(
      aes(x = 39.38,y = -6.31))+


    geom_text(
      aes(x = 0.23,y = 16.8),label="2579-2409 BC",nudge_x=0.8,nudge_y=0,hjust="left")+
      geom_text(
      aes(x = 14.18,y = 1.13),label="392-115 BC",nudge_x=0.8,nudge_y=-1,hjust="left")+
      geom_text(
      aes(x = 7.7,y = 9.62),label="1616-1506 BC",nudge_x=0.8,nudge_y=-1,hjust="left")+
      geom_text(
      aes(x = -0.78,y = 10.53),label="2401-1300 BC",nudge_x=0.8,nudge_y=1,hjust="left")+
      geom_text(
      aes(x = 9.92,y = 2.9),label="756-197 BC",nudge_x=0.8,nudge_y=1,hjust="left")+
      geom_text(
      aes(x = 39.38,y = -6.31),label="605-760 AD", nudge_x=0.8, nudge_y=-1.5, hjust="left")+

  theme_bw() +
    theme(axis.line = element_blank(),
          axis.title = element_blank(),
          legend.key = element_blank(),
          panel.border = element_rect(colour = "black",
                                      size = .5,
                                      linetype = "solid"
          )
    ) + coord_quickmap(xlim = c(-20,55),ylim =c(-37,35))+
  ggtitle("Map of archeological remains retained for calibration")
map_archeo_retained

ggsave(filename = "map_archeo_retained_20180516.pdf",plot = map_archeo_retained,width = 9,height = 6,units = "in")

lst <- list.files("~/Mil_Splatche/FINAL/timing_new/",pattern = "*noncorr")[grep(list.files("~/Mil_Splatche/FINAL/timing_new/",pattern = "*noncorr"),pattern = "*.tar.gz")]
results <- NULL
lst.dir <- list.dirs("./")
for(d in lst.dir){
r.lst <- list.files(path = d)[grep(pattern = "*output.txt", x=list.files(path = d))]
for(i in r.lst){
  tmp <- scan(paste(d,"/",i,sep=""),
     skip=1,what = "character",sep = ":")
tmp.results <- tmp[seq(2,22,2)]
results <- cbind(results,tmp.results)
rownames(results) <- tmp[seq(1,21,2)]
}
  }

#Filtering for non-complete simulations
tst <- as.numeric(apply(results,2,min))==-1
results <- results[,which(tst==FALSE)]
#unlink(list.dirs(path = ".", full.names = TRUE,recursive = FALSE), recursive = TRUE)

# defining calibration dates
real <- rownames(results)
real <- function(x){c(
  rnorm(1,4494,42.5),
  rnorm(1,3848,143),
  rnorm(1,2253.5,69.25),
  1700,
  3800,
  rnorm(1,3561,27.5),
  rnorm(1,3850.5,275.25),
  rnorm(1,2476.5,139.75),
  rnorm(1,3021,85),
  rnorm(1,1122.5,51.25),
  rnorm(1,1317.5,38.75))}

### Linear regression
# We estimated the interception of the real time axis by a linear regression between simulated and real dates. This allowed us to infer a distribution for the timing of the begining of the expansion.  
# We performed several attempts for the linear regression, including or not some of the calibration dates.  
# After the linear regressions were made for each of the simulated datasets, we estimated and plotted the density of the estimated time for the onset of the expansion. We also extracted the mode of the distribution and reported it.

timing <- NULL
for(i in 1:ncol(results)){
  tmp <- real()
  timing <- cbind(timing,(lm(tmp[-c(2,4,5,9,10)]~as.numeric(results[-c(2,4,5,9,10),i])))$coefficients[1])
}
timing <- timing[,timing>1]

#Plotting
density_timing_retained <- ggplot() +
  geom_density(aes(as.numeric(timing)),
                        color="orange",fill="orange",alpha=0.25) +
  geom_vline(xintercept =density(timing)$x[which(density(timing)$y ==max(density(timing)$y))], color="red" ) +
  xlab("Estimated date for the onset of expansion") + xlim(0,NA)

timing.mode <- density(timing)$x[which(density(timing)$y ==max(density(timing)$y))]
#ggsave(filename = "density_timing_retained_20180613.pdf", density_timing_retained,width = 7,height = 4.84)
density_timing_retained

#The obtained mode was in this case 
timing.mode
# with a 95% confidence interval of
quantile(timing,probs=c(0.025,0.975))

#Note that in all cases this timing might be underestimated.