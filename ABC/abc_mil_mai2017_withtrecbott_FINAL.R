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
nat.earth <- stack('/Users/cubry/Documents/GIS/naturalearth/HYP_50M_SR_W.tif')
lakes <- readOGR("/Users/cubry/Documents/GIS/naturalearth",layer = "ne_10m_lakes")
rivers <- readOGR("/Users/cubry/Documents/GIS/naturalearth",layer = "ne_50m_rivers_lake_centerlines")
playas <- readOGR("/Users/cubry/Documents/GIS/naturalearth",layer = "ne_50m_playas")
ocean <- readOGR("/Users/cubry/Documents/GIS/naturalearth",layer = "ne_10m_ocean")
shapefile <- readOGR("/Users/cubry/Documents/GIS/naturalearth",layer = "ne_50m_admin_0_countries")
data <- fortify(shapefile)
data.lakes <- fortify(lakes)
data.playas <- fortify(playas)
data.ocean <- fortify(ocean)
data.rivers <- fortify(rivers)


densities_plot <- function(prior = prior, post = post, title = "") {ggplot()  +
    #Plot the density graphs
    geom_density(aes(prior,fill="priors"),alpha=0.5) +
    geom_density(aes(post,fill="posteriors"),alpha=0.5) +
    #guide
    scale_fill_manual(values=c('#E69F00','#999999')) +
    #labels
    labs(title= title,x="") +
    #Define the theme  
    #theme_bw() +
    theme(panel.grid=element_blank(),legend.title=element_blank())
}

posteriorposition_map <- function(posteriorposition = posteriorposition,title="",n=50,h=25) {
  ggplot()+  
    geom_polygon(aes(x = long, y = lat, group = group),
                 fill="lightblue", data = data.ocean, size = .3) +
    geom_polygon(aes(x = long, y = lat, group = group), data = data,
                 colour ='antiquewhite4', fill = "white", size = .3) +
    geom_polygon(aes(x = long, y = lat, group = group), data = data.lakes,
                 colour ='lightblue',fill="lightblue", size = .3) +
    stat_density2d(data = posteriorposition,
                   aes(x=Long1,y=Lat1,  fill = ..level.., alpha = ..level..),
                   size = 0.1, geom = 'polygon',n = n,h=h) +
    geom_path(aes(x = long, y = lat, group = group), data = data.rivers,
              colour ='lightblue') +
    geom_point(aes(x=0.2333333,y=16.80))+
    scale_fill_gradient(low = "yellow", high = "red", guide = FALSE) +
    scale_alpha(range = c(0, 0.25), guide = FALSE) +
    ggtitle(title) +
    theme_bw() +
    theme(axis.line = element_blank(),
          axis.title = element_blank(),
          legend.key = element_blank(),
          panel.border = element_rect(colour = "black",
                                      size = .5,
                                      linetype = "solid"
          )
    ) + coord_quickmap(xlim = c(-20,50),ylim =c(-0,35))
}



#### Data analysis ####

list <- list.files(path="~/Mil_Splatche/FINAL/",pattern = "noncorr.sum.stat.txt")
setwd("~/Mil_Splatche/FINAL/")
param1pop <- read.table(list[1], header = T)
for(i in 2:length(list)){
  param1pop = rbind(param1pop,
                             read.table(list[i], header = T))}

stat.obs = read.table("~/Mil_Splatche/stat.obs.maf5_6.txt",header = T)


param1pop_OLD <- read.table("C:/Users/cubry/Documents/Mil_Splatche/OLD/retrieve cluster/Splatche_mil_1pop_ne_tau/mil_1pop_ne_tausum.stat.txt", header = T)

stat.obs = read.table("~/Mil_Splatche/stat.obs.maf5_6.txt",header = T)

setwd("~/")

boxplot(param1pop[,(length(param1pop)-20):(length(param1pop))],las=2, main = "Simulated vs Observed sum stat\n 1 pop with admixture and missing data\n non-correlated allele frequencies, Large Geo Priors",xaxt="n")
axis(1, at= 1:21,c(paste("SFS",seq(1:7)),paste("RV",seq(1:14))), las=2)
points(1:21,stat.obs[stat.obs$Ref=="MAF_3",2:22],pch=3,col="red")
legend("topright",c("MAF version 3"),pch=3,col=c("red"))

boxplot(param1pop_OLD[,(length(param1pop_OLD)-20):(length(param1pop_OLD))],las=2, main = "Simulated vs Observed sum stat\n 1 pop",xaxt="n")
axis(1, at= 1:21,c(paste("SFS",seq(1:7)),paste("RV",seq(1:14))), las=2)
points(1:21,stat.obs[stat.obs$Ref=="MAF_3",2:22],pch=3,col="red")
legend("topright",c("MAF version 3"),pch=3,col=c("red"))

library(abc)

 sum.stats.simul <- rbind(
   param1pop[,(length(param1pop)-20):(length(param1pop))],
   param1pop_OLD[,(length(param1pop_OLD)-20):(length(param1pop_OLD))]
 )
 
 sum.stats.simul$Model <- c(rep("Admixture model",nrow(param1pop)),
                            rep("No admixture model",nrow(param1pop_OLD))
 )

 sum.stats.simul <- (na.omit(sum.stats.simul))

#gfitpca(target=stat.obs[1,2:22],sumstat=(sum.stats.simul[,][1:21]),index=sum.stats.simul$Model,xlim = c(-11,6),ylim=c(-5,25))
gfitpca(target=stat.obs[1,2:22],sumstat=(sum.stats.simul[,][1:21]),index=sum.stats.simul$Model,xlim = c(-15,15),ylim=c(-6,15))

pdf(file = "/Users/cubry/Documents/Mil_Splatche/gfit PCA plot for Mil Splatche ABC.pdf")
gfitpca(target=stat.obs[1,2:22],sumstat=(sum.stats.simul[,][1:21]),index=sum.stats.simul$Model,xlim = c(-15,15),ylim=c(-6,15))
dev.off()

# ABC

nnOld <- abc(target = stat.obs[1,2:22],sumstat = param1pop_OLD[,9:29],param = param1pop_OLD[,c(1:5)],method = "neuralnet",tol = 0.1,numnet = 20,sizenet=13)

nn1 <- abc(target = stat.obs[1,2:22],sumstat = param1pop[,16:36],
           param = param1pop[,c(1:8,11,13:15)],method = "neuralnet",
           tol = 0.005,numnet = 20,sizenet=7,hcorr = FALSE#,
           #transf = c("none","none","none","none","none","log","log","log","none","none","none","none")
           )

# map <- as.data.frame(nnOld$adj.values[,1:2])
# maps <- plot_grid(posteriorposition_map(map,n = c(300,100),h = 50)+
#                     geom_point(aes(median(nnOld$adj.values[,1]),median(nnOld$adj.values[,2])),colour="blue")+
#                     geom_vline(xintercept = quantile(nnOld$adj.values[,1],probs=0.1),colour="gray75")+
#                     geom_vline(xintercept = quantile(nnOld$adj.values[,1],probs=0.9),colour="gray75")+
#                     geom_hline(yintercept = quantile(nnOld$adj.values[,2],probs=0.1),colour="gray75")+
#                     geom_hline(yintercept = quantile(nnOld$adj.values[,2],probs=0.9),colour="gray75")+
#                     labs(title="Non admixture model"))
# maps
# 

alpha<- densities_plot(param1pop$alpha,nn1$adj.values[,"alpha"],title = "alpha")
beta<- densities_plot(param1pop$beta,nn1$adj.values[,"beta"],title = "beta")
gamma<- densities_plot(param1pop$gamma,nn1$adj.values[,"gamma"],title = "gamma")
long <- densities_plot(param1pop$Long1,nn1$adj.values[,"Long1"],title = "Longitude")
lat <- densities_plot(param1pop$Lat1,nn1$adj.values[,"Lat1"],title = "Latitude")
gen <- densities_plot(param1pop$generations,nn1$adj.values[,"generations"],title = "generations")
acc <- densities_plot(param1pop$accroissement,nn1$adj.values[,"accroissement"],title = "accroissement")
mig <- densities_plot(param1pop$migration,nn1$adj.values[,"migration"],title = "migration")
sizeBeforeExp <- densities_plot(param1pop$SizeBeforeExpansion,nn1$adj.values[,"SizeBeforeExpansion"],title = "SizeBeforeExpansion")
TimeBott <- densities_plot(param1pop$TimeOfBottleneck,nn1$adj.values[,"TimeOfBottleneck"],title = "TimeBott")
AncSize <- densities_plot(param1pop$AncestralSize,nn1$adj.values[,"AncestralSize"],title = "AncestralSize")
K <- densities_plot(param1pop$MainCarryingCapacity,nn1$adj.values[,"MainCarryingCapacity"],title = "MainCarryingCapacity")

prow<-plot_grid(alpha + theme(legend.position="none"),
                beta + theme(legend.position="none"),
                gamma + theme(legend.position="none"),
                lat + theme(legend.position="none"),
                long + theme(legend.position="none"),
                gen + theme(legend.position="none"),
                mig + theme(legend.position="none"),
                acc + theme(legend.position="none"),
                sizeBeforeExp + theme(legend.position="none"),
                TimeBott + theme(legend.position="none"),
                AncSize + theme(legend.position="none"),
                K + theme(legend.position="none")
)
grobs <- ggplotGrob(lat)$grobs
legend <- grobs[[which(sapply(grobs, function(x) x$name) == "guide-box")]]
plot_grid( prow, legend, rel_widths = c(3, .5))


map <- as.data.frame(nn1$adj.values[,1:2])
maps <- plot_grid(posteriorposition_map(map,n = c(150,50),h = 50)+
                    geom_point(aes(median(nn1$adj.values[,1]),median(nn1$adj.values[,2])),colour="blue")+
                    geom_vline(xintercept = quantile(nn1$adj.values[,1],probs=0.1),colour="gray75")+
                    geom_vline(xintercept = quantile(nn1$adj.values[,1],probs=0.9),colour="gray75")+
                    geom_hline(yintercept = quantile(nn1$adj.values[,2],probs=0.1),colour="gray75")+
                    geom_hline(yintercept = quantile(nn1$adj.values[,2],probs=0.9),colour="gray75")+
                    labs(title="Non-correlated allelic frequencies model"))
maps
