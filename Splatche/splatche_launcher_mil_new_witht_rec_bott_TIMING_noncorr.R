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

####################################################
# Script de lancement de Splatche sur le Cluster   #
# Auteur : Philippe Cubry                          #
# Contributeur : Olivier Francois                  #
# Date : 23 mai 2016                               #
# Version 5                                        #
####################################################

#### ReadMe ####
# script lancant Splatche sur le cluster
# Ecrit par Philippe Cubry sur la base de codes d'Olivier Francois
# Usage : splatche_launcher.R <nombre de runs> <nombre de populations sources>
# Requis pour Splatche: fichiers settings_XXX.txt", "dens_init.txt", contenu du dossier "dataset_1_layer"



#### getting Job_Id ####
job <- paste(Sys.getenv("JOB_ID"),"-",Sys.getenv("SGE_TASK_ID"),sep="")

#### getting call arguments ####
args <- commandArgs(TRUE)

#### source useful file ####
source("./splatche_input_files_creation_mil_new_witht_rec_bott_TIMING.R")

#### Loading posterior distributions ####
load("./nn3_tol005_default")

#### Definition du nombre de runs et de populations sources si specifie dans l'appel du script (defaut a 1 et 1) ####
if (is.na(args[1])==TRUE){
  nrun <- 1 #nombre de runs
  norigin <- 1 #nombre de populations sources
} else {nrun<-as.numeric(args[1]); norigin<-as.numeric(args[2])} #valeurs specifiees dans l'appel du script

#### Boucle sur le nombre de runs ####
for (it in 1:nrun){
  print(paste(it,"over",nrun,sep=" ")) #compteur de run

  ##### tirage des coordonnees des populations sources #####
  coordo <- NULL #Initialisation d'une variable de stockage des coordonnees des populations sources

  for (ori in 1:norigin){
    ox = sample(x = nn3$adj.values[,"Long1"],size = 1) #tirage de la longitude
    oy = sample(x = nn3$adj.values[,"Lat1"],size = 1) #tirage de la latitude

    while(!is.african(cbind(ox,oy))){ #test des coordonnees pour verifier qu'on est pas dans l'eau, sinon on retire
    ox = sample(x = nn3$adj.values[,"Long1"],size = 1) #tirage de la longitude
    oy = sample(x = nn3$adj.values[,"Lat1"],size = 1) #tirage de la latitude
    }
    coordo = rbind(coordo,c(ox,oy)) #stockage des coordonnees tirees
  }

  ##### tirage des variables #####
  ng <- sample(x = nn3$adj.values[,"generations"],size = 1) #tirage du nombre de generations a simuler
  rate <- sample(x = nn3$adj.values[,"accroissement"],size = 1) #tirage du taux d'accroissement des populations
  mig <- sample(x = nn3$adj.values[,"migration"],size = 1) #tirage du taux de migration
  res <- sample(x = nn3$adj.values[,"SizeBeforeExpansion"],size = 1) #tirage de la taille avant expansion
  Tres <- sample(x = nn3$adj.values[,"TimeOfBottleneck"],size = 1) #tirage de la taille avant expansion
  resB <- sample(x = nn3$adj.values[,"AncestralSize"],size = 1) #tirage de la taille avant expansion
  MutRate<- sample(MutRate_prior,1)
  K <- sample(x = nn3$adj.values[,"MainCarryingCapacity"],size = 1)

  TrecBott <- 0

  splatche(input = "settings.txt", ng, rate, mig, norigin,res,Tres,resB,TrecBott, MutRate, K, coord.o = coordo) #appel de la fonction Splatche
  system(paste("splatche2-01 settings",ng,"_",rate,"_",mig,"_",res,"_",
               Tres,"_",resB,"_",TrecBott,"_",MutRate,"_",K, ".txt", sep="")) #lancement de splatche avec les parametres tires

  ##### stockage des parametres uniquement si Splatche a genere des fichiers de sortie #####
  if(file.exists(paste("./datasets_1layer/GeneticsOutput/settings",ng,"_",rate,"_",mig,"_",res,"_",
                       Tres,"_",resB,"_",TrecBott,"_",MutRate,"_",K,"_GeneSamples_2.arp" ,sep=""))){

    coo1 <- NULL ; noms <- NULL
    for (npop in 1:norigin){
      coo1<-cbind(coo1,coordo[npop,1],coordo[npop,2])
      noms <-cbind(noms,paste("Long",npop,sep=""),paste("Lat",npop,sep=""))
    }
    param<- rbind(param,c(coo1,ng,rate,mig,res,Tres,resB,TrecBott,MutRate,K)) #stockage des parametres dans un fichier
    colnames(param) <- c(noms,"generations","accroissement","migration","SizeBeforeExpansion",
                         "TimeOfBottleneck","AncestralSize","TimeForRecentBott","MutationRate","MainCarryingCapacity")
  }
  ##### Sauvegarde des parametres des runs realises dans un fichier #####
  write.table(file = paste("param_",job,".txt",sep=""), param, row.names = F, quote=F)
}
