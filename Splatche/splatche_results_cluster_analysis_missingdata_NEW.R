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

#########################################################
# R script for Splatche outputs ABC analysis on cluster #
# will produce summarized statistics of the simulations #
# including binned SFS and rare variants per group      #
# January 14th 2016                                     #
# Auth: P. Cubry based on code from O. Francois         #
#########################################################

## Caution: only useful for cultivated pearl millet data with 146 individuals,
## to apply to other datasets modify SFS binning and maf functions
## requires file group.txt containing groups definition
## require script "mixture_for_splatche_dev.R"
## require files "freq_wild.sample10000.txt" "cult.sample.miss.freq"
## Usage : run using qsub from a launching script providing project name (e.g. mil_1pop) and output file name

#### Gather arguments from script call ####
args <- commandArgs(TRUE)

#### Source R script implementing functions to make admixture
source("./mixture_for_splatche_dev.R")

#### Define useful functions to calculate summarized statistics ####

##### Function to convert *.arp results in genotype object #####
myprocess1 = function(
  input= string,
  Nbsample = 146,
  Nbloci = 100,
  Nbindiv = rep(2, 146))
{
  ## process DNA file
  X <- scan(file = input, what = character(), sep = "\t", quiet = TRUE, skip = 0, nlines = 0, comment.char = "#")
  Totindiv <- sum(Nbindiv)
  lx <- length(X)
  X <- X[30:lx]
  M <- data.frame(matrix(NA, nrow = Totindiv, ncol = (2 + Nbloci)))
  for(k in 1:Nbsample){
    shift <- 1
    if (k == 1) k.s <- 0 else k.s <- sum(Nbindiv[1:(k-1)])
    for (i in 1:(Nbindiv[k])){
      X[shift + (i-1)*(2+Nbloci)] <- k
      M[k.s + i, ] <- X[ ( shift + (i-1)*(2+Nbloci) ):(shift -1 +i*(2+Nbloci) )]
      shift <- shift + 1
    }
    lx <- length(X)
    X<-X[(10 + shift+i*(2+Nbloci)):lx]
  }
  return(M[,-2])
}

##### Function to calculate minimal allele frequency #####
maf = function(x){ 
  x = x[x != 9]
  f = sum(x) ; (min(f , length(x) - f ))}

##### Function to bin SFS #####
histo.bin = function(x){
  n = length(x)
  y = sapply(1:73, FUN = function(i) sum(x == i) )
  c(y[1], y[2], y[3], sum(y[4:5]), sum(y[6:10]), sum(y[11:20]), sum(y[20:73]) )/n
}

##### Function to compute mean per group #####
z.groupm = function(z){
  sapply(1:14, FUN = function(i) mean(z[group == i]) )
}

#### Load group definition ####
group = scan("./group.txt")

#### Variables definition ####
sum.stat = NULL
param = NULL

#### Loading simulations parameters ####
print("loading parameters")
job_id <- gsub(".tar.gz","",gsub(paste(args[1],"_",sep=""),"",list.files(pattern = args[1])))
print("job_id =")
print(job_id)
for (i in 1:length(job_id)){
  param.simu <- read.table(file = paste("param_",job_id[i],".txt",sep=""),sep="",header=T,colClasses="character")
  param.simu$job_id <- job_id[i]
  param = rbind(param,param.simu)
  rm(param.simu, i)
}

#### Loop for creating summarized statistics from simulated datasets, adding mixture ####
for (n in 1:nrow(param)){
  print(paste("loop",n, "of",nrow(param)))
  ##### construction d'un objet contenant les snps #####
  print("extracting files from archive")
  untar(tarfile=paste(args[1],"_",param["job_id"][n,],".tar.gz",sep=""),
        files = paste(args[1],"_",param["job_id"][n,],"/GeneticsOutput/settings",param$generations[n],"_",
                      param$accroissement[n],"_",param$migration[n],"_",param$SizeBeforeExpansion[n],"_",
                      param$TimeOfBottleneck[n],"_",param$AncestralSize[n],"_",param$TimeForRecentBott[n],"_",
                      param$MutationRate[n],"_",param$MainCarryingCapacity[n],"_GeneSamples_",seq(1,10,1),".arp" ,sep=""),
        list = FALSE, exdir = ".",compressed = "gzip", verbose = FALSE)
  
  
  string.snp = paste("./",args[1],"_",param["job_id"][n,],"/GeneticsOutput/settings",param$generations[n],"_",
                     param$accroissement[n],"_",param$migration[n],"_",param$SizeBeforeExpansion[n],"_",
                     param$TimeOfBottleneck[n],"_",param$AncestralSize[n],"_",param$TimeForRecentBott[n],"_",
                     param$MutationRate[n],"_",param$MainCarryingCapacity[n],"_GeneSamples_1.arp" ,sep="")
  genotype = myprocess1(string.snp, Nbsample = 146, Nbloci = 500, Nbindiv = rep(2, 146) )[,-1]
  
  for (r in 2:10){
    string1 = paste("./",args[1],"_",param["job_id"][n,],"/GeneticsOutput/settings",param$generations[n],"_",
                    param$accroissement[n],"_",param$migration[n],"_",param$SizeBeforeExpansion[n],"_",
                    param$TimeOfBottleneck[n],"_",param$AncestralSize[n],"_",param$TimeForRecentBott[n],"_",
                    param$MutationRate[n],"_",param$MainCarryingCapacity[n],"_GeneSamples_",r,".arp" ,sep="")
    genotype = cbind(genotype, myprocess1(string1, Nbsample = 146, Nbloci = 500, Nbindiv = rep(2, 146) )[,-1])
  }
  genotype1 = genotype[ (1:292)%%2 == 1, ]
  genotype2 = genotype[ (1:292)%%2 == 0, ]
  genotype = cbind(genotype1, genotype2)
  
  lst = which(apply(genotype, MARGIN = 2, FUN = function(x) length(unique(x)) ) == 1) #liste les homozygotes
  if(length(lst)!=0){genotype = genotype[,-lst]}
  
  genotype = apply(genotype, MARGIN = 2, FUN = function (x) {as.numeric(as.factor(x))-1})
  
  print("removing temp archive directory")
  unlink(x = paste("./",args[1],"_",param["job_id"][n,],sep=""),recursive = TRUE)
  
  # Go to the next step in the loop if there is some incongruent data
  if(length(which(is.na(genotype)==TRUE))) next
  
  # Adding wild and missing frequencies to genotype object
  freq.wild <- t(read.table("./freq_wild.sample10000_moins_PE08743.txt",header=TRUE))
  cult.miss.freq <- t(read.table("./cult.sample.miss.freq.txt"))
  row.names(genotype) <- group
  freq.wild <- rbind(freq.wild,cult.miss.freq)
  genotype <- rbind(genotype,freq.wild[,1:ncol(genotype)])
  rownames(genotype)[nrow(genotype)] <- "missing.freq"
  
  #### Loop implementing admixture and computing summary statistics ####
  for(r in 1:10){
    # Implementing admixture
    genotype.rempla <- genotype
    
    # Replacement in Western genotypes
    alpha <- sample_admix_coef(0.07,0.015) ; names(alpha)<-"alpha"
    genotype.rempla[row.names(genotype.rempla)==2,] <-
      apply(genotype.rempla[row.names(genotype.rempla)==2|row.names(genotype.rempla)=="wild.west.sample.freq",],
            MARGIN=2,
            FUN = function(x){
              freq.rempla <- x[length(x)]
              x <- x[-length(x)]
              sapply(x,FUN=function(x){if(rbinom(1,1,alpha)==1){ x = rbinom(1,1,freq.rempla)} else x = x})
            })
    
    # Replacement in Eastern genotypes
    beta <- 0 ; names(beta) <- "beta"
       beta <- sample_admix_coef(0.07,0.015) ; names(beta)<-"beta"
       genotype.rempla[row.names(genotype.rempla)==12|row.names(genotype.rempla)==10,] <-
         apply(genotype.rempla[row.names(genotype.rempla)==12|row.names(genotype.rempla)==10|row.names(genotype.rempla)=="wild.east.sample.freq",],
               MARGIN=2,
               FUN = function(x){
                 freq.rempla <- x[length(x)]
                 x <- x[-length(x)]
                 sapply(x,FUN=function(x){if(rbinom(1,1,beta)==1){ x = rbinom(1,1,freq.rempla)} else x = x})
               })

    # Replacement in Center genotypes
    gamma <- 0 ; names(gamma) <- "gamma"
       gamma <- sample_admix_coef(0.07,0.015) ; names(gamma)<-"gamma"
       genotype.rempla[row.names(genotype.rempla)==7|
                         row.names(genotype.rempla)==11|
                         row.names(genotype.rempla)==6|
                         row.names(genotype.rempla)==8|
                         row.names(genotype.rempla)==4|
                         row.names(genotype.rempla)==5,] <-
         apply(genotype.rempla[row.names(genotype.rempla)==7|
                                 row.names(genotype.rempla)==11|
                                 row.names(genotype.rempla)==6|
                                 row.names(genotype.rempla)==8|
                                 row.names(genotype.rempla)==4|
                                 row.names(genotype.rempla)==5|
                                 row.names(genotype.rempla)=="wild.center.sample.freq",],
               MARGIN=2,
               FUN = function(x){
                 freq.rempla <- x[length(x)]
                 x <- x[-length(x)]
                 sapply(x,FUN=function(x){if(rbinom(1,1,gamma)==1){ x = rbinom(1,1,freq.rempla)} else x = x})
               })
 
       # Adding missing data
       genotype.rempla[1:length(group),] <-
         apply(genotype.rempla[row.names(genotype.rempla)%in%(1:14)|row.names(genotype.rempla)=="missing.freq",],
               MARGIN=2,
               FUN = function(x){
                 freq.rempla <- x[length(x)]
                 x <- x[-length(x)]
                 sapply(x,FUN=function(x){if(rbinom(1,1,freq.rempla)==1){x = 9} else x = x})
               })
       
       
      genotype.rempla <- genotype.rempla[1:146,]
      
      # Computations of summary statistics  
      print("computing summary statistics")
      spect = apply(genotype.rempla, MARGIN = 2, FUN = function(x) {maf(x)}   ) #cree le spectre de frequence
      
      lst1 = which(spect == 1) # liste les singletons dans le spectre
      geno1 = genotype.rempla[,lst1]
      ind1 = apply(geno1, MARGIN = 2, FUN = function(x) {
        y = x[x != 9] ; i1 = which(x != 9)[1] ;
        i = which( (x != y[1]) & (x != 9) ) ;
        if (length(i) > 1) i1 else i 
      })
      
      z = sapply(1:146, FUN = function(x) sum(ind1 == x) )
      
      stat.sim <- c(histo.bin(spect), z.groupm(z)/sum(z.groupm(z)) ) #compile les stats resumees
      names(stat.sim) <- c(paste("SFS",seq(1,7,1), sep=""),paste("RareVariants",seq(1,14,1), sep=""))
      stat.temp <- as.data.frame(c(param[n,],alpha,beta,gamma,stat.sim))
      
      sum.stat <- rbind(sum.stat, stat.temp )
      cond <- c(rep("SFS",7),rep("RareVariants",14))
  }
print("saving computed statistics to file")
write.table(file = paste(args[2],"sum.stat.txt",sep=""), sum.stat, row.names = F, quote=F)
}
print("saving computed statistics to file")
write.table(file = paste(args[2],"sum.stat.txt",sep=""), sum.stat, row.names = F, quote=F)
