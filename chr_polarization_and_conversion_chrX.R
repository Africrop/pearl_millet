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

#### Code pour extraire SNPs avec etat ancestral verifie

#### Define chr under study
args <- commandArgs(TRUE)
i <- args[1]

#### First we extract information about ancestral information on SNP ####
ancestral_state <- read.table(paste("AncestralSNP_state_FINAL_chr",i,".txt",sep=""))
#ancestral_state <- read.table("~/data/mil_polym/AncestralSNP_state.txt")

colnames(ancestral_state) <- c("Chromosome","Position","References","Alternative","Ped1","Ped2","Ancestral")
ancestral_state$Polym_Id <- paste(as.character(ancestral_state$Chromosome),ancestral_state$Position,sep="_")


#### We then our analysis chr by chr ####
  ##### Set a counter #####
  j = 0
  
  ##### definition of links to be used #####
  inputFile <- paste("chr",i,"_sample_215_of_allSNP_chr",i,"_inclRef.txt.gz",sep="")
  #inputFile <- "~/data/mil_polym/chr1_sample_215_of_allSNP_chr1_inclRef.tar.gz"
  outputFile1<- paste("chr",i,"_215samples_polarized_SNPs_Ref_and_Anc.txt",sep="")
  outputFile2<- paste("chr",i,"_215samples_polarized_SNPs_data.txt",sep="")
  # outputFile3<- paste("/scratch/cubry/chr",i,"_215samples_polarized_SNPs_data.geno",sep="")
  #outputFile<- "~/Mil_data/genotypes_215/chr1_sample_215_of_polarizedSNP_chr1.txt"
  #outputFile<- "~/data/mil_polym/chr1_sample_215_of_polarizedSNP_chr1.txt"
  
  ##### open the connection to the consider chromosome #####
  chr  <- file(inputFile, open = "r")
  
  ##### extract the title line and create two files, one which will #####
  ##### contain the data for the 215 genotypes, the other with the  #####
  ##### reference (as ancestral) and derived base                   #####
  title_line <- readLines(chr,n=1,warn=FALSE)
  #    writeLines(title_line[1:3], con = outputFile,sep = "\n")
  title_line.splt <- unlist(strsplit(title_line,split=" "))
  write.table(t(c((title_line.splt[1:3]),"Ancestral")), outputFile1, append = T,col.names = F,row.names = F,quote=F)
  write.table(t(title_line.splt[-3]), outputFile2, append = T,col.names = F,row.names = F,quote=F)
  
  ##### now we read the file line by line and extract SNPs that can be polarized ####
  while (length(oneLine <- readLines(chr, n = 1, warn = FALSE)) > 0) {
    ###### Increment the counter and print it ######
    j = j + 1
    print(paste("Reading line ",j," of chr ",i,sep=""))
    
    ###### splitting the raw data ######
    oneLine<- unlist(strsplit(oneLine, split="\t"))
    oneLine<-unlist(strsplit(oneLine,split=" "))
    
    ###### test that there is no "N" call in the line, that an ancestral state ######
    ###### can be given at the position and that the reference is consistent   ######
    ###### between "ancestral" file and considered data                        ######
    if(length(which(oneLine=="N"))==0 &
       length(which(ancestral_state$Polym_Id==paste(oneLine[1],oneLine[2],sep="_")))!=0){
       if(ancestral_state[which(ancestral_state$Polym_Id==paste(oneLine[1],oneLine[2],sep="_")),]$References == oneLine[3]) {
      
        ####### If these conditions are fullfilled, store ancestral state in a variable #######
        anc <- ancestral_state[which(ancestral_state$Polym_Id==paste(oneLine[1],oneLine[2],sep="_")),]$Ancestral
        
        ####### write the reference and the ancestral state in a file #######
        write.table(t(c((oneLine[1:3]),as.character(anc))), outputFile1, append = T,col.names = F,row.names = F,quote=F)
        
        ####### And proceed to the data transformation #######
          ######## First simplify to homozygote by replacing ambiguous bases ######
          oneLine[oneLine == "R"] = sample(c("A","G"), sum(oneLine[-c(1,2)] == "R"),replace=T)
          oneLine[oneLine == "Y"] = sample(c("C","T"), sum(oneLine[-c(1,2)] == "Y"),replace=T)
          oneLine[oneLine == "W"] = sample(c("A","T"), sum(oneLine[-c(1,2)] == "W"),replace=T)
          oneLine[oneLine == "S"] = sample(c("C","G"), sum(oneLine[-c(1,2)] == "S"),replace=T)
          oneLine[oneLine == "M"] = sample(c("A","C"), sum(oneLine[-c(1,2)] == "M"),replace=T)
          oneLine[oneLine == "K"] = sample(c("T","G"), sum(oneLine[-c(1,2)] == "K"),replace=T)
          oneLine[oneLine == "U"] = "T"
          ######## Convert missing data to 9 ########
          oneLine[oneLine == "-"] = 9
          ######## Defining the reference allele as the ancestral one and convert to O (derived) and 1 (ancestral) ########
          oneLine[oneLine!=9][-c(1,2)]= as.numeric(oneLine[oneLine!=9][-c(1,2)] == anc)
          ######## Write data to file ########
          write.table(t(oneLine[-3]), outputFile2, append = T,col.names = F,row.names = F,quote=F)
          # write.table(paste(oneLine[-c(1:3)],collapse = ""), outputFile3, append = T,col.names = F,row.names = F,quote=F)
          
    }
   } 
  }
  
  close(chr) 
  

