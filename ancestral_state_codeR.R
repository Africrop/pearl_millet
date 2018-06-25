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

#code for verification of the ancestral genotype in the vcf file
#output the chromosme, position, the reference, the ancestral state if different from the reference, the genotype of the two ancestral individuals, and the ancestral state inferred 

inputAnc <- "C:/Users/cubry/Documents/Mil_data/Ancestral/Calling.vcf"
Anc  <- file(inputAnc, open = "r")
outputFileAnc<- "C:/Users/cubry/Documents/Mil_data/Ancestral/AncestralSNP_state.txt"

i=1
j=1


while (length(oneLine <- readLines(Anc, n = 1, warn = FALSE)) > 0) {
  if (length(grep("#",oneLine))==0) {
#    print(oneLine)
	oneLineR<-unlist(strsplit(oneLine, split="\t"))

#	print(substr(oneLineR[11],1,3))
						#genotype should be identifical between individuals and not heterozygotes = ambiguity
if(substr(oneLineR[11],1,3)!="./."|substr(oneLineR[10],1,3)!="./."){		
		if (substr(oneLineR[11],1,3)=="./.") {sum<-(abs(as.numeric(substr(oneLineR[10],1,1))-as.numeric(substr(oneLineR[10],3,3)))) 
			} else {
				if (substr(oneLineR[10],1,3)=="./.") {sum<-(abs(as.numeric(substr(oneLineR[11],1,1))-as.numeric(substr(oneLineR[11],3,3))))
					} else {
					if (substr(oneLineR[10],1,3)==substr(oneLineR[11],1,3)) {sum<-abs(as.numeric(substr(oneLineR[10],1,1))-as.numeric(substr(oneLineR[10],3,3)))} else {sum<-1
						}
					}
			}

#	print(sum)
	if(sum==0) {	
		j<-j+1
		if ( (substr(oneLineR[11],1,3)=="1/1") || (substr(oneLineR[10],1,3)=="1/1") ) {ancestralstate<-oneLineR[5]} else {
																ancestralstate<-oneLineR[4]
																}
		text<-t(c(oneLineR[1],oneLineR[2],oneLineR[4],oneLineR[5],substr(oneLineR[10],1,3), substr(oneLineR[11],1,3),ancestralstate[1]))
		write.table(text, file = outputFileAnc, sep = " ", quote=FALSE, col.names=FALSE, row.names=FALSE, append=TRUE)		


		}
				
   	print(i)
	print(j)	
	i<-i+1
}	
  }
}
close(Anc) 
