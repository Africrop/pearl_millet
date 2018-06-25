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

#!/bin/sh

# Before running, ensure it is correctly formatted by using dos2unix script_name.sh
# Make it executable by chmod 755 script_name.sh

# Name: run_splatche_riz.sh
# commandline: qsub /run_splatche_V4.sh <nb_simulations> <nb_origins> <nom_sortie> <nom_dossier_sortie>



# ecrit les erreurs dans le fichier de sortie standard 
#$ -j y 

# shell que l'on veut utiliser 
#$ -S /bin/bash 

# indiquer son email pour suivre l'execution : 
#$ -M philippe.cubry@ird.fr

# obtenir un message au demarrage (b) , a la fin (e), en cas d'abandon (a) 
#$ -m bea 

# la queue que l'on veut utiliser : 
#$ -q bioinfo.q

#$ -N mil_splatche

###### Definition des variables de chemin

path_to_dir="/home/cubry/scripts/splatche_mil_witht_rec_bott";
path_to_tmp="/scratch/cubry-$JOB_ID-$SGE_TASK_ID" 

###### Creation du repertoire temporaire sur noeud

mkdir $path_to_tmp

scp -rp nas:/$path_to_dir/* $path_to_tmp/
echo "tranfert donnees master -> noeud";
cd $path_to_tmp

###### Execution du programme
if [ $# != 0 ]

	then
		Rscript $path_to_tmp/splatche_launcher_mil_new_witht_rec_bott.R $1 $2
	else	
		Rscript $path_to_tmp/splatche_launcher_mil_new_witht_rec_bott.R
fi

rm settings* datasets_1layer/dens_init* datasets_1layer/dyn* datasets_1layer/GeneS* datasets_1layer/genetic_* datasets_1layer/ppv* datasets_1layer/r* datasets_1layer/v* 
mv datasets_1layer/ $3_$JOB_ID-$SGE_TASK_ID/

##### Transfert des données du noeud vers master apres compression
tar -czvf $3_$JOB_ID-$SGE_TASK_ID.tar.gz $3_$JOB_ID-$SGE_TASK_ID/
rm -r $3_$JOB_ID-$SGE_TASK_ID/
scp -rp $path_to_tmp/* cubry@nas:$4

echo "Transfert donnees node -> master";


#### Suppression du repertoire tmp noeud

rm -rf $path_to_tmp

echo "Suppression des donnees sur le noeud";
