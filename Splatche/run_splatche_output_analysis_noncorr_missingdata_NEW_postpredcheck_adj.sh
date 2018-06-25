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

# Name: run_splatche_output_analysis.sh
# commandline: qsub /run_splatche_output_analysis_VFreq_west_east_center_missingdata.sh <project name> <output_file_name> <data directory>

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

#$ -V

#$ -N splatche

###### Definition des variables de chemin

path_to_dir="/home/cubry/scripts/splatche_mil_witht_rec_bott/";
path_to_data=$3
path_to_tmp="/scratch/splatche_analysis-$JOB_ID" 

###### Creation du repertoire temporaire sur noeud

mkdir $path_to_tmp

scp -rp nas:/$path_to_data/*.tar.gz $path_to_tmp/
scp -rp nas:/$path_to_dir/group.txt $path_to_tmp/
scp -rp nas:/$path_to_dir/freq_wild_moins_PE08743.txt $path_to_tmp/
scp -rp nas:/$path_to_dir/cult.miss.freq.txt $path_to_tmp/
scp -rp nas:/$path_to_dir/sfs_cult.txt $path_to_tmp/
scp -rp nas:/$path_to_dir/splatche_results_cluster_analysis_missingdata_NEW_postpredcheck_adj.R $path_to_tmp/
scp -rp nas:/$path_to_dir/nn3_tol005_default $path_to_tmp/
scp -rp nas:/$path_to_dir/admix_posteriors_corrfreq_FINAL.txt $path_to_tmp/
scp -rp nas:/$path_to_dir/freq_wild.sample10000_moins_PE08743.txt $path_to_tmp/
scp -rp nas:/$path_to_dir/cult.sample.miss.freq.txt $path_to_tmp/
scp -rp nas:/$path_to_data/param* $path_to_tmp/
echo "tranfert donnees master -> noeud";
cd $path_to_tmp

###### Execution du programme
if [ $# != 0 ]

	then
		/usr/local/bin/Rscript $path_to_tmp/splatche_results_cluster_analysis_missingdata_NEW_postpredcheck_adj.R $1 $2
fi

##### Transfert des données du noeud vers master apres effacement des dossiers compresses
rm -rf $path_to_tmp/*.tar.gz
rcp -rp $path_to_tmp/ nas:/$path_to_dir/ppc/


echo "Transfert donnees node -> master";

#### Suppression du repertoire tmp noeud

rm -rf $path_to_tmp/

echo "Suppression des donnees sur le noeud";
