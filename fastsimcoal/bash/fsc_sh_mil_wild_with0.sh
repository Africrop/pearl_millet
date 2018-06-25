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

#!/bin/bash

# Adapted from script by Remi Tournebize
# Usage: qsub fsc_sh_mil_wild.sh <model to test>
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -M philippe.cubry@ird.fr
#$ -N fastsimcoal
#$ -m be
#$ -V
#$ -q bioinfo.q

model=$1
obsPath="/data2/projects/africrop_models/mil_fastsimcoal/scenarii/$model"
## check if option --removeZeroSFS


minL=10
maxL=100

minSim=500000
maxSim=1000000

path="/scratch/cubry/$model"

##############
### SCRIPT ###
##############
mkdir /scratch/cubry
mkdir $path
scp -rp nas:$obsPath/* $path/
cd $path

for run in $(seq $startRun $endRun)
do
	echo "=========================================================================="
	echo "Model $model"
	echo "Parameters MinL $minL MaxL $maxL with at most $maxSim simulations"
	echo "=========================================================================="
	fsc25221 -t multisfs_wild_${model}.tpl -e multisfs_wild_${model}.est -n $minSim -N $maxSim -M 0.01 -l $minL -L $maxL --multiSFS -d
	mv $path/multisfs_wild_${model} $path/f_${model}_${run}
done

rcp -rp $path/ nas:/$obsPath/
rm -r /scratch/cubry
