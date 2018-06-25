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
# Usage: qsub fsc_sh_mil_wild.sh <cultivated group> <model to test>
#$ -cwd
#$ -j y
#$ -S /bin/bash
#$ -M philippe.cubry@ird.fr
#$ -N fastsimcoal
#$ -m be
#$ -V
#$ -q bioinfo.q

group=$1
model=$2
obsPath="/data2/projects/africrop_models/mil_fastsimcoal/scenarii/wild-cult/$group/"
## check if option --removeZeroSFS

minL=40
maxL=150

minSim=250000
maxSim=1000000

path="/scratch/cubry-$JOB_ID/$model"

##############
### SCRIPT ###
##############
mkdir /scratch/cubry-$JOB_ID/
mkdir $path
scp -rp nas:$obsPath/* $path/
cd $path


echo "=========================================================================="
echo "RUN $run"
echo "SPLATCHE parameters"
echo "min/maxL $minL $maxL"
echo "min/maxSim $minSim $maxSim"
echo "Running analysis on model $model"
echo "=========================================================================="

fsc25221 -t multisfs_${model}.tpl -e multisfs_${model}.est -n $minSim -N $maxSim -M 0.001  -l $minL -L $maxL --multiSFS -d -B 1 -C 10
mv $path/multisfs_${model} $path/f_${model}_$JOB_ID

scp -r $path/f_${model}_$JOB_ID nas:/$obsPath/new_ter/
rm -r /scratch/cubry-$JOB_ID/
