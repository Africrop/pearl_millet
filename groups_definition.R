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

wild_w <- read.table("groups/wild_w.txt")
wild_c <- read.table("groups/wild_c.txt")
wild_e <- read.table("groups/wild_e.txt")
cult_w <- read.table("groups/cult_w.txt")
cult_c <- read.table("groups/cult_c.txt")
cult_e <- read.table("groups/cult_e.txt")
cult_i <- read.table("groups/cult_i.txt")
cult_s <- read.table("groups/cult_s.txt")

all_ind <- rbind(wild_e,wild_c,wild_w,cult_s,cult_i,cult_e,cult_c,cult_w)
names(all_ind) <- "Id"
all_ind$Pop <- c(rep("wild_e",nrow(wild_e)),rep("wild_c",nrow(wild_c)),rep("wild_w",nrow(wild_w)),rep("cult_s",nrow(cult_s)),rep("cult_i",nrow(cult_i)),rep("cult_e",nrow(cult_e)),rep("cult_c",nrow(cult_c)),rep("cult_w",nrow(cult_w)))
write.table(all_ind,"groups/all_ind.txt")
