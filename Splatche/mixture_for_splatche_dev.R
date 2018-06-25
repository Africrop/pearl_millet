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

#######################################################
# R script to compute set of parameters for admixture #
# implementation based on Splatche outputs            #
# January 18th 2016 - Version 1                       #
# Written by Philippe Cubry                           #
# Contributors O. Francois, Y. Vigouroux              #
#######################################################

#### Function to calculate parameters for Beta distribution ####
#' Compute alpha and beta parameters for a Beta distribution given mean and
#' variance parameters
#'
#' @param u the mean parameter
#' @param v the variance parameter
#'
#' @return a vector of two values, alpha and beta parameters
#'
#' @examples >beta_parameters(0.15,0.01)
#' [1] 1.7625 9.9875
beta_parameters <- function(u,v){
  a <- u*((u*(1-u)/v)-1)
  b <- (1-u)*((u*(1-u)/v)-1)
  c(a,b)}

#### Function to sample admix coef (probability of gene flow with given mean and variance) ####
#' sample an admixture coefficient from a Beta distribution given a mean and
#' a variance
#'
#' @param admix_coef_mean the mean of the considered admixture coefficient
#' @param admix_coef_var the variance of the considered admixture coefficient
#'
#' @return a value to be used as admixture probability
#'
#' @examples > sample_admix_coef(0.05,0.01)
#' [1] 0.01477277
sample_admix_coef <- function(
  admix_coef_mean = 0.15,
  admix_coef_var = 0.01)
  {
  rbeta(1, beta_parameters(admix_coef_mean,admix_coef_var)[1],beta_parameters(admix_coef_mean,admix_coef_var)[2])
  }

