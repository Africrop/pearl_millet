### Set of functions to use to estimate SFS for mil pops
### Warning : only usefull with O/1 coded lfmm-style matrix
### Adaptation might be needed to use with other kind of data
### Especially regarding the way MAF/DAF is estimated

#' Function to compute the number of missing data in a sample
#'
#' @param x # the sample to investigate
#'
#' @return number of missing data coded as "9" in the sample
#'
missing_data_computation <- function(x){
  length(which(x==9))
}


#' Function to compute MAF wih regards to missing data
#'
#' @param x # the matrix with data to use
#'
#' @return the MAF (estimated frequency) for the considered data
#'
maf4 = function(x){
  x = x[x != 9]
  f = sum(x) ; (min(f , length(x) - f ))/length(x)
  }

#' Function to compute MAF wih regards to missing data.
#' CAUTION this statistics is intended to compute MAF within
#' the sample under investigation without regards to reference
#' allele.
#'
#' @param x # the matrix with data to use
#'
#' @return the MAF as allele counts for the considered data
#'
maf3 = function(x){
  x = x[x != 9]
  f = sum(x) ; min(f , length(x) - f ) }


#' Function to compute MAF based on global reference allele wih regards to missing data
#'
#' @param x # the matrix with data to use
#'
#' @return the MAF estimated from present data as frequencies for the considered data
#'
maf.x = function(x){
  x = x[x != 9]
  f = (length(x) - sum(x))/length(x) ; f}

#' Function to compute SFS based on mafx results
#'
#' @param x # the matrix with data to use
#'
#' @return the SFS with MAF as allele counts for the considered data
#'
sfs.x = function(x,N=nrow(x)){ y = apply(x,2,maf.x);
                              sapply(y,function(y){rbinom(1,as.integer(0.75*N),y)})}

#' Function to compute MAF based on minimal allele freq wih regards to missing data
#'
#' @param x # the matrix with data to use
#'
#' @return the MAF estimated from present data as frequencies for the considered data
#'
maf.xn = function(x){
  x = x[x != 9]
  f = min(sum(x),length(x) - sum(x))/length(x) ; f}

#' Function to compute SFS based on maf.xn results
#'
#' @param x # the matrix with data to use
#'
#' @return the SFS with MAF as allele counts for the considered data
#'
sfs.xn = function(x,N=nrow(x)){ y = apply(x,2,maf.xn);
sapply(y,function(y){rbinom(1,as.integer(0.75*N),y)})}


#' Function to compute DAF wih regards to missing data
#'
#' @param x # the matrix with data to use
#'
#' @return the DAF (observed frequency) for the considered data
#'
daf4 = function(x){
  x = x[x != 9]
  (length(x)-length(x[x==1]))/length(x)
  }

#' Function to compute DAF wih regards to missing data.
#'
#' @param x # the matrix with data to use
#'
#' @return the DAF as allele counts for the considered data
#'
daf3 = function(x){
  x = x[x != 9]
  f = length(x[x==1]) - sum(x)
  }
  
#' Function to compute SFS based on daf4 results
#'
#' @param x # the matrix with data to use
#'
#' @return the SFS with DAF as allele counts for the considered data
#'
sfs.daf4 = function(x,N=nrow(x)){ y = apply(x,2,daf4);
                              sapply(y,function(y){rbinom(1,as.integer(0.75*N),y)})}
