# LOAD DEPENDENCIES
library(tidyverse)
library(pivotprobs)
library(gtools)
library(doParallel)
source("code/utils/new_sv_iter.R")

# Load data
load("output/big_list_2.RData")
vap <- read.csv("data/case_vap.csv", sep = "") # voting age pop.
cat("Data imported. \n")
source("code/prep_cses.R")  # data prep

nn = names(big_list_na_omit)

return_obj = function(case, lambda, s){
  load(paste0("output/files/", lambda, "/", s, "_", case, "_iterout.Rdata"))
  return(inner_list)
}

bind_together = map(nn[1:10], ~ return_obj(.x, 1, 85))