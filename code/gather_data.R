# LOAD DEPENDENCIES
library(tidyverse)
library(pivotprobs)
library(gtools)
library(doParallel)
source("code/utils/new_sv_iter.R")
library(devtools)
source_url("https://raw.githubusercontent.com/tobiasnowacki/RTemplates/master/plottheme.R")

# Load data
load("output/big_list_2.RData")
vap <- read.csv("data/case_vap.csv", sep = "") # voting age pop.
cat("Data imported. \n")
source("code/prep_cses.R")  # data prep

nn = names(big_list_na_omit)

# This should become a variable so I can run the script serverside.
s = 85
lambda = 1

fpath = function(lambda, s, ext){
  paste0("output/files/", lambda, "/", s, "_", ext, ".Rdata")
}


# Function to load all cases for param values
return_obj = function(case, lambda, s){
  load(paste0("output/files/", lambda, "/", s, "_", case, "_iterout.Rdata"))
  return(inner_list)
}

bind_together = map(nn, ~ return_obj(.x, lambda, s))
names(bind_together) = nn


# Function to get summary statistics
get_sum_stats = function(obj){
  w = obj$rcv[[i]]$weights
  map_dfr(c("rcv", "plur"), function(y){
    names(obj[[y]]) = 1:length(obj[[y]])
    map(obj[[y]], ~
          tibble(
            prev = weighted.mean(.x$tau > 0, w),
            mag  = weighted.mean(.x$tau[.x$tau > 0], w[.x$tau > 0]),
            eb = prev * mag)) %>%
      bind_rows(.id = "iter") %>%
      mutate(system = y)
  })
}

sum_df = map_dfr(bind_together, ~ get_sum_stats(.x), .id = "case")
save(sum_df, file = fpath(lambda, s, "sum"))

# -- needs further processing

# Function to save v_vec data
get_vvecs = function(obj){
  w = obj$rcv[[i]]$weights
  outobj = list()
  outobj$rcvvec = map_dfr(obj$rcv, ~ .x$v.vec.before)
  outobj$plurvec = map_dfr(obj$plur, ~ .x$v.vec.before)
  outobj$rcvbr = map_dfr(obj$rcv, ~ .x$best.response.v.vec)
  outobj$plurbr = map_dfr(obj$plur, ~ .x$best.response.v.vec)
  return(outobj)
}

vvecdf = map(bind_together, ~ get_vvecs(.x))
names(vvecdf) = nn[1:10]
save(vvecdf, file = fpath(lambda, s, "vvec"))

# Create distance plots
joint_v_vec_plot(vvecdf, "output/")
# Create vvec plots
plot_v_vec_distance(vvecdf, "output/", n_lag = 5, avg_span = 10)

# vvec plots...
