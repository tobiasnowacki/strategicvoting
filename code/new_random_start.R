cmd_line_args <- commandArgs(trailingOnly = TRUE)

# LOAD DEPENDENCIES
library(tidyverse)
library(pivotprobs)
library(gtools)
library(doParallel)
source("code/utils/new_sv_iter.R")
source("code/utils/sv_theme_template.R")

# Load data
load("output/big_list_2.RData")
vap <- read.csv("data/case_vap.csv", sep = "") # voting age pop.
cat("Data imported. \n")
source("code/prep_cses.R")  # data prep
set.seed(cmd_line_args[1])

names_vec <- names(big_list_na_omit)
v_vec_list <- list()

uniform_ternary <- rdirichlet(10, rep(1, 6))

# Set up lists.
big_rcv_sum_sense <- list()
big_rcv_vec_sense <- list()

# Set up parameters.
lambda <- 0.05
s_val <- 85

# Override previous .txt file
close(file("clusterlog_rand.txt", open="w"))

clno <- detectCores()
cl <- makeCluster(clno, outfile = "clusterlog_rand.txt")
registerDoParallel(cl)

# For loop -- every iteration is a random vvec
out <- foreach(rand_iter = 1:10,
               .packages = c("gtools", "stringr", "tidyverse", "pivotprobs")
) %dopar% {
  prec <- 85
  s_val <- 85
  cat(paste0("\n === starting point = ", rand_iter, " =============== \n"))
  rcv_vec <- list()
  rand_v_vec <- uniform_ternary[rand_iter, ] %>% as.numeric
  for (case in 1:160) {
    cat(paste0(case, ": ", names(big_list_na_omit)[case], "   "))
    out <- sv_iter(U = big_list_na_omit[[case]]$U,
                   weights = big_list_na_omit[[case]]$weights,
                   rule = "AV",
                   s = s_val,
                   lambda = 0.05,
                   starting.v.vec = rand_v_vec,
                   max.iterations = 61)
    names(out) = 1:61
    dat = map_dfr(out, ~ .x$v.vec.before %>% as.matrix) %>% t %>% as.data.frame
    colnames(dat) = c("abc", "acb", "bac", "bca", "cab", "cba")
    rcv_vec[[case]] <- tibble(dat, 
                             case = names(big_list_na_omit)[[case]],
                             system = "IRV", 
                             iter = 1:61, 
                             part = rand_iter)
  }
  rcv_vec
}
stopCluster(cl)

save.image(paste0("output/files/random_", cmd_line_args[1], ".Rdata"))

cat("Done.")
