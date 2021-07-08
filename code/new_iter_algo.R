# LOAD DEPENDENCIES
library(tidyverse)
library(pivotprobs)
library(gtools)
library(doParallel)

# Load functions
source("code/utils/new_sv_iter.R")
# Load data
load("output/big_list_2.RData")
vap <- read.csv("data/case_vap.csv", sep = "") # voting age pop.
cat("Data imported. \n")
source("code/prep_cses.R")  # data prep

# which(names(big_list_na_omit) == "SWE_2014")

# SET PARAMETERS
cmd_line_args <- commandArgs(trailingOnly = TRUE)
cat(cmd_line_args, sep = "n")

s_list <- as.list(c(10, 25, 40, 55, 70, 85)) # precision (s)
lambda_list <- as.list(c(0.05, 0.1, 0.01))   # responsiveness ()
# epsilon_thresh <- 0.001            # threshold (epsilon)
max_iter_val <- 100              # no of iterations
which_cases <- 32          # which cases?

# If command line does not pick s, lambda:
ifelse(length(cmd_line_args >= 1),
       s_choice <- as.numeric(cmd_line_args[1]),
       s_choice <- 6)
ifelse(length(cmd_line_args >= 2),
       l_choice <- as.numeric(cmd_line_args[2]),
       l_choice <- 1)

lambda_val <- lambda_list[[l_choice]]
s_val <- s_list[[s_choice]]

# Override previous .txt file
close(file("clusterlog.txt", open="w"))

# MULTIPLE CORES
clno <- detectCores()
cl <- makeCluster(clno, outfile = "clusterlog.txt")
registerDoParallel(cl)
result_list = list()
# parallelise
cases_converge <- foreach(j = which_cases, 
            .packages = c("gtools", "stringr", "tidyverse", "pivotprobs",
                          "questionr")
            ) %dopar% {
              case = big_list_na_omit[[j]]
	      inner_list = list()
              cat(paste0("Starting case ", j, " out of ", max(which_cases), ".", "\n"))
              inner_list$rcv = sv_iter(
                U = case$U,
                s = s_val,
                rule = "AV",
                weights = case$weights,
                lambda = lambda_val,
                max.iterations = max_iter_val,
                noisy = FALSE
              )
              inner_list$plur = sv_iter(
                U = case$U,
                s = 85,
                rule = "plurality",
                weights = case$weights,
                lambda = lambda_val,
                max.iterations = max_iter_val,
                noisy = FALSE
              )
              
	      filepath = paste0("output/files/", l_choice, "/", s_val, "_", names(big_list_na_omit)[j], "_iterout.Rdata")
              save(inner_list,  file = filepath)
              cat(paste0("Done with case ", j, " out of ", max(which_cases), ".", "\n"))
	      inner_list
            }
stopCluster(cl)
# Save results
cat("Script complete.")
