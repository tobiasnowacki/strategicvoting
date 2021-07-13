# LOAD DEPENDENCIES
library(tidyverse)
library(pivotprobs)
library(gtools)
library(doParallel)

# Load functions
source("code/utils/sv_iter_four.R")
# Load data
load("output/big_list_4_parties.RData")
vap <- read.csv("data/case_vap.csv", sep = "") # voting age pop.
cat("Data imported. \n")
source("code/prep_cses.R") # data prep

# NZL_2014 only has three parties!

# which(names(big_list_na_omit) == "SWE_2014")

# SET PARAMETERS
s_val <- 85
lambda_val <- 0.05
max_iter_val <- 75 # no of iterations
which_cases <- 1:100 # which cases?

# Override previous .txt file
close(file("clusterlog_four.txt", open = "w"))

# MULTIPLE CORES
clno <- detectCores()
cl <- makeCluster(clno, outfile = "clusterlog_four.txt")
registerDoParallel(cl)
result_list <- list()
# parallelise
cases_converge <- foreach(
    j = which_cases,
    .packages = c(
        "gtools", "stringr", "tidyverse", "pivotprobs",
        "questionr"
    )
) %dopar% {
    case <- big_list_na_omit[[j]]
    inner_list <- list()
    cat(paste0("Starting case ", j, " out of ", max(which_cases), ".", "\n"))
    inner_list$rcv <- sv_iter(
        U = case$U,
        s = s_val,
        rule = "AV",
        weights = case$weights,
        lambda = lambda_val,
        max.iterations = max_iter_val,
        noisy = FALSE
    )
    inner_list$plur <- sv_iter(
        U = case$U,
        s = 85,
        rule = "plurality",
        weights = case$weights,
        lambda = lambda_val,
        max.iterations = max_iter_val,
        noisy = FALSE
    )

    filepath <- paste0(
        "output/files/",
        "1",
        "/",
        s_val,
        "_",
        names(big_list_na_omit)[j], "_iterout_fourparties.Rdata"
    )
    save(inner_list, file = filepath)
    cat(paste0("Done with case ", j, " out of ", max(which_cases), ".", "\n"))
    inner_list
}
stopCluster(cl)
# Save results
cat("Script complete.")
