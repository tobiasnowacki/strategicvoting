##################################################
## Project: Strategic Voting in RCV
## Script purpose: Classification
## Date: 17/12/2018
## Author:
##################################################

### DEPENDENCIES
library(here)
library("scatterplot3d")
library(car)
install.packages("rgl")
library(rgl)
library(plotly)

# Load my own functions:
source(here("utils/functions.r"))
# Load pivotal probability functions:
av_piv_path <- here("utils/av_pivotal_probs_analytical_general_v2.r")
source(av_piv_path)
plur_piv_path <- here("utils/plurality_pivotal_probabilities_analytical.r")
source(plur_piv_path)

# To replicate Andy's function(s):
sim_appr2 <- here("utils", "general_iteration_simulation_approach.r")
source(sim_appr2)
sv_file <-  here("utils/sv.r")
source(sv_file)

# Load CSES data:
load(here("../output/cses_big_list_2.RData"))

remove_nas <- function(x){
  mat <- cbind(x$U, x$weights)
  mat <- na.omit(mat)
  return(list(U = mat[, 1:3], weights = as.numeric(mat[, 4])))
}

create_v_vec <- function(x){
  x$U <- x$U + runif(nrow(x$U) * ncol(x$U), min = 0, max = 0.001)
  sin_vote <- as.numeric(apply(x$U, 1, function(x) sin_vote_scalar(x)))
  num_list <- c(1:6)
  sin_df <- (sapply(num_list, function(x) as.numeric(sin_vote == x)))
  sin_df <- sin_df * x$weights
  sin_vec <- colSums(sin_df)
  return(sin_vec / sum(sin_vec))
}

# Create list with v_vecs from CSES utility dfs:
big_list_na_omit <- lapply(big_list, function(x) remove_nas(x))
sin_vote_list <- lapply(big_list_na_omit, function(x) sincere.vote.mat.from.U(x$U, rule = "AV"))
for(i in 1:length(sin_vote_list)){
  weights <- big_list_na_omit[[i]]$weights
  v_vec <- ballot.props.from.vote.mat.and.weights(sin_vote_list[[i]], weights)
  big_list_na_omit[[i]]$v_vec <- as.numeric(v_vec)
}

### PREPARE DATA

# Create dataframe of v_vec cases
v_vec_df <- as.data.frame(t(sapply(big_list_na_omit, function(x) x$v_vec)))
names(v_vec_df) <- c("ABC", "ACB", "BAC", "BCA", "CAB", "CBA")

second_pref_df <- data.frame(mAB = v_vec_df[, 1] / (v_vec_df[, 1] + v_vec_df[, 2]), mBA = v_vec_df[, 3] / (v_vec_df[, 3] + v_vec_df[, 4]), mCB = v_vec_df[, 6] / (v_vec_df[, 5] + v_vec_df[, 6]))

scatterplot3d(second_pref_df$mAB, second_pref_df$mCB, second_pref_df$mBA)
scatter3d(second_pref_df$mAB, second_pref_df$mCB, second_pref_df$mBA)
plot_ly(second_pref_df, x = ~mAB, y = ~mBA, z = ~mCB)
