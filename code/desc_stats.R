##################################################
## Project: Strategic Voting in RCV
## Script purpose: CSES Descriptive Statistics
## Date: 25/11/2018
## Author:
##################################################

### 
### Dependencies
###

# Set WD etc.
library(here)

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

# Load my own functions:
functions <- here("utils/functions.r")
source(functions)

# Load necessary libraries:
library(ggplot2)
library(tidyr)
library(dplyr)
library(purrr)
library(ggtern)
library(plotly)

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
v_vec_list <- list()
for(i in 1:length(sin_vote_list)){
  weights <- big_list_na_omit[[i]]$weights
  v_vec <- ballot.props.from.vote.mat.and.weights(sin_vote_list[[i]], weights)
  big_list_na_omit[[i]]$v_vec <- as.numeric(v_vec)
}

# Drop NA cases
names(big_list)[[39]]
names(big_list)[[145]]
big_list_na_omit[[39]] <- NULL
big_list_na_omit[[144]] <- NULL

# Create ternary plot with distribution of first preferences.
v_vec_df <- matrix(NA, nrow = 160, ncol = 6)
simplex_df <- matrix(NA, nrow = 160, ncol = 3)
for(i in 1:length(big_list_na_omit)){
  v_vec <- big_list_na_omit[[i]]$v_vec
  simplex_vec <- c(v_vec[1] + v_vec[2], v_vec[3] + v_vec[4], v_vec[5] + v_vec[6])
  v_vec_df[i, ] <- as.numeric(v_vec)
  simplex_df[i, ] <- simplex_vec
}
simplex_df <- as.data.frame(simplex_df)
names(simplex_df) <- c("A", "B", "C")

ggtern(simplex_df, aes(A, B, C)) +
  geom_point()

# Create second-preference plot.

neutral_loc <- apply(second_pref_df, 1, function(x) which.min(abs(0.5 - x))))
strongest_loc <- apply(second_pref_df, 1, function(x) which.max(abs(0.5 - x))))
order <- data.frame(unlist(neutral_loc), unlist(strongest_loc))

order$medium_loc <- apply(order, 1, function(x) c(1, 2, 3)[!(c(1, 2, 3) %in% x)])
plot(strongest, neutral)

second_ordered_df <- matrix(NA, nrow = 160, ncol = 3)
for (i in 1:160){
  col_order <- as.numeric(order[i, ])
  if(col_order[1] != 2){second_pref_df[i, 2] <- 1 - second_pref_df[i, 2]}
  second_ordered_df[i, ] <- as.numeric(second_pref_df[i, col_order])
}

second_ordered_df <- as.data.frame(second_ordered_df)
plot_ly(second_ordered_df, x = ~V2, y = ~V3, z = ~V1)
plot(second_ordered_df$V2, second_ordered_df$V3)
plot(second_ordered_df$V1, second_ordered_df$V3)
plot3d(second_ordered_df$V2, second_ordered_df$V3, second_ordered_df$V1)

plane <- lm("V1 ~ V2 + V3", second_ordered_df)
planes3d(plane$coefficients["V2"], plane$coefficients["V3"], -1, plane$coefficients["(Intercept)"], alpha = 0.1, front = "line")

v_vec_df_long <- gather(v_vec_df)

ggplot(v_vec_df_long) +
  geom_density(aes(value)) +
  facet_wrap(~ key)
