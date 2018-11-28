##################################################
## Project: Strategic Voting in RCV
## Script purpose: Replication of Andy's plots
##					 in the CSES case
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
library(reshape2)
library(dplyr)
library(purrr)

# Load CSES data:
load(here("../output/cses_big_list_2.RData"))

###
### ANALYSIS
###

# Create list with v_vecs from CSES utility dfs:
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

big_list_na_omit <- lapply(big_list, function(x) remove_nas(x))
sin_vote_list <- lapply(big_list_na_omit, function(x) sincere.vote.mat.from.U(x$U, rule = "AV"))

v_vec_list <- list()
for(i in 1:length(sin_vote_list)){
  weights <- big_list_na_omit[[i]]$weights
  v_vec <- ballot.props.from.vote.mat.and.weights(sin_vote_list[[i]], weights)
  big_list_na_omit[[i]]$v_vec <- as.numeric(v_vec)
}

s_list <- as.list(seq(from = 10, to = 120, by = 10))

names(big_list)[[39]]
names(big_list)[[145]]
big_list_na_omit[[39]] <- NULL
big_list_na_omit[[144]] <- NULL

sv_list <- list()
for(i in 144:length(big_list_na_omit)){
  print(i)
  sv_list[[i]] <- return_sv_tau(c(big_list_na_omit[[i]]$v_vec, 0, 0, 0), big_list_na_omit[[i]]$U, s_list)
}

save(file = here("../output/sv_list.Rdata"), sv_list)


# Proportion of optimal strategic votes
prop_list <- list()
for(i in 1:length(sv_list)){
  print(i)
  prop_list[[i]] <- sv_prop(sv_list[[i]], big_list_na_omit[[i]]$weights)
  df <- as.data.frame(prop_list[[i]])
  df$case <- as.character(names(big_list_na_omit)[[i]])
  df$s <- as.numeric(s_list)
  n <- apply(df[, 1:3], 1, function(x) sum(x))
  df[, 1:5] <- df[, 1:5] / sum(big_list_na_omit[[i]]$weights)
  prop_list[[i]] <- df
}
prop_df <- do.call(rbind, prop_list)
names(prop_df)[1:5] <- c("rcv_first", "rcv_sec", "rcv_third", "plur_first", "plur_sec")
prop_df_long <- melt(prop_df[, c(2, 3, 5, 6, 7)], id.vars = c("case", "s"))
prop_df_agg <- as.data.frame(prop_df_long %>% 
                               group_by(variable, s) %>% 
                               summarize(mean(value)))
names(prop_df_agg)[3] <- "value"

# Need to include weights to fully replicate Andy's function.

cses_prop <- ggplot(prop_df_long, aes(x = s, y = value)) +
  geom_line(aes(colour = variable, group = interaction(case, variable)), alpha = 0.06) +
  geom_line(data = prop_df_agg, aes(colour = variable, group = variable, x = s, y = value),
            size = 1) + 
  labs(x = "Information (s)", 
       y = "Proportion of voters in CSES (case) casting ballot type",
       colour = "Ballot order") +
  theme_bw() +
  scale_color_manual(values = c("lime green", "blue", "red")) + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0.1, 0)) + 
  theme(legend.position = "bottom")
ggsave(here("../output/figures/cses_prop.pdf"), cses_prop)

# Distribution of incentives 
prop_df$inc_rcv <- prop_df$rcv_sec + prop_df$rcv_third
prop_df$inc_plur <- prop_df$plur_sec

cses_inc <- ggplot(prop_df[prop_df$s == 60, ], aes(x = inc_plur, y = inc_rcv)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) + 
  facet_wrap(~ s) +
  theme_bw() +
  scale_x_continuous(limits = c(0, 0.7), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.7), expand = c(0, 0)) +
  labs(x = "Proportion of CSES respondents with positive SI under Plurality",
       y = "Proportion of CSES respondents with positive SI under RCV")
ggsave(here("../output/figures/cses_inc.pdf"), cses_inc)


# Check alternative method of computing incentive proportions
sv_prop_alt <- function(tau_obj, weights, s = 60){
  tau_obj <- tau_obj[tau_obj$s == s, ]
  inc_rcv <- sum(weights[tau_obj$tau_rcv > 0]) / sum(weights)
  inc_plur <- sum(weights[tau_obj$tau_plur > 0]) / sum(weights)
  return(c(inc_rcv, inc_plur))
}

prop_list_alt <- list()
for(i in 1:length(sv_list)){
  prop_list_alt[[i]] <- sv_prop_alt(sv_list[[i]], big_list_na_omit[[i]]$weights, 60)
}
prop_df_alt <- as.data.frame(do.call(rbind, prop_list_alt))
ggplot(prop_df_alt, aes(V2, V1)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) + 
  theme_bw() +
  scale_x_continuous(limits = c(0, 0.7), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.7), expand = c(0, 0))

# Ok, still not the same as Andy's results. Wonder why? Next step is to compare the precise tau output.


# Voting paradoxes

# Interdependence



