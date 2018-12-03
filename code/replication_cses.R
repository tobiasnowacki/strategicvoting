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

s_list <- as.list(c(15, 30, 45, 60, 75, 85, 90))

names(big_list)[[39]]
names(big_list)[[145]]
big_list_na_omit[[39]] <- NULL
big_list_na_omit[[144]] <- NULL

# Loop that creates the tau objects
# uses my own function -- see below for faster implementation
# sv_list <- list()
# for(i in 1:length(big_list_na_omit)){
#   print(i)
#   sv_list[[i]] <- return_sv_tau(c(big_list_na_omit[[i]]$v_vec, 0, 0, 0), big_list_na_omit[[i]]$U, s_list)
# }

# Use Andy's function to get us started and convert into my dataframe
sv_list <- list()
for(i in 1:length(big_list_na_omit)){
  print(i)
  this_list <- big_list_na_omit[[i]]
  df_list <- lapply(s_list, function(x) convert_andy_to_sv_item(this_list, s = x))
  df <- as.data.frame(do.call(rbind, df_list))
  df$case <- names(big_list_na_omit)[[i]]
  #df <- apply(df, 2, as.numeric)
  sv_list[[i]] <- df
}

####

# Check whether they really do produce the same as Andy's function!
# test_toby <- return_sv_tau(c(big_list_na_omit[[109]]$v_vec, 0, 0, 0), big_list_na_omit[[109]]$U, list(60))
# test_andy <- sv(U = big_list_na_omit[[109]]$U, weights = big_list_na_omit[[109]]$weights, s = 60, rule = "AV")
# # OK, what I need to do is to run the entire loop with Andy's function and check the entire DF for discrepancies.
# 
# andy_list <- list()
# andy_list_prop <- list()
# for(i in 1:length(big_list)){
#   this_list <- big_list[[i]]
#   print(i)
#   sv_obj <- sv(U = this_list$U, weights = this_list$weights, s = 60, rule = "AV")
#   sv_obj_plur <- sv(U = this_list$U, weights = this_list$weights, s = 60)
#   andy_list[[i]] <- sv_obj
#   prop_av <- sum(sv_obj$weights[!is.na(sv_obj$tau) & sv_obj$tau > 0]) / sum(sv_obj$weights[!is.na(sv_obj$tau)])
#   prop_plur <- sum(sv_obj_plur$weights[!is.na(sv_obj_plur$tau) & sv_obj_plur$tau > 0]) / sum(sv_obj_plur$weights[!is.na(sv_obj$tau)]) 
#   andy_list_prop[[i]] <- c(prop_av, prop_plur)
# }
# 
# andy_list_prop_df <- do.call(rbind, andy_list_prop)
# 
# # compare:
# andy_list_prop_df == prop_df[prop_df$s == 60, c(8, 9)]

####

save(file = here("../output/sv_list.Rdata"), sv_list)
# load(here("../output/sv_list.Rdata"))


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
names(prop_df)[1:5] <- c("rcv_first", "rcv_second", "rcv_third", "plur_first", "plur_second")
prop_df_long <- melt(prop_df[, c(2, 3, 5, 6, 7)], id.vars = c("case", "s"))
prop_df_agg <- as.data.frame(prop_df_long %>% 
                               group_by(variable, s) %>% 
                               summarize(mean(value)))
names(prop_df_agg)[3] <- "value"


cses_prop <- ggplot(prop_df_long, aes(x = s, y = value)) +
  geom_line(aes(colour = variable, group = interaction(case, variable)), alpha = 0.05) +
  geom_line(data = prop_df_agg, aes(colour = variable, group = variable, x = s, y = value),
            size = 2) + 
  labs(x = "Information (s)", 
       y = "Proportion of voters in CSES (case) casting ballot type",
       colour = "Ballot order") +
  theme_bw() +
  # scale_color_manual(values = c("lime green", "blue", "red")) + 
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0.1, 0), limits = c(0, 0.5)) + 
  theme(legend.position = "bottom", legend.direction = "vertical")
ggsave(here("../output/figures/cses_freq.pdf"), cses_prop, height = 5, width = 4)

# Distribution of incentives (in scatterplot)
prop_df$inc_rcv <- prop_df$rcv_second + prop_df$rcv_third
prop_df$inc_plur <- prop_df$plur_second

cses_inc <- ggplot(prop_df, aes(x = inc_plur, y = inc_rcv)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) + 
  facet_wrap(~ s) +
  theme_bw() +
  scale_x_continuous(limits = c(0, 0.7), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.7), expand = c(0, 0)) +
  labs(x = "Proportion of CSES respondents with positive SI under Plurality",
       y = "Proportion of CSES respondents with positive SI under RCV")
ggsave(here("../output/figures/cses_prop.pdf"), cses_inc, width = 5, height = 5)


# Check alternative method of computing incentive proportions
# sv_prop_alt <- function(tau_obj, weights, s = 60){
#   tau_obj <- tau_obj[tau_obj$s == s, ]
#   inc_rcv <- sum(weights[tau_obj$tau_rcv > 0]) / sum(weights)
#   inc_plur <- sum(weights[tau_obj$tau_plur > 0]) / sum(weights)
#   return(c(inc_rcv, inc_plur))
# }
# 
# prop_list_alt <- list()
# for(i in 1:length(sv_list)){
#   prop_list_alt[[i]] <- sv_prop_alt(sv_list[[i]], big_list_na_omit[[i]]$weights, 60)
# }
# prop_df_alt <- as.data.frame(do.call(rbind, prop_list_alt))
# ggplot(prop_df_alt, aes(V2, V1)) +
#   geom_point() +
#   geom_abline(slope = 1, intercept = 0) + 
#   theme_bw() +
#   scale_x_continuous(limits = c(0, 0.7), expand = c(0, 0)) +
#   scale_y_continuous(limits = c(0, 0.7), expand = c(0, 0))

## QQ-Plots

qq_mega_list <- list()
for(i in 1:length(sv_list)){
  print(i)
  x <- sv_list[[i]]
  utils <- big_list_na_omit[[i]]$U
  case <- names(big_list_na_omit)[[i]]
  qq <- qq_function_two(x, utils)
  qq$case <- case
  qq_mega_list[[i]] <- qq
}
qq_mega_df <- as.data.frame(do.call(rbind, qq_mega_list))
qq_mega_by_s <- split(qq_mega_df, qq_mega_df$s)
qq_agg <- lapply(as.list(names(qq_mega_by_s)), function(s) {
    z <- qq_mega_by_s[[s]]
    df <- as.data.frame(qqplot(x = z$x, y = z$y, plot.it = FALSE))
    df$s <- rep(s, nrow(df))
    return(df)
  })
qq_agg_df <- as.data.frame(do.call(rbind, qq_agg))
cses_qq <- ggplot(qq_mega_df, aes(x = x, y = y)) +
  geom_line(aes(x = x, y = y, group = case), alpha = 0.1) +
  geom_line(data = qq_agg_df, aes(x = x, y = y), colour = "red", lwd = 2) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", colour = "blue") +
  theme_bw()  + 
  facet_wrap(vars(s)) +
  xlim(-30, 30) + ylim(-30, 30) +
ggsave(here("../output/figures/cses_qq.pdf"), cses_qq)


# Voting paradoxes

# Interdependence



