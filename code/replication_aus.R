##################################################
## Project: Strategic Voting in RCV
## Script purpose: Replication of Andy's plots
##					 in the Australia case
## Date: 30/10/2018
## Author:
##################################################

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


# Import AES and Ballot data
# Note that these are normalised utilities. To obtain like-dislike scores I will need to re-run the original script.
aes_utils <- read.csv(here("../data", "australia", "AES_utility.csv"))[, -1]
aes_utils <- aes_utils[, c("GRN", "LIB", "LAB")] #common ordering

nsw <- read.csv(here("../data/australia/nsw_ballots.csv"))[, -1]
resampling <- read.csv(here("../data/australia/nsw_resampling.csv"))[, -1]

# Create DF with ballotprofiles
const_bp <- data.frame(district = nsw$District,
                       AB = nsw$`GRN.LIB` + nsw$`GRN.LIB.LAB`,
                       AC = nsw$`GRN.LAB` + nsw$`GRN.LAB.LIB`,
                       BA = nsw$`LIB.GRN` + nsw$`LIB.GRN.LAB`,
                       BC = nsw$`LIB.LAB` + nsw$`LIB.LAB.GRN`,
                       CA = nsw$`LAB.GRN` + nsw$`LAB.GRN.LIB`,
                       CB = nsw$`LAB.LIB` + nsw$`LAB.LIB.GRN`,
                       A = nsw$GRN,
                       B = nsw$LIB,
                       C = nsw$LAB)

# For now, let's not use truncated ballots:
const_bp_no_trunc <- const_bp
const_bp_no_trunc[, 8:10] <- 0
const_bp_no_trunc[, 2:10] <- t(apply(const_bp_no_trunc[, 2:10], 1, function(x) x / sum(x)))

# TO-DO: for data *WITH* truncated prefs, make sure sin_vec is evaluated correctly
# 	--> functions.r
# TO-DO (GENERAL): check whether opt. vote results make sense (e.g. very little severe pushover / third pref opt)
# DONE -- SEEMS TO BE WORKING


# TEST IF OPT VOTE YIELDS SAME RESULTS AS ANDY'S FUNCTION
v_vec <- c(as.numeric(const_bp_no_trunc[4, 2:10]))
test_out <- sv(U = aes_utils[, 1:3], v.vec = v_vec[1:6], s = 80, rule = "AV")
sum(table(test_out$opt.votes.strategic, test_out$opt.votes.sincere)
test_out_toby <- return_sv_prop(v_vec, aes_utils[, 1:3], list(80))
# YES!


# LEVELS OF STRATEGIC VOTING.

# Set levels of s at which to evaluate.
s_list <- as.list(seq(from = 10, to = 130, by = 10))

# Run loop over all constituencies.
const_opt_dist <- list()
for(i in 1:nrow(const_bp_no_trunc)){
	print(i)
	prop <- return_sv_prop(const_bp_no_trunc[i, 2:10], aes_utils[, 1:3], s_list)
	prop[, 1:3] <- as.data.frame(t(apply(prop[, 1:3], 1, function(x) x / sum(x))))
	prop$const <- const_bp_no_trunc[i, 1]
	const_opt_dist[[i]] <- prop
}

const_opt_dist_df_wide <- do.call(rbind, const_opt_dist)
const_opt_dist_df <- melt(const_opt_dist_df_wide[, c("second", "third", "plurality_second", "const", "s")], 
							id.vars = c("const", "s"))

# Get means.
const_opt_dist_df_agg <- as.data.frame(const_opt_dist_df %>% 
							group_by(variable, s) %>% 
							summarize(mean(value)))
names(const_opt_dist_df_agg)[3] <- "value"

# Plot. (SI over levels of s, by SI type)
aus_freq <- ggplot(const_opt_dist_df, aes(x = s, y = value)) +
	geom_line(aes(colour = variable, group = interaction(const, variable)), alpha = 0.3) +
	geom_line(data = const_opt_dist_df_agg, aes(colour = variable, group = variable, x = s, y = value),
		size = 3) + 
	labs(x = "Information (s)", 
		y = "Proportion of voters in AES casting ballot type",
		colour = "Sincere pref. as first on ballot") +
	theme_bw() +
	theme(legend.position = "bottom")
gg_path <- here("output/figures/australia_sv_freq.pdf")
ggsave(gg_path)

# Next, let's have a distribution of cases according to whether they are more amenable to SV in RCV or Plurality.
# head(const_opt_dist_df_wide)

## Table: by constituency
const_opt_dist_df_wide$inc_rcv <- const_opt_dist_df_wide[, 2] + const_opt_dist_df_wide[, 3]
const_opt_dist_df_wide$inc_plur <- const_opt_dist_df_wide$plurality_second

# Plot. (Total SI, for a given level of s)
aus_prop <- ggplot(const_opt_dist_df_wide, aes(x = inc_plur, y = inc_rcv)) +
	geom_point() +
	geom_abline(slope = 1, intercept = 0) +
	facet_wrap(~ s) +
	theme_bw() +
	scale_x_continuous(limits = c(0, 0.7), expand = c(0, 0)) +
	scale_y_continuous(limits = c(0, 0.7), expand = c(0, 0)) +
	labs(x = "Proportion of AES respondents with positive SI under Plurality",
		y = "Proportion of AES respondents with positive SI under RCV")
gg_path2 <- here("output/figures/australia_sv_prop.pdf")
ggsave(gg_path2)


## QQ-PLOT

# Create dataframe of qq-plot coordinates, by constituency and by s
qq_list <- lapply(s_list, function(x) qq_function(const_bp_no_trunc, aes_utils[, 1:3], x))
qq_df <- do.call(rbind, qq_list)
qq_df$s <- unlist(rep(as.vector(s_list), each = nrow(aes_utils) * nrow(const_bp_no_trunc)))

# Also get coordinates for qq-plot aggregated over constituencies, by s
big_qq_list <- lapply(qq_list, function(z) as.data.frame(qqplot(x = z$x, y = z$y, plot.it = FALSE)))
big_qq_df <- do.call(rbind, big_qq_list)
big_qq_df$s <- rep(unlist(s_list), each = nrow(aes_utils) * nrow(const_bp_no_trunc))

# Plot.
qq_plot_faceted <- ggplot(qq_df) +
	geom_line(aes(x = x, y = y, group = const), alpha = 0.1) +
	geom_line(data = big_qq_df, aes(x = x, y = y), colour = "red", lwd = 2) + 
	geom_abline(intercept = 0, slope = 1, linetype = "dotted", colour = "blue") +
	theme_bw()  + 
	facet_wrap(vars(s))
ggsave(here("../output/figures/australia_sv_qq.pdf"), qq_plot_faceted)

## OCCURRENCE OF VOTING PARADOXES

# todo: write functions (wasted vote; non-monotonicity)


## STRATEGIC INTERDEPENDENCE
## TO-DO: Fix plurality; compare quantities to "level-1 strat voting"
## Second part of interdependence replication?


# Todo: Pack all of this into a function at the end and shift to functions.r
lambda <- 0.3

# Get the optimal votes for one constituency.
v_vec <- const_bp_no_trunc[5, 2:10]
mega_df <- return_sv_tau(v_vec, aes_utils[, 1:3], s_list)
by_s_df <- split(mega_df, mega_df$s)

# For each s, get 6x6 (and 3x3) mat
vote_matrix <- function(df, type = "rcv"){
	if(type == "rcv"){
		df$opt_rcv <- factor(df$opt_rcv, levels = 1:6)
		tab <- tapply(df$opt_rcv, df$sin_rcv, table)
		tab <- do.call(rbind, tab)
		return(tab)
	}
	if(type == "plur"){
		df$opt_plur <- factor(df$opt_plur, levels = 1:3)
		tab <- tapply(df$opt_plur, df$sin_rcv, table)
		tab <- do.call(rbind, tab)
		return(tab)
	}
}

vote_mat_rcv <- lapply(by_s_df, function(x) vote_matrix(x, type = "rcv"))
vote_mat_plur <- lapply(by_s_df, function(x) vote_matrix(x, type = "plur"))

v_vec_init_weighted <- as.numeric(v_vec[1:6] / sum(v_vec[1:6]))
v_vec_init_weighted_plur <- c(v_vec_init_weighted[1] + v_vec_init_weighted[2],
	v_vec_init_weighted[3] + v_vec_init_weighted[4], 
	v_vec_init_weighted[5] + v_vec_init_weighted[6])

# I can wrap this into a new function.
new_vec <- lapply(vote_mat_rcv, function(x) v_vec_init_weighted %*% x)
new_vec <- lapply(new_vec, function(x) x / sum(x))
new_vec <- lapply(new_vec, function(x) lambda * x + (1 - lambda) * v_vec_init_weighted)

# How will this work for plurality? Will need 6 x 3 matrix, rather than 6 x 6.
new_vec_plur <- lapply(vote_mat_plur, function(x) v_vec_init_weighted %*% x)
new_vec_plur <- lapply(new_vec_plur, function(x) x / sum(x))
new_vec_plur <- lapply(new_vec_plur, function(x) lambda * x + (1 - lambda) * v_vec_init_weighted_plur)

# Get proportion of level-2 strategic voters

# Doesn't make sense right now -- produces output for multiple values for s in each list element even though each new_vec comes from a specific value for s
# work with for loop instead?
# level_2_strat_incent_rcv <- lapply(new_vec, function(x) return_sv_prop(c(x, 0, 0, 0), aes_utils[, 1:3], s_list))
# Given for-loop, obsolete?

# To obtain the proportion of voters in plurality, I will need to do either of the following:
# (a) Take 3-item vec in new_vec_plur and split it into 6 such that I can run return_sv_prop
new_vec_plur_six <- lapply(new_vec_plur, function(x) rep(x, each = 2) / 2)
# (b) Write a new function for plurality specifically (why if can avoid?).

# For loop: For each level of S, compute proportion of voters voting for first, second, third vote in AV and first, second in plurality.

inter_df <- matrix(NA, ncol = 5, nrow = length(s_list))

for (i in 1:length(s_list)){
	out_rcv <- return_sv_prop(c(new_vec[[i]], 0, 0, 0), aes_utils[, 1:3], list(s_list[[i]]))
	# to-do (still): add plurality proportions, see above comment.
	out_plur <- return_sv_prop(c(new_vec_plur_six[[i]], 0, 0, 0), aes_utils[, 1:3], list(s_list[[i]]))
	inter_df[i, 1:3] <- as.matrix(out_rcv[1, 1:3])
	inter_df[i, 4:5] <- as.matrix(out_rcv[1, 4:5])
}

# compare to original
return_sv_prop(c(v_vec_init_weighted, 0, 0, 0), aes_utils[, 1:3], s_list)

# Loop over constituencies and put into mega-dataframe for plotting.
