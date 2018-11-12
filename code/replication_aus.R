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


# Return data.frame with levels of strat voting by s and constituency.

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
head(const_opt_dist_df_wide)

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


## QQ_PLOT

qq_function <- function(const_bp_no_trunc, utils, s){
	const_taus <- list()
	const_taus_qq <- list()
	for(i in 1:nrow(const_bp_no_trunc)){
			# print(const_bp_no_trunc[i, ])
			v_vec <- const_bp_no_trunc[i, 2:10]
			const_taus[[i]] <- return_sv_tau(as.numeric(v_vec), utils, s)
			const_taus_qq[[i]] <- as.data.frame(qqplot(x = unlist(const_taus[[i]]$tau_vec_plur[const_taus[[i]]$V4 == s]), 
				y = unlist(const_taus[[i]]$tau_vec_rcv[const_taus[[i]]$V4 == s]), plot.it = FALSE))
	}
	const_taus_qq_df <- do.call(rbind, const_taus_qq)
	const_taus_qq_df$const <- rep(const_bp_no_trunc$district, each = nrow(aes_utils))
	return(const_taus_qq_df)
}

qq_list <- lapply(s_list, function(x) qq_function(const_bp_no_trunc, aes_utils[, 1:3], x))

qq_df <- do.call(rbind, qq_list)
qq_df$s <- unlist(rep(as.vector(s_list), each = nrow(aes_utils) * nrow(const_bp_no_trunc)))

big_qq_list <- lapply(qq_list, function(z) as.data.frame(qqplot(x = z$x, y = z$y, plot.it = FALSE)))
big_qq_df <- do.call(rbind, big_qq_list)
big_qq_df$s <- rep(unlist(s_list), each = nrow(aes_utils) * nrow(const_bp_no_trunc))

qq_plot_faceted <- ggplot(qq_df) +
	geom_line(aes(x = x, y = y, group = const), alpha = 0.1) +
	geom_line(data = big_qq_df, aes(x = x, y = y), colour = "red", lwd = 2) + 
	geom_abline(intercept = 0, slope = 1, linetype = "dotted", colour = "blue") +
	theme_bw()  + 
	facet_wrap(vars(s))
ggsave(here("output/figures/australia_sv_qq.pdf"))
