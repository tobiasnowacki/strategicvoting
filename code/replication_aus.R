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
library(ggtern)


# Import AES and Ballot data
# Note that these are normalised utilities. To obtain like-dislike scores I will need to re-run the original script.
aes_utils <- read.csv(here("../data", "australia", "AES_utility.csv"))[, -1]
aes_utils <- aes_utils[, c("GRN", "LIB", "LAB")] #common ordering

# Import non-standardised utilities
aes_utils_raw <- read.csv(here("../data", "australia", "aes_nsw_full.csv"))[, -1]

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

### ---------------------------------
### SUMMARY STATS
### ---------------------------------

simplex_df <- data.frame(A = const_bp[, 2] + const_bp[, 3] + const_bp[, 8],
                         B = const_bp[, 4] + const_bp[, 5] + const_bp[, 9],
                         C = const_bp[, 6] + const_bp[, 7] + const_bp[, 10])
ggtern(simplex_df, aes(A, B, C)) +
  geom_point()

### ---------------------------------
### ANALYSIS
### ---------------------------------


# Set levels of s at which to evaluate.
s_list <- as.list(c(15, 30, 45, 60, 75, 85, 90))

### PLAYGROUND ### ---
### END PLAYGROUND ### ----

# Run return_sv_tau loop over all constituencies.
set.seed(23112018)
mega_tau_list <- list()

for(i in 1:nrow(const_bp)){
	print(i)
	v_vec <- as.numeric(const_bp[i, 2:10])
	mega_tau_list[[i]] <- return_sv_tau(v_vec, aes_utils_raw, s_list)
	mega_tau_list[[i]]$const <- const_bp[i, 1]
}

# Do the same as for-loop, just with 
library(pbapply)
mega_tau_list_2 <- pbapply(const_bp[, 2:10], 1, function(x) return_sv_tau(as.numeric(x), aes_utils_raw, s_list))

mega_tau_list <- mega_tau_list_2

for(i in 1:nrow(const_bp)){
  mega_tau_list[[i]]$const <- const_bp[i, 1]
}

# save as separate object to avoid having to run it every time.
save(mega_tau_list, file = here("../output/mega_tau_list_trunc_plur.Rdata"))
load(here("../output/mega_tau_list_trunc_plur.Rdata"))

##
# From resulting loop, run:
##

# (1) levels of strategic voting
prop_list <- lapply(mega_tau_list, function(x) try(sv_prop(x)))
prop_df <- as.data.frame(do.call(rbind, prop_list))
prop_df$const <- rep(c(const_bp$district), each = length(s_list))
prop_df$s <- rep(unlist(s_list), nrow(const_bp))
prop_df[, 1:5] <- prop_df[, 1:5] / 1133
prop_df <- prop_df[, c(2, 3, 5, 6, 7)]
names(prop_df)[1:3] <- c("RCV_second", "RCV_third", "plur_second")

prop_df_long <- melt(prop_df, id.vars = c("const", "s"))
prop_df_agg <- as.data.frame(prop_df_long %>% 
                                         group_by(variable, s) %>% 
                                         summarize(mean(value)))
names(prop_df_agg)[3] <- "value"

aus_freq <- ggplot(prop_df_long, aes(x = s, y = value)) +
  geom_line(aes(colour = variable, group = interaction(const, variable)), alpha = 0.05) +
  geom_line(data = prop_df_agg, aes(colour = variable, group = variable, x = s, y = value),
            size = 2) + 
  labs(x = "Information (s)", 
       y = "Proportion of voters in AES casting ballot type",
       colour = "Sincere pref. as first on ballot") +
  theme_bw() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0.1, 0), limits = c(0, 0.5)) + 
  theme(legend.position = "bottom", legend.direction = "vertical")

gg_path <- here("../output/figures/australia_sv_freq.pdf")
ggsave(gg_path, height = 5, width = 4)

# (2) Total strategic incentives

prop_df$inc_rcv <- prop_df$RCV_second + prop_df$RCV_third
prop_df$inc_plur <- prop_df$plur_second

aus_inc <- ggplot(prop_df, aes(x = inc_plur, y = inc_rcv)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0) + 
  facet_wrap(~ s) +
  theme_bw() +
  scale_x_continuous(limits = c(0, 0.7), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 0.7), expand = c(0, 0)) +
  labs(x = "Proportion of AES respondents with positive SI under Plurality",
       y = "Proportion of AES respondents with positive SI under RCV")
gg_path2 <- here("../output/figures/australia_sv_prop.pdf")
ggsave(gg_path2, aus_inc, width = 5, height = 5)

# This still seems weird. What is happening here?
inc60 <- prop_df[prop_df$s == 60, ]
sum(apply(aes_utils_raw, 1, function(x) max(x) == x[1])) / nrow(aes_utils_raw)
# Basically, all Green voters want to vote strategically in plurality, with a few exceptions.

# Check out who the outliers are
# inc90 <- prop_df[prop_df$s == 90, ]
# outliers <- which(inc90$inc_plur > 0.2)
# ggtern(simplex_df[outliers, ], aes(A, B, C)) +
#   geom_point(data = simplex_df, aes(A, B, C), colour = "light grey") +  
#   geom_point()

# (3) q-q plots

qq_mega_list <- lapply(mega_tau_list, function(x) qq_function_two(x, aes_utils_raw))
qq_mega_df <- as.data.frame(do.call(rbind, qq_mega_list))
qq_mega_df$const <- rep(c(const_bp$district), each = nrow(aes_utils_raw) * length(unlist(s_list)))

qq_mega_by_s <- split(qq_mega_df, qq_mega_df$s)
qq_agg <- lapply(qq_mega_by_s, function(z) as.data.frame(qqplot(x = z$x, y = z$y, plot.it = FALSE)))
qq_agg_df <- as.data.frame(do.call(rbind, qq_agg))
qq_agg_df$s <- rep(unlist(s_list), each = nrow(aes_utils_raw) * nrow(const_bp))

aus_qq <- ggplot(qq_mega_df, aes(x = x, y = y)) +
  geom_line(aes(x = x, y = y, group = const), alpha = 0.1) +
  geom_line(data = qq_agg_df, aes(x = x, y = y), colour = "red", lwd = 2) + 
  geom_abline(intercept = 0, slope = 1, linetype = "dotted", colour = "blue") +
  theme_bw()  + 
  facet_wrap(vars(s)) +
  xlim(-30, 30) + ylim(-30, 30)
ggsave(here("../output/figures/australia_sv_qq.pdf"), aus_qq, width = 6, height = 6)

# (4) occurence of voting paradoxes
### Note that I could do this much faster if I incorporated it into the main loop --
### that way I wouldn't have to calculate pprobs twice.

# Non-monotonicity under AV
nonmon_list <- c()
s <- 15
for(i in 1:nrow(const_bp)){
  print(i)
  v_vec <- as.numeric(const_bp[i, 2:10] / sum(const_bp[i, 2:10]))
  pprobs <- av.pivotal.event.probs.general(v_vec, rep(s, 4))
  nonmon_list[[i]] <- non_monoton(mega_tau_list[[i]][mega_tau_list[[i]]$s == s, ], pprobs)
}

# Wasted vote under Plurality
# need to decide what to do with trunc in plurality.
wasted_list <- c()
for(i in 1:nrow(const_bp)){
  print(i)
  v_vec <- as.numeric(const_bp[i, 2:10] / sum(const_bp[i, 2:10]))
  v_vec_three <- c(v_vec[1] + v_vec[2] + v_vec[7], v_vec[3] + v_vec[4] + v_vec[8], v_vec[5] + v_vec[6] + v_vec[9])
  pprobs <- plurality.pivotal.probabilities(v_vec_three, s)
  wasted_list[[i]] <- wasted_vote(mega_tau_list[[i]][mega_tau_list[[i]]$s == s, ], pprobs)
}

paradox_df <- matrix(NA, ncol = 3, nrow = nrow(const_bp))
for(i in 1:nrow(const_bp)){
  no_show <- nonmon_list[[i]]$no_show / nonmon_list[[i]]$total
  nonmon <- (nonmon_list[[i]]$nonmon1 + nonmon_list[[i]]$nonmon2) / nonmon_list[[i]]$total
  wasted <- wasted_list[[i]]$wasted / wasted_list[[i]]$total
  paradox_df[i, ] <- c(no_show, nonmon, wasted)
}
paradox_df <- as.data.frame(paradox_df)
names(paradox_df) <- c("no_show", "nonmon", "wasted")

ggplot(paradox_df, aes(x = wasted)) +
  geom_point(aes(y = no_show, colour = "No-show"), size = 1.5, alpha = 0.2) +
  geom_point(aes(y = nonmon, colour = "Non-mon"), size = 1.5, alpha = 0.2) +
  geom_abline(intercept = 0, slope = 1, lty = "dotted") +
  geom_smooth(method = "loess", aes(y = no_show), colour = "blue") +
  geom_smooth(method = "loess", aes(y = nonmon), colour = "red") +
  xlim(0, 0.4) + ylim(0, 0.4) +
  labs(x = "Pr(Wasted Vote, Plurality)", y = "Pr(Voting Paradox, RCV)", colour = "Paradox type") +
  scale_colour_manual(breaks = c("No-show", "Non-mon"), values = c("No-show" = "blue", "Non-mon" = "red")) +
  theme_bw() +
  theme(legend.position = "bottom", legend.direction = "horizontal")
ggsave(here("../output/figures/paradoxes_aus.pdf"), width = 4, height = 4)

# Check if Andy's function yields the same result...
# v_vec <- const_bp_no_trunc[4, 2:10] / sum(const_bp_no_trunc[1, 4:10])
# names(aes_utils_raw) <- c("A", "B", "C")
# test <- sv(U = aes_utils_raw, s = 85, v.vec = as.numeric(v_vec)[1:6], rule = "AV")
# 
# test_toby_p <- av.pivotal.event.probs.general(as.numeric(v_vec), rep(85, 4))
# non_monoton(mega_tau_list[[4]], test_toby_p)
# no_show_non_mon_from_sv_object(test)
# Yes! (without truncated ballots, that is.)


# (5) interdependence

# Set s = 80
mega_tau_fixed_s <- lapply(mega_tau_list, function(x) x[x$s == 85, ])

# Run loop over constituencies. This will take a long time...
inter_df <- list()
lambda_list <- as.list(seq(0, 0.5, 0.05))
# This is going to run for a long time!
for(i in 1:nrow(const_bp)){
  print(i)
  v_vec <- as.numeric(const_bp[i, 2:10]) / sum(as.numeric(const_bp[i, 2:10]))
  tau <- mega_tau_fixed_s[[i]]
  const_props <- level_two_props(v_vec, lambda_list, aes_utils_raw, tau, list(80))
  #const_props$const <- const_bp[i, 1]
  inter_df[[i]] <- const_props
}

save(inter_df, file = here("../output/interdependence.Rdata"))

inter_df_full <- as.data.frame(do.call(rbind, inter_df))
inter_df_full$const <- rep(const_bp$district, each = length(lambda_list))

l1_plot <- ggplot(inter_df_full, aes(x = lambda, group = const)) +
  geom_line(aes(y = L1RCV), colour = "blue", alpha = 0.25) +
  geom_line(aes(y = L1PLUR), colour = "orange", alpha = 0.25) +
  theme_bw()
ggsave(here("../output/figures/level1_diff.pdf"))

l0_plot <- ggplot(inter_df_full, aes(x = lambda, group = const)) +
  geom_line(aes(y = L0RCV), colour = "blue", alpha = 0.25) +
  geom_line(aes(y = L0PLUR), colour = "orange", alpha = 0.25) +
  theme_bw()
ggsave(here("../output/figures/level0_diff.pdf"))

