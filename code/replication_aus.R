##################################################
## Project: Strategic Voting in RCV
## Script purpose: Replication of Andy's plots
##					 in the Australia case
## Date: 30/10/2018
## Author:
##################################################

# Set WD etc.
library(here)

av_piv_path <- here("code/utils/av_pivotal_probs_analytical_general_v2.r")
source(av_piv_path)
functions <- here("code/utils/functions.r")
source(functions)

library(ggplot2)
library(reshape2)

# Import AES and Ballot data
# Note that these are normalised utilities. To obtain like-dislike scores I will need to re-run the original script.
aes_utils <- read.csv(here("data", "australia", "AES_utility.csv"))[, -1]
aes_utils <- aes_utils[, c("GRN", "LIB", "LAB")] #common ordering

nsw <- read.csv(here("data/australia/nsw_ballots.csv"))[, -1]
resampling <- read.csv(here("data/australia/nsw_resampling.csv"))[, -1]

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

const_opt_dist_df <- do.call(rbind, const_opt_dist)
const_opt_dist_df <- melt(const_opt_dist_df, id.vars = c("const", "s"))

# Plot.

hist(const_opt_dist_df$value[const_opt_dist_df$variable == "third" ])

ggplot(const_opt_dist_df, aes(x = s, y = value)) +
	geom_line(aes(colour = variable, group = interaction(const, variable)), alpha = 0.3) +
	labs(x = "Information (s)", 
		y = "Proportion of voters in AES casting ballot type",
		colour = "Sincere pref. as first on ballot") +
	theme_bw() +
	theme(legend.position = "bottom")