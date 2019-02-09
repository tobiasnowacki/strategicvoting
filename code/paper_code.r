##################################################
## Project: Strategic Voting in IRV
## Script purpose: Analysing CSES data
## Date: 08/02/2019
## Author: 
##################################################

##
## DEPENDENCIES
##

# Set WD etc.
library(here)

# Load packages
requiredPackages <- c("ggplot2", "ggtern", "dplyr", "purrr", "tidyr", "lmtest", "sandwich", "plm", "extrafont", "RColorBrewer")
ipak <- function(pkg){
        new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
        if (length(new.pkg))
                install.packages(new.pkg, dependencies = TRUE)
        sapply(pkg, require, character.only = TRUE)
}
ipak(requiredPackages)

# Load functions
source(here("utils/functions.r"))
source(here("utils/av_pivotal_probs_analytical_general_v2.r"))
source(here("utils/plurality_pivotal_probabilities_analytical.r"))
source(here("utils", "general_iteration_simulation_approach.r"))
source(here("utils/sv.r"))

# Load data
load(here("../output/cses_big_list_2.RData"))
vap <- read.csv(here("../data/case_vap.csv"), sep = "") # voting age pop.

# Load fonts
font_import()

# Define ggplot theme
theme_sv <- function(){
  theme_bw(base_size=11, base_family="Roboto Light") %+replace%
  theme(
    panel.grid.major =  element_line(
      colour = "grey50", 
      size = 0.2,
      linetype = "dotted"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "grey97"),
    plot.margin = unit(c(0.2, 1, 0.2, 1), "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.title = element_text(size = 10, family = "Roboto Medium", face = "bold"),
    strip.background = element_rect(fill= NULL, colour = "white", linetype = NULL), 
    strip.text = element_text(colour = 'grey50', size = 9, vjust = 0.5, family = "Roboto Medium")
  )
}

##
## DATA CLEANING AND PREP
##

# Create list with sincere v_vecs from CSES utility dfs:
big_list_na_omit <- lapply(big_list, function(x) remove_nas(x))
sin_vote_list <- lapply(big_list_na_omit, function(x) sincere.vote.mat.from.U(x$U, rule = "AV"))
v_vec_list <- list()
# Append to big list
for(i in 1:length(sin_vote_list)){
  weights <- big_list_na_omit[[i]]$weights
  v_vec <- ballot.props.from.vote.mat.and.weights(sin_vote_list[[i]], weights)
  big_list_na_omit[[i]]$v_vec <- as.numeric(v_vec)
}

# Set list of s values for future analysis
s_list <- as.list(c(10, 25, 40, 55, 70, 85))

# Drop NA cases (LTU, BEL)
names(big_list)[[39]]
names(big_list)[[145]]
big_list_na_omit[[39]] <- NULL
big_list_na_omit[[144]] <- NULL

# Create main sv_list object (this could be done much more elegantly with a function)
sv_list <- list()
for(i in 1:length(big_list_na_omit)){
  print(i)
  this_list <- big_list_na_omit[[i]]
  df_list <- lapply(s_list, function(x) convert_andy_to_sv_item_two(this_list$U, this_list$weights, x, this_list$v_vec))
  df <- as.data.frame(do.call(rbind, df_list))
  df$case <- names(big_list_na_omit)[[i]]
  df$weight <- big_list_na_omit[[i]]$weights
  df$country <- substr(df$case, 1, 3)
  df$weight_sum <- sum(big_list_na_omit[[i]]$weights)
  df$VAP <- vap$VAP[vap$cntry == df$country[1]]
  df$m <- vap$Freq[vap$cntry == df$country[1]]
  df$weight_rep <- df$weight * (df$VAP / (df$weight_sum * df$m))
  #df <- apply(df, 2, as.numeric)
  sv_list[[i]] <- df
}

# Put everything together into one big DF
big_sv_df <- do.call(rbind, sv_list)

##
## DESCRIPTIVE
##

simplex_df <- matrix(NA, nrow = 160, ncol = 3) # FP (3-item)
v_vec_df <- matrix(NA, nrow = 160, ncol = 6) # 6-item

# Check first preferences
for(i in 1:length(big_list_na_omit)){
  v_vec <- big_list_na_omit[[i]]$v_vec
  simplex_vec <- c(v_vec[1] + v_vec[2], v_vec[3] + v_vec[4], v_vec[5] + v_vec[6])
  v_vec_df[i, ] <- as.numeric(v_vec)
  simplex_df[i, ] <- simplex_vec
}
simplex_df <- as.data.frame(simplex_df)
names(simplex_df) <- c("A", "B", "C")

# Check second preferences	
second_prefs <- data.frame(mAB = v_vec_df[, 1] / (v_vec_df[, 1] + v_vec_df[, 2]), mBA = v_vec_df[, 3] / (v_vec_df[, 3] + v_vec_df[, 4]), mCB = v_vec_df[, 6] / (v_vec_df[, 5] + v_vec_df[, 6]))

# Get classification
class_vec <- apply(v_vec_df, 1, classify.vec)

##
## PREVALENCE
##

# Indicator whether positive SI or not
big_sv_df$pos_rcv <- big_sv_df$tau_rcv > 0
big_sv_df$pos_plur <- big_sv_df$tau_plur > 0

# Create tidy dataframe for regression
prev_df <- big_sv_df[, c("s", "case", "weight_rep", "pos_rcv", "pos_plur")]
prev_df$respondent <- 1:nrow(prev_df)
prev_df <- gather(prev_df, key = "type", value = "positive", 4:5)

# Run regression(s) and store results
prev_estimates <- function(s){
	prev_mod <- lm("positive ~ as.factor(type) - 1", weights = weight_rep, prev_df[prev_df$s == s, ])
	prev_coef <- coeftest(prev_mod, vcov = vcovHC(prev_mod, "HC0", cluster = "respondent"))
	return(data.frame(
		coef = c(prev_coef[1, 1], prev_coef[2, 1]), 
		se = c(prev_coef[1, 2], prev_coef[2, 2]),
		type = c("Plurality", "IRV"), 
		s = s
		))
}
prev_results <- lapply(s_list, function(x) prev_estimates(x))
prev_results <- do.call(rbind, prev_results)

# Plot this thing
ggplot(prev_results, aes(x = s)) + 
	geom_point(aes(y = coef, color = type)) +
	geom_linerange(aes(
		ymin = coef - 1.96 * se, 
		ymax = coef + 1.96 * se, 
		color = type), 
		lwd = 2, alpha = 0.3) +
	theme_sv() + 
	ylim(0, 0.28) +
	scale_x_continuous(breaks = unlist(s_list)) +
	labs(x = "s", y = expression(paste("Pr(", tau > 0, " | s)")))
ggsave(here("../output/figures_paper/prevalence_reg.pdf"), width = 6, height = 4, device = cairo_pdf)

# Scatterplot?

##
## MAGNITUDE
##

# Create data frame from results
mag_df <- big_sv_df[big_sv_df$s %in% c(25, 85), c("tau_rcv", "tau_plur", "weight_rep", "case", "s")]
mag_df$respondent <- 1:nrow(mag_df)
mag_df <- gather(mag_df, type, tau, 1:2)

# Run regression(s) and store results
# Need to change function such that df is seen as variable!
rcv_diff <- function(df, epsilon, s){
  df[, "above_epsilon"] <- df[, "tau"] > epsilon 
  weighting <- df[, "weight_rep"]
  # return(length(df[, "weight_rep"]))
  model <- lm("above_epsilon ~ as.factor(type) - 1", data = df, weights = mag_df[mag_df$s == s, "weight_rep"])
  result <- coeftest(model, vcov = vcovHC(model, type = "HC0", cluster = "respondent"))
  return(data.frame(coef = result[1:2, 1], se = result[1:2, 2], type = c("Plurality", "IRV"), s = s, epsilon = epsilon))
}

epsilon_list <- list(0.00000001, 0.0000001, 0.000001, 0.00001, 0.0001, 0.0005, 0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.5, 2)
epsilon_results_1 <- lapply(epsilon_list, function(x) rcv_diff(mag_df[mag_df$s == 25, ], x, 25))
epsilon_results_1 <- do.call(rbind, epsilon_results_1)
epsilon_results_2 <- lapply(epsilon_list, function(x) rcv_diff(mag_df[mag_df$s == 85, ], x, 85))
epsilon_results_2 <- do.call(rbind, epsilon_results_2)
epsilon_results <- rbind(epsilon_results_1, epsilon_results_2)

# Plot results
epsilon_true_scale <- ggplot(epsilon_results, aes(x = as.factor(epsilon), colour = type)) +
	geom_point(aes(y = coef)) +
	geom_linerange(aes(
		ymin = coef - 1.96 * se, 
		ymax = coef + 1.96 * se), 
		lwd = 2, alpha = 0.3) +
	theme_sv() +
	facet_wrap(. ~ s) +
	labs(x = "Cost threshold, c", y = "Proportion")