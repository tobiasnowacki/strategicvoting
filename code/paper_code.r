##################################################
## Project: Strategic Voting in IRV
## Script purpose: Analysing CSES data
## Date: 08/02/2019
## Author:
##################################################

##
## DEPENDENCIES
##

# Load packages
requiredPackages <- c("here", "ggplot2", "ggtern", "dplyr", "purrr", "tidyr", "lmtest", "sandwich", "plm", "extrafont", "RColorBrewer", "boot")
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
load(here("../output/big_list_2.RData"))
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
country_weight <- c()
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
  country_weight[i] <- df$VAP[1] / df$m[1]
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
neutral <- sapply(class_vec, function(x) grepl("N", x))
dm <- sapply(class_vec, function(x) grepl("DM", x))
sp <- sapply(class_vec, function(x) grepl("SP", x))

big_sv_df$neutral <- 0
big_sv_df$neutral[big_sv_df$case %in% names(big_list_na_omit)[neutral]] <- 1

big_sv_df$dm <- as.numeric(big_sv_df$case %in% names(big_list_na_omit)[dm])
big_sv_df$sp <- as.numeric(big_sv_df$case %in% names(big_list_na_omit)[sp])
big_sv_df$default <- 1

##
## PREVALENCE
##

# Indicator whether positive SI or not
big_sv_df$pos_rcv <- big_sv_df$tau_rcv > 0
big_sv_df$pos_plur <- big_sv_df$tau_plur > 0

# Create tidy dataframe for regression
prev_df <- big_sv_df[, c("s", "case", "weight_rep", "pos_rcv", "pos_plur", "neutral", "sp", "dm", "default")]
prev_df$respondent <- 1:nrow(prev_df)
prev_df <- gather(prev_df, key = "type", value = "positive", 4:5)

# Run regression(s) and store results
prev_estimates <- function(s, cond = "default", val = 1){
	prev_mod <- lm("positive ~ as.factor(type) - 1", weights = weight_rep, prev_df[prev_df$s == s & prev_df[, cond] == val, ])
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
	labs(x = "Precision of Beliefs", y = expression(paste("Pr(", tau > 0, " | s)")))
ggsave(here("../output/figures_paper/prevalence_reg.pdf"), width = 6, height = 4, device = cairo_pdf)

# Do the same thing for neutral v. non-neutral
prev_results_n <- lapply(s_list, function(x) prev_estimates(x, "neutral", 1)) %>% do.call(rbind, .)
prev_results_n$class <- "neutral"
prev_results_not_n <- lapply(s_list, function(x) prev_estimates(x, "neutral", 0)) %>% do.call(rbind, .)
prev_results_not_n$class <- "not neutral"
prev_results_by_n <- rbind(prev_results_n, prev_results_not_n)

ggplot(prev_results_by_n, aes(x = s)) +
  geom_point(aes(y = coef, color = type)) +
  geom_linerange(aes(
    ymin = coef - 1.96 * se,
    ymax = coef + 1.96 * se,
    color = type),
    lwd = 2, alpha = 0.3) +
  theme_sv() +
  ylim(0, 0.28) +
  facet_wrap(. ~ class) +
  scale_x_continuous(breaks = unlist(s_list)) +
  labs(x = "Precision of Beliefs", y = expression(paste("Pr(", tau > 0, " | s)")))
ggsave(here("../output/figures_paper/prevalence_by_n.pdf"), width = 6, height = 4, device = cairo_pdf)

# Scatterplot

# ...

# Difference between prevalences by case and s
scatter <- prev_df %>% group_by(s, case, type) %>% summarise(prop = sum(as.numeric(positive) * weight_rep) / sum(weight_rep))
diffs <- scatter %>% group_by(s, case) %>% summarise(diff = prop[2] - prop[1])
diffs$weights <- rep(country_weight / sum(country_weight), 6)

diff_quants <- diffs %>% group_by(s) %>% summarise(mu = weighted.mean(diff, weights))

prev_box <- ggplot(diffs, aes(x = as.factor(s), y = diff)) +
  geom_boxplot(aes(y = diff, weight = weights), width = 0.4) +
  geom_hline(yintercept = 0, lty = "dotted") +
  theme_sv() +
  labs(x = "s", y = "Pr(tau_irv > 0) - Pr(tau_plur > 0)")
ggsave("../output/figures_paper/prev_box.pdf", width = 5, height = 4, device = cairo_pdf)

prev_dens <- ggplot(diffs, aes(x = diff)) +
  geom_density(aes(weight = weights), fill = "grey50") +
  geom_vline(data = diff_quants, aes(xintercept = mu), lty = "dashed") +
  facet_wrap(. ~ s, nrow = 6) +
  theme_sv() +
  labs(x = "s", y = "Pr(tau_irv > 0) - Pr(tau_plur > 0)")
ggsave("../output/figures_paper/prev_dens.pdf", width = 5, height = 6, device = cairo_pdf)

##
## EXPECTED BENEFIT (I can write a function to make this smoother)
## (Right now it's just a lot of copy-paste from prevalence section)

big_sv_df$exben_rcv <- 0
big_sv_df$exben_rcv[big_sv_df$tau_rcv > 0] <- big_sv_df$tau_rcv[big_sv_df$tau_rcv > 0]
big_sv_df$exben_plur <- 0
big_sv_df$exben_plur[big_sv_df$tau_plur > 0] <- big_sv_df$tau_plur[big_sv_df$tau_plur > 0]

exben_df <- big_sv_df[, c("s", "case", "weight_rep", "exben_rcv", "exben_plur", "neutral", "dm", "sp", "default")]
exben_df$respondent <- 1:nrow(exben_df)
exben_df <- gather(exben_df, key = "type", value = "exben", 4:5)

exben_estimates <- function(s, cond = "default", val = 1){
  exben_mod <- lm("exben ~ as.factor(type) - 1", weights = weight_rep, exben_df[exben_df$s == s & exben_df[, cond] == val, ])
  exben_coef <- coeftest(exben_mod, vcov = vcovHC(exben_mod, "HC0", cluster = "respondent"))
  return(data.frame(
    coef = c(exben_coef[1, 1], exben_coef[2, 1]),
    se = c(exben_coef[1, 2], exben_coef[2, 2]),
    type = c("Plurality", "IRV"),
    s = s
    ))
}
exben_results <- lapply(s_list, function(x) exben_estimates(x))
exben_results <- do.call(rbind, exben_results)

# Plot this thing
ggplot(exben_results, aes(x = s)) +
  geom_point(aes(y = coef, color = type)) +
  geom_linerange(aes(
    ymin = coef - 1.96 * se,
    ymax = coef + 1.96 * se,
    color = type),
    lwd = 2, alpha = 0.3) +
  theme_sv() +
  ylim(0, 0.28) +
  scale_x_continuous(breaks = unlist(s_list)) +
  labs(x = "Precision of Beliefs", y = "Expected Benefit")
ggsave(here("../output/figures_paper/ex_ben.pdf"), width = 6, height = 4, device = cairo_pdf)

exben_results_n <- lapply(s_list, function(x) exben_estimates(x, "neutral", 1)) %>% do.call(rbind, .)
exben_results_n$class <- "neutral"
exben_results_not_n <- lapply(s_list, function(x) exben_estimates(x, "neutral", 0)) %>% do.call(rbind, .)
exben_results_not_n$class <- "not neutral"

exben_results_by_n <- rbind(exben_results_n, exben_results_not_n)

ggplot(exben_results_by_n, aes(x = s)) +
  geom_point(aes(y = coef, color = type)) +
  geom_linerange(aes(
    ymin = coef - 1.96 * se,
    ymax = coef + 1.96 * se,
    color = type),
    lwd = 2, alpha = 0.3) +
  theme_sv() +
  ylim(0, 0.28) +
  facet_wrap(. ~ class) +
  scale_x_continuous(breaks = unlist(s_list)) +
  labs(x = "Precision of Beliefs", y = "Expected Benefit")
ggsave(here("../output/figures_paper/ex_ben_by_neutral.pdf"), width = 6, height = 4, device = cairo_pdf)

# Do the same thing for neutral v. non-neutral preferences

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

epsilon_list <- list(0, 0.00000001, 0.0000001, 0.000001, 0.00001, 0.0001, 0.0005, 0.001, 0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 1.5, 2)
epsilon_results_1 <- lapply(epsilon_list, function(x) rcv_diff(mag_df[mag_df$s == 25, ], x, 25))
epsilon_results_1 <- do.call(rbind, epsilon_results_1)
epsilon_results_2 <- lapply(epsilon_list, function(x) rcv_diff(mag_df[mag_df$s == 85, ], x, 85))
epsilon_results_2 <- do.call(rbind, epsilon_results_2)
epsilon_results <- rbind(epsilon_results_1, epsilon_results_2)

# Plot results
epsilon_factor_scale <- ggplot(epsilon_results, aes(x = as.factor(epsilon), colour = type)) +
	geom_point(aes(y = coef)) +
	geom_linerange(aes(
		ymin = coef - 1.96 * se,
		ymax = coef + 1.96 * se),
		lwd = 2, alpha = 0.3) +
	theme_sv() +
  theme(legend.position = "bottom", legend.direction = "horizontal", axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
	facet_wrap(. ~ s) +
	labs(x = "Cost threshold", y = "Pr(Incentive > Threshold)")
ggsave("../output/figures_paper/epsilon_factor_scale_props.pdf", width = 7, height = 4, device = cairo_pdf)

##
## CONDORCET WINNER
##

cw_win_rcv <- list()
cw_win_plur <- list()

for(i in 1:length(sv_list)){
  print(i)
  v_zero <- gen_v_zero(sv_list[[i]]$sin_rcv[sv_list[[i]]$s == 85])
  v_opt_rcv <- gen_v_zero(sv_list[[i]]$opt_rcv[sv_list[[i]]$s == 85])[[1]]
  v_opt_plur <- gen_v_zero_plur(sv_list[[i]]$opt_plur[sv_list[[i]]$s == 85])
  cw_win_rcv[[i]] <- cbind(evaluate_success_of_CW_given_U_and_V.mat(U = big_list_na_omit[[i]]$U, V.mat = v_opt_rcv, V0 = v_zero[[1]], lambdas = c(.1, .2, .3, .4, .5), big_list_na_omit[[i]]$weights, rule = "AV", m = 500, M = 1000), names(big_list_na_omit)[[i]])
  cw_win_plur[[i]] <- cbind(evaluate_success_of_CW_given_U_and_V.mat(U = big_list_na_omit[[i]]$U, V.mat = v_opt_plur, V0 = v_zero[[2]], lambdas = c(.1, .2, .3, .4, .5), big_list_na_omit[[i]]$weights, rule = "plurality", m = 500, M = 1000), names(big_list_na_omit)[[i]])
}

cw_df <- return_cw_df(cw_win_rcv, cw_win_plur, c(0.1, 0.2, 0.3, 0.4, 0.5), country_weight)

cw_win_rcv_25 <- list()
cw_win_plur_25 <- list()

for(i in 1:length(sv_list)){
  print(i)
  v_zero <- gen_v_zero(sv_list[[i]]$sin_rcv[sv_list[[i]]$s == 25])
  v_opt_rcv <- gen_v_zero(sv_list[[i]]$opt_rcv[sv_list[[i]]$s == 25])[[1]]
  v_opt_plur <- gen_v_zero_plur(sv_list[[i]]$opt_plur[sv_list[[i]]$s == 25])
  cw_win_rcv_25[[i]] <- cbind(evaluate_success_of_CW_given_U_and_V.mat(U = big_list_na_omit[[i]]$U, V.mat = v_opt_rcv, V0 = v_zero[[1]], lambdas = c(.1, .2, .3, .4, .5), big_list_na_omit[[i]]$weights, rule = "AV", m = 500, M = 1000), names(big_list_na_omit)[[i]])
  cw_win_plur_25[[i]] <- cbind(evaluate_success_of_CW_given_U_and_V.mat(U = big_list_na_omit[[i]]$U, V.mat = v_opt_plur, V0 = v_zero[[2]], lambdas = c(.1, .2, .3, .4, .5), big_list_na_omit[[i]]$weights, rule = "plurality", m = 500, M = 1000), names(big_list_na_omit)[[i]])
}

cw_df_25 <- return_cw_df(cw_win_rcv_25, cw_win_plur_25, c(0.1, 0.2, 0.3, 0.4, 0.5), country_weight)

cw_joint <- rbind(cw_df, cw_df_25)
cw_joint$s <- rep(c(85, 25), each = 10)

ggplot(cw_joint, aes(x = as.factor(lambda), color = type)) +
  geom_point(aes(y = mu), position = position_dodge(0.2)) +
  geom_linerange(aes(ymin = lower, ymax = upper), lwd = 2, alpha = 0.3, position = position_dodge(0.2)) +
  theme_sv() +
  ylim(0.5, 1) +
  labs(x = "Proportion of strategic voters", y = "Pr(Condorcet Winner elected)") +
  facet_wrap(.~ s) +
  theme(legend.position = "bottom", legend.direction = "horizontal")
ggsave("../output/figures_paper/cw_agg_two_s.pdf", width = 5, height = 4, device = cairo_pdf)


# Plot
ggplot(cw_df, aes(x = as.factor(lambda), color = type)) +
  geom_point(aes(y = mu), position = position_dodge(0.2)) +
  geom_linerange(aes(ymin = lower, ymax = upper), lwd = 2, alpha = 0.3, position = position_dodge(0.2)) +
  theme_sv() +
  ylim(0.5, 1) +
  labs(x = "Proportion of strategic voters", y = "Pr(Condorcet Winner elected)") +
  theme(legend.position = "bottom", legend.direction = "horizontal")
ggsave("../output/figures_paper/cw_agg.pdf", width = 5, height = 4, device = cairo_pdf)

##
## INTERDEPENDENCE
##

s <- 85 # Set level at which to evaluate
sv_list_fixed_s <- lapply(sv_list, function(x) x[x$s == s, ])

lambda_list <- as.list(seq(0, 0.5, 0.05))

inter_list <- list()
for(i in 1:length(big_list_na_omit)){
  print(i)
  df <- big_list_na_omit[[i]]
  sv_item <- sv_list_fixed_s[[i]]
  out <- level_two_props_cses(c(df$v_vec, 0, 0, 0), lambda_list, df$U, sv_item, s, df$weights)
  out$case <- names(big_list_na_omit)[[i]]
  inter_list[[i]] <- out
}

inter_df <- do.call(rbind, inter_list)
inter_df$cntryweight <- rep(country_weight, each = 11)

boot_wmean_2 <- function(x, weight = weights, d){
  weighted.mean(x[d], weight[d])
}

get_boot_ci <- function(var){
  out <- tapply(as.numeric(inter_df[[var]]), inter_df$lambda, function(x) {
  z <- boot(x, boot_wmean_2, 1000, weight = country_weight)
  z_mu <- mean(z$t)
  ci <- boot.ci(z, type = "perc")
  return(c(z_mu, ci[[4]][4:5]))
})
  out <- as.data.frame(do.call(rbind, out))
  names(out) <- c("mu", "lower", "upper")
  out$lambda <- seq(0, 0.5, 0.05)
  return(out)
}

inter_agg <- lapply(list("L1RCV", "L0RCV", "L1PLUR", "L0PLUR"), get_boot_ci)
inter_agg <- do.call(rbind, inter_agg)
inter_agg$system <- rep(c("IRV", "Plurality"), each = 22)
inter_agg$type <- rep(c("L2 versus L1", "L2 versus L0"), each = 11, times = 2)

inter_plot <- ggplot(inter_agg, aes(x = lambda, y = mu, color = system)) +
  geom_point() +
  geom_linerange(aes(ymin = lower, ymax = upper), lwd = 2, alpha = 0.5) +
  facet_wrap(. ~ type) +
  theme_sv() +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  labs(x = "Proportion of strategic voters", y = "Proportion with different optimal vote")
ggsave("../output/figures_paper/inter.pdf", width = 7, height = 5, inter_plot, device = cairo_pdf)

##
## MULTIPLE ITERATIONS
##

# Write general function that takes original big_list item, no. of iterations, lambda, s
# returns: list of v_vec_strat after every iteration
# DF of strategic incentives after k iterations

# return_iterations <- function(object, lambda, k, s){
#   obj <- object
#   v_vec <- obj$v_vec
#   v_vec_rcv <- list(v_vec)
#   v_vec_plur <- list(v_vec)
#   for(i in 1:(k+1)){
#     item_rcv <- convert_andy_to_sv_item_two(obj$U, obj$weights, 85, v_vec_rcv[[i]])
#     v_vec_rcv[[i + 1]] <- lambda * as.numeric(table(factor(item$opt_rcv, 1:6))/nrow(item_rcv)) + (1 - lambda) * v_vec_rcv[[i]]
#   }
#   # Do the same for plurality
#   for(i in 1:(k+1)){
#     item_plur <- convert_andy_to_sv_item_two(obj$U, obj$weights, 85, v_vec_plur[[i]])
#     vec_temp <- rep(as.numeric(table(factor(item_plur$opt_plur, 1:3))/nrow(item_plur)), each = 2) /2
#     v_vec_plur[[i + 1]] <- lambda * vec_temp + (1 - lambda) * v_vec_plur[[i]]
#   }
#   return(list(v_vec_rcv, v_vec_plur, item_rcv, item_plur))
# }

# # Loop the same thing over all cases. (I'm sure there's a more elegant way...)
# sv_iter_rcv <- list()
# sv_iter_plur <- list()
# for(j in 1:length(big_list_na_omit)){
#     print(j)
#     out <- return_iterations(big_list_na_omit[[j]], 0.1, 20, 85)

#     df_rcv <- out[[3]]
#     df_rcv$case <- names(big_list_na_omit)[[j]]
#     df_rcv$weight <- big_list_na_omit[[j]]$weights
#     df_rcv$country <- substr(df_rcv$case, 1, 3)
#     df_rcv$weight_sum <- sum(big_list_na_omit[[j]]$weights)
#     df_rcv$VAP <- vap$VAP[vap$cntry == df_rcv$country[1]]
#     df_rcv$m <- vap$Freq[vap$cntry == df_rcv$country[1]]
#     df_rcv$weight_rep <- df_rcv$weight * (df_rcv$VAP / (df_rcv$weight_sum * df_rcv$m))
#     sv_iter_rcv[[j]] <- df_rcv

#     # Compute summary statistics

#     df_plur <- out[[4]]
#     df_plur$case <- names(big_list_na_omit)[[j]]
#     df_plur$weight <- big_list_na_omit[[j]]$weights
#     df_plur$country <- substr(df_plur$case, 1, 3)
#     df_plur$weight_sum <- sum(big_list_na_omit[[j]]$weights)
#     df_plur$VAP <- vap$VAP[vap$cntry == df_plur$country[1]]
#     df_plur$m <- vap$Freq[vap$cntry == df_plur$country[1]]
#     df_plur$weight_rep <- df_plur$weight * (df_plur$VAP / (df_plur$weight_sum * df_plur$m))
#     sv_iter_plur[[j]] <- df_plur
# }
# sv_iter_rcv_df <- do.call(rbind, sv_iter_rcv)
# sv_iter_plur_df <- do.call(rbind, sv_iter_plur)

# # Analysis of the resulting L20 incentives
# ## Difference between RCV and Plurality
# head(sv_iter_plur_df)
# sv_iter_rcv_df$pos <- sv_iter_rcv_df$tau_rcv > 0
# sv_iter_plur_df$pos <- sv_iter_plur_df$tau_plur > 0

# # Difference in prevalence, pooled and weighted
# weighted.mean(sv_iter_rcv_df$tau_rcv > 0, weights = sv_iter_rcv_df$weight_rep)
# weighted.mean(sv_iter_plur_df$tau_plur > 0, weights = sv_iter_plur_df$weight_rep)

# # Plot distribution of magnitudes
# pdf("iterated_magnitude.pdf")
# plot(density(log(sv_iter_rcv_df$tau_rcv[sv_iter_rcv_df$pos == TRUE]), weights = sv_iter_rcv_df$weight_rep[sv_iter_rcv_df$pos == TRUE] / sum(sv_iter_rcv_df$weight_rep[sv_iter_rcv_df$pos == TRUE])),
# xlim = c(-35, 20))

# lines(density(log(sv_iter_plur_df$tau_plur[sv_iter_plur_df$pos == TRUE]), weights = sv_iter_plur_df$weight_rep[sv_iter_plur_df$pos == TRUE] / sum(sv_iter_rcv_df$weight_rep[sv_iter_rcv_df$pos == TRUE])), col = "blue")

# title("Distribution of ln(E[tau | tau > 0])")
# dev.off()

# What I ideally want is the following output:
# for each s, prevalence, magnitude and expected benefit after each iteration (with associated standard errors)

# Ugly code: improvve functions and write documentation.

### EQUILIBRIUM ANALYSIS
iteration_return_summary <- function(object, v_vec, lambda, s){
  obj <- object
  tab <- convert_andy_to_sv_item_two(obj$U, obj$weights, s, v_vec)

  prev_rcv <- weighted.mean(tab$tau_rcv > 0, weights = obj$weights)
  mag_rcv <- weighted.mean(tab$tau_tilde_rcv[tab$tau_tilde_rcv > 0], weights = obj$weights[tab$tau_tilde_rcv > 0])
  eb_rcv <- prev_rcv * mag_rcv
  strat_vec_rcv <- as.numeric(table(factor(tab$opt_rcv, 1:6))/nrow(tab))
  new_vec_rcv <- lambda * strat_vec_rcv + (1 - lambda) * v_vec

  prev_plur <- weighted.mean(tab$tau_plur > 0, weights = obj$weights)
  weight_mag_plur <- sum(obj$weights[tab$tau_plur > 0])
  # print(weight_mag_plur)
  mag_plur <- weighted.mean(tab$tau_tilde_plur[tab$tau_tilde_plur > 0], weights = obj$weights[tab$tau_tilde_plur > 0]) / weight_mag_plur
  eb_plur <- prev_plur * mag_plur

  strat_vec_plur <- rep(as.numeric(table(factor(tab$opt_plur, 1:3))/nrow(tab)), each = 2) /2
  new_vec_plur <- lambda * strat_vec_plur + (1 - lambda) * v_vec

  return(list(rcv_sum = c(prev_rcv, mag_rcv, eb_rcv),
              rcv_vec = new_vec_rcv,
              plur_sum = c(prev_plur, mag_plur, eb_plur),
              plur_vec = new_vec_plur))
}

iteration_wrapper <- function(object, v_vec, lambda, s, k){
    rcv_sum_list <- list()
    rcv_v_vec_list <- list(v_vec)
    plur_sum_list <- list()
    plur_v_vec_list <- list(v_vec)

    # RCV loop
    for (j in 1:k) {
       out <- iteration_return_summary(object, rcv_v_vec_list[[j]], lambda, s)
       rcv_sum_list[[j]] <- out$rcv_sum
       rcv_v_vec_list[[j + 1]] <- out$rcv_vec
    }
    rcv_sum <- as.data.frame(do.call(rbind, rcv_sum_list))
    names(rcv_sum) <- c("Prevalence", "Magnitude", "ExpBenefit")
    rcv_sum$k <- 1:k
    rcv_v_vec <- as.data.frame(do.call(rbind, rcv_v_vec_list))

    # Plurality loop
    for (j in 1:k) {
       out <- iteration_return_summary(object, plur_v_vec_list[[j]], lambda, s)
       plur_sum_list[[j]] <- out$plur_sum
       plur_v_vec_list[[j + 1]] <- out$plur_vec
    }
    plur_sum <- as.data.frame(do.call(rbind, plur_sum_list))
    names(plur_sum) <- c("Prevalence", "Magnitude", "ExpBenefit")
    plur_sum$k <- 1:k
    plur_v_vec <- as.data.frame(do.call(rbind, plur_v_vec_list))

    return(list(rcv_sum, rcv_v_vec, plur_sum, plur_v_vec))
}


big_rcv_sum <- list()
big_plur_sum <- list()
big_rcv_vec <- list()
big_plur_vec <- list()

lambda <- 0.05
k <- 60

# Repeat loop for multiple values of s (list structure)
# Warning, takes a long imte!
for(prec in 6:6){
  s_val <- s_list[[prec]]
  print(s_val)
  rcv_sum <- list()
  plur_sum <- list()
  rcv_vec <- list()
  plur_vec <- list()
  for (case in 1:1) {
    print(case)
    out <- iteration_wrapper(big_list_na_omit[[case]], big_list_na_omit[[case]]$v_vec, lambda, s_val, k)
    rcv_sum[[case]] <- out[[1]]
    rcv_sum[[case]]$case <- names(big_list_na_omit)[[case]]
    plur_sum[[case]] <- out[[3]]
    plur_sum[[case]]$case <- names(big_list_na_omit)[[case]]
    rcv_vec[[case]] <- out[[2]]
    plur_vec[[case]] <- out[[4]]
  }
  rcv_sum <- do.call(rbind, rcv_sum)
  rcv_sum$s <- s_val
  plur_sum <- do.call(rbind, plur_sum)
  plur_sum$s <- s_val
  big_rcv_sum[[prec]] <- rcv_sum
  big_rcv_vec[[prec]] <- rcv_vec
  big_plur_sum[[prec]] <- plur_sum
  big_plur_vec[[prec]] <- plur_vec
}

big_rcv_sum <- do.call(rbind, big_rcv_sum)
big_plur_sum <- do.call(rbind, big_plur_sum)

# Test and check diff. with Andy's code
test <- as.data.frame(big_plur_vec[[6]])
test[, c(1, 3, 5)] + test[, c(2, 4, 6)]
big_list_na_omit[[1]]$v_vec[c(1, 3, 5)] + big_list_na_omit[[1]]$v_vec[c(2, 4, 6)]

# Aggregate

iteration_return_summary(big_list_na_omit[[124]], as.numeric(big_rcv_vec[[124]][13, ]), 0.2, 85)

eb_rcv_agg <- as.data.frame(big_rcv_sum %>% group_by(k) %>% summarise(avg = weighted.mean(ExpBenefit, weights = country_weight, na.rm = TRUE)))
eb_plur_agg <- as.data.frame(big_plur_sum %>% group_by(k) %>% summarise(avg = weighted.mean(ExpBenefit, weights = country_weight, na.rm = TRUE)))

ggplot(big_rcv_sum, aes(x = k, y = Prevalence)) + 
  geom_line(data = big_rcv_sum, aes(x = k, y = Prevalence, group = case), alpha = 0.1, colour = "blue") + 
  geom_line(data = big_plur_sum, aes(x = k, y = Prevalence, group = case), alpha = 0.1, colour = "orange") + 
  geom_line(data = rcv_agg, aes(x = k, y = prev), colour = "blue", lwd = 2) + 
  geom_line(data = plur_agg, aes(x = k, y = prev), colour = "orange", lwd = 2) + 
  facet_wrap(. ~ s) +
  theme_sv() +
  labs(x = "Iteration", y = "Prevalence")
ggsave(here("../output/figures/iterated_prevalence.pdf"), device = cairo_pdf)

ggplot(big_rcv_sum, aes(x = k, y = Magnitude)) + 
  geom_line(data = big_rcv_sum, aes(x = k, y = Magnitude, group = case), alpha = 0.1, colour = "blue") + 
  geom_line(data = big_plur_sum, aes(x = k, y = Magnitude, group = case), alpha = 0.1, colour = "orange") + 
  geom_line(data = rcv_agg, aes(x = k, y = mag), colour = "blue", lwd = 2) + 
  geom_line(data = plur_agg, aes(x = k, y = mag), colour = "orange", lwd = 2) + 
  facet_wrap(. ~ s) +
  theme_sv() +
  labs(x = "Iteration", y = "Magnitude")
ggsave(here("../output/figures/iterated_magnitude.pdf"), device = cairo_pdf)

ebplot <- ggplot(big_rcv_sum, aes(x = k, y = ExpBenefit)) + 
  geom_line(data = big_rcv_sum, aes(x = k, y = ExpBenefit, group = case), alpha = 0.1, colour = "blue") + 
  geom_line(data = big_plur_sum, aes(x = k, y = ExpBenefit, group = case), alpha = 0.1, colour = "orange") + 
  geom_line(data = rcv_agg, aes(x = k, y = eb), colour = "blue", lwd = 2) + 
  geom_line(data = plur_agg, aes(x = k, y = eb), colour = "orange", lwd = 2) + 
  facet_wrap(. ~ s) +
  theme_sv() +
  labs(x = "Iteration", y = "Expected Benefit")
ggsave(here("../output/figures/iterated_expbenefit.pdf"), ebplot, device = cairo_pdf)

# same but with zoomed-in y-axis
ebplot + ylim(0, 0.25)
ggsave(here("../output/figures/iterated_expbenefit_zoomed.pdf"), device = cairo_pdf)

# v-vec euclidean distance
bind_vecs_rcv <- function(vec_list){
  vec_list_mod <- lapply(vec_list, function(x) {
    orig <- x[1, ]
    x$dist <- c(NA, sqrt(rowSums((x[-1, 1:6] - orig[rep(1, 60), ])^2)))
    x$diff <- c(NA, sqrt(rowSums((x[-1, 1:6] - x[-61, 1:6])^2)))
    x$k <- 1:61
    return(x)
    })
  out <- do.call(rbind, vec_list_mod) 
  out$case <- rep(names(big_list_na_omit), each = 61)
  return(out)   
 } 

 bind_vecs_plur <- function(vec_list){
  vec_list_mod <- lapply(vec_list, function(x) {
    x_three <- x[, c(1, 3, 5)] + x[, c(2, 4, 6)]
    orig <- x_three[1, ]
    x_three$dist <- c(NA, sqrt(rowSums((x_three[-1, 1:3] - orig[rep(1, 60), ])^2)))
    x_three$diff <- c(NA, sqrt(rowSums((x_three[-1, 1:3] - x_three[-61, 1:3])^2)))
    x_three$k <- 1:61
    return(x_three)
    })
  out <- do.call(rbind, vec_list_mod) 
  out$case <- rep(names(big_list_na_omit), each = 61)
  return(out)   
 } 

vec_list_rcv <- lapply(big_rcv_vec, bind_vecs_rcv)
vec_rcv <- do.call(rbind, vec_list_rcv)
vec_rcv$s <- rep(as.numeric(s_list), each = 9760)
vec_rcv$type <- "IRV"
vec_list_plur <- lapply(big_plur_vec, bind_vecs_plur)
vec_plur <- do.call(rbind, vec_list_plur)
vec_plur$s <- rep(as.numeric(s_list), each = 9760)
vec_plur$type <- "Plurality"

vec_df <- rbind(vec_rcv[, c("dist", "diff", "k", "case", "s", "type")], 
                vec_plur[, c("dist", "diff", "k", "case", "s", "type")])

vec_df_agg <- as.data.frame(vec_df %>% group_by(k, s, type) %>% summarise(wtd_diff = as.numeric(weighted.mean(log(diff), weights = country_weight)), wtd_dist = as.numeric(weighted.mean(log(dist), weights = country_weight))))

ggplot(vec_df, aes(x = k, y = log(diff))) +
  geom_line(aes(group = interaction(case, type), colour = type), alpha = 0.05) +
  geom_line(data = vec_df_agg, aes(x = k, y = as.numeric(wtd_diff), group = type, colour = type), lwd = 2) +
  facet_wrap(. ~ s) +
  theme_sv() +
  labs(x = "kth Iteration", y = "log(Distance to k - 1)") +
  scale_color_manual(values = c("blue", "orange"))
ggsave(here("../output/figures/iteration_euclid_diff.pdf"), device = cairo_pdf)
  geom_line(data = eb_rcv_agg, aes(x = k, y = avg), colour = "blue", lwd = 2) + 
  geom_line(data = eb_plur_agg, aes(x = k, y = avg), colour = "orange", lwd = 2) + 
  theme_sv()
ggsave(here("../output/figures/iteration_expected_benefit.pdf"), height = 4, width = 4, device = cairo_pdf)


# list structure (how output should look like):
  # s = 25
  ## cases (1 - 160)
  ### DF of summary results for each iteration:
  ###   prevalence, magnitude, expected benefit, v_vec
  ### full item of first iteration? (to get at distribution of magnitudes)
  ### full item of last iteration? (to get at distribution of magnitudes)
  # s = 40
  ## ...

  ggplot(vec_df, aes(x = k, y = log(dist))) +
  geom_line(aes(group = interaction(case, type), colour = type), alpha = 0.05) +
  geom_line(data = vec_df_agg, aes(x = k, y = as.numeric(wtd_dist), group = type, colour = type), lwd = 2) +
  facet_wrap(. ~ s) +
  theme_sv() +
  labs(x = "kth Iteration", y = "log(Distance to first iteration)") +
  scale_color_manual(values = c("blue", "orange"))
  ggsave(here("../output/figures/iteration_euclid_dist.pdf"), device = cairo_pdf)

inspect <- cbind(vec_df[vec_df$s == 55 & vec_df$k == 60 & vec_df$type == "IRV", ], country_weight)
inspect[, 1:2] <- log(inspect[, 1:2])
inspect[order(inspect$diff), ]
weighted
