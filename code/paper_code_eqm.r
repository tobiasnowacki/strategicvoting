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
requiredPackages <- c("here", "ggplot2", "ggtern", "dplyr", "purrr", "tidyr", "lmtest", "sandwich", "plm", "extrafont", "RColorBrewer", "boot", "svMisc", "ggtern")
ipak <- function(pkg){
        new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
        if (length(new.pkg))
                install.packages(new.pkg, dependencies = TRUE)
        sapply(pkg, require, character.only = TRUE)
}
ipak(requiredPackages)

# Load functions
source(here("code/utils/functions.r"))
source(here("code/utils/av_pivotal_probs_analytical_general_v2.r"))
source(here("code/utils/plurality_pivotal_probabilities_analytical.r"))
source(here("code/utils/general_iteration_simulation_approach.r"))
source(here("code/utils/sv.r"))

# Load data
load(here("output/big_list_2.RData"))
vap <- read.csv(here("data/case_vap.csv"), sep = "") # voting age pop.

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

# Colourblind palette for plots etc.
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


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
big_list_na_omit[[39]] <- NULL
big_list_na_omit[[144]] <- NULL

# Create main sv_list object (this could be done much more elegantly with a function)
sv_list <- list()
country_weight <- c()
n <- length(big_list_na_omit)
for(i in 1:n){
  progress(i)
  if(i == n) cat("Done! \n")
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

###
### EQUILIBRIUM ANALYSIS
###

rearrange <- function(x){
  voter <- as.numeric(x[1])
  order_utils <- unlist(c(cmat_u[voter, ]))
  utils <- x[26:28]
  utils <- utils[order_utils]
  order_pprobs <- unlist(c(cmat_plur[voter, ]))
  order_pprobs_rcv <- unlist(c(cmat_rcv[voter, ]))
  pprobs_plur <- x[22:24]
  pprobs_plur <- unlist(c(pprobs_plur[order_pprobs]))
  pprobs_rcv <- x[10:21]
  pprobs_rcv <- unlist(c(pprobs_rcv[order_pprobs_rcv]))
  out <- unlist(c(utils, pprobs_plur, pprobs_rcv))
  names(out) <- c("uA", "uB", "uC", "ABpo", "ACpo", "BCpo", "ABo", "ACo", "BCo", "AB_ABo", "AB_ACo", "AB_CBo", "AC_ACo", "AC_ABo", "AC_BCo", "BC_BCo", "BC_BAo", "BC_ACo")
  return(out)
}

big_rcv_sum <- list()
big_plur_sum <- list()
big_rcv_vec <- list()
big_plur_vec <- list()
piv_ratio_rcv <- list()
piv_ratio_plur <- list()

# Choice parameters
lambda <- 0.05
k <- 60

# Repeat loop for multiple values of s (list structure)
# Warning, takes a long time!
for(prec in c(1, 4, 6)){
  s_val <- s_list[[prec]]
  cat(paste0("=== s = ", s_val, "=============== \n"))
  rcv_sum <- list()
  plur_sum <- list()
  rcv_vec <- list()
  plur_vec <- list()
  rcv_piv <- list()
  plur_piv <- list()
  for (case in 1:160) {
    cat(paste0(case, ": ", names(big_list_na_omit)[case], "   "))
    out <- many_iterations(big_list_na_omit[[case]], big_list_na_omit[[case]]$v_vec, lambda, s_val, k)
    rcv_sum[[case]] <- out[[1]] %>% mutate(case = names(big_list_na_omit)[[case]])
    plur_sum[[case]] <- out[[3]] %>% mutate(case = names(big_list_na_omit)[[case]])
    rcv_vec[[case]] <- out[[2]]
    plur_vec[[case]] <- out[[4]]
    rcv_piv[[case]] <- piv_ratio(out[[1]])
    plur_piv[[case]] <- piv_ratio(out[[3]])
  }
  rcv_sum <- do.call(rbind, rcv_sum)
  plur_sum <- do.call(rbind, plur_sum)
  big_rcv_sum[[prec]] <- rcv_sum
  big_rcv_vec[[prec]] <- rcv_vec
  big_plur_sum[[prec]] <- plur_sum
  big_plur_vec[[prec]] <- plur_vec
}

save.image(here("output/intermediate2.Rdata"))
# load(here("output/intermediate2.RData"))

big_rcv_sum <- do.call(rbind, big_rcv_sum)
big_plur_sum <- do.call(rbind, big_plur_sum)


## QUANTITIES OF INTEREST

## CONJECTURE 1
# Summed pivotal probabilities
rcv_summary <- big_rcv_sum %>% group_by(case, s, iter) %>% summarise_at(vars(AB:BCp), first) 
rcv_summary$rcv_sum <- rowSums(rcv_summary[, 4:16])
rcv_summary$plur_sum <- rowSums(rcv_summary[, 7:19])
rcv_summary <- rcv_summary %>% mutate(diff = plur_sum - rcv_sum)

ggplot(rcv_summary, aes(iter, diff)) +
  geom_line(aes(group = interaction(case, s)), alpha = 0.05) +
  geom_hline(yintercept = 0, lty = "dotted") +
  facet_wrap(. ~ s) +
  theme_sv() +
  labs(x = "Iteration", y = "Pr(Pivotal | Plur) - Pr(Pivotal | RCV)")
ggsave(here("output/figures/conj1_rcv.pdf"), device = cairo_pdf, height = 3, width = 6)

# Do the same thing for plurality
plur_summary <- big_plur_sum %>% group_by(case, s, iter) %>% summarise_at(vars(AB:BCp), first) 
plur_summary$rcv_sum <- rowSums(plur_summary[, 4:16])
plur_summary$plur_sum <- rowSums(plur_summary[, 7:19])
plur_summary <- plur_summary %>% mutate(diff = plur_sum - rcv_sum)

ggplot(plur_summary, aes(iter, diff)) +
  geom_line(aes(group = interaction(case, s)), alpha = 0.05) +
  geom_hline(yintercept = 0, lty = "dotted") +
  facet_wrap(. ~ s) +
  theme_sv() +
  labs(x = "Iteration", y = "Pr(Pivotal | Plur) - Pr(Pivotal | RCV)")
ggsave(here("output/figures/conj1_plur.pdf"), device = cairo_pdf, height = 3, width = 6)

## CONJECTURE 2

cmat_plur <- matrix(c(1, 2, 3, 2, 1, 3, 1, 3, 2, 3, 1, 2, 2, 3, 1, 3, 2, 1), byrow = T, ncol = 3)

cmat_rcv <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 2, 1, 3, 7, 8, 9, 4, 5, 6, 10, 12, 11, 1, 3, 2, 4, 6, 5, 10, 11, 12, 7, 8, 9, 3, 1, 2, 10, 11, 12, 4, 6, 5, 7, 9, 8, 2, 3, 1, 7, 9, 8, 10, 12, 11, 4, 5, 6, 3, 2, 1, 10, 12, 11, 7, 9, 8, 4, 6, 5), byrow = T, ncol = 12)

cmat_u <- matrix(c(1, 2, 3, 1, 3, 2, 2, 1, 3, 2, 3, 1, 3, 1, 2, 3, 2, 1), byrow = T, ncol = 3)

testdf <- big_rcv_sum %>% filter(s == 85)
ordered_u_probs <- t(apply(testdf, 1, rearrange))
ordered_u_probs <- apply(ordered_u_probs, 2, as.numeric)

testdf2 <- as.data.frame(ordered_u_probs)
testdf2 <- testdf2 %>% mutate(ben_p = BCpo * ((uB - uC) / 2), cost_p = ABpo * (uA - uB) + ACpo * ((uB - uC)/2),
ben_rcv = AB_CBo * (uB - uC) + BC_BCo * ((uB - uC)/2) + BC_ACo * ((uA - uC)/2), cost_rcv = ABo * (uA - uB) + AB_ABo * (uA - uB) + AB_ACo * (uA - uC) + AC_ACo * ((uA - uC) / 2) + AC_ABo * ((uA - uB) / 2) + AC_BCo * ((uB - uC)/2) + BC_BAo * ((uA - uB)/2))

# compute correlations
testdf3 <- cbind(testdf, testdf2)
testdf3 <- testdf3 %>% group_by(case, iter) %>% summarise(corr_plur = cor(ben_p, cost_p), corr_rcv = cor(ben_rcv, cost_rcv))

ggplot(testdf3, aes(iter, corr_rcv)) +
  geom_line(aes(group = case), alpha = 0.05) +
  geom_line(aes(y = corr_plur, group = case), colour = "blue", alpha = 0.05) +
  geom_hline(yintercept = 0, lty = "dotted") +
  theme_sv() +
  labs(x = "Iteration", y = "Correlation")
ggsave(here("output/figures/conj2_rcv.pdf"), device = cairo_pdf, height = 3, width = 6)


## SUMMARY OF INCENTIVES
# Prevalence, Magnitude, EB
big_rcv_sum %>% group_by(case, iter) %>% summarise(prev = mean(tau_rcv > 0), mag = mean(tau_rcv[tau_rcv > 0]), eb = prev * mag) 





###
### EUCLIDEAN DISTANCES
###

euclid <- function(x) {
  # Euclidean distance between one iteration and next
  if (isdf <- is.data.frame(x)) {
    x <- data.matrix(x)
  }
  dij <- c(NA, sqrt(rowSums((tail(x, -1) - head(x, -1))^2)))
  x <- cbind(x, diff = dij)
  if (isdf) {
    x <- as.data.frame(x)
  }
  x
}

euclid_first <- function(df){
  # Euclidean distance between one iteration and first
  df_follow <- df[2:nrow(df), ]
  dist <- apply(df_follow, 1, function(x) sqrt(sum((x - df[1, ])^2)))
  return(dist)
}

euclid_together <- function(df){
  a <- euclid(df)
  b <- c(NA, euclid_first(df))
  return(cbind(a, b))
}

rcv_vec <- lapply(big_rcv_vec[c(1, 4, 6)], function(x) lapply(x, function(y){
  euclid_together(y)
}))
rcv_vec <- lapply(rcv_vec, function(x) {do.call(rbind, x) %>% mutate(iter = rep(1:61, 160), case = rep(names(big_list_na_omit), each = 61))})
rcv_vec <- do.call(rbind, rcv_vec)
rcv_vec <- rcv_vec %>% mutate(s = rep(c(10, 55, 85), each = 9760), system = "RCV")


plur_vec <- lapply(big_plur_vec[c(1, 4, 6)], function(x) lapply(x, function(y){
  euclid_together(y)
}))
plur_vec <- lapply(plur_vec, function(x) {do.call(rbind, x) %>% mutate(iter = rep(1:61, 160), case = rep(names(big_list_na_omit), each = 61))})
plur_vec <- do.call(rbind, plur_vec)
plur_vec <- plur_vec %>% mutate(s = rep(c(10, 55, 85), each = 9760), system = "Plurality")

dist_df <- rbind(rcv_vec, plur_vec)

# Figure: Euclidean distances from one vector to the next iteration
ggplot(dist_df, aes(x = iter, y = diff)) +
  geom_line(aes(group = interaction(case, s, system), colour = system), alpha = 0.05) + 
  facet_wrap(. ~ s) +
  ylim(0, 0.05) + 
  labs(x = "Iteration", y = "Euclidean Distance") + 
  scale_color_manual(values = c("orange", "blue")) + 
  theme_sv() +
  theme(legend.position = "bottom", legend.direction = "horizontal")
ggsave(here("output/figures/euclidean.pdf"), device = cairo_pdf, height = 3, width = 6)

# Figure: Euclidean distances from one vector to the first iteration
ggplot(dist_df, aes(x = iter, y = b)) +
  geom_line(aes(group = interaction(case, s, system), colour = system), alpha = 0.05) + 
  facet_wrap(. ~ s) +
  labs(x = "Iteration", y = "Euclidean Distance") + 
  scale_color_manual(values = c("orange", "blue")) + 
  theme_sv() +
  theme(legend.position = "bottom", legend.direction = "horizontal")
ggsave(here("output/figures/euclidean.pdf"), device = cairo_pdf, height = 3, width = 6)

# Statistic: Avg change at first and last iteration
dist_df %>% filter(iter %in% c(2, 61)) %>% group_by(system, s, iter) %>% summarise(mean(diff))


# Simplex paths
rcv_vec_df <- big_rcv_vec[[6]] %>% do.call(rbind, .) %>% mutate(case = rep(1:160, each = 61), iter = rep(1:61, 160), A = V1 + V2, B = V3 + V4, C = V5 + V6)
plur_vec_df <- big_plur_vec[[6]] %>% do.call(rbind, .) %>% mutate(case = rep(1:160, each = 61), iter = rep(1:61, 160), A = V1 + V2, B = V3 + V4, C = V5 + V6)
line_df <- data.frame(A = c(0.5, 1/3), B = c(0.5, 1/3), C = c(0, 1/3))

ggtern(plur_vec_df, aes(A, B, C)) + 
  geom_line(aes(group = case), alpha = 0.1) + 
  geom_line(data = line_df, lty = "dashed", colour = "red") +
  geom_point(data = plur_vec_df[plur_vec_df$iter == 1, ], alpha = 0.5) + 
  geom_point(data = plur_vec_df[plur_vec_df$iter == 61, ], size = 0.5, colour = "blue") + 
  theme_sv()
ggsave(here("output/figures/tatonnement_plur.pdf"), device = cairo_pdf, width = 3, height = 3)

ggtern(rcv_vec_df, aes(A, B, C)) + 
  geom_line(aes(group = case), alpha = 0.1) + 
  geom_line(data = line_df, lty = "dashed", colour = "red") +
  geom_point(data = rcv_vec_df[plur_vec_df$iter == 1, ], alpha = 0.5) + 
  geom_point(data = rcv_vec_df[plur_vec_df$iter == 61, ], size = 0.5, colour = "blue") + 
  theme_sv()
ggsave(here("output/figures/tatonnement_rcv.pdf"), device = cairo_pdf, width = 3, height = 3)

###
### EXPECTED BENEFIT AND OTHERS
### 

# Aggregate means
rcv_agg <- as.data.frame(big_rcv_sum %>% group_by(k, s) %>% summarise(
  Prevalence = weighted.mean(Prevalence, weights = country_weight, na.rm = TRUE),
  Magnitude = weighted.mean(Magnitude, weights = country_weight, na.rm = TRUE),
  ExpBenefit = weighted.mean(ExpBenefit, weights = country_weight, na.rm = TRUE)
)) %>% gather(., key = "Type", value = "Statistic", 3:5) %>% mutate(System = "IRV")

plur_agg <- as.data.frame(big_plur_sum %>% group_by(k, s) %>% summarise(  
  Prevalence = weighted.mean(Prevalence, weights = country_weight, na.rm = TRUE),
  Magnitude = weighted.mean(Magnitude, weights = country_weight, na.rm = TRUE),
  ExpBenefit = weighted.mean(ExpBenefit, weights = country_weight, na.rm = TRUE)
)) %>% gather(., key = "Type", value = "Statistic", 3:5) %>% mutate(System = "Plurality")
agg_df <- rbind(rcv_agg, plur_agg)


# Tidying data
big_rcv_sum <- gather(big_rcv_sum, key = "Type", value = "Statistic", 1:3) %>% mutate(System = "IRV")

big_plur_sum <- gather(big_plur_sum, key = "Type", value = "Statistic", 1:3) %>% mutate(System = "Plurality")

big_df <- rbind(big_rcv_sum, big_plur_sum)


ggplot(big_df, aes(x = k, y = Statistic)) +
    geom_line(aes(color = System, group = interaction(System, case)), alpha = 0.1) +
    geom_line(data = agg_df, aes(color = System), lwd = 2) + 
    facet_wrap(s ~ Type, scales = "free") +
    theme_sv() +
    labs(x = "Iteration", y = "Statistic") +
    scale_color_manual(values = cbbPalette[c(3, 2)]) +
    theme(legend.position = "bottom", legend.direction = "horizontal")
ggsave(here("output/figures/iterated_complete.pdf"), device = cairo_pdf, width = 6, height = 5)


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
ggsave(here("../output/figures/iterated_magnitude2.pdf"), device = cairo_pdf)

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
