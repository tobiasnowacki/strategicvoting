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

### EQUILIBRIUM ANALYSIS
big_rcv_sum <- list()
big_plur_sum <- list()
big_rcv_vec <- list()
big_plur_vec <- list()

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
  for (case in 1:160) {
    cat(paste0(case, ": ", names(big_list_na_omit)[case], "   "))
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

load(here("../output/intermediate.RData"))

###
### EUCLIDEAN DISTANCES
###

# v_vec convergence
euclid <- function(x) {
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

rcv_vec <- lapply(big_rcv_vec[c(1, 4, 6)], function(x) {do.call(rbind, x) %>% euclid %>% mutate(iter = rep(1:61, 160), case = rep(names(big_list_na_omit), each = 61))})
rcv_vec <- do.call(rbind, rcv_vec)
rcv_vec <- rcv_vec %>% mutate(s = rep(c(10, 55, 85), each = 9760), system = "RCV")

plur_vec <- lapply(big_plur_vec[c(1, 4, 6)], function(x) {do.call(rbind, x) %>% euclid %>% mutate(iter = rep(1:61, 160), case = rep(names(big_list_na_omit), each = 61))})
plur_vec <- do.call(rbind, plur_vec)
plur_vec <- plur_vec %>% mutate(s = rep(c(10, 55, 85), each = 9760), system = "Plurality")

dist_df <- rbind(rcv_vec, plur_vec)

ggplot(dist_df, aes(x = iter, y = diff)) +
  geom_line(aes(group = interaction(case, s, system), colour = system), alpha = 0.05) + 
  facet_wrap(. ~ s) +
  ylim(0, 0.05) + 
  labs(x = "Iteration", y = "Euclidean Distance") + 
  scale_color_manual(values = c("orange", "blue")) + 
  theme_sv()
ggsave(here("../output/figures/euclidean.pdf"), device = cairo_pdf)

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
ggsave(here("../output/figures/tatonnement_plur.pdf"), device = cairo_pdf)

ggtern(rcv_vec_df, aes(A, B, C)) + 
  geom_line(aes(group = case), alpha = 0.1) + 
  geom_line(data = line_df, lty = "dashed", colour = "red") +
  geom_point(data = rcv_vec_df[plur_vec_df$iter == 1, ], alpha = 0.5) + 
  geom_point(data = rcv_vec_df[plur_vec_df$iter == 61, ], size = 0.5, colour = "blue") + 
  theme_sv()
ggsave(here("../output/figures/tatonnement_rcv.pdf"), device = cairo_pdf)

###
### EXPECTED BENEFIT AND OTHERS
### 

# Aggregate means
rcv_agg <- as.data.frame(big_rcv_sum %>% group_by(k) %>% summarise(
  Prevalence = weighted.mean(Prevalence, weights = country_weight, na.rm = TRUE),
  Magnitude = weighted.mean(Magnitude, weights = country_weight, na.rm = TRUE),
  ExpBenefit = weighted.mean(ExpBenefit, weights = country_weight, na.rm = TRUE)
)) %>% gather(., key = "Type", value = "Statistic", 2:4) %>% mutate(System = "IRV")

plur_agg <- as.data.frame(big_plur_sum %>% group_by(k) %>% summarise(  
  Prevalence = weighted.mean(Prevalence, weights = country_weight, na.rm = TRUE),
  Magnitude = weighted.mean(Magnitude, weights = country_weight, na.rm = TRUE),
  ExpBenefit = weighted.mean(ExpBenefit, weights = country_weight, na.rm = TRUE)
)) %>% gather(., key = "Type", value = "Statistic", 2:4) %>% mutate(System = "Plurality")
agg_df <- rbind(rcv_agg, plur_agg)


# Tidying data
big_rcv_sum <- gather(big_rcv_sum, key = "Type", value = "Statistic", 1:3) %>% mutate(System = "IRV")
big_plur_sum <- gather(big_plur_sum, key = "Type", value = "Statistic", 1:3) %>% mutate(System = "Plurality")

big_df <- rbind(big_rcv_sum, big_plur_sum)


ggplot(big_df %>% filter(s == 85), aes(x = k, y = Statistic)) +
    geom_line(aes(color = System, group = interaction(System, case)), alpha = 0.1) +
    geom_line(data = agg_df, aes(color = System), lwd = 2) + 
    facet_wrap(. ~ Type) +
    theme_sv() +
    labs(x = "Iteration", y = "Statistic") +
    scale_color_manual(values = cbbPalette) +
    theme(legend.position = "bottom", legend.direction = "horizontal")





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
