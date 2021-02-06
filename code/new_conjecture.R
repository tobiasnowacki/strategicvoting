# PREPARE EVERYTHING -------------------------------------------
# Load dependencies
library(tidyverse)
library(pivotprobs)
library(gtools)
library(doParallel)
library(Hmisc)
source("code/utils/new_sv_iter.R")
source("code/utils/sv_theme_template.R")
source("code/utils/new_dist_helpers.R")
source("code/utils/new_plot_helpers.R")

# Load data
load("output/big_list_2.RData")
vap <- read.csv("data/case_vap.csv", sep = "") # voting age pop.
cat("Data imported. \n")
source("code/prep_cses.R")  # data prep

which_cases = 1:160 

# Get case names
nn = names(big_list_na_omit)

# Get command arguments / values
cmd_line_args <- commandArgs(trailingOnly = TRUE)
cat(cmd_line_args, sep = "n")
ifelse(length(cmd_line_args >= 1),
       s <- as.numeric(cmd_line_args[1]),
       s <- 85)
ifelse(length(cmd_line_args >= 2),
       lambda<- as.numeric(cmd_line_args[2]),
       lambda <- 1)

# Helper functions
fpath = function(lambda, s, ext){
  paste0("output/files/", lambda, "/", s, "_", ext, ".Rdata")
}
ppath = function(lambda, s, ext){
  paste0("output/figures/", lambda, "/", s, "_", ext, ".pdf")
}

# Function to load all cases for param values
return_obj = function(case, lambda, s){
  load(paste0("output/files/", lambda, "/", s, "_", case, "_iterout.Rdata"))
  return(inner_list)
}
# Bind all cases together into one DF
bind_together = map(nn[1:160], ~ return_obj(.x, lambda, s))
names(bind_together) = nn[1:160]

# Load vvec data
load("output/files/1/85_vvec.Rdata")

# Input: DF of utilities
# Output: Vector of sincere preferences (as scalars) -- under AV
sin_vote_scalar <- function(util){
  # Input: DF of utilities
  # Output: Vector of sincere preferences (as scalars) -- under AV
  stopifnot(length(util) == 3)
  max <- which(util == max(util))
  min <- which(util == min(util))
  if(length(max) > 1){
    print(util)
  }
  if (max == 1){
    if (min == 3){
      out <- 1
    }
    if (min == 2){
      out <- 2
    }
  }
  if (max == 2){
    if (min == 3){
      out <- 3
    }
    if (min == 1){
      out <- 4
    }
  }
  if (max == 3){
    if (min == 2){
      out <- 5
    }
    if (min == 1){
      out <- 6
    }
  }
  return(out)
}

# Function to grab relevant info (pivot probs, utilities) from each RCV case
grab_stuff = function(x, this.is.U, is.plur = FALSE){
  names(this.is.U) = c("uA", "uB", "uC")
  pp = matrix(unlist(x$piv.probs), 
              nrow = nrow(this.is.U), 
              ncol = length(x$piv.probs), 
              byrow = T, 
              dimnames = list(NULL, names(unlist(x$piv.probs))))
  voter = apply(this.is.U, 1, sin_vote_scalar) # adapt so it yields range 1:6 in plurality, too.
  mat = cbind(pp, this.is.U, x$V0, voter, x$weights) %>% as.data.frame 
  
  if(is.plur == TRUE){
    mat = mat %>%
      select(c(-b_a, -c_a, -c_b))
    # mat = mat[, c(2, 3, 1, 4:10)] # rearrange to fit later code
  }
  
  return(mat)
}

# LOAD ITERATIONS DATA -------------------------------------------------------
# Get data from iterations
big_rcv_sum <- list()
big_plur_sum <- list()
for(i in which_cases){
  # RCV cases
  mini_list1 <- map_dfr(bind_together[[i]]$rcv[1:61], 
                        grab_stuff, 
                        this.is.U = big_list_na_omit[[i]]$U, 
                        .id = "iter") 
  mini_list1$case <- names(big_list_na_omit)[i]
  mini_list1$system <- "IRV"
  big_rcv_sum[[i]] <- mini_list1
  
  # Plurality cases
  mini_list2 <- map_dfr(bind_together[[i]]$plur[1:61], 
                        grab_stuff,
                        is.plur = TRUE,
                        this.is.U = big_list_na_omit[[i]]$U, 
                        .id = "iter")
  mini_list2$case <- names(big_list_na_omit)[i]
  mini_list2$system <- "Plurality"
  big_plur_sum[[i]] <- mini_list2
}

big_rcv_sum <- do.call(rbind, big_rcv_sum)
big_plur_sum <- do.call(rbind, big_plur_sum)

# DATA OBJECTS FOR ANALYSIS -------------------------------------------------------
#   ## CONJECTURE 1a
cmat_plur <- matrix(c(1, 2, 3, 
                      2, 1, 3, 
                      1, 3, 2, 
                      3, 1, 2, 
                      2, 3, 1, 
                      3, 2, 1), byrow = T, ncol = 3)

cmat_rcv <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 
                     2, 1, 3, 7, 8, 9, 4, 5, 6, 10, 12, 11, 
                     1, 3, 2, 4, 6, 5, 10, 11, 12, 7, 8, 9, 
                     3, 1, 2, 10, 11, 12, 4, 6, 5, 7, 9, 8, 
                     2, 3, 1, 7, 9, 8, 10, 12, 11, 4, 5, 6, 
                     3, 2, 1, 10, 12, 11, 7, 9, 8, 4, 6, 5), 
                   byrow = T, ncol = 12)

cmat_u <- matrix(c(1, 2, 3, 
                   1, 3, 2, 
                   2, 1, 3, 
                   2, 3, 1, 
                   3, 1, 2, 
                   3, 2, 1), byrow = T, ncol = 3)

# Function to rearrange pivot probabilities and utilities according to voter type 
rearrange = function(x, type){
  out = list()
  for(v in 1:6){
    # need to rearrange utils as well
    subdf = x %>% filter(voter == v)
    subdf[, c("uA", "uB", "uC")] = subdf[, c("uA", "uB", "uC")][, cmat_u[v, ]]
    if(type == "IRV"){
      subdf[, 2:13] = subdf[, 2:13][, cmat_rcv[v, ]]
      names(subdf)[2:13] = c("ABo", "ACo", "BCo", "AB_ABo", "AB_ACo", "AB_CBo", "AC_ACo", "AC_ABo", "AC_BCo", "BC_BCo", "BC_BAo", "BC_ACo")
    }
    if(type == "Plurality"){
      subdf[, 2:4] = subdf[, 2:4][, cmat_plur[v, ]]
      names(subdf)[2:4] = c("ABpo", "ACpo", "BCpo")
    }
    out[[v]] = subdf
  }
  return(do.call(rbind, out))
}

# RCV DATA ---------------------------------------
rcvdf85 = rearrange(big_rcv_sum, type = "IRV")

# Add cost and benefit analyses
rcvdf852 <- rcvdf85 %>% mutate(
  ben_rcv2 = AB_CBo * (uB - uC) + BC_BCo * ((uB - uC)/2) + BC_ACo * ((uA - uC)/2),
  cost_rcv2 = ABo * (uA - uB) + AB_ABo * (uA - uB) + AB_ACo * (uA - uC) + AC_ACo * ((uA - uC) / 2) + AC_ABo * ((uA - uB) / 2) + AC_BCo * ((uB - uC)/2) + BC_BAo * ((uA - uB)/2),
  ben_rcv3 = AB_CBo * ((uB - uC) / 2) + BC_BAo * ((uA - uB) / 2),
  cost_rcv3 = ACo * (uA - uC) + BCo * (uB - uC) + AB_ABo * ((uA - uB) / 2) + AB_ACo * (uA - uC)/2 + AC_ACo * (uA - uC)/2 + AC_ABo * (uA - uB) + BC_BCo * (uC - uB)/2 + BC_ACo * (uA - uC)/2,
  pprob2 = AB_CBo + BC_BCo + BC_ACo,
  pprob3 = AB_CBo + BC_BAo)

# Compute correlations
rcvdf853 <- rcvdf852 %>% group_by(case, iter) %>% 
  summarise(corr_rcv2 = cor(ben_rcv2, cost_rcv2), 
            corr_rcv3 = cor(ben_rcv3, cost_rcv3),
            pprob2 = weighted.mean(pprob2, weight = weights), 
            pprob3 = weighted.mean(pprob3, weight = weights), 
            path = "IRV")

# Prepare DF for plotting
rcvdf854 = rcvdf853 %>%
  pivot_longer(corr_rcv2:pprob3, names_to = "quantity") %>%
  mutate(type = ifelse(quantity %in% c("corr_rcv2", "pprob2"), "IRV_second", "IRV_third"), 
         quantity = ifelse(quantity %in% c("corr_rcv2", "corr_rcv3"), "correlation", "pprob")) %>%
  pivot_wider(names_from = "quantity", values_from = "value")

# PLURALITY DATA -------------------------------------------
plurdf85 = rearrange(big_plur_sum, type = "Plurality") %>% 
  mutate(ben_p = BCpo * ((uB - uC) / 2), 
         cost_p = ABpo * (uA - uB) + ACpo * ((uB - uC)/2), 
         pprob_plur = BCpo)

# Compute correlation etc.
plurdf853 = plurdf85 %>%
  group_by(case, iter) %>% 
  summarise(correlation = cor(ben_p, cost_p), 
            pprob = weighted.mean(pprob_plur, weight = weights), 
            type = "Plurality", 
            path = "Plurality")

# get country weights
sv_list <- list()
n <- length(big_list_na_omit)
country_weight <- matrix(nrow = n, ncol = 2)
for(i in 1:n){
  print(i)
  if(i == n) cat("Done! \n")
  this_list <- big_list_na_omit[[i]]
  dfcase <- names(big_list_na_omit)[[i]]
  dfcountry <- substr(dfcase, 1, 3)
  dfweight_sum <- sum(big_list_na_omit[[i]]$weights)
  dfVAP <- vap$VAP[vap$cntry == dfcountry]
  dfm <- vap$Freq[vap$cntry == dfcountry]
  country_weight[i, 1] <- names(big_list_na_omit)[[i]]
  country_weight[i, 2] <- dfVAP / dfm
}
weightdf <- as.data.frame(country_weight)
names(weightdf) <- c("case", "ctryweight")
weightdf$ctryweight <- as.numeric(as.character(weightdf$ctryweight))

# Bind RCV and Plurality together
conjdf <- rbind(plurdf853, rcvdf854) %>% left_join(weightdf) %>%
  mutate(iter = as.numeric(iter))
conjdf_quant <- conjdf %>% 
  group_by(iter, type) %>% 
  summarise(
    corr_q025 = wtd.quantile(correlation, probs = 0.025, weight = ctryweight),
    corr_q25 = wtd.quantile(correlation, probs = 0.25, weight = ctryweight),
    corr_q50 = wtd.quantile(correlation, probs = 0.5, weight = ctryweight),
    corr_q75 = wtd.quantile(correlation, probs = 0.75, weight = ctryweight),
    corr_q975 = wtd.quantile(correlation, probs = 0.975, weight = ctryweight),
    corr_mean = wtd.mean(correlation, weight = ctryweight),
    pprob_q025 = wtd.quantile(pprob, probs = 0.025, weight = ctryweight),
    pprob_q25 = wtd.quantile(pprob, probs = 0.25, weight = ctryweight),
    pprob_q50 = wtd.quantile(pprob, probs = 0.5, weight = ctryweight),
    pprob_q75 = wtd.quantile(pprob, probs = 0.75, weight = ctryweight),
    pprob_q975 = wtd.quantile(pprob, probs = 0.975, weight = ctryweight),
    pprob_mean = wtd.mean(pprob, weight = ctryweight)
  )
  
# Path for saving figures 
save(conjdf, file = fpath(lambda, s, "conj"))
  
#   # Pprobs plot (raw)
ggplot(conjdf, aes(x = as.numeric(iter))) +
geom_line(aes(y = pprob, group = interaction(case, type), colour = type), alpha = 0.05) +
geom_line(data = conjdf_quant, aes(x = iter, y = pprob_mean, group = type, colour = type), lwd = 2) +
geom_hline(yintercept = 0, lty = "dashed") +
labs(x = "Iteration", y = "Probability vote is beneficial * electorate size") +
theme_tn() +
coord_cartesian(ylim = c(0, 0.5))
ggsave(ppath(lambda, s, "conj1"),
       width = 4,
       height = 4)

#   # Correlations plot
ggplot(conjdf, aes(x = as.numeric(iter))) +
geom_line(aes(y = correlation, group = interaction(case, type), colour = type), alpha = 0.05) +
geom_line(data = conjdf_quant, aes(x = iter, y = corr_mean, group = type, colour = type), lwd = 2) +
geom_hline(yintercept = 0, lty = "dashed") +
theme_tn() +
labs(x = "Iteration", y = "Correlation") +
theme(legend.position = "bottom",
      legend.direction = "horizontal")
ggsave(ppath(lambda, s, "conj2"),
       width = 4,
       height= 4)


