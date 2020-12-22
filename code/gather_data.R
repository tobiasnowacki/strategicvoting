# LOAD DEPENDENCIES
library(tidyverse)
library(pivotprobs)
library(gtools)
library(doParallel)
source("code/utils/new_sv_iter.R")
source("code/utils/sv_theme_template.R")

# Load data
load("output/big_list_2.RData")
vap <- read.csv("data/case_vap.csv", sep = "") # voting age pop.
cat("Data imported. \n")
source("code/prep_cses.R")  # data prep

source("code/utils/new_dist_helpers.R")
source("code/utils/new_plot_helpers.R")

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

bind_together = map(nn, ~ return_obj(.x, lambda, s))
names(bind_together) = nn

# Function to get summary statistics
get_sum_stats = function(obj){
  w = obj$rcv[[i]]$weights
  map_dfr(c("rcv", "plur"), function(y){
    names(obj[[y]]) = 1:length(obj[[y]])
    map(obj[[y]], ~
          tibble(
            prev = weighted.mean(.x$tau > 0, w),
            mag  = weighted.mean(.x$tau[.x$tau > 0], w[.x$tau > 0]),
            eb = prev * mag)) %>%
      bind_rows(.id = "iter") %>%
      mutate(system = y)
  })
}

sum_df = map_dfr(bind_together, ~ get_sum_stats(.x), .id = "case")
save(sum_df, file = fpath(lambda, s, "sum"))

# Get summary df and plot statistics
summary_stats_wide <- sum_df %>% 
  gather(., key = "stat", 
         value = "value", "prev":"eb") %>%
  mutate(iter = as.numeric(iter)) %>%
  rename(ExpBenefit = eb,
         Magnitude = mag,
         Prevalence = prev)

# weight by case weight
summary_agg <- summary_stats_wide %>% group_by(iter, stat, system) %>%
  summarise(value = weighted.mean(value, case_weight_tbl$case_weight))

# Plot summary statistics for parameter combination
ggplot(summary_stats_wide, aes(iter, value)) +
  geom_line(aes(group = interaction(system, case, stat),
                colour = system), alpha = 0.3) +
  geom_line(data = summary_agg %>% filter(system == "plur"),
            aes(group = interaction(system, stat)), alpha = 1,
            color = "#CC6600", lwd = 1.1) +
  geom_line(data = summary_agg %>% filter(system == "irv"),
            aes(group = interaction(system, stat)), alpha = 1,
            color = "#004C99", lwd = 1.1) +
  facet_wrap(. ~ stat, scales = "free_y") +
  theme_tn() +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  labs(x = "Degree of Strategicness (Iterations)")
ggsave(ppath(lambda, s, "summary"))

# Function to save v_vec data
get_vvecs = function(obj){
  w = obj$rcv[[i]]$weights
  outobj = list()
  outobj$rcvvec = map_dfr(obj$rcv, ~ .x$v.vec.before)
  outobj$plurvec = map_dfr(obj$plur, ~ .x$v.vec.before)
  outobj$rcvbr = map_dfr(obj$rcv, ~ .x$best.response.v.vec)
  outobj$plurbr = map_dfr(obj$plur, ~ .x$best.response.v.vec)
  return(outobj)
}

vvecdf = map(bind_together, ~ get_vvecs(.x))
names(vvecdf) = nn
# Save vvec data (for later distance plots)
save(vvecdf, file = fpath(lambda, s, "vvec"))

# Create distance plots
joint_v_vec_plot(vvecdf, paste0("output/figures/", lambda, "/", s))
# Create vvec plots
plot_v_vec_distance(vvecdf, paste0("output/figures/", lambda, "/", s), 
                    n_lag = 20, avg_span = 10)


