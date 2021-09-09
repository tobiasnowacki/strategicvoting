# PREPARE EVERYTHING -------------------------------------------
# Load dependencies
library(tidyverse)
library(pivotprobs)
library(gtools)
library(doParallel)
source("code/utils/new_sv_iter.R")
source("code/utils/sv_theme_template.R")
source("code/utils/new_dist_helpers.R")
source("code/utils/new_plot_helpers.R")
source("code/utils/new_summary_helpers.R")

# Load data
load("output/big_list_2.RData")
vap <- read.csv("data/case_vap.csv", sep = "") # voting age pop.
cat("Data imported. \n")
source("code/prep_cses.R")  # data prep

# Grab case names
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

# Helping functions for saving files / figures
fpath = function(lambda, s, ext){
  paste0("output/files/", lambda, "/", s, "_", ext, ".Rdata")
}
ppath = function(lambda, s, ext){
  paste0("output/figures/", lambda, "/", s, "_", ext, ".pdf")
}

# LOAD SIMULATION DATA -----------------------------------------
# Function to load all cases for param values
return_obj = function(case, lambda, s){
  load(paste0("output/files/", lambda, "/", s, "_", case, "_iterout.Rdata"))
  return(inner_list)
}

# Grab all 160 simulations
bind_together = map(nn, ~ return_obj(.x, lambda, s))
names(bind_together) = nn

# SUMMARY STATISTICS -------------------------------------------
sum_df = map_dfr(bind_together, ~ get_sum_stats(.x), .id = "case")
save(sum_df, file = fpath(lambda, s, "sum")) # save data

# Prep summary df and plot statistics
summary_stats_wide <- sum_df %>% 
  rename(ExpBenefit = eb,
         Magnitude = mag,
         Prevalence = prev) %>%
  pivot_longer(Prevalence:ExpBenefit) %>%
  mutate(iter = as.numeric(iter)) %>%
  mutate(system = case_when(
    system == "plur" ~ "Plurality",
    system == "rcv" ~ "RCV"
  ))

<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
# weight by case weight
summary_agg <- summary_stats_wide %>% 
  group_by(iter, name, system) %>%
  summarise(value = weighted.mean(value, case_weight_tbl$case_weight)) %>%
  mutate(iter = as.numeric(iter))

=======
=======
>>>>>>> 95425d540f4efbe579af50da3c43f21f50cc9383
=======
>>>>>>> 95425d540f4efbe579af50da3c43f21f50cc9383
summary_agg <- summary_stats_wide %>%
  left_join(case_weight_tbl, by = c("case" = "case")) %>%
  group_by(iter, system, name) %>%
  summarise(value = weighted.mean(value, case_weight)) %>%
  mutate(iter = as.numeric(iter))

# weight by case weight
# summary_agg <- summary_stats_wide %>% 
#   group_by(iter, name, system) %>%
#   summarise(value = weighted.mean(value, case_weight_tbl$case_weight)) %>%
#   mutate(iter = as.numeric(iter))

<<<<<<< HEAD
<<<<<<< HEAD
>>>>>>> 95425d540f4efbe579af50da3c43f21f50cc9383
=======
>>>>>>> 95425d540f4efbe579af50da3c43f21f50cc9383
=======
>>>>>>> 95425d540f4efbe579af50da3c43f21f50cc9383
# Plot summary statistics for parameter combination and save
ggplot(summary_stats_wide, aes(iter, value)) +
  geom_line(aes(group = interaction(system, case, name),
                colour = system), alpha = 0.05) +
  geom_line(data = summary_agg %>% filter(system == "Plurality"),
            aes(group = interaction(system, name)), alpha = 1,
            color = "#CC6600", lwd = 1.1) +
  geom_line(data = summary_agg %>% filter(system == "RCV"),
            aes(group = interaction(system, name)), alpha = 1,
            color = "#004C99", lwd = 1.1) +
  facet_wrap(. ~ name, scales = "free_y") +
  theme_tn() +
  scale_color_manual(values = c("#CC6600", "#004C99")) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  labs(x = "Degree of Strategicness (Iterations)")
<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
ggsave(ppath(lambda, s, "summary"), width = 8, height = 3.5)

=======
=======
>>>>>>> 95425d540f4efbe579af50da3c43f21f50cc9383
=======
>>>>>>> 95425d540f4efbe579af50da3c43f21f50cc9383

ggsave(ppath(lambda, s, "summary"), width = 8, height = 3.5)

# WINNER DATA ------------------------------------------------

# Function to prepare winner data
get_expected_utility <- function(obj){
  rcv_df <- map_dbl(obj$rcv[1:(length(obj$rcv) - 1)], 
    ~ .x$exp_win_mean
   ) %>% as.data.frame()
  names(rcv_df)[1] <- "value"

  rcv_df$iter <- 1:(length(obj$rcv) - 1)
  rcv_df$system <- "rcv"

  plur_df <- map_dbl(obj$plur[1:(length(obj$rcv) - 1)], 
    ~ .x$exp_win_mean
  ) %>% as.data.frame()
  names(plur_df)[1] <- "value"

  plur_df$iter <- 1:(length(obj$plur) - 1)
  plur_df$system <- "plurality"

  return(rbind(rcv_df, plur_df))
}

win_df <- map(bind_together, ~ get_expected_utility(.x), .id = "case")
names(win_df) <- nn
win_df <- bind_rows(win_df, .id = "case")
save(win_df, file = fpath(lambda, s, "winners"))

# Do the same but with winners
get_winners <- function(obj){
  rcv_df <- map(obj$rcv[1:(length(obj$rcv) - 1)], 
    ~ .x$wintbl
   ) 

  rcv_df <- do.call(rbind, rcv_df) %>% as.data.frame()
  rcv_df$iter <- 1:(length(obj$rcv) - 1)
  rcv_df$system <- "rcv"

  plur_df <- map(obj$plur[1:(length(obj$rcv) - 1)], 
    ~ .x$wintbl
  ) 

  plur_df <- do.call(rbind, plur_df) %>% as.data.frame()
  plur_df$iter <- 1:(length(obj$plur) - 1)
  plur_df$system <- "plurality"

  return(rbind(rcv_df, plur_df))
}

winners_df <- map(bind_together, ~ get_winners(.x), .id = "case")
names(winners_df) <- nn
winners_df <- bind_rows(winners_df, .id = "case")
save(winners_df, file = fpath(lambda, s, "winners_tbl"))

<<<<<<< HEAD
<<<<<<< HEAD
>>>>>>> 95425d540f4efbe579af50da3c43f21f50cc9383
=======
>>>>>>> 95425d540f4efbe579af50da3c43f21f50cc9383
=======
>>>>>>> 95425d540f4efbe579af50da3c43f21f50cc9383
# VVEC DATA --------------------------------------------------
# Function to prepare v_vec data
get_vvecs = function(obj){
  w = obj$rcv[[i]]$weights
  outobj = list()
  outobj$rcvvec = map_dfr(obj$rcv, ~ .x$v.vec.before)
  outobj$plurvec = map_dfr(obj$plur, ~ .x$v.vec.before)
  outobj$rcvbr = map_dfr(obj$rcv, ~ .x$best.response.v.vec)
  outobj$plurbr = map_dfr(obj$plur, ~ .x$best.response.v.vec)
  return(outobj)
}

# Run and save vvec data
vvecdf = map(bind_together, ~ get_vvecs(.x))
names(vvecdf) = nn
save(vvecdf, file = fpath(lambda, s, "vvec"))

# Create joint vvec plot (ternary)
joint_v_vec_plot(vvecdf, paste0("output/figures/", lambda, "/", s))
# Unload ggtern because it screws with normal ggplots
devtools::unload("ggtern")
R.methodsS3::setMethodS3("print", "ggplot", ggplot2:::print.ggplot)
R.methodsS3::setMethodS3("plot", "ggplot", ggplot2:::plot.ggplot)
R.methodsS3::setMethodS3("grid.draw", "ggplot", ggplot2:::grid.draw.ggplot)
# Create distance plots
plot_v_vec_distance(vvecdf, paste0("output/figures/", lambda, "/", s), 
                    n_lag = 20, avg_span = 10) # values for paper

cat("new_gather_data.R done! \n")

