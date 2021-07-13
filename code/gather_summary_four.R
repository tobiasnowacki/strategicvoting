# Script to gather four party results
library(tidyverse)
library(gtools)

# Load data
load("output/big_list_4_parties.RData")
vap <- read.csv("data/case_vap.csv", sep = "") # voting age pop.
source("code/utils/sv_iter_four.R")
cat("Data imported. \n")
source("code/prep_cses.R")  # data prep

# Grab case names
nn = names(big_list_na_omit)

# Get command arguments / values
# cmd_line_args <- commandArgs(trailingOnly = TRUE)
# cat(cmd_line_args, sep = "n")
# ifelse(length(cmd_line_args >= 1),
#    s <- as.numeric(cmd_line_args[1]),
#    s <- 85)
# ifelse(length(cmd_line_args >= 2),
#    lambda<- as.numeric(cmd_line_args[2]),
#    lambda <- 1)
s <- 85
lambda <- 1

# Helping functions for saving files / figures
fpath = function(lambda, s, ext){
  paste0("output/files/", lambda, "/", s, "_", ext, ".Rdata")
}
ppath = function(lambda, s, ext){
  paste0("output/figures/", lambda, "/", s, "_", ext, ".pdf")
}

# LOAD SIMULATION DATA -----------------------------------------
# Function to load all cases for param values
return_obj = function(fname){
  load(paste0("output/files/", 1, "/", fname))
  return(inner_list)
}

# Grab all 160 simulations

filelist <- list.files(path = "output/files/1/")
filelist <- filelist[grepl("fourparties", filelist)]
caselist <- substr(filelist, 4, 11)
bind_together = map(filelist, ~ return_obj(.x))
# names(bind_together) = nn

# SUMMARY STATISTICS -------------------------------------------
sum_df = map_dfr(bind_together, ~ get_sum_stats(.x), .id = "case")
save(sum_df, file = fpath(lambda, s, "sum_four")) # save data

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

# weight by case weight
case_weight_tbl <- case_weight_tbl[case_weight_tbl$case %in% caselist, ] %>%
    arrange(case)

summary_agg <- summary_stats_wide %>% 
  group_by(iter, name, system) %>%
  summarise(value = weighted.mean(value, case_weight_tbl$case_weight)) %>%
  mutate(iter = as.numeric(iter))

# Plot summary statistics for parameter combination and save
t <- ggplot(summary_stats_wide, aes(iter, value)) +
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

ggsave(t, filename = ppath(lambda, s, "summary_four"), width = 8, height = 3.5)
