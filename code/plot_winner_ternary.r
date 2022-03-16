# Script to plot the proportion of simulations in each iteration that
# return the same winner as sincere voting

# Load dependencies
library(tidyverse)
library(rio)
library(gtools)
library(devtools)


# Load data
load("output/files/1/85_winners_tbl.Rdata")
load("output/files/1/85_winners.Rdata")

source("code/utils/sv_theme_template.R")
source("code/utils/new_sv_iter.R")
load("output/big_list_2.RData")
vap <- read.csv("data/case_vap.csv", sep = "") # voting age pop.
cat("Data imported. \n")
source("code/prep_cses.R")  # data prep
lambda <- 1
s <- 85

# Function to save plot in folder
ppath = function(lambda, s, ext){
  paste0("output/figures/", lambda, "/", s, "_", ext, ".pdf")
}

# Plot win probabilities in ternary...
# library(ggtern)
# winners_df_clean <- winners_df %>% 
#     filter(iter == 1 | iter == 60) %>%
#     rename(A = V1, B = V2, C = V3)

# p1 <- ggtern(winners_df_clean, aes(B, A, C)) +
#     # geom_line(aes(group = case), alpha = 0.2) +
#     geom_point(aes(colour = system), alpha = 0.1) +
#     facet_wrap(. ~ iter) +
#     theme_tn()

# ggsave(p1,
#   filename = ppath(lambda, s, "win_ternary"),
#   width = 10, height = 10
# )

# Merge together
ddf <- winners_df %>%
    left_join(win_df)

diff_df <- ddf %>%
    group_by(case, iter) %>%
    summarise(
        diff_V1 = V1[system == "plurality"] - V1[system == "rcv"], 
        diff_V2 = V2[system == "plurality"] - V2[system == "rcv"],
        diff_V3 = V3[system == "plurality"] - V3[system == "rcv"],
        diff_eu = value[system == "plurality"] - value[system == "rcv"]) %>%
    mutate(vec_diff = sqrt(diff_V1^2 + diff_V2^2 + diff_V3^2))

diff_df_pruned <- diff_df %>%
    filter(iter == 1 | iter == 60)

p2 <- ggplot(diff_df_pruned, aes(diff_V1, diff_eu)) +
    geom_line(aes(group = case), alpha = 0.1) +
    geom_point(aes(colour = as.factor(iter)), alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1) +
    theme_tn()

ggsave(p2,
  filename = ppath(lambda, s, "win_ternary_scatter"),
  width = 5, height = 5
)

p3 <- ggplot(diff_df_pruned, aes(vec_diff, diff_eu)) +
    geom_line(aes(group = case), alpha = 0.1) +
    geom_point(aes(colour = as.factor(iter)), alpha = 0.7) +
    geom_abline(intercept = 0, slope = 1) +
    theme_tn()

ggsave(p3,
  filename = ppath(lambda, s, "win_ternary_scatter_mse"),
  width = 5, height = 5
)

diff_df_pruned %>% filter(diff_eu < -0.4)

ddf %>% filter(case == "CHE_1999", iter %in% c(1, 60))

ddf %>% filter(case == "ROU_2004", iter %in% c(1, 60))
