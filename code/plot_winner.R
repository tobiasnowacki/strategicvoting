# Script to plot the proportion of simulations in each iteration that return the same winner as sincere voting

# Load dependencies
library(tidyverse)
library(rio)
library(gtools)

# Load data
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

# Compute aggregates
diff_df <- win_df %>%
  group_by(case, iter) %>%
  summarise(diff = value[system == "rcv"] - value[system == "plurality"])

diff_agg <- diff_df %>%
  left_join(case_weight_tbl, by = c("case" = "case")) %>%
  group_by(iter) %>%
  summarise(diff = weighted.mean(diff, case_weight)) %>%
  mutate(iter = as.numeric(iter))

p <- ggplot(diff_df, aes(iter, diff)) +
  geom_line(aes(group = case), alpha = 0.05) +
  geom_line(data = diff_agg,
            alpha = 1,
            color = "blue", lwd = 1.1) +
  scale_color_manual(values = c("#CC6600", "#004C99")) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  scale_x_continuous(expand = c(0, 0)) +
  theme(
    legend.position = "bottom", 
    legend.direction = "horizontal") +
  labs(
    x = "Degree of Strategicness (Iterations)", 
    y = "Delta EU RCV - Plur") +
  theme_tn()

ggsave(p, 
  filename = ppath(lambda, s, "win_util"), 
  width = 8, height = 3.5
)

# Plot by simulated system
win_df_agg <- win_df %>%
  left_join(case_weight_tbl, by = c("case" = "case")) %>%
  group_by(iter, system) %>%
  summarise(value = weighted.mean(value, case_weight)) %>%
  mutate(iter = as.numeric(iter))

p2 <- ggplot(win_df, aes(iter, value)) +
  geom_line(aes(group = case), alpha = 0.05) +
  geom_line(data = win_df_agg,
            alpha = 1, 
            lwd = 1.1) +
  scale_color_manual(values = c("#CC6600", "#004C99")) +
  facet_wrap(~ system) +
  theme_tn()

ggsave(p2,
  filename = ppath(lambda, s, "win_util_by_system"),
  width = 8, height = 3.5
)

# Scatterplots
win_scatter <- win_df %>%
  filter(iter %in% c(1, 60)) %>%
  pivot_wider(
    values_from = value, names_from = iter, names_prefix = "iter"
  ) %>%
  mutate(decrease = iter60 < iter1)

p3 <- ggplot(win_scatter, aes(iter1, iter60)) +
  geom_point(aes(color = system, shape = decrease), alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0) +
  facet_wrap(. ~ system) +
  scale_shape_manual(values = c(16, 21)) +
  labs(x = "Avg EU at 1st Iteration", y = "Avg EU at 60th Iteration") +
  theme_tn()

ggsave(p3,
  filename = ppath(lambda, s, "win_util_scatter"),
  width = 8, height = 3.5
)

win_scatter2 <- win_df %>%
  filter(iter %in% c(1, 60)) %>%
  pivot_wider(
    values_from = value, names_from = system
  ) %>%
  mutate(decrease = rcv < plurality)

p4 <- ggplot(win_scatter2, aes(plurality, rcv)) +
  geom_point(aes(color = as.factor(iter), shape = decrease), alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0) +
  facet_wrap(. ~ iter) +
  scale_shape_manual(values = c(16, 21)) +
  labs(x = "Avg EU in Plurality", y = "Avg EU in RCV") +
  theme_tn()

ggsave(p4,
  filename = ppath(lambda, s, "win_util_scatter_system"),
  width = 8, height = 3.5
)
