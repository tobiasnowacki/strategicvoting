# Script to produce figure
# September 2021
# Creates plot to compare baseline
# against different electoral systems

# Load dependencies
## -------------------------------------------

library(tidyverse)
library(rio)
library(questionr)
library(RColorBrewer)
options(tibble.width = Inf)

lambda <- 1

case_df <- import("data/systems.csv")
vap <- read.csv("data/case_vap.csv", sep = "")
source("code/utils/sv_theme_template.R")

ppath = function(lambda, s, ext){
  paste0("output/figures/", lambda, "/", s, "_", ext, ".pdf")
}

## Load data
## -------------------------------------------

all_sum <- map_dfr(
    c(10, 55, 85),
    ~ import(paste("output/files/1/", .x, "_sum.Rdata", sep = "")) %>%
        mutate(s = .x)
)

# Construct aggregates
# Create table with case weights
case_weight_tbl <- tibble(
    case = unique(all_sum$case),
    cntry = substr(case, 1, 3)
) %>%
    inner_join(vap) %>%
    mutate(case_weight = VAP / Freq)

# Merge case weights in
sum_mg <- all_sum %>%
    rename(
        System = "system",
        Prevalence = "prev",
        Magnitude = "mag",
        ExpBenefit = "eb"
    ) %>%
    left_join(case_df) %>%
    left_join(case_weight_tbl) %>%
    mutate(system_group = case_when(
        system %in% c("pr", "mixed") ~ "pr",
        system %in% c("fptp", "tworound", "parallel") ~ "plur",
        system %in% c("rcv") ~ "rcv",
        TRUE ~ "other"
    )) %>%
    pivot_longer(Prevalence:ExpBenefit)

# Aggregate overall (baseline)
sum_agg <- sum_mg %>%
    group_by(iter, System, s, name) %>%
    summarise(
        value = wtd.mean(value, case_weight)
    ) %>%
    mutate(
        system_group = "baseline",
        iter = as.numeric(iter), 
        s = as.factor(s)
    ) %>%
    filter(iter < 60)

# Plot (by s)
out <- ggplot(sum_agg, aes(x = iter, y = value)) +
    geom_line(aes(
        group = interaction(System, s), 
        colour = s,
        lty = System
        )
    ) +
    facet_wrap(. ~ name, scales = "free") +
    scale_linetype_manual(values = c("solid", "dashed"), labels = c("Plurality", "RCV")) +
    scale_colour_brewer(palette = "Dark2") +
    labs(x = "Iteration", y = "Value", lty = "System", colour = "Precision") +
    theme_tn()

ggsave(out, filename = ppath(lambda, 85, "precision_baseline"), width = 8, height = 3.5)
