# Creates plot to compare baseline against different electoral systems

library(tidyverse)
library(rio)
library(questionr)
library(RColorBrewer)
options(tibble.width = Inf)

sum_df <- import("output/files/1/85_sum.Rdata")
case_df <- import("data/systems.csv")
vap <- read.csv("data/case_vap.csv", sep = "")
source("code/utils/sv_theme_template.R")

ppath = function(lambda, s, ext){
  paste0("output/figures/", lambda, "/", s, "_", ext, ".pdf")
}

# Create table with case weights
case_weight_tbl <- tibble(case = unique(sum_df$case),
	cntry = substr(case, 1, 3)) %>%
	inner_join(vap) %>% 
	mutate(case_weight = VAP / Freq)

# Merge case weights in
sum_mg <- sum_df %>%
    rename(
        System = 'system',
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

# Aggregate by system
sum_agg_system <- sum_mg %>%
    group_by(iter, System, system_group, name) %>%
    summarise(
       value = wtd.mean(value, case_weight)
    ) %>%
    mutate(iter = as.numeric(iter)) %>%
    filter(iter < 60, !is.na(system_group), system_group != "other")

# Aggregate overall (baseline)
sum_agg <- sum_mg %>%
    group_by(iter, System, name) %>%
    summarise(
        value = wtd.mean(value, case_weight)
    ) %>%
    mutate(system_group = "baseline",
    iter = as.numeric(iter)) %>%
    filter(iter < 60)

# Double check the RCV cases really have quite low weight
# only Australia and Ireland
case_weight_tbl %>%
    mutate(pct = case_weight / sum(case_weight)) %>%
    left_join(case_df) %>%
    filter(system == "rcv" | system == "stv")

# Plot (by system)
out <- ggplot(sum_agg_system, aes(x = iter, y = value)) +
    geom_line(
        data = sum_agg, aes(
            group = System, 
            lty = System
        ), 
        colour = "#999999", 
        alpha = 0.6
    ) + 
    geom_line(aes(
        group = interaction(System, system_group), 
        colour = system_group,
        lty = System
        )
    ) +
    facet_wrap(. ~ name, scales = "free") +
    scale_linetype_manual(values = c("solid", "dashed"), labels = c("Plurality", "RCV")) +
    scale_colour_brewer(palette = "Dark2", labels = c("Plurality", "PR", "Ranked")) +
    labs(x = "Iteration", y = "Value", lty = "Simulated System", colour = "Case System") +
    theme_tn()

lambda <- 1; s <- 85
ggsave(out, filename = ppath(lambda, s, "system_baseline"), width = 8, height = 3.5)
