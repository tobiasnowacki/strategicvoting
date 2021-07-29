library(tidyverse)
library(rio)
library(questionr)
library(RColorBrewer)
options(tibble.width = Inf)

sum_df <- import("output/files/1/85summary_1_85.Rdata")
case_df <- import("data/systems.csv")
vap <- read.csv("data/case_vap.csv", sep = "")
source("code/utils/sv_theme_template.R")


case_weight_tbl <- tibble(case = unique(sum_df$case),
	cntry = substr(case, 1, 3)) %>%
	inner_join(vap) %>% 
	mutate(case_weight = VAP / Freq)

vap
head(sum_df)

sum_mg <- sum_df %>%
    left_join(case_df) %>%
    left_join(case_weight_tbl) %>%
    mutate(system_group = case_when(
        system %in% c("pr", "mixed") ~ "pr",
        system %in% c("fptp", "tworound", "parallel") ~ "plur",
        system %in% "rcv" ~ "rcv", 
        TRUE ~ "other"
    )) %>%
    pivot_longer(Prevalence:ExpBenefit)

# Aggregate by system
sum_agg <- sum_mg %>%
    group_by(iter, System, system_group, name) %>%
    summarise(
       value = wtd.mean(value, case_weight)
    )

# Plot (by system)
ggplot(sum_agg, aes(x = iter, y = value)) +
    geom_line(aes(
        group = interaction(System, system_group), 
        colour = system_group,
        lty = System
        )
    ) +
    facet_wrap(. ~ name, scales = "free") +
    scale_colour_brewer(palette = "Dark2", labels = c("Plurality", "PR", "Ranked")) +
    labs(x = "Iteration", y = "Value", lty = "Simulated System", colour = "Case System") +
    theme_tn()
