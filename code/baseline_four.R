library(tidyverse)
library(rio)
library(questionr)
library(RColorBrewer)
options(tibble.width = Inf)

sum_df <- import("output/files/1/85_sum.Rdata")
sum_four_df <- import("output/files/1/85_sum_four.Rdata")
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
    left_join(case_weight_tbl) %>%
    pivot_longer(Prevalence:ExpBenefit) %>%
    mutate(setting = "3 Parties") 

sum_four_mg <- sum_four_df %>%
    select(-case) %>%
    rename(
        System = 'system',
        Prevalence = "prev",
        Magnitude = "mag",
        ExpBenefit = "eb",
        case = "casename"
    ) %>%
    left_join(case_weight_tbl) %>%
    pivot_longer(Prevalence:ExpBenefit) %>%
    mutate(setting = "4 Parties") 

all_df <- rbind(sum_mg, sum_four_mg)

# Aggregate up
sum_all <- all_df %>%
    group_by(iter, System, setting, name) %>%
    summarise(
        value = wtd.mean(value, case_weight)
    ) %>% 
    mutate(iter = as.numeric(iter)) %>%
    filter(iter < 61)

# Plot (by system)
out <- ggplot(sum_all, aes(x = iter, y = value)) +
    geom_line(aes(
        group = interaction(System, setting), 
        colour = setting,
        alpha = setting,
        lty = System
        )
    ) +
    facet_wrap(. ~ name, scales = "free") +
    scale_linetype_manual(values = c("solid", "dashed"), 
        labels = c("Plurality", "RCV")) +
    scale_colour_manual(values = c("#999999", "#0000FF")) +
    scale_alpha_manual(values = c(0.6, 1)) +
    labs(x = "Iteration", y = "Value", lty = "Simulated System", colour = "No of Parties", alpha = "No of Parties") +
    theme_tn()

lambda <- 1; s <- 85
ggsave(out, filename = ppath(lambda, s, "parties_baseline"), width = 8, height = 3.5)

