# calculate average distance between top and bottom preference(s)

# Load packages and data
library(tidyverse)
library(gtools)
library(fields)
library(ggsci)
library(ggtern)
library(devtools)
source_url('https://raw.githubusercontent.com/tobiasnowacki/RTemplates/master/plottheme.R')

# Load functions
source("code/utils/new_sv_iter.R")
# Load data
load("output/big_list_2.RData")
vap <- read.csv("data/case_vap.csv", sep = "") # voting age pop.
cat("Data imported. \n")
source("code/prep_cses.R")  # data prep
systems <- read.csv("data/systems.csv")

return_dat <- function(X){
  diff <- mean(apply(X$U, 1, function(y) max(y) - min(y)))
  top <- mean(apply(X$U, 1, function(y) max(y)))
  bottom <- mean(apply(X$U, 1, function(y) min(y)))
  return(data.frame(Difference = diff, Top = top, Bottom = bottom))
}

easy_summary <- map_dfr(big_list_na_omit, ~ return_dat(.x), .id = "case")
easy_summary_mg <- easy_summary |>
  pivot_longer(!case, names_to = "stat", values_to = "value") |>
  left_join(systems, by = c("case")) |>
  filter(system %in% c("pr", "fptp", "rcv"))

easy_summary_mg$stat <- fct_relevel(easy_summary_mg$stat, c("Top", "Bottom", "Difference"))

# RCV is a little bit of an outlier - but this could just because of the much smaller sample...
ggplot(easy_summary_mg, aes(x = value)) +
  geom_density(aes(fill = system), alpha = 0.5) +
  facet_grid(~ stat) +
  scale_fill_nejm() +
  theme_tn()

ggsave("output/figures/profile_comp.pdf", height = 3, width = 8, device = cairo_pdf)
