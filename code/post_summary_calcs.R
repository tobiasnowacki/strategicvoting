### LOAD DEPENDENCIES -------

# install.packages("here")
# source("code/full_header.R")  # fn's and data
library(tidyverse)
library(gtools)
source("code/utils/new_sv_iter.R")
source("code/utils/sv_theme_template.R")
load("output/big_list_2.RData")
vap <- read.csv("data/case_vap.csv", sep = "") # voting age pop.
cat("Data imported. \n")
source("code/prep_cses.R")  # data prep

# 
s_param <- c("10", "55", "85")
lambda_param <- c("1")

summary_list <- list()
for(i in s_param){
  for(j in lambda_param){
    # load file
    datapath <- paste0("output/files/", j, "/", i, "_sum", ".Rdata")
    load(datapath)
    # append cases w/ relevant data
    summary_list[[paste0(j, "_", i)]] <- sum_df %>% 
      mutate(s = i, 
             lambda = j)
  }
} 

hidden_scale_adj <- data.frame(iter = 3, value = 0.69, name = "Prevalence", s = "s = 10")

summary_df <- do.call(rbind, summary_list) %>% 
  mutate(cntry = substr(case, 1, 3)) %>% 
  inner_join(case_weight_tbl) %>% 
  pivot_longer(prev:eb) %>% 
  mutate(s = recode(s, `10` = "s = 10", `55` = "s = 55", `85` = "s = 85"),
         name = recode(name, prev = "Prevalence", mag = "Magnitude", eb = "ExpBenefit"),
         system = recode(system, plur = "Plurality", rcv = "IRV"),
         iter = as.numeric(iter), 
         wrapcat = paste0(name, ", ", s))

# Change factor order for grid order in plot
level.byrow <- function(vec, nc){
  fac <- factor(vec) #if it is not a factor
  mlev <- matrix(levels(fac), nrow=nc, byrow=T)
  factor(fac, levels= c(mlev))
}
summary_df$wrapcat = level.byrow(summary_df$wrapcat, 3)

# Adjust NaN that occurred in US_2004 in magnitude when s = 10 (basically no-one would want to vote strategically)
summary_df = summary_df %>%
  mutate(value = case_when(
    name %in% c("Magnitude", "ExpBenefit") & is.na(value) ~ 0,
    TRUE ~ value
  ))

# Calculate weighted mean for each iteration
agg_df <- summary_df %>% 
  group_by(iter,  lambda, wrapcat, system) %>% 
  summarise(value = weighted.mean(value, case_weight)) %>%
  mutate(iter = as.numeric(iter))

head(agg_df)

# Get ratio
ratio = agg_df %>%
  group_by(iter, lambda, wrapcat) %>%
  summarise(ratio = value[system == "Plurality"] / value[system == "IRV"] ) %>%
  filter(grepl("ExpBenefit", wrapcat))

ratio %>% filter(iter %in% c(1, 60))

ggplot(ratio, aes(x = iter, y = ratio)) +
  geom_point(aes(colour = wrapcat)) + 
  geom_line(aes(colour = wrapcat)) +
  theme_tn()

# Get scatterplot

filterdf = summary_df %>% 
  filter(iter %in% c(1, 60), s == "s = 85") %>%
  select(-wrapcat) %>%
  pivot_wider(names_from = system, values_from = value)

p = list()
for(var in c("ExpBenefit", "Magnitude", "Prevalence")){
  maxdf = filterdf %>% filter(name == var)
  maxval = max(c(max(maxdf$IRV), max(maxdf$Plurality)))
  p[[paste0(var, ", 1st Iteration")]] = ggplot(maxdf %>% filter(iter == 1), 
                  aes(x = IRV, y = Plurality)) +
                geom_point(shape = 1) +
                geom_abline(slope = 1,
                  lty = "dashed") +
                theme_tn() +
                xlim(c(0, maxval)) +
                ylim(c(0, maxval)) +
                coord_fixed(ratio = 1) +
                ggtitle(paste0(var, ", 1st Iteration")) +
                theme(plot.title = element_text(size = 6))
  p[[paste0(var, ", 60st Iteration")]] = ggplot(maxdf %>% filter(iter == 60), 
                  aes(x = IRV, y = Plurality)) +
                geom_point(shape = 1) +
                geom_abline(slope = 1,
                  lty = "dashed") +
                theme_tn() +
                xlim(c(0, maxval)) +
                ylim(c(0, maxval)) +
                coord_fixed(ratio = 1) +
                ggtitle(paste0(var, ", 60th Iteration")) +
                theme(plot.title = element_text(size = 6))
}


library(patchwork)

(p[[1]] + p[[3]] + p[[5]]) / (p[[2]] + p[[4]] + p[[6]])

ggsave("output/figures/1/crosssection_s85.pdf", 
  device = cairo_pdf, width = 8, height = 5)


# Same thing for s = 10
filterdf = summary_df %>% 
  filter(iter %in% c(1, 60), s == "s = 10") %>%
  select(-wrapcat) %>%
  pivot_wider(names_from = system, values_from = value)

p = list()
for(var in c("ExpBenefit", "Magnitude", "Prevalence")){
  maxdf = filterdf %>% filter(name == var)
  maxval = max(c(max(maxdf$IRV), max(maxdf$Plurality)))
  p[[paste0(var, ", 1st Iteration")]] = ggplot(maxdf %>% filter(iter == 1), 
                  aes(x = IRV, y = Plurality)) +
                geom_point(shape = 1) +
                geom_abline(slope = 1,
                  lty = "dashed") +
                theme_tn() +
                xlim(c(0, maxval)) +
                ylim(c(0, maxval)) +
                coord_fixed(ratio = 1) +
                ggtitle(paste0(var, ", 1st Iteration")) +
                theme(plot.title = element_text(size = 6))
  p[[paste0(var, ", 60st Iteration")]] = ggplot(maxdf %>% filter(iter == 60), 
                  aes(x = IRV, y = Plurality)) +
                geom_point(shape = 1) +
                geom_abline(slope = 1,
                  lty = "dashed") +
                theme_tn() +
                xlim(c(0, maxval)) +
                ylim(c(0, maxval)) +
                coord_fixed(ratio = 1) +
                ggtitle(paste0(var, ", 60th Iteration")) +
                theme(plot.title = element_text(size = 6))
}


library(patchwork)

(p[[1]] + p[[3]] + p[[5]]) / (p[[2]] + p[[4]] + p[[6]])

ggsave("output/figures/1/crosssection_s10.pdf", 
  device = cairo_pdf, width = 8, height = 5)


# Same thing for s = 55
filterdf = summary_df %>% 
  filter(iter %in% c(1, 60), s == "s = 55") %>%
  select(-wrapcat) %>%
  pivot_wider(names_from = system, values_from = value)

p = list()
for(var in c("ExpBenefit", "Magnitude", "Prevalence")){
  maxdf = filterdf %>% filter(name == var)
  maxval = max(c(max(maxdf$IRV), max(maxdf$Plurality)))
  p[[paste0(var, ", 1st Iteration")]] = ggplot(maxdf %>% filter(iter == 1), 
                  aes(x = IRV, y = Plurality)) +
                geom_point(shape = 1) +
                geom_abline(slope = 1,
                  lty = "dashed") +
                theme_tn() +
                xlim(c(0, maxval)) +
                ylim(c(0, maxval)) +
                coord_fixed(ratio = 1) +
                ggtitle(paste0(var, ", 1st Iteration")) +
                theme(plot.title = element_text(size = 6))
  p[[paste0(var, ", 60st Iteration")]] = ggplot(maxdf %>% filter(iter == 60), 
                  aes(x = IRV, y = Plurality)) +
                geom_point(shape = 1) +
                geom_abline(slope = 1,
                  lty = "dashed") +
                theme_tn() +
                xlim(c(0, maxval)) +
                ylim(c(0, maxval)) +
                coord_fixed(ratio = 1) +
                ggtitle(paste0(var, ", 60th Iteration")) +
                theme(plot.title = element_text(size = 6))
}


library(patchwork)

(p[[1]] + p[[3]] + p[[5]]) / (p[[2]] + p[[4]] + p[[6]])

ggsave("output/figures/1/crosssection_s55.pdf", 
  device = cairo_pdf, width = 8, height = 5)