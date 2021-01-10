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
lambda_param <- c("1", "2", "3")

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

agg_df <- summary_df %>% 
  group_by(iter,  lambda, wrapcat, system) %>% 
  summarise(value = weighted.mean(value, case_weight)) %>%
  mutate(iter = as.numeric(iter))

main_plot <- ggplot(summary_df %>% filter(lambda == 1 & iter < 61), aes(as.numeric(iter), value)) +
  geom_line(aes(group = interaction(case, system), colour = system), alpha = 0.1) +
  geom_point(data = hidden_scale_adj, alpha = 0) +
  scale_color_manual(values = c("#004C99", "#CC6600")) +
  geom_line(data = agg_df %>% filter(system == "IRV" & lambda == 1 & iter < 61), 
            color = "#004C99", lwd = 1.1) + 
  geom_line(data = agg_df %>% filter(system == "Plurality" & lambda == 1 & iter < 61), 
            color = "#CC6600", lwd = 1.1) + 
  facet_wrap(wrapcat ~ ., scales = "free", shrink = TRUE) +
  guides(colour = guide_legend(override.aes = list(alpha = 1))) +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  labs(x = "Degree of Strategicness (Iterations)", y = "Value", colour = "System") +
  theme_tn()
ggsave("output/figures/iterated_complete.pdf", main_plot, 
       width = 6, height = 6.5)
