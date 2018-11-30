##################################################
## Project: Strategic Voting in RCV
## Script purpose: Radar plot for piv.probs
##					 in the CSES case
## Date: 30/11/2018
## Author:
##################################################

### 
### Dependencies
###

# Set WD etc.
library(here)

# Load pivotal probability functions:
av_piv_path <- here("utils/av_pivotal_probs_analytical_general_v2.r")
source(av_piv_path)
plur_piv_path <- here("utils/plurality_pivotal_probabilities_analytical.r")
source(plur_piv_path)

# To replicate Andy's function(s):
sim_appr2 <- here("utils", "general_iteration_simulation_approach.r")
source(sim_appr2)
sv_file <-  here("utils/sv.r")
source(sv_file)

# Load my own functions:
functions <- here("utils/functions.r")
source(functions)

# Load necessary libraries:
library(ggplot2)
library(reshape2)
library(dplyr)
library(purrr)

# Load CSES data:
load(here("../output/cses_big_list_2.RData"))

###
### ANALYSIS
###

# get DF of piv_probs
pp_list <- list()
for(i in 1:length(big_list)){
  print(i)
  x <- big_list[[i]]
  sv_obj <- sv(x$U, weights = x$weights, s = 60, rule = "AV")
  pp_list[[i]] <- sv_obj$piv.probs
}
pp_df <- as.data.frame(do.call(rbind, pp_list))
pp_df <- as.data.frame(apply(pp_df, 2, function(x) as.numeric(x)))
pp_df$group <- as.factor(rownames(pp_df))

pp_df_long <- melt(pp_df[, 4:13], id = c("group"))
names(pp_df_long) <- c("group", "event", "prob")

pp_plot <- ggplot(pp_df_long, aes(x = event, y = prob)) +
  geom_path(aes(group = group), alpha = 0.1) +
  coord_polar() + 
  theme_bw()
ggsave(here("../output/figures/cses_piv_probs.pdf"), pp_plot)
