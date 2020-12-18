# LOAD DEPENDENCIES
library(tidyverse)
library(pivotprobs)
library(gtools)
library(doParallel)
source("code/utils/new_sv_iter.R")
library(devtools)
source_url("https://raw.githubusercontent.com/tobiasnowacki/RTemplates/master/plottheme.R")

# Load data
load("output/big_list_2.RData")
vap <- read.csv("data/case_vap.csv", sep = "") # voting age pop.
cat("Data imported. \n")
source("code/prep_cses.R")  # data prep

nn = names(big_list_na_omit)

return_obj = function(case, lambda, s){
  load(paste0("output/files/", lambda, "/", s, "_", case, "_iterout.Rdata"))
  return(inner_list)
}

bind_together = map(nn, ~ return_obj(.x, 1, 85))
names(bind_together) = nn

get_sum_stats = function(obj){
  w = obj$rcv[[i]]$weights
  map_dfr(c("rcv", "plur"), function(y){
    names(obj[[y]]) = 1:length(obj[[y]])
    map(obj[[y]], ~
            tibble(
              prev = weighted.mean(.x$tau > 0, w),
              mag  = weighted.mean(.x$tau[.x$tau > 0], w[.x$tau > 0]),
              eb = prev * mag)) %>%
    bind_rows(.id = "iter") %>%
    mutate(system = y)
 })
}

sum_df = map_dfr(bind_together, ~ get_prevalence(.x), .id = "case")

ggplot(sum_df, aes(as.numeric(iter), mag)) +
  geom_line(aes(group = case, colour = system), alpha = 0.3) +
  facet_wrap(~ system)
  theme_tn()

