library(tidyverse)
library(pivotprobs)
library(devtools)
source_url("https://raw.githubusercontent.com/tobiasnowacki/RTemplates/master/plottheme.R")

library(here)                       # to get dir
source(here("code/full_header.R"))  # fn's and data
source(here("code/prep_cses.R"))    # data prep

# Create list of vvecs
start = rep(c(1/3, 1/3, 1/3)/2, each = 2) 
end = rep(c(1/2, 1/2, 0), each = 2) + runif(n = 6, 0, 0.00000001)
steps = seq(0, 1, length.out = 100)
vecs = lapply(steps, function(x) x * start + (1 - x) * end)

# Calculate pivotal probabilities
piv_old = map(vecs, ~ av.pivotal.event.probs.general(v.vec = c(.x, 0,0,0), s.vec = rep(85, 4)))
piv_new = map(vecs, ~ irv_election(n = 10000) %>%
            election_event_probs(method = "en", alpha = (.x * 85)))
piv_mc = map(vecs, ~ irv_election(n = 10000) %>%
            election_event_probs(method = "mc", 
                                 alpha = (.x * 85),
                                 num_sims = 100000))

# Extract AB piv probs
ab_old = map_dbl(piv_old, ~ .x$AB) / 10000
ab_new = map_dbl(piv_new, ~ .x$a_b$integral)
ab_mc = map_dbl(piv_mc, ~ .x$a_b$integral) 

# Get into table
df = tibble(
    c_share = map_dbl(vecs, ~ .x[5]),
    iter = 1:100,
    old = ab_old,
    new = ab_new,
    mc = ab_mc
    )  %>%
    pivot_longer(old:mc)

ggplot(df, aes(c_share, value)) +
    geom_line(aes(colour = name)) +
    labs(x = "Smallest Share", 
         y = "Pr(AB event)", 
         colour = "Method") +
    theme_tn()
ggsave("output/figures/pivot_trouble.pdf",
    device = cairo_pdf)