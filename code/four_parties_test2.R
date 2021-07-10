# Load libraries
library(tidyverse)
library(devtools)
devtools::install_github("tobiasnowacki/pivotprobs")
library(pivotprobs)
library(gtools)
library(stringr)

# Load data
source("code/utils/sv_iter_four.R")
load("output/big_list_4_parties.RData")
vap <- read.csv("data/case_vap.csv", sep = "") # voting age pop.
cat("Data imported. \n")

# not sure if I still need this, but we'll see...
four_party_ballots <- function() {
    apply(permutations(
        n = length(candidates),
        r = length(candidates),
        v = candidates, repeats.allowed = FALSE
    ),
    1,
    paste,
    collapse = ""
    )
}

source("code/prep_cses.R")  # data prep -- should be ind of no of parties

test_util <- big_list_na_omit[[102]]
names(big_list_na_omit)[[102]] # ITA_2006

# Load 4-candidate data
out <- sv_four(test_util$U, test_util$weights, s = 85, rule = "AV")

glimpse(out)

# plurality
out <- sv_iter(
    test_util$U,
    85,
    rule = "plurality",
    max.iterations = 80
)

# Let's get prevalence...
prevalence <- map_dbl(out, 
    ~ mean(.x$tau > 0, na.rm = TRUE))
pdat <- tibble(
    x = 1:length(prevalence), 
    y = prevalence,
    type = "Plurality"
)

ggplot(pdat, aes(x = x, y = y)) +
    geom_line()

# IRV
out_av <- sv_iter(
    test_util$U,
    85,
    rule = "AV",
    max.iterations = 80
)

# Let's get prevalence...
prevalence_av <- map_dbl(
    out_av,
    ~ mean(.x$tau > 0, na.rm = TRUE)
)
pdat_av <- tibble(
    x = 1:length(prevalence_av), 
    y = prevalence_av,
    type = "AV"
)

ggplot(pdat_av, aes(x = x, y = y)) +
    geom_line()

joint_df <- rbind(pdat, pdat_av)

ggplot(joint_df, aes(x, y)) +
    geom_line(aes(colour = type)) +
    labs(x = "Iteration", y = "Prevalence")
ggsave("output/figures/irv4_test.pdf")

out[[1]]$v.vec.before

out[[4]]$opt.votes.strategic %>% head
out[[4]]$opt.votes.sincere %>% head
test_util$U %>% head

prevalence

# This stuff is super noisy!!!


 [1] 0.4688826 0.4115983 0.2885431 0.3797737 0.3330976 0.2871287 0.3373409 0.1449788 0.2623762 0.2475248 0.3118812 0.2157001 0.3147100 0.3628006 0.2581330 0.2616690 0.4603960 0.2086280 0.2835926 0.3182461 0.1541726 0.2241867 0.2567185 0.3338048
[25] 0.2694484 0.3437058 0.1626591 0.2241867 0.3528996 0.5360679