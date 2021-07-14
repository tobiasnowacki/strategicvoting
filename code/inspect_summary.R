library(tidyverse)

# Load data
load("output/files/1/85_sum_four.Rdata")


summary_stats_wide <- sum_df %>% 
  rename(ExpBenefit = eb,
         Magnitude = mag,
         Prevalence = prev) %>%
  pivot_longer(Prevalence:ExpBenefit) %>%
  mutate(iter = as.numeric(iter)) %>%
  mutate(system = case_when(
    system == "plur" ~ "Plurality",
    system == "rcv" ~ "RCV"
  )) %>%
  filter(system == "RCV") ## temporary zoom -- likely want to fix


# How many RCV cases have a magnitude of zero?

hist(sum_df$mag)
sum(sum_df$mag < 0, na.rm = TRUE)

nrow(sum_df)
sum_df %>% filter(mag < 0)

which(names(big_list_na_omit) == "BGR_2001")
big_list_na_omit$BGR_2001$U %>% head

names(bind_together)
obj <- bind_together[[1]]$rcv[[1]]
sum(obj$tau < 0)
which(tolower(obj$opt.votes.strategic) != obj$opt.votes.sincere)

obj$eu.mat[445, ]
obj$tau[445]
obj$opt.votes.sincere[445]
obj$opt.votes.strategic[445]

eu_test <- obj$eu.mat[445:446, ]
vzero_test <- obj$V0[445:446, ]

V.mat_test <- ballot_mat_from_eu_mat(
    eu_test,
    break_ties_with_sincerity = TRUE,
    sincere_mat = vzero_test, 
    normalize_eu_mat = TRUE
)
V.mat_test

length(obj)

sum_df %>%
    filter(mag == 0)

t <- ggplot(summary_stats_wide, aes(iter, value)) +
  geom_line(aes(group = interaction(system, case, name),
                colour = system), alpha = 0.05) +
  facet_wrap(. ~ name, scales = "free_y") +
  theme(legend.position = "bottom", legend.direction = "horizontal") +
  labs(x = "Degree of Strategicness (Iterations)")

ggsave(t, file = "temp1.pdf")
t
