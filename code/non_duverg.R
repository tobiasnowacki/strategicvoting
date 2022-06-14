# Load dependencies
library(tidyverse)
library(rio)
library(gtools)
library(pivotprobs)

# Load data
load("output/files/1/85_winners_tbl.Rdata")
load("output/files/1/85_winners.Rdata")
load("output/big_list_2.RData")
vap <- read.csv("data/case_vap.csv", sep = "") # voting age pop.

# Load scripts and templates
source("code/utils/sv_theme_template.R")
source("code/utils/new_sv_iter.R")
source("code/prep_cses.R") # data prep


# Load vvec data
load("output/files/1/85_vvec.Rdata")

c_share_60 <- map_dbl(vvecdf, ~ .x$plurvec$c[60])
c_share_250 <- map_dbl(vvecdf, ~ .x$plurvec$c[250])
br_share_60 <- map_dbl(vvecdf, ~ .x$plurbr$C[2])
br_share_250 <- map_dbl(vvecdf, ~ .x$plurbr$C[250])

hist_data <- data.frame(
  c_share_60, 
  c_share_250,
  br_share_60,
  br_share_250
) %>%
  pivot_longer(c_share_60:br_share_250)

ggplot(hist_data %>% filter(value < 0.1), aes(x = value)) +
  geom_density() +
  facet_wrap(. ~ name, scales = "free") +
  theme_tn()


# toy_u <- data.frame(
#   A = c(7, 7, 7, 7, 5, 5, 5, 5, 1.999, 2.001),
#   B = c(4, 4, 4, 4, 7, 7, 7, 7, 2.001, 1.999),
#   C = c(2, 2, 2, 2, 2, 2, 2, 2, 8, 8)
# )

delta <- c(2:10 %o% 10^(-2:-5)) %>% sort()
u_vals <- seq(2.1, 9, by = 0.1)

gen_toy_example <- function(delta, u_c){
  
  min <- 2 - delta
  max <- 2 + delta
  toy_u <- data.frame(
    A = c(7, 7, 7, 7, 5, 5, 5, 5, min, max),
    B = c(4, 4, 4, 4, 7, 7, 7, 7, max, min),
    C = c(2, 2, 2, 2, 2, 2, 2, 2, u_c, u_c)
  )

  x <- sv(
    U = toy_u,
    s = 85,
    # v.vec = c(0.4, 0.4, 0.2),
    rule = "plurality"
  )

  return(x$best.response.v.vec[3])

}

param_grid <- expand.grid(delta = delta, u_c = u_vals)

param_grid$share_br_c <- apply(param_grid, 1, function(x) gen_toy_example(x[1], x[2]))

tbl <- param_grid %>% 
  group_by(delta) %>%
  summarise(
    min_c = min(u_c[share_br_c == 0.2]),
    max_c = max(u_c[share_br_c == 0.2]),
    min_x = min(u_c[share_br_c == 0]),
    max_x = max(u_c[share_br_c == 0])
  ) %>% 
  filter(delta > 0)

tbl$min_c[tbl$min_c == Inf] <- 9

ggplot(tbl, aes(x = delta, y = min_c)) +
  geom_point() +
  scale_x_log10() +
  labs(x = "Difference btw u(A) and u(B)", y = "min u(C) to produce eqm") +
  theme_tn()

ggsave("grids.pdf", device = cairo_pdf)





# On client, load the data as follows:
load("output/files/1/85_workbench.Rdata")

# For SV statistics, convert iteration to numeric
sum_df <- sum_df %>% mutate(iter = as.numeric(iter))

# Function to find condorcet winner from 6-item vvec
return_condorcet <- function(vvec){
    a_b <- sum(vvec[c(1, 2, 5)]) - sum(vvec[c(3, 4, 6)])
    a_c <- sum(vvec[c(1, 2, 3)]) - sum(vvec[c(4, 5, 6)])
    b_c <- sum(vvec[c(1, 3, 4)]) - sum(vvec[c(2, 5, 6)])

    if (a_b > 0 & a_c > 0) {
       return("A")
    } else if (a_c < 0 & b_c < 0) {
       return("C")
    } else if (a_b < 0 & b_c > 0) {
       return("B")
    } else {
       return("NA")
    }
}

# Apply to vvecs
c_winner <- map_chr(
    big_list_na_omit, ~
    return_condorcet(.x$v_vec)
    ) %>% 
    enframe() %>%
    rename(case = name, cwinner = value)

# Merge winners, EU and summary stats together
ddf <- winners_df %>%
    left_join(win_df) %>%
    left_join(sum_df) %>%
    left_join(case_df) %>%
    left_join(c_winner) %>%
    rename(
        winpr_a = V1,
        winpr_b = V2,
        winpr_c = V3,
        avg_exp_eu = value
    ) %>%
    mutate(
        system_group = case_when(
            system %in% c("pr", "mixed") ~ "pr",
            system %in% c("fptp", "tworound", "parallel") ~ "plur",
            system %in% c("rcv") ~ "rcv", 
        TRUE ~ "other"
        ),
        cwin_prob = case_when(
            cwinner == "A" ~ winpr_a,
            cwinner == "B" ~ winpr_b,
            cwinner == "C" ~ winpr_c
        )
    )


head(ddf)

winpr_a_by_case <- ddf %>%
  group_by(iter, case, cwinner) %>%
  summarise(
    diff_winpr_a = winpr_a[system == "rcv"] - 
      winpr_a[system == "plurality"]
    ) %>%
  group_by(cwinner, iter) %>%
  summarise(
    mean_diff = mean(diff_winpr_a)
  ) %>%
  filter(cwinner != "NA")

ggplot(winpr_a_by_case, aes(x = iter, y = mean_diff)) +
  geom_hline(yintercept = 0, lty = "dotted") +
  geom_line(aes(colour = cwinner)) +
  labs(
    x = "Iteration", 
    y = "Pr(A Wins | RCV) - Pr(A Wins | Plur)",
    colour = "Condorcet Winner"  
  ) +
  theme_tn()

ggsave("a_win.pdf", device = cairo_pdf)


head(winpr_a_by_case)


