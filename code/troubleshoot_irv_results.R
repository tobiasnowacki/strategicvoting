# Load dependencies
library(tidyverse)
library(here)
source(here("code/full_header.R"))  # fn's and data
source(here("code/prep_cses.R"))    # data prep
source(here("code/utils/ternary_functions.R"))    # data prep


# Load summary stats to identify odd cases
load("output/files/1/85summary_1_85.Rdata")

summary_stats %>%
    filter(iter > 100 & System == "IRV") %>%
    arrange(-Prevalence)

ggplot(summary_stats %>% filter(System == "IRV" & case == "SWE_2014"),
  aes(x = iter, y = Prevalence)) +
    geom_path(aes(group = case))


# Load v_vecs to plot odd cases
load("output/files/1/85v_vecs_1_85.Rdata")

swe_vvec = out$SWE_2014

glimpse(swe_vvec)

library(ggtern)

v_vec_swe = out$SWE_2014[[1]] %>%
  as.data.frame %>%
  mutate(iter = 1:251)

v_vec_df_tern = v_vec_swe %>%
    mutate(A = V1 + V2,
           B = V3 + V4,
           C = V5 + V6,
           # ternary transformation
           C = C + .5 *B,
           B = sqrt(3/4)*B)

ggplot(v_vec_df_tern, aes(x = C, y = B)) +
  geom_point(aes(colour = iter), alpha = 0.7) +
  geom_path(alpha = 0.3) +
  scale_colour_gradient(low = "red", high = "blue") +
  theme_minimal()



ggtern(out$SWE_2014[[1]], aes(V1 + V2, V3 + V4, V5 + V6)) +
    geom_path()

swe_vvec[[1]][110:140, ]

out = many_iterations_until_convergence(
                        big_list_na_omit[[32]], 
                        big_list_na_omit[[32]]$v_vec, 
                        0.05, 
                        85, 
                        0.001, 
                        200, ae = TRUE)

out$rcv_df %>% filter(iter == 181) %>% head

with(out$rcv_df %>% filter(iter == 181), mean(sin_rcv != opt_rcv))
with(out$rcv_df %>% filter(iter == 181), mean(tau_rcv > 0))

with(out$rcv_df %>% filter(iter == 180), CAB[2] - ABC[2])

with(out$rcv_df %>% filter(iter == 180), table(sin_rcv))

names(out)
out$rcv_br[179:182]
out$rcv_v_vec[179:183, ]

test = sv(big_list_na_omit[[32]]$U, big_list_na_omit[[32]]$weights, 
   v.vec = out$rcv_v_vec[181, ] %>% unlist, 85, "AV")

v1 = out$rcv_v_vec[180, ] %>% unlist
v2 = out$rcv_v_vec[181, ] %>% unlist

linweights = seq(0, 1, by = 0.02)
lin_combs = map(linweights, ~ .x * v1 + (1 - .x) * v2)

tt = map(lin_combs, ~ table(
    sv(big_list_na_omit[[32]]$U, 
       big_list_na_omit[[32]]$weights, 
   v.vec = .x, 85, "AV")$opt.votes.strategic) %>% as.data.frame)

ttdf = tt %>% bind_rows(.id = "lincomb") %>%
  complete(lincomb, Var1, fill = list(Freq = 0)) %>%
  group_by(lincomb) %>%
  mutate(prop = Freq / sum(Freq)) 

ggplot(ttdf, aes(lincomb, prop)) +
  geom_line(aes(colour = Var1, group = Var1)) +
  geom_point(aes(colour = Var1))


tt[[1]][[1]] %>% pivot_longer

    sv(big_list_na_omit[[32]]$U, 
       big_list_na_omit[[32]]$weights, 
   v.vec = lin_combs[[1]], 85, "AV")$opt.votes.strategic

names(test)
head(test[[8]])