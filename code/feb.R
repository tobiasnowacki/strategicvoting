library(tidyverse)

load("output/85_SWE_2014_iterout-2.Rdata")

tau_vec = inner_list$rcv[[177]]$tau
tau_vec_2 = cases_converge[[1]]$rcv[[88]]$tau

mean(tau_vec > 0)

map_dbl(inner_list$rcv, ~ mean(.x$tau > 0))
map_dbl(inner_list$rcv, ~ mean(tolower(.x$opt.votes.strategic) != .x$opt.votes.sincere))

inner_list$rcv[[1]]$weights

diff_vec = tolower(inner_list$rcv[[177]]$opt.votes.strategic) != inner_list$rcv[[177]]$opt.votes.sincere

diff_vec %>% head

tau_vec %>% head
tau_vec[2] - 5e-09
mean(diff_vec)

map_dbl(cases_converge[[1]]$rcv, ~ mean(.x$tau > 0))

w = cases_converge[[1]]$rcv[[1]]$w
map_dfr(cases_converge[[1]]$rcv, ~
      tibble(
        prev = weighted.mean(tolower(.x$opt.votes.strategic) != .x$opt.votes.sincere, w),
        mag  = weighted.mean(.x$tau[tolower(.x$opt.votes.strategic) != .x$opt.votes.sincere], w[tolower(.x$opt.votes.strategic) != .x$opt.votes.sincere]),
        eb = prev * mag)) 