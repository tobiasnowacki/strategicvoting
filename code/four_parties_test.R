# Testing four party

library(tidyverse)
library(pivotprobs)
library(gtools)



#' Compute pivot probabilities and generate P-matrices for four-candidate IRV elections
#'
#' @param sims A tibble of simulations, with an `id` column and
#' a column for each distinct ballot
#' @param n The number of voters.
#' @param reporting The level of reporting detail (0: none; 1:
#' type of pivot event; 2: name of pivot event)
#' @return A list of pivot events, each of which is a list with elements
#' `integral` and `P`
#' @examples
#' sims <- simulate_ordinal_results_from_dirichlet(k = 4, n = 10000)
#' out <- irv_pivot_probs_four_cands(sims)
#' out %>% combine_P_matrices()
#' @export
irv_pivot_probs_four_cands <- function(sims, n = 1000, reporting = 1) {
    PP_LIBRARY <<- list() ## warning: writes to global variable
    if (reporting >= 1) {
        cat("Round 0: ")
    }
    round_0 <- sims %>% round_0_pivot_probs(reporting = reporting)
    if (reporting >= 1) {
        cat("done.\nRound 1: ")
    }
    round_1 <- sims %>% round_1_pivot_probs(reporting = reporting)
    if (reporting >= 1) {
        cat("done.\nRound 2: ")
    }
    round_2 <- sims %>% last_round_pivot_probs(reporting = reporting)
    if (reporting >= 1) {
        cat("done.\n")
    }
    c(round_0, round_1, round_2)
}

condense_mat <- function(mat) {
    # replace NA with zero
    mat[is.na(mat)] <- 0
    cns <- colnames(mat)
    ucns <- unique(colnames(mat))
    new_mat <- matrix(0, nrow = nrow(mat), ncol = length(ucns))
    colnames(new_mat) <- ucns
    for (j in 1:ncol(mat)) {
        new_mat[, cns[j]] <- new_mat[, cns[j]] + mat[, j]
    }
    new_mat
}

drop_candidate_and_condense_matrix <- function(sims, cand_to_drop = "A") {
    # drop
    colnames(sims) <- str_replace(colnames(sims), cand_to_drop, "")
    # condense
    condense_mat(sims)
}

drop_candidate_and_condense <- function(tb, cand_to_drop = "A") {
    ids <- tb$id

    tb %>%
        select(-id) %>%
        as.matrix() %>%
        drop_candidate_and_condense_matrix(cand_to_drop) %>%
        as_tibble() %>%
        mutate(id = ids) %>%
        relocate(id)
}

rank_mat <- function(frs) {
    out <- matrix(1, nrow = nrow(frs), ncol = ncol(frs))
    colnames(out) <- colnames(frs)
    for (j in 1:(ncol(out) - 1)) { # j is the column of the candidate we are ranking
        for (k in 1:ncol(frs)) { # k is the column of the candidate we are comparing to
            if (j == k) {
                next
            }
            out[, j] <- out[, j] + as.integer(frs[, j] < frs[, k])
        }
    }
    out[, ncol(out)] <- sum(1:ncol(out)) - apply(out[, 1:(ncol(out) - 1)], 1, sum)
    out
}

# test: much faster than ranking, and gets it right.
# frs <- gtools::rdirichlet(100000, alpha = rep(5, 3))
# system.time(rank_slow <- apply(frs, 1, rank))
# system.time(rank_fast <- rank_mat(frs))


get_loser2 <- function(sims) {
    # optimized I think
    loser <- rep("", nrow(sims))
    the_min <- rep(1, nrow(sims))
    for (j in 1:ncol(sims)) {
        this_is_lower <- sims[, j] > 0 & sims[, j] < the_min
        the_min[this_is_lower] <- sims[this_is_lower, j]
        loser[this_is_lower] <- colnames(sims)[j]
    }
    loser
}

get_first_rank_shares <- function(sims) {
    colnames(sims) <- str_extract(colnames(sims), "^.")
    condense_mat(sims)
}

get_loser_from_tibble_of_sims <- function(df) {
    df %>%
        select(-id) %>%
        as.matrix() %>%
        get_first_rank_shares() %>%
        get_loser2()
}

get_winner_for_each_elimination <- function(sims) {
    cands <- names(sims)[2] %>%
        str_split("") %>%
        .[[1]] %>%
        sort()

    tibble(dropped_cand = cands, sims = list(sims)) %>%
        mutate(sims = map2(sims, dropped_cand, drop_candidate_and_condense)) %>%
        mutate(winners = map(sims, irv_winners)) %>%
        select(-sims) %>%
        unnest(cols = c(winners)) %>%
        pivot_wider(names_from = dropped_cand, values_from = winner)
}

parse_event_name <- function(event) {
    cands <<- event %>%
        str_replace("_", "") %>%
        str_split("") %>%
        .[[1]]
    w <<- cands[1]
    x <<- cands[2]
    y <<- cands[3]
    z <<- cands[4]
}

round_0_irv_pivot_prob <- function(sims, n = 1000, event = "AB_CD", wfee = NULL, noisy = F) {
    # wfee stands for "winner for each elimination"
    # assign candidate names to variables
    parse_event_name(event)
    # get first rank shares and first rank ranks
    frs <- PP_LIBRARY[["frs"]]
    if (is.null(frs)) {
        frs <- get_first_rank_shares(sims %>% select(-id) %>% as.matrix())
        PP_LIBRARY[["frs"]] <<- frs
    }
    ranks <- PP_LIBRARY[["ranks"]]
    if (is.null(ranks)) {
        ranks <- rank_mat(frs)
        PP_LIBRARY[["ranks"]] <<- ranks
    }
    # condition 1: w and x in last two places
    cond1 <- ranks[, x] + ranks[, w] == 2 * ncol(frs) - 1
    if (noisy) {
        cat("cond1 passed by ", sum(cond1, na.rm = T), " cases.\n")
    }
    if (is.null(wfee)) {
        # who would win given each
        wfee <- get_winner_for_each_elimination(sims)
    }
    cond2 <- wfee[, x] == y
    if (noisy) {
        cat("cond2 passed by ", sum(cond2, na.rm = T), " cases.\n")
    }
    cond3 <- wfee[, w] == z
    if (noisy) {
        cat("cond3 passed by ", sum(cond3, na.rm = T), " cases.\n")
    }
    cond <- cond1 & cond2 & cond3
    if (sum(cond, na.rm = T) == 0) {
        return(0)
    }
    # now compute density
    w_vs_x <- frs[, w] - frs[, x]
    the_density <- try(density_estimate(x = w_vs_x[cond], eval.points = c(0)), silent = T)
    if (class(the_density) == "try-error") {
        msg <- geterrmessage()
        cat("Error: ", msg, "-- setting density to zero.\n")
        the_density <- 0
    }
    mean(cond) * the_density * (1 / n)
}

make_P_from_event_name <- function(event, cands = c("A", "B", "C", "D")) {
    parse_event_name(event)
    P <- make_empty_P(cands)
    ballots <- colnames(P)
    z_wins <- str_detect(ballots, str_c("^", x))
    P[y, !z_wins] <- 1
    P[z, z_wins] <- 1
    P
}

round_0_pivot_probs <- function(sims, n = 1000, reporting = 1) {
    out <- list()
    PP_LIBRARY <<- list() # clearing out by default
    wfee <- get_winner_for_each_elimination(sims)
    cands <- names(sims)[2] %>%
        str_split("") %>%
        .[[1]]
    for (w in cands) {
        for (x in cands) {
            if (w >= x) {
                next
            }
            for (y in cands) {
                if (x == y) {
                    next
                }
                for (z in cands) {
                    if (y == z | w == z) {
                        next
                    }
                    event <- str_c(w, x, "_", y, z)
                    if (reporting >= 2) {
                        cat(event)
                    }
                    pp <- sims %>% round_0_irv_pivot_prob(n = n, event = event, wfee = wfee)
                    this_P <- make_P_from_event_name(event, cands)
                    out[[event]] <- list(integral = pp, P = this_P)
                    # mirror event
                    event <- str_c(x, w, "_", z, y)
                    if (reporting >= 2) {
                        cat(" | ", event, "\n")
                    }
                    this_P <- make_P_from_event_name(event, cands)
                    out[[event]] <- list(integral = pp, P = this_P)
                }
            }
        }
    }
    out
}

PP_LIBRARY <- list()

# initially all events had zero probability because cond1 could not be met -- was operating on frs and not ranks.

parse_1r_event_name <- function(event) {
    cands <<- event %>%
        str_replace("\\.", "") %>%
        str_replace("_", "") %>%
        str_split("") %>%
        .[[1]]
    v <<- cands[1]
    w <<- cands[2]
    x <<- cands[3]
    y <<- cands[4]
    z <<- cands[5]
}

round_1_irv_pivot_prob <- function(sims, n = 1000, event = "D.AB_AC", wfee = NULL, noisy = F) {
    # wfee stands for "winner for each elimination"
    # assign candidate names to variables
    parse_1r_event_name(event)
    # get first rank shares and first rank ranks
    frs <- PP_LIBRARY[["frs"]]
    if (is.null(frs)) {
        frs <- get_first_rank_shares(sims %>% select(-id) %>% as.matrix())
        PP_LIBRARY[["frs"]] <<- frs
    }
    ranks <- PP_LIBRARY[["ranks"]]
    if (is.null(ranks)) {
        ranks <- rank_mat(frs)
        PP_LIBRARY[["ranks"]] <<- ranks
    }
    # condition 1: v is last in 0th round
    cond1 <- ranks[, v] == ncol(frs)
    if (noisy) {
        cat("cond1 passed by ", sum(cond1, na.rm = T), " cases.\n")
    }
    # now we shift to the next round
    # condition 2: neither w nor x wins the first round
    sims1_v <- PP_LIBRARY[[str_c("sims_drop_", v)]]
    if (is.null(sims1_v)) {
        sims1_v <- drop_candidate_and_condense(sims, cand_to_drop = v)
        PP_LIBRARY[[str_c("sims_drop_", v)]] <<- sims1_v
    }
    frs1_v <- PP_LIBRARY[[str_c("frs_drop_", v)]]
    if (is.null(frs1_v)) {
        frs1_v <- get_first_rank_shares(sims1_v %>% select(-id) %>% as.matrix())
        PP_LIBRARY[[str_c("frs_drop_", v)]] <<- frs1_v
    }
    ranks1_v <- PP_LIBRARY[[str_c("ranks_drop_", v)]]
    if (is.null(ranks1_v)) {
        ranks1_v <- rank_mat(frs1_v)
        PP_LIBRARY[[str_c("ranks_drop_", v)]] <<- ranks1_v
    }
    cond2 <- ranks1_v[, w] != 1 & ranks1_v[, x] != 1
    if (noisy) {
        cat("cond2 passed by ", sum(cond2, na.rm = T), " cases.\n")
    }
    if (is.null(wfee)) {
        # who would win given each elimination
        wfee <- get_winner_for_each_elimination(sims1)
    }
    # cond3: eliminate x and get y
    cond3 <- wfee[, x] == y
    if (noisy) {
        cat("cond3 passed by ", sum(cond3, na.rm = T), " cases.\n")
    }
    # cond4: eliminate w and get z
    cond4 <- wfee[, w] == z
    if (noisy) {
        cat("cond4 passed by ", sum(cond4, na.rm = T), " cases.\n")
    }
    cond <- cond1 & cond2 & cond3 & cond4
    if (sum(cond, na.rm = T) == 0) {
        return(0)
    }
    # now compute density
    w_vs_x <- frs1_v[, w] - frs1_v[, x]
    the_density <- try(density_estimate(x = w_vs_x[cond], eval.points = c(0)), silent = T)
    if (class(the_density) == "try-error") {
        msg <- geterrmessage()
        cat("Error: ", msg, "-- setting density to zero.\n")
        the_density <- 0
    }
    mean(cond) * the_density * (1 / n)
}

# this could be expanded to deal with multiple previous rounds
round_1_pivot_probs <- function(sims, n = 1000, reporting = 1) {
    out <- list()
    PP_LIBRARY <<- list() # clear out by default -- slightly inefficient
    cands <- names(sims)[2] %>%
        str_split("") %>%
        .[[1]]
    # v.wx_yz restrictions: v is not w, x, y, or z. w is not z, x is not y, y is not z.
    for (v in cands) {
        drop_v_sims <- drop_candidate_and_condense(sims, cand_to_drop = v)
        wfee <- get_winner_for_each_elimination(drop_v_sims)
        for (w in cands) {
            for (x in cands) {
                if (w >= x) {
                    next
                }
                for (y in cands) {
                    if (x == y) {
                        next
                    }
                    for (z in cands) {
                        if (y == z | w == z | v %in% c(w, x, y, z)) {
                            next
                        }
                        event <- str_c(v, ".", w, x, "_", y, z)
                        if (reporting >= 2) {
                            cat(event)
                        }
                        pp <- sims %>% round_1_irv_pivot_prob(n = n, event = event, wfee = wfee)
                        this_P <- make_P_from_1r_event_name(event, cands)
                        out[[event]] <- list(integral = pp, P = this_P)
                        # mirror event
                        event <- str_c(v, ".", x, w, "_", z, y)
                        if (reporting >= 2) {
                            cat(" | ", event, "\n")
                        }
                        this_P <- make_P_from_1r_event_name(event, cands)
                        out[[event]] <- list(integral = pp, P = this_P)
                    }
                }
            }
        }
    }
    out
}

make_empty_P <- function(cands) {
    k <- length(cands)
    P <- matrix(0, nrow = k, ncol = factorial(k))
    rownames(P) <- cands
    ballots <- gtools::permutations(n = k, r = k, v = cands) %>%
        apply(1, paste, collapse = "")
    colnames(P) <- ballots
    P
}

# v.wx_yz P: y wins unless you rank x first or v first and x second.
make_P_from_1r_event_name <- function(event, cands = c("A", "B", "C", "D")) {
    parse_1r_event_name(event)
    P <- make_empty_P(cands)
    ballots <- colnames(P)
    z_wins <- str_detect(ballots, str_c("^", x)) | str_detect(ballots, str_c("^", v, x))
    P[y, !z_wins] <- 1
    P[z, z_wins] <- 1
    P
}

# last round pivot events: we just crank until we have only two
get_to_last_round <- function(sims) {
    while (str_length(colnames(sims)[2]) > 2) {
        sims <- one_round_of_irv(sims)
    }
    sims
}

make_P_from_last_round_event_name <- function(event, cands = c("A", "B", "C", "D")) {
    event_split <- event %>%
        str_split("") %>%
        .[[1]]
    y <- event_split[1]
    z <- event_split[2]

    P <- make_empty_P(cands)
    ballots <- colnames(P)
    y_before_z <- str_detect(ballots, str_c(y, ".*", z))
    P[y, y_before_z] <- 1
    P[z, !y_before_z] <- 1
    P
}

last_round_pivot_probs <- function(sims, n = 1000, reporting = 1) {
    cands <- names(sims)[2] %>%
        str_split("") %>%
        .[[1]]
    sims_last <- get_to_last_round(sims)
    sims_last %>%
        pivot_longer(cols = -id) %>%
        filter(!is.na(value) & value > 0) %>%
        group_by(id) %>%
        slice(1) %>%
        separate(col = name, into = c("y", "z"), sep = 1) %>%
        mutate(y_vs_z = 2 * (value - .5)) -> margins
    out <- list()
    for (y in cands) {
        for (z in cands[cands > y]) {
            event <- str_c(y, z)
            if (reporting >= 2) {
                cat(event)
            }
            cond <- margins$y == y & margins$z == z
            the_density <- density_estimate(x = margins$y_vs_z[cond], eval.points = c(0))
            pp <- mean(cond) * the_density * (1 / n)
            this_P <- make_P_from_last_round_event_name(event, cands)
            out[[event]] <- list(integral = pp, P = this_P)
            event <- str_c(z, y)
            if (reporting >= 2) {
                cat(" | ", event, "\n")
            }
            this_P <- make_P_from_last_round_event_name(event, cands)
            out[[event]] <- list(integral = pp, P = this_P)
        }
    }
    out
}

#' Compute IRV winners for simulated elections
#'
#' Takes a data frame of simulated ordinal elections (one election per row,
#' one ballot ordering per column) and returns a data fram of IRV winners.
#'
#' @param sims A tibble of simulations, with an `id` column and
#' a (named) column for each distinct ballot (e.g. ABCD, ABDC,...)
#' @return A tibble with columns `id` and `winner`
#' @examples
#' sims <- simulate_ordinal_results_from_dirichlet(k = 3, n = 10)
#' irv_winners(sims)
irv_winners <- function(sims) {
    stopifnot(length(unique(colnames(sims))) == ncol(sims))
    while (!near(max(as.matrix(sims)[1, -which(names(sims) == "id")], na.rm = T), 1)) {
        sims <- one_round_of_irv(sims)
    }
    sims %>%
        pivot_longer(cols = -id, names_to = "winner", values_to = "won") %>%
        filter(near(won, 1)) %>%
        select(id, winner)
}

one_round_of_irv <- function(sims) {
    sims %>%
        mutate(loser = get_loser_from_tibble_of_sims(.)) %>%
        # now we group so that we only have to drop and condense once per losing candidate (rather than once per row)
        group_by(loser) %>%
        nest() %>%
        mutate(reduced = map(data, drop_candidate_and_condense, loser)) %>%
        ungroup() %>%
        select(reduced) %>%
        unnest(cols = c(reduced)) -> out

    out %>%
        colSums(na.rm = T) -> col_sums

    col_sums[col_sums > 0] %>%
        names() %>%
        sort() -> non_zero_cols

    out %>%
        select(all_of(non_zero_cols)) %>%
        arrange(id) %>%
        relocate(id)
}

#' Faster Monte Carlo estimation of pivot probabilities
#'
#' These methods take simulated election results and return a list of
#' pivot events, each with a pivot probability (\code{integral}) and
#' P matrix (\code{P}). They produce identical results to
#' \code{election_event_probs()} (or nearly so) but they are faster because
#' they are written in a less general way.
#'
#' For plurality elections, any number of candidates can be specified.
#' The other voting systems can only handle three candidates.
#'
#' Results have been validated against \code{election_event_probs()}.
#'
#' @param sims A matrix of simulated election results, with one column per
#' ballot type. Must be 6 columns for the ordinal methods (IRV, Kemeny-Young, positional).
#' @param n Size of electorate
#' @param window Window within which two candidates are considered to be tied
#' when \code{method="rectangular"}. Wider window means lower variance but
#' more bias.
#' @param s Score allocated to a second-place ranking in positional
#' and IRV elections.
#' @param cand_names Names of the candidates.
#' @param sep Separation between candidate names.
#' @param merge Merge adjacent pivot events?
#' @param drop Drop a dimension?
#' @param skip_non_pivot_events Skip non pivot events?
#' @param raw Return counts rather than integrals/proportions?
#' @param kemeny For Condorcet method, compute Kemeny-Young pivot event probabilities? If \code{F}, only handles event in which the Condorcet winner
#' is decided.
#'
#' @examples
#' sims <- gtools::rdirichlet(100000, alpha = c(9,7,3,4,4,6))
#' plurality_pivot_probs_from_sims(sims)
#' positional_pivot_probs_from_sims(sims, s=.5)
#' irv_pivot_probs_from_sims(sims)
#' condorcet_pivot_probs_from_sims(sims)
#'
#'
#'@name standalone_monte_carlo_methods
NULL

#' @export
brute_force_mc_eep <- function(num_sims = 10000000, batch_size = 500000, alpha = c(10, 9, 4), election_method = "plurality", n = 1000, s = NULL){
  if(election_method == "positional" & is.null(s)){stop("You must provide a value for `s` if `election_method=positional`.")}
  if(!is.null(s) && !(s >=0 & s<=1)){stop("s must be between 0 and 1.")}
  sims_so_far <- 0
  count_df_list <- list()
  i <- 1
  while(sims_so_far < num_sims){
    cat(".")
    this_batch_size <- min(batch_size, num_sims - sims_so_far)
    sims <- gtools::rdirichlet(this_batch_size, alpha)
    f <- str_c(election_method, "_event_probs_from_sims") %>% get()
    these_counts <- f(sims, method = "rectangular", skip_non_pivot_events = F, window = 1/n, s = s, raw = T) %>%
      map("integral")
    this_count_df <- tibble(event = names(these_counts), tally = these_counts %>% unlist() %>% as.numeric())
    count_df_list[[i]] <- this_count_df
    sims_so_far <- sims_so_far + this_batch_size
    i <- i + 1
  }
  cat("\n")
  bind_rows(count_df_list) %>%
    group_by(event) %>%
    summarize(tally = sum(tally)) %>%
    ungroup() %>%
    mutate(prob = tally/sum(tally))
}

density_estimate <- function(x, bw_divisor = 1, eval.points = c(0)){
  if(length(x) <= 1){return(rep(0, length(eval.points)))} # can't get a bandwidth with only 1 point
  bw <- ks::hpi(x, binned = T)/bw_divisor
  ks::kde(x = x, h = bw, eval.points = eval.points)$estimate
}

# I think this is no longer used - was applied only to Condorcet
pivot_prob_from_delta_and_condition <- function(delta, condition, method = "density", bw_divisor = 1, window = .01, n = 1000, eval.points = c(0,0)){

  if(method == "rectangular"){
    rep(mean(condition & abs(delta) < window/2)/(n*window), 2)
  }else if(method %in% c("naive_density", "density")){
    mean(condition)*density_estimate(delta[condition], bw_divisor = bw_divisor, eval.points = eval.points)*(1/n)
  }else{
    stop("Unknown method for estimating density.")
  }
}

# making more general
pivot_prob_from_delta_and_condition2 <- function(delta, condition, method = "density", bw_divisor = 1, window = .01, n = 1000, drop = T, merge = F, raw = F, increments = 10){
  if(method == "rectangular"){
    offset <- ifelse(merge, 0, 1/(2*n))
    cond1 <- condition & abs(delta - offset) < window/2
    if(merge){
      cond2 <- cond1
    }else{
      cond2 <- condition & abs(delta + offset) < window/2
    }
    if(raw){c(sum(cond1), sum(cond2))}else{c(mean(cond1), mean(cond2))/(window*n)}
  }else if(method %in% c("naive_density", "density")){
    limits1 <- c(0, 1/n)
    if(merge){limits1 <- limits1 - 1/(2*n)}
    if(drop){increments <- 1}
    som1 <- sequence_of_midpoints(limits1[1], limits1[2], increments = increments)
    if(drop){grid_width <- limits1[2] - limits1[1]}else{grid_width <- som1[2] - som1[1]}
    prob1 <- mean(condition)*sum(density_estimate(delta[condition], bw_divisor = bw_divisor, eval.points = som1))*grid_width
    if(merge){
      prob2 <- prob1
    }else{
      limits2 <- c(-1/n, 0)
      som2 <- sequence_of_midpoints(limits2[1], limits2[2], increments = increments)
      if(drop){grid_width <- limits2[2] - limits2[1]}else{grid_width <- som2[2] - som2[1]}
      prob2 <- mean(condition)*sum(density_estimate(delta[condition], bw_divisor = bw_divisor, eval.points = som2))*grid_width
    }
    c(prob1, prob2)
  }else{
    stop("Unknown method for estimating density.")
  }
}


#' @rdname standalone_monte_carlo_methods
#' @export
plurality_event_probs_from_sims <- function(sims = NULL, n = 1000, window = .01, cand_names = NULL, sep = "_", method = "density", merge = F, drop = F, bw_divisor = 1, skip_non_pivot_events = T, s = NULL, raw = F){

  out <- list()
  if(is.null(cand_names)){cand_names <- letters[1:ncol(sims)]}
  if(is.null(sims)){stop("You need to provide `sims`.")}
  for(i in 1:(ncol(sims)-1)){
    for(j in (i+1):ncol(sims)){
      pp <- ab_plurality_tie_for_first_from_sims(cbind(sims[,c(i,j), drop = F], sims[,-c(i,j), drop = F]), method = method, n = n, merge = merge, drop = drop, window = window, bw_divisor = bw_divisor, raw = raw)
      out[[paste0(cand_names[i], sep, cand_names[j])]] <- list(
        integral = pp[1],
        P = plurality_P_matrix_from_indices(i,j,ncol(sims))
      )
      out[[paste0(cand_names[j], sep, cand_names[i])]] <- list(
        integral = pp[2],
        P = plurality_P_matrix_from_indices(j,i,ncol(sims))
      )
    }
  }

  if(!skip_non_pivot_events){
    empty_P <- matrix(0, nrow = length(cand_names), ncol = length(cand_names) + 1)
    for(i in 1:ncol(sims)){
      this_P <- empty_P
      this_P[i,] <- 1
      out[[paste0(cand_names[i], "_")]] <- list(
        integral = a_wins_from_sims_prob(sims = cbind(sims[,i,drop=F], sims[,-i,drop=F]), n = n, raw = raw),
        P = this_P)
    }
  }

  out
}

plurality_P_matrix_from_indices <- function(i,j,k){
  out <- matrix(0, nrow = k, ncol = k+1)
  out[i,] <- 1
  out[c(i,j),j] <- c(0,1)
  out
}

a_wins_from_sims_prob <- function(sims, n = 1000, raw = F){
  cond <- rep(T, nrow(sims))
  for(j in 2:ncol(sims)){
    cond <- cond & sims[,1] - sims[,j] - 1/n > 0
  }
  if(raw){sum(cond)}else{mean(cond)}
}

sequence_of_midpoints <- function(from = NULL, to = NULL, increments = 10){
  boundary_cuts <- seq(from, to, length = increments + 1)
  boundary_cuts[-(increments + 1)] + (boundary_cuts[2] - boundary_cuts[1])/2
}

ab_plurality_tie_for_first_from_sims <- function(sims, method = "density", n = 1000, merge = F, drop = F, increments = 10, window = .01, bw_divisor = 1, raw = F){

  delta <- sims[,1] - sims[,2]
  cond <- rep(T, nrow(sims))
  if(method == "density"){
    sum_12 <- sims[,1] + sims[,2]
    for(j in 3:ncol(sims)){
      cond <- cond & (sims[,j] < sum_12/2 - 1/n)
    }
  }else if(method == "naive_density"){
    for(j in 3:ncol(sims)){
      cond <- cond & (sims[,j] < sims[,1] - 1/n & sims[,j] < sims[,2] - 1/n)
    }
  }else if(method == "rectangular"){
    for(j in 3:ncol(sims)){
      # everyone else is at least 1/n behind either a or b.
      cond <- cond & (sims[,j] < sims[,1] - 1/n | sims[,j] < sims[,2] - 1/n)
    }
  }else{
    stop("Unknown method for plurality pivot prob estimation: ", method, "\n")
  }

  pivot_prob_from_delta_and_condition2(delta, cond, method = method, bw_divisor = bw_divisor, window = window, n = n, drop = drop, merge = merge, raw = raw, increments = increments)

}


winner_name <- function(v, cand_names){
  cand_names[which(v == max(v))]
}

runner_up_names <- function(v, cand_names, n){
  cand_names[which(max(v) > v & max(v) - v < 1/n)] %>%
    sort() %>%
    paste(collapse = "")
}

# this method takes a very long time and is plurality-specific -- not used
plurality_mc_accounting <- function(alpha = c(10, 9, 6), num_sims = 10000000, sims_increment = 500000, n = 1000, cand_names = NULL, raw = F){
  if(is.null(cand_names)){cand_names <- letters[1:length(alpha)]}
  increments <- ceiling(num_sims/sims_increment)
  cat(increments, " to go through.\n")
  all_types <- c()
  for(i in 1:increments){
    cat(".")
    sims <- gtools::rdirichlet(sims_increment, alpha)
    winner_names <- apply(sims, 1, winner_name, cand_names = cand_names)
    runner_up_names <- apply(sims, 1, runner_up_names, cand_names = cand_names, n = n)
   all_types <- c(all_types, paste0(winner_names, "_", runner_up_names))
  }
  cat("done.\n")
  if(raw){return(all_types)}
  out <- list()
  for(type in unique(all_types)){
    out[[type]] <- sum(all_types == type)/length(all_types)
  }
  out
}

#' @rdname standalone_monte_carlo_methods
#' @export
positional_event_probs_from_sims <- function(sims, window = .01, n = 1000, s = .5, cand_names = NULL, sep = "_", method = "density", merge = F, drop = F, increments = 10, bw_divisor = 1, skip_non_pivot_events = T, raw = F){

  if(ncol(sims) != 6){
    stop("sims must have 6 columns.")
  }

  if(is.null(cand_names)){
    cand_names <- names(sims)
    if(is.null(cand_names)){
      cand_names <- letters[1:3]
    }
  }

  # sims assumed to have columns abc, acb, bac, bca, cab, cba
  # each row is a simulation, i.e. a set of ballot shares

  # assemble the scores
  score_a <- sims[,1] + sims[,2] + s*(sims[,3] + sims[,5])
  score_b <- sims[,3] + sims[,4] + s*(sims[,1] + sims[,6])
  score_c <- sims[,5] + sims[,6] + s*(sims[,2] + sims[,4])

  # go through the events
  out <- list()
  ab_P <- rbind(c(1,1,s,0,1,1-s),
        c(0,0,1-s,1,0,s),
        0)

  ab <- ab_plurality_tie_for_first_from_sims(sims = cbind(score_a, score_b, score_c), method = method, n = n, merge = merge, drop = drop, increments = increments, window = window, bw_divisor = bw_divisor, raw = raw)

  out[[paste0(cand_names[1], sep, cand_names[2])]] <- list(
    integral = ab[1],
    P = cbind(ab_P, c(1,0,0))
  )

  out[[paste0(cand_names[2], sep, cand_names[1])]] <-
    list(integral = ab[2],
         P = cbind(ab_P[c(2,1,3), c(3,4,1,2,6,5)], c(0,1,0)))

  ac <- ab_plurality_tie_for_first_from_sims(sims = cbind(score_a, score_c, score_b), method = method, n = n, merge = merge, drop = drop, increments = increments, window = window, bw_divisor = bw_divisor, raw = raw)

  out[[paste0(cand_names[1], sep, cand_names[3])]] <- list(
    integral = ac[1],
    P = cbind(ab_P[c(1,3,2), c(2,1,5,6,3,4)], c(1,0,0))
  )

  out[[paste0(cand_names[3], sep, cand_names[1])]] <-
    list(integral = ac[2],
         P = cbind(ab_P[c(2,3,1), c(4,3, 6,5, 1,2)], c(0,0,1)))

  bc <- ab_plurality_tie_for_first_from_sims(sims = cbind(score_b, score_c, score_a), method = method, n = n, merge = merge, drop = drop, increments = increments, window = window, bw_divisor = bw_divisor, raw = raw)

  out[[paste0(cand_names[2], sep, cand_names[3])]] <-
    list(integral = bc[1],
         P = cbind(ab_P[c(3,1,2), c(5,6,2,1,4,3)], c(0,1,0)))

  out[[paste0(cand_names[3], sep, cand_names[2])]] <-
    list(integral = bc[2],
         P = cbind(ab_P[c(3,2,1), c(6,5,4,3,2,1)], c(0,0,1)))

  if(!skip_non_pivot_events){
    score_sims <- cbind(score_a, score_b, score_c)
    empty_P <- matrix(0, nrow = 3, ncol = 7)
    for(i in 1:3){
      this_P <- empty_P
      this_P[i,] <- 1
      out[[paste0(cand_names[i], "_")]] <- list(
        integral = a_wins_from_sims_prob(sims = cbind(score_sims[,i,drop=F], score_sims[,-i,drop=F]), n = n, raw = raw),
        P = this_P)
    }
  }

  out

}

### TODO: add drop method
#' @rdname standalone_monte_carlo_methods
#' @export
irv_event_probs_from_sims <- function(sims, window = .01, n = 1000, s = 0, cand_names = NULL, sep = "_", method = "density", merge = F, bw_divisor = 1, skip_non_pivot_events = T){

  if(ncol(sims) != 6){
    stop("sims must have 6 columns.")
  }

  if(is.null(cand_names)){
    cand_names <- names(sims)
    if(is.null(cand_names)){
      cand_names <- letters[1:3]
    }
  }

  # sims assumed to have columns abc, acb, bac, bca, cab, cba
  # each row is a simulation, i.e. a set of ballot shares

  # positional scores
  score_a <- sims[,1] + sims[,2] + s*(sims[,3] + sims[,5])
  score_b <- sims[,3] + sims[,4] + s*(sims[,1] + sims[,6])
  score_c <- sims[,5] + sims[,6] + s*(sims[,2] + sims[,4])

  # pairwise margins
  a_vs_b <- 2*(sims[,1] + sims[,2] + sims[,5] - .5)
  a_vs_c <- 2*(sims[,1] + sims[,2] + sims[,3] - .5)
  b_vs_c <- 2*(sims[,3] + sims[,4] + sims[,1] - .5)

  out <- list()

  # second round pivot events
  # a_b and b_a
  ab_P <- rbind(c(1,1,0,0,1,0,1), c(0,0,1,1,0,1,0),0)
  ab_pp <- irv_second_round_pivot_prob_ab(score_a, score_b, score_c, a_vs_b, s, method = method, n = n, merge = merge, window = window, bw_divisor = bw_divisor)
  out[[paste0(cand_names[1], sep, cand_names[2])]] <- list(
    integral = ab_pp[1],
    P = ab_P
  )
  ba_P <- ab_P
  ba_P[,7] <- c(0,1,0)
  out[[paste0(cand_names[2], sep, cand_names[1])]] <- list(
    integral = ab_pp[2],
    P = ba_P
  )

  # a_c and c_a
  ac_P <- rbind(c(1,1,1,0,0,0,1), 0, c(0,0,0,1,1,1,0))
  ac_pp <- irv_second_round_pivot_prob_ab(score_a, score_c, score_b, a_vs_c, s, method = method, n = n, merge = merge, window = window, bw_divisor = bw_divisor)
  out[[paste0(cand_names[1], sep, cand_names[3])]] <- list(
    integral = ac_pp[1],
    P = ac_P
  )
  ca_P <- ac_P
  ca_P[,7] <- c(0,0,1)
  out[[paste0(cand_names[3], sep, cand_names[1])]] <- list(
    integral = ac_pp[2],
    P = ca_P
  )

  # b_c and c_b
  bc_P <- rbind(0, c(1,0,1,1,0,0,1), c(0,1,0,0,1,1,0))
  bc_pp <- irv_second_round_pivot_prob_ab(score_b, score_c, score_a, b_vs_c, s, method = method, n = n, merge = merge, window = window, bw_divisor = bw_divisor)
  out[[paste0(cand_names[2], sep, cand_names[3])]] <- list(
    integral = bc_pp[1],
    P = bc_P
  )
  out[[paste0(cand_names[3], sep, cand_names[2])]] <- list(
    integral = bc_pp[2],
    P = bc_P
  )

  # first-round pivot events -- one set for each pair of candidates

  # a and b
  ab_pps <- irv_first_round_pivot_probs_ab(score_a, score_b, score_c, a_vs_c, b_vs_c, s, method = method, n = n, merge = merge, window = window, bw_divisor = bw_divisor)
  ab_P <- rbind(c(1,1,0,0,1,1,1), c(0,0,1,1,0,0,0), 0)
  ba_P <- rbind(c(1,1,0,0,0,0,0), c(0,0,1,1,1,1,1), 0)

  # a_b|ab
  out[[paste0(cand_names[1], sep, cand_names[2], "|",  cand_names[1], cand_names[2])]] <- list(
    integral = ab_pps[["i_j|ij"]][1],
    P = ab_P
  )
  # b_a|ba
  out[[paste0(cand_names[2], sep, cand_names[1], "|",  cand_names[2], cand_names[1])]] <- list(
    integral = ab_pps[["i_j|ij"]][2],
    P = ba_P
  )

  # a_b|ac
  out[[paste0(cand_names[1], sep, cand_names[2], "|",  cand_names[1], cand_names[3])]] <- list(
    integral = ab_pps[["i_j|ik"]][1],
    P = ab_P[c(1,3,2),]
  )
  # b_a|ca
  out[[paste0(cand_names[2], sep, cand_names[1], "|",  cand_names[3], cand_names[1])]] <- list(
    integral = ab_pps[["i_j|ik"]][2],
    P = ba_P[c(1,3,2),]
  )

  # a_b|cb
  out[[paste0(cand_names[1], sep, cand_names[2], "|",  cand_names[3], cand_names[2])]] <- list(
    integral = ab_pps[["i_j|kj"]][1],
    P = ab_P[c(3,2,1),]
  )
  # b_a|bc
  out[[paste0(cand_names[2], sep, cand_names[1], "|",  cand_names[2], cand_names[3])]] <- list(
    integral = ab_pps[["i_j|kj"]][2],
    P = ba_P[c(3,2,1),]
  )


  # a and c
  ac_pps <- irv_first_round_pivot_probs_ab(score_a, score_c, score_b, a_vs_b, -b_vs_c, s, method = method, n = n, merge = merge, window = window, bw_divisor = bw_divisor)
  ac_P <- rbind(c(1,1,1,1,0,0,1), 0, c(0,0,0,0,1,1,0))
  ca_P <- rbind(c(1,1,0,0,0,0,0), 0, c(0,0,1,1,1,1,1))

  # a_c|ac
  out[[paste0(cand_names[1], sep, cand_names[3], "|",  cand_names[1], cand_names[3])]] <- list(
    integral = ac_pps[["i_j|ij"]][1],
    P = ac_P
  )
  # b_a|ba
  out[[paste0(cand_names[3], sep, cand_names[1], "|",  cand_names[3], cand_names[1])]] <- list(
    integral = ac_pps[["i_j|ij"]][2],
    P = ca_P
  )

  # a_c|ab
  out[[paste0(cand_names[1], sep, cand_names[3], "|",  cand_names[1], cand_names[2])]] <- list(
    integral = ac_pps[["i_j|ik"]][1],
    P = ac_P[c(1,3,2),]
  )
  # c_a|ba
  out[[paste0(cand_names[3], sep, cand_names[1], "|",  cand_names[2], cand_names[1])]] <- list(
    integral = ac_pps[["i_j|ik"]][2],
    P = ca_P[c(1,3,2),]
  )

  # a_c|bc
  out[[paste0(cand_names[1], sep, cand_names[3], "|",  cand_names[2], cand_names[3])]] <- list(
    integral = ac_pps[["i_j|kj"]][1],
    P = ac_P[c(2,1,3),]
  )
  # c_a|cb
  out[[paste0(cand_names[3], sep, cand_names[1], "|",  cand_names[3], cand_names[2])]] <- list(
    integral = ac_pps[["i_j|kj"]][2],
    P = ca_P[c(2,1,3),]
  )


  # b and c
  bc_pps <- irv_first_round_pivot_probs_ab(score_b, score_c, score_a, -a_vs_b, -a_vs_c, s, method = method, n = n, merge = merge, window = window, bw_divisor = bw_divisor)
  bc_P <- rbind(0, c(1,1,1,1,0,0,1), c(0,0,0,0,1,1,0))
  cb_P <- rbind(0, c(0,0,1,1,0,0,0), c(1,1,0,0,1,1,1))

  # b_c|bc
  out[[paste0(cand_names[2], sep, cand_names[3], "|",  cand_names[2], cand_names[3])]] <- list(
    integral = bc_pps[["i_j|ij"]][1],
    P = bc_P
  )
  # c_b|cb
  out[[paste0(cand_names[3], sep, cand_names[2], "|",  cand_names[3], cand_names[2])]] <- list(
    integral = bc_pps[["i_j|ij"]][2],
    P = cb_P
  )

  # b_c|ba
  out[[paste0(cand_names[2], sep, cand_names[3], "|",  cand_names[2], cand_names[1])]] <- list(
    integral = bc_pps[["i_j|ik"]][1],
    P = bc_P[c(3,2,1),]
  )
  # c_b|ab
  out[[paste0(cand_names[3], sep, cand_names[2], "|",  cand_names[1], cand_names[2])]] <- list(
    integral = bc_pps[["i_j|ik"]][2],
    P = cb_P[c(3,2,1),]
  )

  # b_c|ac
  out[[paste0(cand_names[2], sep, cand_names[3], "|",  cand_names[1], cand_names[3])]] <- list(
    integral = bc_pps[["i_j|kj"]][1],
    P = bc_P[c(2,1,3),]
  )
  # c_b|ca
  out[[paste0(cand_names[3], sep, cand_names[2], "|",  cand_names[3], cand_names[1])]] <- list(
    integral = bc_pps[["i_j|kj"]][2],
    P = cb_P[c(2,1,3),]
  )

  if(!skip_non_pivot_events){
    out[[paste0(cand_names[1],"__", cand_names[2])]] <- list(
      integral = mean(score_a - score_c > 1/n & score_b - score_c > 1/n & a_vs_b > 1/n),
      P = rbind(rep(1, 7), 0, 0)
    )

    out[[paste0(cand_names[2],"__", cand_names[1])]] <- list(
      integral = mean(score_a - score_c > 1/n & score_b - score_c > 1/n & a_vs_b < -1/n),
      P = rbind(0, rep(1, 7), 0)
    )

    out[[paste0(cand_names[1],"__", cand_names[3])]] <- list(
      integral = mean(score_a - score_b > 1/n & score_c - score_b > 1/n & a_vs_c > 1/n),
      P = rbind(rep(1, 7), 0, 0)
    )

    out[[paste0(cand_names[3],"__", cand_names[1])]] <- list(
      integral = mean(score_a - score_b > 1/n & score_c - score_b > 1/n & a_vs_c < -1/n),
      P = rbind(0, 0, rep(1, 7))
    )

    out[[paste0(cand_names[2],"__", cand_names[3])]] <- list(
      integral = mean(score_b - score_a > 1/n & score_c - score_a > 1/n & b_vs_c > 1/n),
      P = rbind(0, rep(1, 7), 0)
    )

    out[[paste0(cand_names[3],"__", cand_names[2])]] <- list(
      integral = mean(score_b - score_a > 1/n & score_c - score_a > 1/n & b_vs_c < -1/n),
      P = rbind(0, 0, rep(1, 7))
    )

  }

  out
}

irv_second_round_pivot_prob_ab <- function(score_a, score_b, score_c, a_vs_b, s, method = "density", n = 1000, merge = F, window = .01, bw_divisor = 1){
  if(method %in% c("density", "naive_density")){
    if(merge){
      limits <- c(0,0)
    }else{
      limits <- c(1, -1)/(2*n)
    }
    if(method == "density"){
      cond_1 <- score_a - (a_vs_b/2)*(1 - s/2) - score_c > 1/n
      cond_2 <- score_b + (a_vs_b/2)*(1 - s/2) - score_c > 1/n
    }else if(method == "naive_density"){
      cond_1 <- score_a - score_c > 1/n
      cond_2 <- score_b - score_c > 1/n
    }
    cond <- cond_1 & cond_2
    the_density <- density_estimate(x = a_vs_b[cond], eval.points = limits, bw_divisor = bw_divisor)
    mean(cond)*the_density*(1/n)
  }else if(method == "rectangular"){
    pp <- mean(score_a - score_c > 1/n & score_b - score_c > 1/n & abs(a_vs_b) < window/2)/(n*window)
    c(pp, pp) # merge by default
  }else{
    stop("Unknown method for positional pivot prob estimation: ", method, "\n")
  }
}

irv_first_round_pivot_probs_ab <- function(score_a, score_b, score_c, a_vs_c, b_vs_c, s, method = "density", n = 1000, merge = F, window = .01, bw_divisor = 1){
  out <- list()
  if(method %in% c("density", "naive_density")){
    if(merge){
      limits <- c(0,0)
    }else{
      limits <- c(1, -1)/(2*n)
    }
    delta <- score_a - score_b
    if(method == "density"){
      cond_1 <- score_c - score_a > 1/n - (delta/2)*(1 - s*(1-s)/2)
      cond_2 <- score_c - score_b > 1/n + (delta/2)*(1 - s*(1-s)/2)
      diff_4 <- a_vs_c - (delta/2)*(1 + s)
      diff_5 <- b_vs_c + (delta/2)*(1 + s)
      a_b_ab_cond <- cond_1 & cond_2 & diff_4 > 1/n & diff_5 > 1/n
      a_b_cb_cond <- cond_1 & cond_2 & diff_4 < -1/n & diff_5 > 1/n
      a_b_ac_cond <- cond_1 & cond_2 & diff_4 > 1/n & diff_5 < -1/n
    }else if(method == "naive_density"){
      c_first_cond <- score_c - score_a > 1/n & score_c - score_b > 1/n
      a_b_ab_cond <- c_first_cond & a_vs_c > 1/n & b_vs_c > 1/n
      a_b_cb_cond <- c_first_cond & a_vs_c < -1/n & b_vs_c > 1/n
      a_b_ac_cond <- c_first_cond & a_vs_c > 1/n & b_vs_c < -1/n
    }
    out[["i_j|ij"]] <- mean(a_b_ab_cond)*density_estimate(x = delta[a_b_ab_cond], eval.points = limits, bw_divisor = bw_divisor)*(1/n)
    out[["i_j|kj"]] <- mean(a_b_cb_cond)*density_estimate(x = delta[a_b_cb_cond], eval.points = limits, bw_divisor = bw_divisor)*(1/n)
    out[["i_j|ik"]] <- mean(a_b_ac_cond)*density_estimate(x = delta[a_b_ac_cond], eval.points = limits, bw_divisor = bw_divisor)*(1/n)
  }else if(method == "rectangular"){
    c_first_cond <- score_c - score_a > 1/n & score_c - score_b > 1/n
    a_b_ab_cond <- c_first_cond & a_vs_c > 1/n & b_vs_c > 1/n
    a_b_cb_cond <- c_first_cond & a_vs_c < -1/n & b_vs_c > 1/n
    a_b_ac_cond <- c_first_cond & a_vs_c > 1/n & b_vs_c < -1/n
    ab_tie_cond <- abs(score_a - score_b) < window/2
    out[["i_j|ij"]] <- rep(mean(a_b_ab_cond & ab_tie_cond)/(n*window), 2)
    out[["i_j|kj"]] <- rep(mean(a_b_cb_cond & ab_tie_cond)/(n*window), 2)
    out[["i_j|ik"]] <- rep(mean(a_b_ac_cond & ab_tie_cond)/(n*window), 2)
  }else{
    stop("Unknown method for positional pivot prob estimation: ", method, "\n")
  }
  out
}

# TODO: add drop method
#' @rdname standalone_monte_carlo_methods
#' @export
condorcet_event_probs_from_sims <- function(sims, n = 1000, window = .01, cand_names = NULL, sep = "_", kemeny = T, method = "density", merge = F, bw_divisor = 1, s = NULL, skip_non_pivot_events = T){

  if(is.null(cand_names)){
    if(kemeny & !is.null(cand_names) & (cand_names %>% sort() %>% paste(collapse = "") != "abc")){
      warning("Kemeny-Young pivot event names will use a, b, c, not the supplied candidate names.")
    }
    cand_names <- letters[1:3]
  }

  if(ncol(sims) != 6){
    stop("sims must have 6 columns.")
  }

  if(method == "density"){method <- "naive_density"; bw_divisor = 1} # formerly undersmoothing this (bw_divisor = 2). trying without  # implementation of density is not yet correct, but naive_density works pretty well.

  # pairwise tallies: what each gets against each other
  # note this is different from meaning above
  a_vs_b <- sims[,1] + sims[,2] + sims[,5]
  a_vs_c <- sims[,1] + sims[,2] + sims[,3]
  b_vs_a <- sims[,3] + sims[,4] + sims[,6]
  b_vs_c <- sims[,3] + sims[,4] + sims[,1]
  c_vs_a <- sims[,5] + sims[,6] + sims[,4]
  c_vs_b <- sims[,5] + sims[,6] + sims[,2]


  out <- list()

  # Condorcet winner event a_b
  delta <- a_vs_b - b_vs_a
  delta_adjust <- ifelse(method == "density", delta, 0)
  cond_1 <- a_vs_c - c_vs_a > 1/n + delta_adjust/2
  cond_2 <- b_vs_c - c_vs_b > 1/n - delta_adjust/2
  pps <- pivot_prob_from_delta_and_condition2(delta, cond_1 & cond_2, method = method, bw_divisor = bw_divisor, window = window, n = n, drop = T, merge = merge)
  a_b_P <- rbind(c(1,1,0,0,1,0), c(0,0,1,1,0,1),0)
  out[[paste0(cand_names[1], sep, cand_names[2])]] <- list(
    integral = pps[1],
    P = cbind(a_b_P, c(1,0,0))
  )
  out[[paste0(cand_names[2], sep, cand_names[1])]] <- list(
    integral = pps[2],
    P = cbind(a_b_P, c(0,1,0))
  )

  # Condorcet winner event a_c
  delta <- a_vs_c - c_vs_a
  delta_adjust <- ifelse(method == "density", delta, 0)
  cond_1 <- a_vs_b - b_vs_a > 1/n + delta_adjust/2
  cond_2 <- c_vs_b - b_vs_c > 1/n - delta_adjust/2
  pps <- pivot_prob_from_delta_and_condition2(delta, cond_1 & cond_2, method = method, bw_divisor = bw_divisor, window = window, n = n, drop = T, merge = merge)
  a_c_P <- rbind(c(1,1,1,0,0,0), 0, c(0,0,0,1,1,1))
  out[[paste0(cand_names[1], sep, cand_names[3])]] <- list(
    integral = pps[1],
    P = cbind(a_c_P, c(1,0,0))
  )
  out[[paste0(cand_names[3], sep, cand_names[1])]] <- list(
    integral = pps[2],
    P = cbind(a_c_P, c(0,0,1))
  )

  # Condorcet winner event b_c
  delta <- b_vs_c - c_vs_b
  delta_adjust <- ifelse(method == "density", delta, 0)
  cond_1 <- b_vs_a - a_vs_b > 1/n + delta_adjust/2
  cond_2 <- c_vs_a - a_vs_c > 1/n - delta_adjust/2
  pps <- pivot_prob_from_delta_and_condition2(delta, cond_1 & cond_2, method = method, bw_divisor = bw_divisor, window = window, n = n, drop = T, merge = merge)
  b_c_P <- rbind(0, c(1,0,1,1,0,0), c(0,1,0,0,1,1))
  out[[paste0(cand_names[2], sep, cand_names[3])]] <- list(
    integral = pps[1],
    P = cbind(b_c_P, c(0,1,0))
  )
  out[[paste0(cand_names[3], sep, cand_names[2])]] <- list(
    integral = pps[2],
    P = cbind(b_c_P, c(0,0,1))
  )

  if(kemeny){

    forward_cycle <- a_vs_b - b_vs_a > 1/n & b_vs_c - c_vs_b > 1/n & c_vs_a - a_vs_c > 1/n
    reverse_cycle <- a_vs_b - b_vs_a < 1/n & b_vs_c - c_vs_b < 1/n & c_vs_a - a_vs_c < 1/n

    # ab_forward
    delta <- a_vs_c - b_vs_a # gap between what a gets in its losing race and what b gets in its losing race
    delta_adjust <- ifelse(method == "density", delta, 0)
    cond_1 <- a_vs_c - c_vs_b > 1/n + delta_adjust/2
    cond_2 <- b_vs_a - c_vs_b > 1/n - delta_adjust/2
    pps <- pivot_prob_from_delta_and_condition2(delta, forward_cycle & cond_1 & cond_2, method = method, bw_divisor = bw_divisor, window = window, n = n, drop = T, merge = merge)
    out[["ac_ba|abca"]] <- list(
      integral = pps[1],
      P = cbind(c(1,0,0), c(1,0,0), c(1,0,0), c(0,1,0), c(1,0,0), c(0,1,0), c(1,0,0))
    )
    out[["ba_ac|abca"]] <- list(
      integral = pps[2],
      P = cbind(c(1,0,0), c(1,0,0), c(0,1,0), c(0,1,0), c(0,1,0), c(0,1,0), c(0,1,0))
    )

    # ab_reverse
    delta <- a_vs_b - b_vs_c
    delta_adjust <- ifelse(method == "density", delta, 0)
    cond_1 <- a_vs_b - c_vs_a > 1/n + delta_adjust/2
    cond_2 <- b_vs_c - c_vs_a > 1/n - delta_adjust/2
    pps <- pivot_prob_from_delta_and_condition2(delta, reverse_cycle & cond_1 & cond_2, method = method, bw_divisor = bw_divisor, window = window, n = n, drop = T, merge = merge)
    out[["ab_bc|bacb"]] <- list(
      integral = pps[1],
      P = cbind(c(1,0,0), c(1,0,0), c(0,1,0), c(0,1,0), c(1,0,0), c(1,0,0), c(1,0,0))
    )
    out[["bc_ab|bacb"]] <- list(
      integral = pps[2],
      P = cbind(c(0,1,0), c(1,0,0), c(0,1,0), c(0,1,0), c(1,0,0), c(0,1,0), c(0,1,0))
    )

    # ac_forward
    delta <- (a_vs_c - c_vs_b) # gap between what a gets in its losing race and what c gets in its losing race
    delta_adjust <- ifelse(method == "density", delta, 0)
    # this one is "altered" because I realized the conditions should incorporate delta, but I abandoned at this point.
    altered_forward_cycle <- a_vs_b - b_vs_a > 1/n & b_vs_c - c_vs_b > 1/n + delta_adjust & c_vs_a - a_vs_c > 1/n - delta_adjust
    cond_1 <- a_vs_c - b_vs_a > 1/n + delta_adjust/2
    cond_2 <- c_vs_b - b_vs_a > 1/n - delta_adjust/2
    pps <- pivot_prob_from_delta_and_condition2(delta, altered_forward_cycle & cond_1 & cond_2, method = method, bw_divisor = bw_divisor, window = window, n = n, drop = T, merge = merge)
    out[["ac_cb|cabc"]] <- list(
      integral = pps[1],
      P = cbind(c(1,0,0), c(1,0,0), c(1,0,0), c(1,0,0), c(0,0,1), c(0,0,1), c(1,0,0))
    )
    out[["cb_ac|cabc"]] <- list(
      integral = pps[2],
      P = cbind(c(1,0,0), c(0,0,1), c(1,0,0), c(0,0,1), c(0,0,1), c(0,0,1), c(0,0,1))
    )

    # ac_reverse
    delta <- a_vs_b - c_vs_a
    delta_adjust <- ifelse(method == "density", delta, 0)
    cond_1 <- a_vs_b - b_vs_c > 1/n + delta_adjust/2
    cond_2 <- c_vs_a - b_vs_c > 1/n - delta_adjust/2
    pps <- pivot_prob_from_delta_and_condition2(delta, reverse_cycle & cond_1 & cond_2, method = method, bw_divisor = bw_divisor, window = window, n = n, drop = T, merge = merge)
    out[["ab_ca|acba"]] <- list(
      integral = pps[1],
      P = cbind(c(1,0,0), c(1,0,0), c(1,0,0), c(0,0,1), c(1,0,0), c(0,0,1), c(1,0,0))
    )
    out[["ca_ab|acba"]] <- list(
      integral = pps[2],
      P = cbind(c(1,0,0), c(1,0,0), c(0,0,1), c(0,0,1), c(0,0,1), c(0,0,1), c(0,0,1))
    )

    # bc_forward
    delta <- b_vs_a  - c_vs_b
    delta_adjust <- ifelse(method == "density", delta, 0)
    cond_1 <- b_vs_a - a_vs_c > 1/n + delta_adjust/2
    cond_2 <- c_vs_b - a_vs_c > 1/n - delta_adjust/2
    pps <- pivot_prob_from_delta_and_condition2(delta, forward_cycle & cond_1 & cond_2, method = method, bw_divisor = bw_divisor, window = window, n = n, drop = T, merge = merge)
    out[["ba_cb|bcab"]] <- list(
      integral = pps[1],
      P = cbind(c(0,1,0), c(0,0,1), c(0,1,0), c(0,1,0), c(0,0,1), c(0,1,0), c(0,1,0))
    )
    out[["cb_ba|bcab"]] <- list(
      integral = pps[2],
      P = cbind(c(0,0,1), c(0,0,1), c(0,1,0), c(0,1,0), c(0,0,1), c(0,0,1), c(0,0,1))
    )

    # bc_reverse
    delta <- b_vs_c - c_vs_a
    delta_adjust <- ifelse(method == "density", delta, 0)
    cond_1 <- b_vs_c - a_vs_b > 1/n + delta_adjust/2
    cond_2 <- c_vs_a - a_vs_b > 1/n - delta_adjust/2
    pps <- pivot_prob_from_delta_and_condition2(delta, reverse_cycle & cond_1 & cond_2, method = method, bw_divisor = bw_divisor, window = window, n = n, drop = T, merge = merge)
#     pps <- pivot_prob_from_delta_and_condition(delta, reverse_cycle & cond_1 & cond_2, method = method, bw_divisor = bw_divisor, window = window, eval.points = limits, n = n)
    out[["bc_ca|cbac"]] <- list(
      integral = pps[1],
      P = cbind(c(0,1,0), c(0,1,0), c(0,1,0), c(0,1,0), c(0,0,1), c(0,0,1), c(0,1,0))
    )
    out[["ca_bc|cbac"]] <- list(
      integral = pps[2],
      P = cbind(c(0,1,0), c(0,0,1), c(0,1,0), c(0,0,1), c(0,0,1), c(0,0,1), c(0,0,1))
    )
  }

  if(!skip_non_pivot_events){

    # a is condorcet winner
    out[[paste0(cand_names[[1]], "_")]] <- list(
      integral = mean(a_vs_b - b_vs_a > 1/n & a_vs_c - c_vs_a > 1/n),
      P = rbind(rep(1, 7), 0, 0)
    )

    # b is condorcet winner
    out[[paste0(cand_names[[2]], "_")]] <- list(
      integral = mean(b_vs_a - a_vs_b > 1/n & b_vs_c - c_vs_b > 1/n),
      P = rbind(0, rep(1, 7), 0)
    )

    # c is condorcet winner
    out[[paste0(cand_names[[3]], "_")]] <- list(
      integral = mean(c_vs_a - a_vs_c > 1/n & c_vs_b - b_vs_c > 1/n),
      P = rbind(0, 0, rep(1, 7))
    )

    forward_cycle <- a_vs_b - b_vs_a > 1/n & b_vs_c - c_vs_b > 1/n & c_vs_a - a_vs_c > 1/n
    reverse_cycle <- b_vs_a - a_vs_b > 1/n & c_vs_b - b_vs_c > 1/n & a_vs_c - c_vs_a > 1/n
    # a wins in cycle
    out[[paste0(cand_names[[1]], "_CYCLE")]] <- list(
      integral = mean((forward_cycle & a_vs_c > b_vs_a & a_vs_c > c_vs_b) | (reverse_cycle & a_vs_b > b_vs_c & a_vs_b > c_vs_a)),
      P = rbind(rep(1, 7), 0, 0)
    )
    # b wins in cycle
    out[[paste0(cand_names[[2]], "_CYCLE")]] <- list(
      integral = mean((forward_cycle & b_vs_a > c_vs_b & b_vs_a > a_vs_c) | (reverse_cycle & b_vs_c > a_vs_b & b_vs_c > c_vs_a)),
      P = rbind(0, rep(1, 7), 0)
    )
    # c wins in cycle
    out[[paste0(cand_names[[3]], "_CYCLE")]] <- list(
      integral = mean((forward_cycle & b_vs_a > a_vs_c & b_vs_a > c_vs_b) | (reverse_cycle & c_vs_a > a_vs_b & c_vs_a > b_vs_c)),
      P = rbind(0, 0, rep(1, 7))
    )

  }


  out

}
