##################################################
## Project: Strategic Voting in AV
## Script purpose: Functions for Calculating
##    Optimal Ballots
## Date: 04/08/2018
## Author:
##################################################

## I could just morph the AV and plurality functions into one...
##  after all, they mostly use the same components.

library(pbapply)

id_pref_order <- function(util){
  # Function that gets the "number" of pref. ordering type
  stopifnot(length(util) == 3)
  max <- which(util == max(util))
  min <- which(util == min(util))
  
  if (max == 1){
    if (min == 3){
      out <- 1
    }
    if (min == 2){
      out <- 2
    }
  }
  if (max == 2){
    if (min == 3){
      out <- 3
    }
    if (min == 1){
      out <- 4
    }
  }
  if (max == 3){
    if (min == 2){
      out <- 5
    }
    if (min == 1){
      out <- 6
    }
  }
  return(out)
}


# ----------------
# Alternative Vote
# ----------------

# Dataframe that has the probabilities of winning based on each ballot
A <- c(1, 0, 0)
B <- c(0, 1, 0)
C <- c(0, 0, 1)
AB <- c(0.5, 0.5, 0)
AC <- c(0.5, 0, 0.5)
BC <- c(0, 0.5, 0.5)
w_av <- as.data.frame(matrix(data = c(A, A, B, B, A, B,                    # AB
                                      A, A, A, C, C, C,                    # AC
                                      B, C, B, B, C, C,                    # BC
                                      A, A, B, B, AB, AB,                  # AB.AB
                                      C, C, B, B, BC, BC,                  # AB.CB
                                      A, A, C, C, AC, AC,                  # AB.AC
                                      A, A, AB, AB, B, B,                  # AC.AB
                                      B, B, BC, BC, C, C,                  # AC.BC
                                      A, A, AC, AC, C, C,                  # AC.AC
                                      AC, AC, A, A, C, C,                  # BC.AC
                                      AB, AB, B, B, A, A,                  # BC.BA
                                      BC, BC, B, B, C, C),                 # BC.BC
                          byrow = TRUE, ncol = 6*3, nrow = 12))

# Dataframe with "default" voter types
default_util <- data.frame(A = c(1, 1, 0.5, 0, 0.5, 0),
                           B = c(0.5, 0, 1, 1, 0, 0.5),
                           C = c(0, 0.5, 0, 0.5, 1, 1))


# Function to calculate optimal vote given a list of pivotal probabilities and
# a matrix of utilities. Does NOT resolve ties.
# sincere: vector of position of sincere preference
av_opt_vote <- function(p_list, utility_df, sincere){
  p_list <- c(p_list$AB, p_list$AC, p_list$BC,
              p_list$AB.AB, p_list$AB.CB, p_list$AB.AC,
              p_list$AC.AB, p_list$AC.BC, p_list$AC.AC,
              p_list$BC.AC, p_list$BC.BA, p_list$BC.BC)
  
  w_test <- w_av[, rep(1:18, nrow(utility_df))]
  u_test <- utility_df[, rep(1:3, 6)]
  u_test <- as.vector(t(u_test))
  
  test <- matrix((mapply(`*`,w_test,u_test,SIMPLIFY=TRUE)), nrow = 12)
  
  cond_utils_sum <- test[, c(T, F, F)] + test[, c(F, T, F)] +
    test[, c(F, F, T)]
  cond_utils_sum_event <- cond_utils_sum * p_list
  vote_utils <- colSums(cond_utils_sum_event)
  
  opt <- as.numeric()
  diff <- as.numeric()
  for (i in 1:(nrow(utility_df))){
    eval <- vote_utils[(i * 6 - 5):(i*6)]
    
    eval_max <- which(eval == max(eval))
    if (length(eval_max) > 1) {eval_max <- 99}
    diff[i] <- eval[eval_max] - eval[sincere[i]]
    ## add the following condition: if eval == sincere,
    ## take second-highest eval and return the difference.
    # print(eval[eval_max], eval[sincere[i]])
    if ((eval[eval_max] == eval[sincere[i]]) == TRUE & eval_max != 99){
      # grab second highest value
      eval_sec_value <- max(eval[-eval_max])
      diff[i] <- eval_sec_value - eval[sincere[i]]
    }
    opt[i] <- eval_max
  }
  return(list(opt = opt, diff = diff))
}

# Function to calculate optimal vote given a list of pivotal probabilities and
# a matrix of utilities. DOES resolve ties.
av_opt_vote_no_ties <- function(p_list, utility_df, sincere){
  p_list <- c(p_list$AB, p_list$AC, p_list$BC,
              p_list$AB.AB, p_list$AB.CB, p_list$AB.AC,
              p_list$AC.AB, p_list$AC.BC, p_list$AC.AC,
              p_list$BC.AC, p_list$BC.BA, p_list$BC.BC)
  
  w_test <- w_av[, rep(1:18, nrow(utility_df))]
  u_test <- utility_df[, rep(1:3, 6)]
  u_test <- as.vector(t(u_test))
  
  test <- matrix((mapply(`*`, w_test, u_test, SIMPLIFY=TRUE)), nrow = 12)
  
  cond_utils_sum <- test[, c(T, F, F)] + test[, c(F, T, F)] +
    test[, c(F, F, T)]
  cond_utils_sum_event <- cond_utils_sum * mpfr(p_list, 256)
  vote_utils <- cond_utils_sum_event[1, ] +
    cond_utils_sum_event[2, ] + cond_utils_sum_event[3, ] + cond_utils_sum_event[4, ] +
    cond_utils_sum_event[5, ] + cond_utils_sum_event[6, ] + cond_utils_sum_event[7, ] +
    cond_utils_sum_event[8, ] + cond_utils_sum_event[9, ] + cond_utils_sum_event[10, ] +
    cond_utils_sum_event[11, ] + cond_utils_sum_event[12, ]
  
  opt <- as.numeric()
  diff <- mpfr(NA)
  for (i in 1:(nrow(utility_df))){
    eval <- vote_utils[(i * 6 - 5):(i*6)]

    eval_max <- which(eval == max(eval))
    if (length(eval_max) > 1) {eval_max <- 99}
    diff[i] <- eval[eval_max] - eval[sincere[i]]
    ## add the following condition: if eval == sincere,
    ## take second-highest eval and return the difference.
    if ((eval[eval_max] == eval[sincere[i]]) == TRUE){
      # grab second highest value
      eval_sec_value <- max(eval[-eval_max])
      diff[i] <- eval_sec_value - eval[sincere[i]]
    }
    opt[i] <- eval_max
  }
  diff <- as.numeric(diff)
  return(list(opt = opt, diff = diff))
}

# Takes a list of pivotal probability lists and a utility dataframe;
# returns matrix of optimal votes. Avoids ties.
# Wraps av_opt_vote and av_opt_vote_no_ties together for multiple isntances of p_list.
av_opt_vote_mat <- function(sim_list, utility_df){
  sincere <- apply(aes_utils[, 1:3], 1, id_pref_order)
  out <- pblapply(sim_list, function(x) av_opt_vote(x, utility_df, sincere))
  
  out_opt <- do.call(rbind, out)[, 1]

  out_opt <- do.call(rbind, out_opt)
  out_diff <- do.call(rbind, out)[, 2]
  out_diff <- do.call(rbind, out_diff)
  
  ties <- unique(which(out_opt == 99, arr.ind = TRUE)[, 1])
  
  if (length(ties) > 0){
    out_ties <- lapply(sim_list[ties], 
                       function(x) av_opt_vote_no_ties(x, utility_df, sincere))
    out_ties_opt <- do.call(rbind, out_ties)[, 1]
    out_ties_opt <- do.call(rbind, out_ties_opt)
    out_ties_diff <- do.call(rbind, out_ties)[, 2]
    out_ties_diff <- do.call(rbind, out_ties_diff)
    out_opt[ties, ] <- out_ties_opt
    out_diff[ties, ] <- as.numeric(out_ties_diff)
  }
  return(list(opt = out_opt, diff = out_diff))
}

# Test to then insert into the new function
# s_85 <- av_opt_vote_mat(s_p_list[[4]], aes_utils[, 1:3])
# ties <- unique(which(s_85$opt == 99, arr.ind = T)[ ,1])
# 
# out_ties <- lapply(s_p_list[[4]][ties],
#                    function(x) av_opt_vote_no_ties(x, aes_utils[, 1:3], sincere))  
# out_ties_opt <- do.call(rbind, out_ties)[ ,1]
# out_ties_opt <- do.call(rbind, out_ties_opt)  
# out_ties_diff <- do.call(rbind, out_ties)[, 2]
# out_ties_diff <- do.call(rbind, out_ties_diff)
# s_85$opt[ties, ] <- out_ties_opt
# s_85_opt <- s_85$opt
# s_85_diff <- s_85$diff
# s_85_diff[ties, ] <- as.numeric(out_ties_diff) # this works although the problem is the magnitude..

# Todo:
# Re-insert into function
# By type of voter, show proportion of strategic voting for different values of \epsilon.

# ---------
# Plurality
# ---------

w_plurality <- as.data.frame(matrix(data = c(1, 0, 0, 0, 1, 0, .5, .5, 0,
                               1, 0, 0, .5, 0, .5, 0, 0, 1,
                               0, .5, .5, 0, 1, 0, 0, 0, 1),
                      byrow = T, nrow = 3))

plur_opt_vote <- function(p_list, utility_df, sincere){
  # Function that calculates the optimal vote given a list of piv. probabilities.
  # Input: 
  #     p_list: list of pivotal probabilities for plurality
  #     utility_df: matrix of utilities for (A, B, C); each row is a respondent
  # Output: v
  p_list <- c(p_list$AB, p_list$AC, p_list$BC)
  w_test <- w_plurality[ , rep(1:9, nrow(utility_df))]
  u_test <- utility_df[, rep(1:3, 3)]
  u_test <- as.vector(t(u_test))
  
  test <- matrix((mapply(`*`,w_test, u_test, SIMPLIFY=TRUE)), nrow = 3)
  cond_utils_sum <- test[, c(T, F, F)] + test[, c(F, T, F)] +
    test[, c(F, F, T)]
  cond_utils_sum_event <- cond_utils_sum * p_list
  
  vote_utils <- colSums(cond_utils_sum_event)
  opt <- as.numeric()
  diff <- as.numeric()
  for (i in 1:(nrow(utility_df))){
    eval <- vote_utils[(i * 3 - 2):(i*3)]
    
    eval_max <- which(eval == max(eval))
    if (length(eval_max) > 1) {eval_max <- 99}
    diff[i] <- eval[eval_max] - eval[sincere[i]]
    if ((eval[eval_max] == eval[sincere[i]]) == TRUE & eval_max != 99){
      # grab second highest value
      eval_sec_value <- max(eval[-eval_max])
      diff[i] <- eval_sec_value - eval[sincere[i]]
    }
    opt[i] <- eval_max
  }
  return(list(opt = opt, diff = diff))
}

plur_opt_vote_mat <- function(sim_list, utility_df){
  # Function that runs plur_opt_vote for a number of situations
  
  sincere <- apply(aes_utils[, 1:3], 1, id_pref_order)
  # Recoding necessary because only three voting options:
    sincere[sincere %in% c(1, 2)] <- 1
    sincere[sincere %in% c(3, 4)] <- 2
    sincere[sincere %in% c(5, 6)] <- 3 
  out <- pblapply(sim_list, function(x) plur_opt_vote(x, utility_df, sincere))
  out_opt <- do.call(rbind, out)[, 1]
  out_opt <- do.call(rbind, out_opt)
  
  out_diff <- do.call(rbind, out)[, 2]
  out_diff <- do.call(rbind, out_diff)
  
  return(list(opt = out_opt, diff = out_diff))
}

# -----------
# General use
# -----------

list_to_summary <- function(mat, levels, labels, type){
  ## takes list of optimal votes and aggregates into 6 x 6 matrix, for multiple
  ## levels of s
  ##    Inputs
  ##        mat:    matrix with optimal votes
  ##        levels: vector of unique optimal votes (1, 2, 3, ...)
  ##        labels: vector of unique voter types ("GRN.COAL.LAB",...)
  ##        type:   vector of voter types in utility dataframe
  
  # Consider whether it's useful to get levels / labels from "mat", rather than specifying
  # (potentially less scope for errors)
  
  # summarise by voter
  count_by_voter <- vapply(levels, function(x) colSums(mat == x), numeric(ncol(mat)))

  # summarise by voter type
  agg_summary <- t(vapply(labels, 
             function(x) colSums(count_by_voter[type == x, ]),
             numeric(ncol(count_by_voter))))
  
  return(agg_summary)
}

# Next step: convert summary into first/second/third voting option.

replace_epsilon <- function(epsilon, x, sincere = sincere, type = aes_utils$type){
  # For a given set of \epsilons, replace optimal vote with sincere vote if U_opt - U_sin < \epsilon.
  
  sincere_mat <- matrix(rep(sincere, 10000), ncol = 736, byrow = T)
  
  costly_list <- list()
  for (e in epsilon){
    location <- x$diff < e
    costly_mat <- x$opt
    costly_mat[location] <- sincere_mat[location]
    name <- as.character(e)
    costly_list[[name]] <- costly_mat
  }
  costly_summary <- lapply(costly_list,
                           function (x) 
                             list_to_summary(x, 1:6, labels, aes_utils$type))
  return(costly_summary)
}

reorder_summary_av<- function(x){
  stopifnot(dim(x) == c(6, 6))
  out <- matrix(c(x[1, 1], x[1, 3], x[1, 5],
                  x[2, 2], x[2, 5], x[2, 3],
                  x[3, 3], x[3, 1], x[3, 6],
                  x[4, 4], x[4, 6], x[4, 1],
                  x[5, 5], x[5, 2], x[5, 4],
                  x[6, 6], x[6, 4], x[6, 2]), 
                ncol = 3, byrow = T)
  return(out)
}

reorder_summary_plur <- function(x){
  stopifnot(dim(x) == c(6, 6))
  out <- matrix(c(x[1, 1], x[1, 2], x[1, 3],
                  x[2, 1], x[2, 3], x[2, 2],
                  x[3, 2], x[3, 1], x[3, 3],
                  x[4, 2], x[4, 3], x[4, 1],
                  x[5, 3], x[5, 1], x[5, 2],
                  x[6, 3], x[6, 2], x[6, 1]), 
                ncol = 3, byrow = T)
  return(out)
}
