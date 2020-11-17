library(tidyverse)
library(devtools)
library(pivotprobs)
library(gtools)


### --- FUNCTIONS ---
ballot.mat.from.eu.mat = function(eu.mat, break.ties.with.sincerity = F, sincere.vote.mat = NULL){
  # this does not need weights. 
  
  if(break.ties.with.sincerity){
    if(is.null(sincere.vote.mat)){
      cat("you must pass a sincere.vote.mat!\n")
    }
    eu.mat <- eu.mat + sincere.vote.mat*(10^(-10))

   }

  
  max.eus = apply(eu.mat, 1, max, na.rm = T)
  ballot.mat = matrix(NA, ncol = ncol(eu.mat), nrow = nrow(eu.mat))
  for(j in 1:ncol(eu.mat)){
    ballot.mat[,j] = as.integer(eu.mat[,j] == max.eus)
  }

  if(break.ties.with.sincerity){
    if(is.null(sincere.vote.mat)){
      cat("you must pass a sincere.vote.mat!\n")
    }
    # get rows where ballot.mat has multiple entries
    multi_ballot <- apply(ballot.mat, 1, sum) > 1
    # replace with sincere vote
    ballot.mat[multi_ballot, ] <- sincere.vote.mat[multi_ballot, ]
    # problem here: no guarantee that ties include sincere option...
  }

  colnames(ballot.mat) = colnames(eu.mat)
  ballot.mat       
}

# Best thing to do is to have a n x 6 matrix that is 1 if vote is "natural" and 0 otherwise (if AV, at least)
optimal.vote.from.V.mat = function(V.mat){
  out = rep(NA, nrow(V.mat))
  for(j in 1:ncol(V.mat)){
    out[V.mat[,j] == 1] = colnames(V.mat)[j]
  }
  out
}

optimal.vote.from.eu.mat = function(eu.mat){
  optimal.vote.from.V.mat(ballot.mat.from.eu.mat(eu.mat))
}

sincere.vote.mat.from.U = function(U, rule = "plurality", candidates = c("A", "B", "C")){
  if(rule %in% c("AV")){
    ballots = apply(permutations(n = length(candidates), r = length(candidates), v = candidates, repeats.allowed = F), 1, paste, collapse = "")
    sincere.P = matrix(NA, nrow = length(ballots), ncol = length(candidates))
    for(i in 1:length(ballots)){
      for(j in 1:length(candidates)){
        sincere.P[i,j] = length(candidates) - which(grepl(candidates[j], str_split(ballots[i], "")[[1]])) 
      }
    }
  }else{
    ballots = candidates
    sincere.P = diag(length(candidates)) # generally 3
  }
  
  sincere.eu.by.ballot = t(sincere.P%*%t(U))
  colnames(sincere.eu.by.ballot) = ballots
  ballot.mat.from.eu.mat(sincere.eu.by.ballot) 

}

### --- MWE SCRIPT ---

# Define starting data
# Same as THA_2011 case
v_vec_start = c(0.2736, 0.3521, 0.1385, 0.1928, 0.0265, 0.0161)
s = 85

util_mat = tribble(~A, ~B, ~C,
                8.57, 3.33, 4.40,
                0.1, 0.14, 0.13,
                3.1, 9.13, 4.25,
                1.72, 8.81, 3.50,
                9.969, 4.96, 4.76,
                10, 6.5, 3.7,
                10.4, 5.8, 5.1)
V0_temp = sincere.vote.mat.from.U(util_mat, rule = "AV")

# Define first election (w/ 10,000 electorate) --------------------
# Sincerity adjustment breaks things bc it's bigger than the differences in EU
e_big = irv_election(n = 10000) %>%
            election_event_probs(method = "en",  
                                 alpha = (v_vec_start * s)) %>% 
            combine_P_matrices()
e_big = e_big[, -7] %>% t()
eu_big = t(e_big%*%t(util_mat))

# Get optimal vote
big_without = ballot.mat.from.eu.mat(eu_big) %>% colSums
names(big_without) = c("ABC", "ACB", "BAC", "BCA", "CAB", "CBA")

big_with = ballot.mat.from.eu.mat(eu_big, 
                       break.ties.with.sincerity = TRUE,
                       sincere.vote.mat = V0_temp) %>% colSums
print(list(big_without = big_without, big_with = big_with))

# Define second election (w/ 1 electorate) ----------------------
# Sincerity adjustment does not break things
e_small = irv_election(n = 1) %>%
            election_event_probs(method = "en",  
                                 alpha = (v_vec_start * s)) %>%
            combine_P_matrices()
e_small = e_small[, -7] %>% t()
eu_small = t(e_small%*%t(util_mat))

# Get optimal vote
small_without = ballot.mat.from.eu.mat(eu_small) %>% colSums
names(small_without) = c("ABC", "ACB", "BAC", "BCA", "CAB", "CBA")

small_with = ballot.mat.from.eu.mat(eu_small, 
                       break.ties.with.sincerity = TRUE,
                       sincere.vote.mat = V0_temp) %>% colSums
print(list(small_without = small_without, small_with = small_with))


