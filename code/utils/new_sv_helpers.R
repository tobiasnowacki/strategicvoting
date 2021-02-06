

# Function to get a 0/1 matrix for voters' sincere preferences from their utilities
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
  ballot.mat.from.eu.mat(sincere.eu.by.ballot, break.ties.with.sincerity = FALSE) 

}

# Function to get sincere preferences from their utilities
sincere_pref_mat_from_U <- function(U, rule = "plurality", candidates = c("A", "B", "C")){
  if(rule %in% c("AV")){
    ballots = apply(gtools::permutations(n = length(candidates), r = length(candidates), v = candidates, repeats.allowed = F), 1, paste, collapse = "")
    sincere.P = matrix(NA, nrow = length(candidates), ncol = length(ballots))
    for(i in 1:length(candidates)){
      for(j in 1:length(ballots)){
        sincere.P[i,j] = length(candidates) - which(grepl(candidates[i], str_split(ballots[j], "")[[1]]))
      }
    }
  }else if(rule == "plurality"){
    ballots = candidates
    sincere.P = diag(length(candidates)) # generally 3
  }else{stop("We only know how to deal rule=plurality and rule=AV.")}

  sincere.eu.by.ballot = as.matrix(U)%*%sincere.P
  
  out <- t(apply(sincere.eu.by.ballot, 1, rank))
  colnames(out) <- ballots
  out
  
}

# Function to get new v_vec from vote matrix and weights
ballot.props.from.vote.mat.and.weights = function(V, weights){
  weight.mat = matrix(weights, 
                      nrow = nrow(V), 
                      ncol = ncol(V), 
                      byrow = FALSE)
  weight.mat[is.na(V)] = 0
  out = apply(V*weight.mat, 2, sum, na.rm = TRUE)/apply(weight.mat, 2, sum, na.rm = TRUE)
  out/sum(out) # not sure why I have to normalize, but there it is 
}

# Function to produce optimal vote matrix from EU matrix
# with sincerity adjustments as of November 2020.
ballot_mat_from_eu_mat <- function(
  eu_mat, 
  break_ties_with_sincerity = TRUE, 
  sincere_mat = NULL, 
  weight = 1e-10, 
  normalize_eu_mat = TRUE){
  
  if(normalize_eu_mat){
    max_eus = apply(eu_mat, 1, max, na.rm = T)
    eu_mat <- eu_mat/matrix(max_eus, nrow = nrow(eu_mat), ncol = ncol(eu_mat), byrow = F)
  }

  if(break_ties_with_sincerity){
    if(is.null(sincere_mat)){stop("you must pass a sincere.vote.mat!\n")}
    eu_mat <- eu_mat + sincere_mat*weight
  }
  # i wonder if there's a better way to do this...
  max_eus = apply(eu_mat, 1, max, na.rm = T)
  
  ballot_mat = matrix(NA, ncol = ncol(eu_mat), nrow = nrow(eu_mat))
  for(j in 1:ncol(eu_mat)){
    ballot_mat[,j] = as.integer(eu_mat[,j] == max_eus)
  }
  
  colnames(ballot_mat) = colnames(eu_mat)
  ballot_mat
  
}

# Reduces optimal vote matrix to a vector of opt votes
optimal.vote.from.V.mat = function(V.mat){
  out = rep(NA, nrow(V.mat))
  for(j in 1:ncol(V.mat)){
    out[V.mat[,j] == 1] = colnames(V.mat)[j]
  }
  out
}

# Reduces EU matrix to difference between best non-sincere and sincere vote
get.tau.from.eu.by.ballot.and.V0 = function(eu.by.ballot, V0){
  eu.without.sincere.fave = eu.with.only.sincere.fave = eu.by.ballot
  eu.without.sincere.fave[V0 == 1] = NA
  max.not.fave = apply(eu.without.sincere.fave, 1, max, na.rm = T) 
  eu.with.only.sincere.fave[V0 == 0] = NA
  tau = max.not.fave - apply(eu.with.only.sincere.fave, 1, sum, na.rm = T)
  tau[is.infinite(tau)] = NA
  tau
}

# Get sincere matrix of preferences
sincere_P = function(rule, 
                     K = 3, 
                     candidates = c("a", "b", "c"), 
                     ballots = c("abc", "acb", "bac", "bca", "cab", "cba")){
  if(rule == "plurality"){
    out = diag(K)
    rownames(out) = colnames(out) = candidates 
  }else if(rule == "AV"){
    out = t(rbind(
      c(2,1,0),
      c(2,0,1),
      c(1,2,0),
      c(0,2,1),
      c(1,0,2),
      c(0,1,2)))
    rownames(out) = candidates 
    colnames(out) = ballots
  }
  out
}

# Legacy function
ballot.mat.from.eu.mat = function(eu.mat, break.ties.with.sincerity = T, sincere.vote.mat = NULL){
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