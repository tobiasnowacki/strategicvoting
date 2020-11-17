nat_vote_mat <- function(sincere.vote.mat){
  if(ncol(sincere.vote.mat) != 6){
    cat("Not AV vote mat!")
  }
  nat_mat <- matrix(NA, nrow = nrow(sincere.vote.mat), 
                    ncol = ncol(sincere.vote.mat))
  ABC_vote <- sincere.vote.mat[, 1] == 1
  nat_mat[ABC_vote] <- rep(c(0, -100, 0, -100, 0, -100), each = sum(ABC_vote))
  ACB_vote <- sincere.vote.mat[, 2] == 1
  nat_mat[ACB_vote] <- rep(c(-100, 0, 0, -100, 0, -100), each = sum(ACB_vote))
  BAC_vote <- sincere.vote.mat[, 3] == 1
  nat_mat[BAC_vote] <- rep(c(0, -100, 0, -100, -100, 0), each = sum(BAC_vote))
  BCA_vote <- sincere.vote.mat[, 4] == 1
  nat_mat[BCA_vote] <- rep(c(0, -100, -100, 0, -100, 0), each = sum(BCA_vote))
  CAB_vote <- sincere.vote.mat[, 5] == 1
  nat_mat[CAB_vote] <- rep(c(-100, 0, -100, 0, 0, -100), each = sum(CAB_vote))
  CBA_vote <- sincere.vote.mat[, 6] == 1
  nat_mat[CBA_vote] <- rep(c(-100, 0, -100, 0, -100, 0), each = sum(CBA_vote))
  return(nat_mat)
}


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

ballot_mat_from_eu_mat <- function(eu_mat, break_ties_with_sincerity = TRUE, sincere_mat = NULL, weight = 1e-15, normalize_eu_mat = TRUE){
  
  if(normalize_eu_mat){
    max_eus = apply(eu_mat, 1, max, na.rm = T)
    eu_mat <- eu_mat/matrix(max_eus, nrow = nrow(eu_mat), ncol = ncol(eu_mat), byrow = F)
  }

  if(break_ties_with_sincerity){
    if(is.null(sincere_mat)){stop("you must pass a sincere.vote.mat!\n")}
    eu_mat <- eu_mat + sincere_mat*weight
  }

  max_eus = apply(eu_mat, 1, max, na.rm = T)
  
  ballot_mat = matrix(NA, ncol = ncol(eu_mat), nrow = nrow(eu_mat))
  for(j in 1:ncol(eu_mat)){
    ballot_mat[,j] = as.integer(eu_mat[,j] == max_eus)
  }
  
  colnames(ballot_mat) = colnames(eu_mat)
  ballot_mat
  
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

### refactored code: to make simple and add stuff later 

## what do we need? a function to compute expected utility, best responses tc  given some v.vec and s. 

sincere_P = function(rule, K = 3, candidates = c("a", "b", "c"), ballots = c("abc", "acb", "bac", "bca", "cab", "cba")){
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

# this is the function used as of e.g. Toulouse talk 
best_response_etc_given_U_v_vec_s = function(U, v.vec, s, weights = NULL, rule = "plurality", the.floor = .01, candidates = c("a", "b", "c"), sincere.vote.mat = NULL, normalize.P.mat = F){
  
  # normalize.P.mat makes the analysis (e.g. expected benefit of strategic voting) be conditional on being pivotal. Not normalizing means using \tilde{P},   
  
  # get P matrix
  if(rule == "AV"){
    stopifnot(length(v.vec) == 6)
    this.v.vec = c(v.vec + rep(the.floor, 6), 0,0,0)
    pps = av.pivotal.event.probs.general(v.vec = this.v.vec/sum(this.v.vec), s.vec = rep(s, 4)) 
    ballots = c("abc", "acb", "bac", "bca", "cab", "cba")
  }else if(rule == "plurality"){
    stopifnot(length(v.vec) == 3)
    this.v.vec = v.vec + rep(the.floor, 3)
    pps = plurality.pivotal.probabilities(v.vec = this.v.vec/sum(this.v.vec), s = s)  
    ballots = candidates
  }
  P.mat = P_mat_at_pivotal_events(pps, rule = rule, ballots = ballots, normalize = normalize.P.mat)
  probability.pivotal = sum(P.mat[,1]) # unique(apply(P.mat, 2, sum))
  # stopifnot(length(probability.pivotal) == 1)
  normalized.P.mat = P.mat/probability.pivotal

  rownames(P.mat) = candidates
  colnames(P.mat) = ballots
  
  # get expected utility matrix
  eu.by.ballot = U%*%normalized.P.mat #  conditional on being pivotal
  colnames(eu.by.ballot) = ballots

  # get optimal votes 
  if(is.null(sincere.vote.mat)){
    sincere.P = sincere_P(rule)
    sincere.eu.by.ballot = U%*%sincere.P
    colnames(sincere.eu.by.ballot) = ballots
    sincere.vote.mat = ballot.mat.from.eu.mat(sincere.eu.by.ballot) 
  }
  
  V.mat = ballot.mat.from.eu.mat(eu.by.ballot, break.ties.with.sincerity = T, sincere.vote.mat = sincere.vote.mat)
  if(length(table(apply(V.mat, 1, sum))) > 1){cat("Ties in V.mat!!\n")} # better error handling. 
  
  if(is.null(weights)){weights = rep(1, nrow(V.mat))}
  
  best.response.v.vec = ballot.props.from.vote.mat.and.weights(V.mat, weights)
  
  list(best.response.v.vec = best.response.v.vec, sincere.vote.mat = sincere.vote.mat, eu.by.ballot = eu.by.ballot, weights = weights, V.mat = V.mat, P.mat = P.mat, v.vec.before = v.vec, probability.pivotal = probability.pivotal)
  
}


iterated_best_response_sequence = function(U, s, weights = NULL, rule = "plurality", lambda = .1, epsilon = .001, max.iterations = 200, until.convergence = T, sincere.proportion = 0, candidates = c("a", "b", "c"), ballots = c("abc", "acb", "bac", "bca", "cab", "cba"), the.floor = .01, noisy = F){

  if(rule == "plurality"){ballots = candidates}
  if(is.null(weights)){weights = rep(1, nrow(U))}
  
  # get sincere profile 
  sincere.P = sincere_P(rule)
  sincere.eu.by.ballot = U%*%sincere.P
  colnames(sincere.eu.by.ballot) = ballots
  sincere.vote.mat = ballot.mat.from.eu.mat(sincere.eu.by.ballot)
  sincere.v.vec = ballot.props.from.vote.mat.and.weights(sincere.vote.mat, weights)
  
  out = list()
  out[[1]] = best_response_etc_given_U_v_vec_s(U = U, v.vec = sincere.v.vec, s = s, weights = weights, rule = rule, the.floor = the.floor, candidates = candidates, sincere.vote.mat = sincere.vote.mat)

  for(i in 2:max.iterations){
    if(noisy){cat("Iteration ", i, "\n", sep ="")}
    strategic.v.vec.i = lambda*out[[i-1]]$best.response.v.vec + (1 - lambda)*out[[i-1]]$v.vec.before
    overall.v.vec.i = sincere.proportion*sincere.v.vec + (1-sincere.proportion)*strategic.v.vec.i
    out[[i]] = best_response_etc_given_U_v_vec_s(U = U, v.vec = overall.v.vec.i, s = s, weights = weights, rule = rule, the.floor = the.floor, candidates = candidates, sincere.vote.mat = sincere.vote.mat)
    # convergence when strategic part is just like this best response
    vd = vec.distance(out[[i]]$best.response.v.vec, strategic.v.vec.i)
    out[[i]]$distance.from.last = vd
    if(is.na(vd)){cat("Can't compute distance.\n"); return(NULL)}
    if(until.convergence & vd < epsilon){
      if(noisy)(cat("Converged after ", k, " iterations!\n", sep = ""))
      return(out)
    }
  } 
  if(noisy & until.convergence){cat("Did not converge!!!!!!!!!\n")}
  out
}