
#### See refactored code far below
### uses some of the functions at the top -- could still use some cleaning out. 

library(stringr)
library(here)

# source(here("utils", "winner_vec_from_ballot_proportions_for_various_rules.r"))  -- this is used in iteration_simulation, which is now not used 
source(here("code", "utils", "EU_given_piv_probs_and_utility.r"))
source(here("code", "utils", "av_pivotal_probs_analytical_general_v2.r"))
source(here("code", "utils", "plurality_pivotal_probabilities_analytical.r"))

## various functions below that are useful for CSES analysis 


# a matrix with n rows and a 1 in the column corresponding to the optimal ballot
# ballot.mat.from.eu.mat = function(eu.mat, break.ties.with.sincerity = F, sincere.vote.mat = NULL){
#   # this does not need weights. 
  
#   if(break.ties.with.sincerity){
#     if(is.null(sincere.vote.mat)){
#       cat("you must pass a sincere.vote.mat!\n")
#       sincere.vote.mat - 1
#     }
#     eu.mat = eu.mat + sincere.vote.mat*(10^(-10)) # this breaks ties in favor of the sincere vote. but actually we need to break all ties -- why is this happening? PICK UP HERE. 
#   }
  
#   max.eus = apply(eu.mat, 1, max, na.rm = T)
#   ballot.mat = matrix(NA, ncol = ncol(eu.mat), nrow = nrow(eu.mat))
#   for(j in 1:ncol(eu.mat)){
#     ballot.mat[,j] = as.integer(eu.mat[,j] == max.eus)
#   }
#   colnames(ballot.mat) = colnames(eu.mat)
#   ballot.mat       
# }

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

get_two_levels_of_strategicness = function(U, lambda = .2, m = 300, M = 5000, weights = NULL, min.step.denom = 1, rule = "plurality", candidates = NULL, dm = NULL, threshold = NULL){
  out1 = iteration_simulation(U, lambda = lambda, m = m, M = M, weights = weights, min.step.denom = min.step.denom, rule = rule, candidates = candidates, dm = dm, thrsehold = NULL)
  out2 = iteration_simulation(U, V0 = out1$V0, V1 = out1$V1, lambda = lambda, m = m, M = M, weights = weights, min.step.denom = min.step.denom, rule = rule, candidates = candidates, dm = dm, thrsehold = NULL)
  list(out1, out2)
}

best.nonsincere.vote.from.eu.by.ballot.and.V0 = function(eu.by.ballot, V0){
  eu.without.sincere.fave = eu.by.ballot
  eu.without.sincere.fave[V0 == 1] = NA
  optimal.vote.from.eu.mat(eu.without.sincere.fave)
}

get.tau.from.eu.by.ballot.and.V0 = function(eu.by.ballot, V0){
  eu.without.sincere.fave = eu.with.only.sincere.fave = eu.by.ballot
  eu.without.sincere.fave[V0 == 1] = NA
  max.not.fave = apply(eu.without.sincere.fave, 1, max, na.rm = T) 
  eu.with.only.sincere.fave[V0 == 0] = NA
  tau = max.not.fave - apply(eu.with.only.sincere.fave, 1, sum, na.rm = T)
  tau[is.infinite(tau)] = NA
  tau
}

ballot.props.from.vote.mat.and.weights = function(V, weights){
  weight.mat = matrix(weights, nrow = nrow(V), ncol = ncol(V), byrow = F)
  weight.mat[is.na(V)] = 0
  out = apply(V*weight.mat, 2, sum, na.rm = T)/apply(weight.mat, 2, sum, na.rm = T)
  out/sum(out) # not sure why I have to normalize, but there it is 
}

sincere.ballot.props.from.U.and.weights = function(U, weights, rule = "plurality", candidates = NULL){
  if(is.null(candidates)){candidates = colnames(U)}
  svm = sincere.vote.mat.from.U(U, rule = rule, candidates = candidates)
  ballot.props.from.vote.mat.and.weights(svm, weights)
}

vec.distance = function(v1, v2){
  sqrt(sum((v1 - v2)^2))
}

extract.from.list = function(l, name){l[[name]]}

repeated_myopic_strategic_vote = function(U, lambda = .2, weights = NULL, rule = "plurality", s.vec = NULL, iterations = 5, noisy = F, the.floor = .01, until.convergence = F, convergence.tolerance = .001, convergence.iteration.limit = 100, sincere.proportion = 0){
  
  if(is.null(s.vec)){s.vec = rep(50, 4)}
  
  # if(noisy){cat("Iteration 1 ")}
  out = list()
  # so we want to do a first sincere poll, and then lambda responding in level-1 way each time.
  out[[1]] = iteration_simulation(U, lambda = lambda, weights = weights, rule = rule, s.vec = s.vec, piv.probs.analytical = T, the.floor = the.floor)
  if(until.convergence){iterations = convergence.iteration.limit}
  for(k in 2:iterations){
    # if(noisy){cat(k, " ", sep = "")}
    last.v.vec = out[[k-1]]$v.vec.after
    out[[k]] = iteration_simulation(U, lambda = lambda, V0 = out[[k-1]]$V0, V1 = out[[k-1]]$V.mat, weights = weights, rule = rule, s.vec = s.vec, piv.probs.analytical = T, v.vec.before = last.v.vec, update.based.on = "previous", the.floor = the.floor, sincere.proportion = sincere.proportion)
    this.v.vec = out[[k]]$fully.strategic.v.vec # this is the proportions of ballots if everyone updated. if this is the same as last.v.vec, then we have converged. this gives us a lambda-independent measure of convergence. 
    if(until.convergence){
      vd = vec.distance(last.v.vec, this.v.vec)
      out[[k]]$distance.from.last = vd
      if(is.na(vd)){cat("Can't compute distance.\n"); return(NULL)}
      # if(noisy){cat("Distance from last: ", round(vd, 4), ".\n", sep = "")}
      if(vd < convergence.tolerance){
        if(noisy)(cat("Converged after ", k, " iterations!\n", sep = ""))
        return(out)
      }
    }
  }
  if(noisy){cat("Did not converge!!!!!!!!!\n")}
  out
  
}

v.vec.mat.from.repeated.polls = function(out){
  this.out = 
  rbind(
    out[[1]]$v.vec.before,
    out[[1]]$v.vec.after
  )
  for(k in 2:length(out)){
    this.out = rbind(this.out, out[[k]]$v.vec.after)
  }
  this.out
}

victory.probs.from.sims = function(sims, rule = "plurality", return.matrix = F){
  if(rule == "AV"){
    stopifnot(ncol(sims) == 6) # don't know how to handle yet
    fp.mat = cbind(sims[,1] + sims[,2], sims[,3] + sims[,4], sims[,5] + sims[,6])
    row.mins = apply(fp.mat, 1, min)
    no.1 = which(fp.mat[,1] == row.mins); no.2 = which(fp.mat[,2] == row.mins); no.3 = which(fp.mat[,3] == row.mins)
    fp.mat[no.1, 2] = fp.mat[no.1, 2] + sims[no.1, 1]
    fp.mat[no.1, 3] = fp.mat[no.1, 3] + sims[no.1, 2]
    fp.mat[no.2, 1] = fp.mat[no.2, 1] + sims[no.2, 3]
    fp.mat[no.2, 3] = fp.mat[no.2, 3] + sims[no.2, 4]
    fp.mat[no.3, 1] = fp.mat[no.3, 1] + sims[no.3, 5]
    fp.mat[no.3, 2] = fp.mat[no.3, 2] + sims[no.3, 6]
    sims = fp.mat   # now it's just like plurality 
  }
  row.maxes = apply(sims, 1, max) # shouldn't have NAs.
  vps = as.integer(sims[,1] == row.maxes)
  for(j in 2:ncol(sims)){
    vps = cbind(vps, as.integer(sims[,j] == row.maxes))
  }
  if(return.matrix){return(vps)}
  apply(vps, 2, mean)
}

victory.prob.mat.from.v.vec.mat = function(v.vec.mat, s, rule = "plurality", M = 1000){
  vpm = matrix(NA, nrow = nrow(v.vec.mat), ncol = 3) # three candidates only
  for(i in 1:nrow(vpm)){
    sims = rdirichlet(M, alpha = s*v.vec.mat[i,]) # not allowing s.vec here
    vpm[i,] = victory.probs.from.sims(sims, rule = rule)
  }
  vpm
}

victory.prob.mat.from.repeated.polls.object = function(obj, s, rule = "plurality", M = 1000){
  victory.prob.mat.from.v.vec.mat(v.vec.mat.from.repeated.polls(obj), s = s, rule = rule, M = M)  
}

victory.prob.mat.from.U.and.lambda = function(U, s.vec, lambda = .2, weights = NULL, rule = "plurality", iterations = 5){
  victory.prob.mat.from.repeated.polls.object(repeated_myopic_strategic_vote(U = U, lambda = lambda, weights= weights, rule = rule, s.vec = s.vec, iterations = iterations), s = s.vec[1], rule = rule)
}

condorcet.winner = function(U, weights = NULL){
  # returns a logical 3-vec identifying the CW if there is one 
  if(is.null(weights)){weights = rep(1, nrow(U))}
  stopifnot(ncol(U) == 3)
  win.mat = matrix(0, nrow = 3, ncol = 3)
  for(i in 1:nrow(win.mat)){
    for(j in 1:ncol(win.mat)){
      if(i == j){next}      
      win.mat[i,j] = as.integer(sum(weights*(as.integer(U[,i] > U[,j])), na.rm = T)/sum(weights[!is.na(U[,i]) & !is.na(U[,j])], na.rm = T) > .5)
    }
  }
  apply(win.mat, 1, sum) == 2
}


# this is the old function -- superseded 
iteration_simulation = function(U, V0 = NULL, V1 = NULL, lambda = .2, m = 300, M = 5000, weights = NULL, min.step.denom = 1, rule = "plurality", dm = NULL, threshold = NULL, break.ties = F, delta = .025, second.round.s = 85, noisy = F, piv.probs.analytical = F, s.vec = NULL, v.vec.before = NULL, update.based.on = "sincere", the.floor = .01, sincere.proportion = 0){

  
  ## U is the matrix of utilities: n rows and k columns, where n is the number of survey respondents and k is the number of candidates. colnames must include the names of the candidates
  ## V0 is the (optional) matrix of sincere votes. it is n x K, where n is the number of survey respondents and K is the number of possible ballots
  ## V1 is the (optional) matrix of level-1 strategic votes. same dimensions. 
  ## This function calculates the optimal votes of voters who believe that the expected outcome is V0 mixed with V1 at a mixture rate of 1 - lambda, lambda
  
  ## lambda is the rate at which V1 should be mixed with V0 in determining the voter's expected result. With lambda = .2, for example, the voter thinks that a proportion .2 of voters is level-1 strategic
  ## m is the number of (weighted) draws from the survey. this depends on/reflects the precision of the voter's beliefs based on the survey: for large m, the voter thinks the survey is correct; for small m, anything is possible. 
  ## M is the number of simulation runs. 
  ## weights in the survey. this determines the probability of selection in resampling.
  ## min.step.denom affects how much noise goes into U if break.ties = T
  ## rule: plurality, dhondt, saintelague, 2R, AV. 
  ## dm: district magnitude, required for dhondt and saintelague
  ## threshold: optional threshold for pr systems (leaving null is same as assuming 0)
  ## break.ties: should we add some noise to the U matrix so that no one has ties? 
  ## delta: voting weight of hypothetical voter. bias-variance tradeoff.
  ## second.round.s: precision of beliefs about second round (in beta distribution that characterizes probability of victory based on sincere like-dislike between the two who advance)
  ## piv.probs.analytical: do you want to get your pivotal probs from shuffling the survey, or from a Dirichlet?
  ## s.vec: for the analytical AV and plurality version 
  
  ## v.vec.before: what is the expected result at which we want to calculate the optimal ballots? this is used for iteration. 
  ## update.based.on can be sincere or previous: should the v.vec.after be the weighted avg of the strategic voting result and the sincere voting result, or the previous result (v.vec.before)? 
  ## the.floor -- minimum value of v.vec. stopgap to deal with weird behavior near edges of Dirichlet
  
  if(break.ties){
    # break ties in utility matrix 
    U.vec = na.omit(sort(unique(as.vector(as.matrix(U)))))
    min.step = min((U.vec[2:length(U.vec)] - U.vec[1:(length(U.vec) - 1)])/2, na.rm = T) # the smallest distance between any unique utilities. 
    U = U + matrix(runif(ncol(U)*nrow(U), min = -min.step/min.step.denom, max = min.step/min.step.denom), nrow(U), ncol(U)) # this allows us to break ties randomly ex ante. larger min.step.denom means we keep more indifference. 
  }
  
  stopifnot(!is.null(colnames(U)))
  stopifnot(length(unique(colnames(U))) == length(colnames(U)))
  # the columnnames of U should be in alphabetical order for AV, because this is how the permutations get made.
  candidates = sort(colnames(U))
  U = U[, candidates]
  
  if(rule %in% c("AV")){
    K = factorial(ncol(U))
    ballots = apply(permutations(n = length(candidates), r = length(candidates), v = candidates, repeats.allowed = F), 1, paste, collapse = "")
    sincere.P = matrix(NA, nrow = length(ballots), ncol = length(candidates))
    for(i in 1:length(ballots)){
      for(j in 1:length(candidates)){
        sincere.P[i,j] = length(candidates) - which(grepl(candidates[j], str_split(ballots[i], "")[[1]])) 
      }
    }
    rownames(sincere.P) = ballots
    colnames(sincere.P) = candidates
  }else{
    K = ncol(U)
    ballots = candidates
    sincere.P = diag(K)
  }
  
  if(is.null(V0)){ # this is referred as "sincere vote mat" below, but note that it need not be sincere -- could be intended votes from survey, though this will complicate things a little because when we update we might get more sincere (if deviations from sincerity were not strategic)  
    sincere.eu.by.ballot = t(sincere.P%*%t(U))
    colnames(sincere.eu.by.ballot) = ballots
    V0 = ballot.mat.from.eu.mat(sincere.eu.by.ballot) 
  }
  
  if(is.null(V1)){  # must be calculating the L1 incentives. 
    V1 = V0 
  }
  
  if(is.null(weights)){cat("(No weights passed.)\n"); weights = rep(1, nrow(U))}
  
  if(rule %in% c("dhondt", "saintelague")){stopifnot(!is.null(dm))}

  if(piv.probs.analytical){
    stopifnot(rule %in% c("plurality", "AV"))
    stopifnot(ncol(U) == 3) # we can't do more than three with this technique -- expanding would require generalizing the P_mat_at_pivotal_events step below 
    stopifnot(!is.null(s.vec))
    # need to get get pivotal probabilities
    # first get expected result -- does not depend on the rule
    if(!is.null(v.vec.before)){
      v.vec = v.vec.before
    }else{
      sincere.v.vec = ballot.props.from.vote.mat.and.weights(V0, weights)
      strategic.v.vec = ballot.props.from.vote.mat.and.weights(V1, weights) # if we passed a matrix of strategic votes, which is what the level-2 voter thinks others will do based on the first poll, then we use that here.
      v.vec = lambda*strategic.v.vec + (1 - lambda)*sincere.v.vec # expected result
    }
    # now get pivotal probabilities -- depends on the rule
    if(rule == "AV"){
      this.v.vec = c(v.vec + rep(the.floor, 6), 0,0,0)
      pps = av.pivotal.event.probs.general(v.vec = this.v.vec/sum(this.v.vec), s.vec = s.vec)  # touch of mass
    }else if(rule == "plurality"){
      this.v.vec = v.vec + rep(the.floor, 3)
      pps = plurality.pivotal.probabilities(v.vec = this.v.vec/sum(this.v.vec), s = s.vec[1]) # a touch of mass to prevent any alpha component from being 0 
    }
    P.mat = t(P_mat_at_pivotal_events(pps, rule = rule, ballots = ballots))
    rownames(P.mat) = ballots
    colnames(P.mat) = candidates
    ballot.prop.mat = V1.ballot.prop.mat = V0.ballot.prop.mat = NULL
  }else{
    P.mat = matrix(0, nrow = length(ballots), ncol = length(candidates))
    ballot.prop.mat = matrix(NA, nrow = M, ncol = ncol(V0))
    colnames(ballot.prop.mat) = ballots
    V1.ballot.prop.mat = V0.ballot.prop.mat = ballot.prop.mat
    v.vec = NULL
    
    # Calculate P.mat: K rows (one for each possible ballot) and k columns (one for each possible winner), saying the proportion of M draws in which each candidate wins, given a delta push for each ballot. 
    for(i in 1:M){
      #if(noisy){ # } & i%%100 == 0){
      #  cat("Iteration ", i, ", ")
      #}
      indices = sample(1:nrow(U), size = m, replace = T, prob = weights)  # this is the only place the weighte should matter
      # get assumed votes for these voters
      V1.ballot.props = apply(V1[indices, ], 2, mean, na.rm = T)
      V0.ballot.props = apply(V0[indices, ], 2, mean, na.rm = T)
      V1.ballot.prop.mat[i,] = V1.ballot.props  # TODO: this is diagnostic. remove when no longer needed.  
      V0.ballot.prop.mat[i,] = V0.ballot.props
      # assumption is that result for this draw of voters is a weighted average of the sincere ballot props and the level 1 strategic ballot props. (if I don't pass a V1, they are the same thing.) 
      ballot.props = lambda*V1.ballot.props + (1 - lambda)*V0.ballot.props  
      ballot.props = ballot.props + runif(length(ballot.props), min = -.001, max = .001) # break ties in vote counts
      ballot.props = ballot.props/sum(ballot.props) # renormalizing
      names(ballot.props) = ballots
      ballot.prop.mat[i,] = ballot.props
      # now we perturb for each possible ballot
      for(j in 1:length(ballot.props)){
        tbp = ballot.props
        tbp[j] = tbp[j] + delta
        this.vec = P.mat[j, ] + winner.vec.from.ballot.props(tbp, rule = rule, dm = dm, threshold = threshold, candidates = candidates, second.round.s = second.round.s)
        P.mat[j, ] = this.vec
      }
    }
    if(noisy){cat("Iterations done.\n")}
    # now we have the P matrix.
    # we mix it with the sincere matrix so that ties are broken in favor of sincerity
    # how much weight should be put on sincerity? any increment could in principle reverse the "true" order, but these are simulations so differences that small are kind of noise.  
    sincere.fraction = 1/(200*M)
    P.mat = (1 - sincere.fraction)*P.mat + sincere.fraction*sincere.P
  }
  
  # we get the optimal vote given this P for everyone in the survey
  eu.by.ballot = t(P.mat%*%t(U))/M
  colnames(eu.by.ballot) = ballots
  V.mat = ballot.mat.from.eu.mat(eu.by.ballot, break.ties.with.sincerity = T, sincere.vote.mat = V0)
  if(length(table(apply(V.mat, 1, sum))) > 1){cat("Ties in V.mat!!\n")} # better error handling. 
  
  # this is probably redundant -- I'm sure this calc has been made elsewhere. 
  fully.strategic.v.vec = ballot.props.from.vote.mat.and.weights(V.mat, weights)
  sincere.v.vec = ballot.props.from.vote.mat.and.weights(V0, weights)


  if(update.based.on == "sincere"){
    v.vec.after =  lambda*fully.strategic.v.vec + (1 -lambda)*ballot.props.from.vote.mat.and.weights(V0, weights) # this is the expected result if a proportion lambda votes strategically and the rest sincerely
  }else if(update.based.on == "previous"){
    v.vec.after =  sincere.proportion*sincere.v.vec + (1-sincere.proportion)*(lambda*fully.strategic.v.vec + (1 -lambda)*v.vec.before)
  }
  # normalizing -- should not be necessary. 
  v.vec.after = v.vec.after/sum(v.vec.after)
  
  # optimal vote
  opt.votes.strategic = optimal.vote.from.V.mat(V.mat)
  opt.votes.sincere = optimal.vote.from.V.mat(V0)
  
  # calculate tau 
  tau = get.tau.from.eu.by.ballot.and.V0(eu.by.ballot, V0)
  
  # now output 
  list(P.mat = P.mat/M, eu.by.ballot = eu.by.ballot, tau = tau, V.mat = V.mat, V0 = V0, V1 = V1, weights = weights, ballot.prop.mat = ballot.prop.mat, V1.ballot.prop.mat = V1.ballot.prop.mat, V0.ballot.prop.mat = V0.ballot.prop.mat, opt.votes.strategic = opt.votes.strategic, opt.votes.sincere = opt.votes.sincere, v.vec.before = v.vec, v.vec.after = v.vec.after, fully.strategic.v.vec = fully.strategic.v.vec)
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

## now for the iteration part 

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

distance.vec.from.ibrs = function(ibrs){
  unlist(lapply(ibrs, extract.from.list, "distance.from.last"))
}

v.vec.mat.from.ibrs = function(ibrs){
  do.call(rbind, lapply(ibrs, extract.from.list, "v.vec.before"))
}

expected.benefit.vec.given.eu.mat.and.sincere.vote.mat = function(eu.mat, sincere.vote.mat){
  max.eus = apply(eu.mat, 1, max, na.rm = T)
  eu.of.sincere = eu.mat*sincere.vote.mat
  eu.of.sincere[sincere.vote.mat == 0] = NA
  sincere.eus = apply(eu.of.sincere, 1, max, na.rm = T)
  max.eus - sincere.eus
}

expected.benefit.vec.given.iteration.output = function(iteration.output){
  eu.mat = iteration.output[["eu.by.ballot"]]
  max.eus = apply(eu.mat, 1, max, na.rm = T)
  sincere.vote.mat = iteration.output[["sincere.vote.mat"]]
  eu.of.sincere = eu.mat
  eu.of.sincere[sincere.vote.mat == 0] = NA
  sincere.eus = apply(eu.of.sincere, 1, max, na.rm = T)
  as.numeric(max.eus - sincere.eus)
}

expected.benefit.mat.from.ibrs = function(ibrs){
  do.call(cbind, lapply(ibrs, expected.benefit.vec.given.iteration.output))
}

pivotal.prob.vec.from.ibrs = function(ibrs){
  unlist(lapply(ibrs, function(l) l[["probability.pivotal"]]))
}


plot.ternary.basic = function(fp.vec, vertex.labels = c("a", "b", "c"), label.offset = .05, fp.result.cex = 1, fp.result.col = "black", secondary.line.col = "gray", secondary.line.lwd = 1, secondary.line.lty = 3, main = NULL, space = .1, add.av.lines = F, add.pl.lines = F){
  
  stopifnot(length(fp.vec) == 3)
  
  xs = c(0, 1) + c(-space, space); ys = c(0, sqrt(3/4)) + sqrt(3/4)*c(-space, space)
  plot(xs, ys, type = "n", xlab = "", ylab = "", axes = F, main = main)
  add.ternary.boundary()
  
  # add sincere point 
  add.ternary.point(fp.vec[c(1,3,2)], pch = 19, cex = fp.result.cex, col = fp.result.col)
  
  # label vertices
  add.ternary.text(c(1,0,0), vertex.labels[1], x.offset = -label.offset, y.offset = -sqrt(3/4)*label.offset)
  add.ternary.text(c(0, 0,1), vertex.labels[2], x.offset = 0, y.offset = sqrt(3/4)*label.offset)
  add.ternary.text(c(0,1,0), vertex.labels[3], x.offset = label.offset, y.offset = -sqrt(3/4)*label.offset)
  
  if(add.av.lines){
    # majority lines 
    add.ternary.lines(c(1/2, 1/2, 0), c(0, 1/2, 1/2), col = secondary.line.col, lwd = secondary.line.lwd, lty = secondary.line.lty)
    add.ternary.lines(c(0, 1/2, 1/2), c(1/2, 0, 1/2), col = secondary.line.col, lwd = secondary.line.lwd, lty = secondary.line.lty)
    add.ternary.lines(c(1/2, 0, 1/2), c(1/2, 1/2, 0), col = secondary.line.col, lwd = secondary.line.lwd, lty = secondary.line.lty)
    
    # elimination lines 
    add.ternary.lines(c(1/2, 1/4, 1/4), c(1/3, 1/3, 1/3), col = secondary.line.col, lwd = secondary.line.lwd, lty = 2)
    add.ternary.lines(c(1/4, 1/2, 1/4), c(1/3, 1/3, 1/3), col = secondary.line.col, lwd = secondary.line.lwd, lty = 2)
    add.ternary.lines(c(1/4, 1/4, 1/2), c(1/3, 1/3, 1/3), col = secondary.line.col, lwd = secondary.line.lwd, lty = 2)		
    
  }
  
  if(add.pl.lines){
    # plurality lines 
    add.ternary.lines(c(1/2, 1/2, 0), c(1/3, 1/3, 1/3), col = secondary.line.col, lwd = secondary.line.lwd, lty = 1)
    add.ternary.lines(c(1/2, 0, 1/2), c(1/3, 1/3, 1/3), col = secondary.line.col, lwd = secondary.line.lwd, lty = 1)
    add.ternary.lines(c(0, 1/2, 1/2), c(1/3, 1/3, 1/3), col = secondary.line.col, lwd = secondary.line.lwd, lty = 1)		
  }
  
}



# what about utility? 
# once I have the expected result, can do expected utility calc. 

winning.prob.mat.given.ibrs = function(ibrs, s, M = 10000, rule = "plurality"){
  v.vec.mat = v.vec.mat.from.ibrs(ibrs) # m x B
  out = matrix(NA, nrow = nrow(v.vec.mat), ncol = 3)
  for(i in 1:nrow(v.vec.mat)){
    results = rdirichlet(M, v.vec.mat[i,]*s)
    out[i,] = victory.probs.from.sims(results, rule) 
  }
  out
}


expected.utility.mat.and.winning.prob.mat.given.ibrs.and.U = function(ibrs, U, M = 10000, s, rule = "plurality"){
  winning.prob.mat = winning.prob.mat.given.ibrs(ibrs, s, M, rule) # m x k
  eum = matrix(NA, nrow = nrow(U), ncol = nrow(winning.prob.mat))
  for(j in 1:ncol(eum)){
    eum[,j] = U%*%winning.prob.mat[j,]  # n x k   k x 1
  }
  list("winning.prob.mat" = winning.prob.mat, "expected.utility.mat" = eum)
}
## pick up here: compare the winning probably matrix. 

expected.utility.vec.given.U.and.v.vec = function(U, v.vec, s, M = 10000, rule = "plurality"){
  results = rdirichlet(M, v.vec*s)
  winning.probs = as.matrix(victory.probs.from.sims(results, rule), ncol = 1)
  U%*%winning.probs
}

expected.utility.mat.given.ibrs.and.U = function(ibrs, U, M = 10000, s, rule = "plurality"){
  v.vec.mat = v.vec.mat.from.ibrs(ibrs) # m x B
  out = matrix(NA, nrow = nrow(U), ncol = nrow(v.vec.mat))
  for(j in 1:nrow(v.vec.mat)){
    out[,j] = expected.utility.vec.given.U.and.v.vec(U, v.vec.mat[j,], s = s, M = M, rule = rule)
  }
  out # n x m
}

weighted.average.utilities.given.bl.entry = function(l){
  U = l$U
  weights = l$weights
  U.is.na = matrix(0, nrow = nrow(U), ncol = ncol(U)); U.is.na[is.na(U)] = 1 
  missing.U = apply(U.is.na, 1, sum, na.rm = T)
  for(j in 1:3){U[missing.U > 0,j] = NA}
  weights[missing.U > 0] = 0
  weights = weights/sum(weights, na.rm = T)
  weight.mat = cbind(weights,weights,weights)
  apply(U*weight.mat, 2, sum, na.rm = T)
} 



  
