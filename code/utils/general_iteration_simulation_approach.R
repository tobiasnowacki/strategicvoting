
library(stringr)
library(here)

# source(here("utils", "winner_vec_from_ballot_proportions_for_various_rules.r"))  -- this is used in iteration_simulation, which is now not used 
source(here("utils", "EU_given_piv_probs_and_utility.r"))
source(here("utils", "av_pivotal_probs_analytical_general_v2.r"))
source(here("utils", "plurality_pivotal_probabilities_analytical.r"))

## various functions below that are useful for CSES analysis 


# a matrix with n rows and a 1 in the column corresponding to the optimal ballot
ballot.mat.from.eu.mat = function(eu.mat){
  # this does not need weights. 
  max.eus = apply(eu.mat, 1, max, na.rm = T)
  ballot.mat = matrix(NA, ncol = 0, nrow = nrow(eu.mat))
  for(j in 1:ncol(eu.mat)){
    ballot.mat = cbind(ballot.mat, as.integer(eu.mat[,j] == max.eus))
  }
  colnames(ballot.mat) = colnames(eu.mat)
  ballot.mat       
}

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
  apply(V*weight.mat, 2, sum, na.rm = T)/apply(weight.mat, 2, sum, na.rm = T)
}

sincere.ballot.props.from.U.and.weights = function(U, weights, rule = "plurality", candidates = NULL){
  svm = sincere.vote.mat.from.U(U, rule = rule, candidates = candidates)
  ballot.props.from.vote.mat.and.weights(svm, weights)
}


repeated_myopic_strategic_vote = function(U, lambda = .2, weights = NULL, rule = "plurality", s.vec = NULL, iterations = 5, noisy = F){
  
  if(is.null(s.vec)){s.vec = rep(50, 4)}
  
  if(noisy){cat("Iteration 1 ")}
  out = list()
  # so we want to do a first sincere poll, and then lambda responding in level-1 way each time.
  out[[1]] = iteration_simulation(U, lambda = lambda, weights = weights, rule = rule, s.vec = s.vec, piv.probs.analytical = T)
  for(k in 2:iterations){
    if(noisy){cat(k, " ", sep = "")}
    out[[k]] = iteration_simulation(U, lambda = lambda, V0 = out[[k-1]]$V0, V1 = out[[k-1]]$V.mat, weights = weights, rule = rule, s.vec = s.vec, piv.probs.analytical = T, v.vec.before = out[[k-1]]$v.vec.after)
  }
  if(noisy){cat("Done.\n")}
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


iteration_simulation = function(U, V0 = NULL, V1 = NULL, lambda = .2, m = 300, M = 5000, weights = NULL, min.step.denom = 1, rule = "plurality", dm = NULL, threshold = NULL, break.ties = F, delta = .025, second.round.s = 85, noisy = F, piv.probs.analytical = F, s.vec = NULL, v.vec.before = NULL){
  
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
      pps = av.pivotal.event.probs.general(v.vec = c(v.vec, 0,0,0), s.vec = s.vec)  # adding truncated ballots
    }else if(rule == "plurality"){
      pps = plurality.pivotal.probabilities(v.vec = v.vec + runif(3, 0, .000001), s = s.vec[1]) # a touch of noise to prevent any alpha component from being 0 
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
  V.mat = ballot.mat.from.eu.mat(eu.by.ballot)
  v.vec.after =  lambda*ballot.props.from.vote.mat.and.weights(V.mat, weights) + (1 -lambda)*ballot.props.from.vote.mat.and.weights(V0, weights) # this is the expected result if a proportion lambda votes strategically and the rest sincerely
  
  # optimal vote
  opt.votes.strategic = optimal.vote.from.V.mat(V.mat)
  opt.votes.sincere = optimal.vote.from.V.mat(V0)
  
  # calculate tau 
  tau = get.tau.from.eu.by.ballot.and.V0(eu.by.ballot, V0)
  
  # now output 
  list(P.mat = P.mat/M, eu.by.ballot = eu.by.ballot, tau = tau, V.mat = V.mat, V0 = V0, V1 = V1, weights = weights, ballot.prop.mat = ballot.prop.mat, V1.ballot.prop.mat = V1.ballot.prop.mat, V0.ballot.prop.mat = V0.ballot.prop.mat, opt.votes.strategic = opt.votes.strategic, opt.votes.sincere = opt.votes.sincere, v.vec.before = v.vec, v.vec.after = v.vec.after)
}

  
