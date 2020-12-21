# Function to run iterations

source("code/utils/new_sv.R")

sv_iter = function(
  U, 
  s, 
  starting.v.vec = NULL,
  weights = NULL, 
  rule = "plurality", 
  lambda = .1, 
  epsilon = .001, 
  max.iterations = 200, 
  until.convergence = TRUE, 
  sincere.proportion = 0, 
  candidates = c("a", "b", "c"), 
  ballots = c("abc", "acb", "bac", "bca", "cab", "cba"), 
  the.floor = .01, 
  noisy = FALSE){

  if(rule == "plurality"){ballots = candidates}
  if(is.null(weights)){weights = rep(1, nrow(U))}
  
  # get sincere profile(s) 
  sincere.P = sincere_P(rule)
  sincere.eu.by.ballot = as.matrix(U)%*%sincere.P
  colnames(sincere.eu.by.ballot) = ballots
  sincere.vote.mat = ballot.mat.from.eu.mat(sincere.eu.by.ballot, break.ties.with.sincerity = FALSE)
  sincere.pref.mat = sincere_pref_mat_from_U(U, rule = rule)
  if(is.null(starting.v.vec)){
    sincere.v.vec = ballot.props.from.vote.mat.and.weights(sincere.vote.mat, weights)
  } else{
    sincere.v.vec = starting.v.vec
  }
  # Run first iteration 
  out = list()
  out[[1]] = sv(U = U,
    weights = weights,
    v.vec = sincere.v.vec,
    s = s,
    rule = rule,
    V0 = sincere.vote.mat,
    sin_pref_mat = sincere.pref.mat)
  # cat(sincere.v.vec)
  # Run subsequent iterations
  for(i in 2:max.iterations){
    if(noisy){cat("Iteration ", i, "\n", sep ="")}
    strategic.v.vec.i = lambda*out[[i-1]]$best.response.v.vec + (1 - lambda)*out[[i-1]]$v.vec.before
    overall.v.vec.i = sincere.proportion*sincere.v.vec + (1-sincere.proportion)*strategic.v.vec.i

    # Replace this below with sv function 
    out[[i]] = sv(U = U,
      weights = weights,
      v.vec = overall.v.vec.i,
      s = s,
      rule = rule,
      V0 = sincere.vote.mat,
      sin_pref_mat = sincere.pref.mat
    )


  #   # convergence when strategic part is just like this best response
  #   vd = vec.distance(out[[i]]$best.response.v.vec, strategic.v.vec.i)
  #   out[[i]]$distance.from.last = vd
  #   if(is.na(vd)){cat("Can't compute distance.\n"); return(NULL)}
  #   if(until.convergence & vd < epsilon){
  #     if(noisy)(cat("Converged after ", k, " iterations!\n", sep = ""))
  #     return(out)
  #   }
    } 
  # if(noisy & until.convergence){cat("Did not converge!!!!!!!!!\n")}
  return(out)
}