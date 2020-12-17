library(here)
cses.location = here("..", "data", "cses", "data/")  # alternative without using here(): setwd("~/Dropbox" etc etc to the location of this file. ) 

c4 = read.csv(paste0(cses.location, "cses4.csv"))
c3 = read.csv(paste0(cses.location, "cses3.csv"))
c2 = read.csv(paste0(cses.location, "cses2.csv"))
c1 = read.csv(paste0(cses.location, "cses1.csv"))

# get utilities for top 3 

prefix = c("A", "B", "C", "D") # prefixes for variables in the CSES waves
varnames = c("3020", "3037", "3009", "3011") # the varname for like-dislike in each survey
big_list = list()
bigger_list = list()

#we add some noise below to break ties. initially it was a very small amount, but I think we want to assume that integer CSES scores are rounded from an underlying scale that is approximately uniform. 
ubound = .5

# putting utility matrices in a list
# perhaps I can organize here so that the expected order of finish is ABC given a sincere plurality contest? 
for(j in 4:1){
  d = get(paste0("c", j))
  p = prefix[j]
  case_name = paste0(p, "1004")
  varname = varnames[j]
  cases = unique(as.character(d[[case_name]]))
  for(i in 1:length(cases)){
    weights = d[d[[case_name]] == cases[i], paste0(p, "1010_3")]
    # I think I did this to have a list with more parties in the Australian case. 
    if(grepl("AUS", cases[i])){
      X = d[d[[case_name]] == cases[i], paste0(p, varname, "_", LETTERS[1:9])]
      X[X > 10] = NA    # replacing 99 with NA
      # add some noise
      X = X + runif(n = nrow(X)*ncol(X), min = -ubound, max = ubound)
      colnames(X) = LETTERS[1:9]
      bigger_list[[cases[i]]] = list("U" = X, "weights" = weights)
    }
    X = d[d[[case_name]] == cases[i], paste0(p, varname, "_", c("A", "B", "C"))]
    X[X > 10] = NA    # replacing 99 with NA
    weights = d[d[[case_name]] == cases[i], paste0(p, "1010_3")]
    # add some noise
    X = X + runif(n = nrow(X)*ncol(X), min = -ubound, max = ubound)
    # now we rename such that the expected plurality (first-preference) order of finish is ABC.
    # determine order of finish
    V0 = sincere.vote.mat.from.U(X, rule = "plurality", candidates = c("A", "B", "C"))
    v.vec = ballot.props.from.vote.mat.and.weights(V0, weights = weights) + runif(n = 3, min = 0, max = .0001) # tie always possible
    v.vec = v.vec/sum(v.vec)
    # reorganize so columns are in order of finish
    X = X[, order(v.vec, decreasing = T)]
    # rename so A is the expected winner
    colnames(X) = c("A", "B", "C")
    if(length(unique(X[, 3])) == 1){next}
    cat(" == ", cases[i], " == ")
    big_list[[cases[i]]] = list("U" = X, "weights" = weights)
  }
}
cat("\n")

cases = sort(names(big_list))
save(big_list, cases, file = here("..", "data", "cses", "data", "big_list_2.RData"))

