#### plurality pivotal probabilities 

# three-party version, paralleling the approach for AV 

## utility copied from AV...v2 used below 

library(gtools)
# library(MCMCpack) identical to gtools
# library(mixtools) slightly different 

get.increment.midpoints = function(start, end, increments, increments.at.end = NULL){
	inc.length = (end - start)/increments
	out = seq(start + inc.length/2, end - inc.length/2, length = increments)
	if(!is.null(increments.at.end)){
		to.add.at.end = get.increment.midpoints(mean(out[c(length(out) - 1, length(out))]), end, increments.at.end) 
		out = c(out[-length(out)], to.add.at.end)
	}
	out
}

pr.tie.for.first.1.and.2.analytical.beta = function(v.vec, s){
  new.v.vec = normalize.to.1(v.vec[1:2])
  alpha.vec = new.v.vec*s
  .5*dbeta(.5, alpha.vec[1], alpha.vec[2])
}

pr.tie.for.first.1.and.2.numerical.dirichlet = function(v.vec, s, increments = 50){
	ys = get.increment.midpoints(1/3, 1/2, increments)
	incr.length = ys[2] - ys[1]
	incr.length*sum(apply(cbind(ys, ys, 1-2*ys), 1, gtools::ddirichlet, s*v.vec))
}

# we want to weight the beta almost nothing until the third component gets pretty small, i.e. sum of first two gets big
pr.tie.for.first.1.and.2 = function(v.vec, s, increments = 50, mix.with.beta = T, mix.start = 2){
  beta.weight = 0
  beta.component = 0
  third.alpha = as.numeric(v.vec[3]*s)
  if(mix.with.beta){
    if(third.alpha < mix.start){
      beta.component = pr.tie.for.first.1.and.2.analytical.beta(v.vec, s)
      beta.weight = ifelse(third.alpha < 1, 1, 1 - (third.alpha - 1)/(mix.start - 1))
    }
  }
  dirichlet.component = pr.tie.for.first.1.and.2.numerical.dirichlet(v.vec, s, increments)
  (1-beta.weight)*dirichlet.component + beta.weight*beta.component 
#  old approach
#  beta.component = pr.tie.for.first.1.and.2.analytical.beta(v.vec, s)
#  beta.weight = (1 - v.vec[3])^k
#  beta.weight*the.beta + (1 - beta.weight)*dirichlet.component
}

# specifying the function from scratch -- gives identical answer to pr.tie.for.first.1.and.2.numerical.dirichlet 
dirichlet.integration.tie.for.first.1.and.2 = function(v.vec, s, increments = 100){
  D.term = gamma(sum(v.vec*s))/exp(sum(lgamma(v.vec*s)))
  zs = get.increment.midpoints(1/3, 1/2, increments = increments)
  incr.length = zs[2] - zs[1]
  ys = D.term*zs^(sum(v.vec[1:2])*s - 2)*(1 - 2*zs)^(v.vec[3]*s - 1)
  incr.length*sum(ys)
}

simpson_rule_tie_for_first_1_and_2 = function(v.vec, s, increments, eps = .000001){
  stopifnot(increments %% 2 == 0)  # n must be odd
  ys = seq(1/3, 1/2 - eps, length = increments + 1)
  simpson_weights = c(1, 4, rep(c(2,4), (increments - 2)/2), 1)
  ((ys[2] - ys[1])/3)*sum(simpson_weights*apply(cbind(ys, ys, 1-2*ys), 1, gtools::ddirichlet, s*v.vec))
}

incomplete.beta = function(x,a,b){
  pbeta(x,a,b)*beta(a,b)
}

# OK that's our solution
pr.tie.for.first.1.and.2.definite.integral = function(v.vec, s){
  alpha = v.vec*s
  D.term = gamma(sum(alpha))/exp(sum(lgamma(alpha)))
  denom = 2^(sum(alpha[1:2]) - 1) # only thing unclear is exactly why this is -1 and not -2: I think in the change of variables from 2z to t I need to put another 2 in the denominator, but have not completely worked it out. 
  (D.term/denom)*incomplete.beta(1/3, alpha[3], sum(alpha[1:2]) - 1)
}


normalize.to.1 = function(vec){vec/sum(vec, na.rm = T)}

plurality.pivotal.probabilities = function(v.vec, s){ # , increments = 50, mix.with.beta = T, mix.start = 2){
	list("AB" = pr.tie.for.first.1.and.2.definite.integral(v.vec, s), #, increments, mix.with.beta, mix.start), 
	"AC" = pr.tie.for.first.1.and.2.definite.integral(v.vec[c(1,3,2)], s), #, increments, mix.with.beta, mix.start), 
	"BC" = pr.tie.for.first.1.and.2.definite.integral(v.vec[c(2,3,1)], s)
	) # , increments, mix.with.beta, mix.start) 
	# )	
}

# could add myatt-fisher approach -- code copied from Nick project below. 

posterior.pivotal.probability = function(pie.vec, s, which.ones = c(1,1,0), prior = 1){
  pie = pie.vec[which(which.ones == 0)]
  inc.beta.part = incomplete.beta(1/3, pie*s + prior, (1 - pie)*s + prior)
  2^(pie*s)*inc.beta.part/2^(s + prior)  # the 2016 wp version. to ask Myatt: is the +1 in the denominator the prior? is 2^(pie*s) actually 2^(pie*s + 1 - prior)?
  # (2^-((1 - pie)*s + 1))*inc.beta.part  # their old specification -- they are the same.
}

# but this is only proportional -- the absolute difference does matter. 


### does mixtools and gtools give the same density? 

#v.vec = c(20, 18, .01)
#v.vec = v.vec/sum(v.vec)
#s = 50
#mixtools::ddirichlet(s*v.vec)


