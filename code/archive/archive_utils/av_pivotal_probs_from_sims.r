# code for outputting pivotal probabilities from simulations given Dirichlet belief parameters 

av.pivotal.probs.from.sims = function(alpha.vec, v.vec = NULL, s.vec = NULL, M = 100000, tol = .005, type = "dirichlet"){
		
	# dirichlet case can take truncated ballots, listed as AB, AC, BA, BC, CA, CB, AX, BX, CX	
		
	# draw M simulated results
	if(type == "dirichlet"){
		if(length(alpha.vec) == 6){alpha.vec = c(alpha.vec, 0, 0, 0)}
		stopifnot(length(alpha.vec) == 9)
		sims = rdirichlet(n = M, alpha = alpha.vec)		
	}else if(type == "general"){
		stopifnot(length(v.vec) == 6)
		# first-preference vote distributions
		fp.alpha = s.vec[1]*c(sum(v.vec[1:2]), sum(v.vec[3:4]), sum(v.vec[5:6]))
		fp.sims = rdirichlet(n = M, alpha = fp.alpha)
		# second-preference shares 
		sp.shares.A = rbeta(n = M, s.vec[2]*v.vec[1], s.vec[2]*v.vec[2])
		sp.shares.B = rbeta(n = M, s.vec[3]*v.vec[3], s.vec[3]*v.vec[4])
		sp.shares.C = rbeta(n = M, s.vec[4]*v.vec[5], s.vec[4]*v.vec[6])
		sims = cbind(fp.sims[,1]*sp.shares.A, fp.sims[,1]*(1 - sp.shares.A), fp.sims[,2]*sp.shares.B, fp.sims[,2]*(1 - sp.shares.B), fp.sims[,3]*sp.shares.C, fp.sims[,3]*(1 - sp.shares.C), 0, 0, 0)  # we add zeroes here because there can't be truncated ballots in the general case.  
	}

	
	# get pivotal probabilities from the simulations
	
	# first preference vote shares for each candidate
	A.sum.1 = apply(sims[,c(1,2,7)], 1, sum)
	B.sum.1 = apply(sims[,c(3,4,8)], 1, sum)
	C.sum.1 = apply(sims[,c(5,6,9)], 1, sum)

	# identify cases where there is a second round -- if not, you can't be pivotal (assuming all ballot orders get some support).  
	FP.mat = cbind(A.sum.1, B.sum.1, C.sum.1)
	second.round = apply(FP.mat, 1, max) < .5    # used to say "- tol" here, but this seems wrong.  
			
	# FP margins
	AB.diff.1 = A.sum.1 - B.sum.1
	AC.diff.1 = A.sum.1 - C.sum.1
	BC.diff.1 = B.sum.1 - C.sum.1
	
	# second-round matchups 
	second.round.AB = AC.diff.1 > 0 & BC.diff.1 > 0 
	second.round.BC = AB.diff.1 < 0 & AC.diff.1 < 0 
	second.round.AC = AB.diff.1 > 0 & BC.diff.1 < 0 

	# majority margins for each pairing
	AB.diff.2 =  A.sum.1 + sims[,5] - (B.sum.1 + sims[,6]) 
	AC.diff.2 =  A.sum.1 + sims[,3] - (C.sum.1 + sims[,4]) 
	BC.diff.2 =  B.sum.1 + sims[,1] - (C.sum.1 + sims[,2]) 
	
	# first round ties for second 
	first.round.ab = AC.diff.1 < 0 & BC.diff.1 < 0 & abs(AB.diff.1) < tol/2 # C in first; A beats B (or vice versa) by less than tol/2. 
	first.round.ac = AB.diff.1 < 0 & BC.diff.1 > 0 & abs(AC.diff.1) < tol/2 # B in first; A beats C (or vice versa) by less than tol/2.
	first.round.bc = AB.diff.1 > 0 & AC.diff.1 > 0 & abs(BC.diff.1) < tol/2 # A in first; B beats C (or vice versa) by less than tol/2.
		
	piv.prob.list = list(
		"ab.ab" = mean(second.round & first.round.ab & AC.diff.2 > 0 & BC.diff.2 > 0),
		"ab.cb" = mean(second.round & first.round.ab & AC.diff.2 < 0 & BC.diff.2 > 0),
		"ab.ac" = mean(second.round & first.round.ab & AC.diff.2 > 0 & BC.diff.2 < 0),

		"ac.ac" = mean(second.round & first.round.ac & AB.diff.2 > 0 & BC.diff.2 < 0),
		"ac.bc" = mean(second.round & first.round.ac & AB.diff.2 < 0 & BC.diff.2 < 0),		
		"ac.ab" = mean(second.round & first.round.ac & AB.diff.2 > 0 & BC.diff.2 > 0),

		"bc.bc" = mean(second.round & first.round.bc & AB.diff.2 < 0 & AC.diff.2 < 0),
		"bc.ac" = mean(second.round & first.round.bc & AB.diff.2 > 0 & AC.diff.2 < 0),
		"bc.ba" = mean(second.round & first.round.bc & AB.diff.2 < 0 & AC.diff.2 > 0),

		"ab" = mean(second.round & second.round.AB & abs(AB.diff.2) < tol/2),  # why was this not divided by 2 before? 
		"bc" = mean(second.round & second.round.BC & abs(BC.diff.2) < tol/2),
		"ac" = mean(second.round & second.round.AC & abs(AC.diff.2) < tol/2),
		
		"ab.1" = mean(second.round & first.round.ab),
		"ac.1" = mean(second.round & first.round.ac),
		"bc.1" = mean(second.round & first.round.bc),
		
		"a.last" = mean(second.round & second.round.BC)*tol,
		"b.last" = mean(second.round & second.round.AC)*tol,
		"c.last" = mean(second.round & second.round.AB)*tol
						
		)
	
	as.list(unlist(piv.prob.list)/tol)  # scale by tol. This yields approximately n times the pivotal probability, so the real pivotal probability is this divided by n.
	
}


# plan: using simulation, check our ability to integrate along a line in a 3-way dirichlet 
prob.being.along.line.sim = function(alpha.vec, x, y, M = 100000, tol = .005, plot = F, plot.sample.rate = .01){
		
	stopifnot(length(alpha.vec) == 3)
	sims.01 = rdirichlet(n = M, alpha = alpha.vec)	# this is on the 0-1 scale.	
	sims = sims.01*(1 - x - y) # this is on the actual scale
	abs.diff.from.line = abs(sims[,2] - (x - y + sims[,1])) # this is on the actual scale
	
	if(plot){
		plot(c(0,1), c(0,1), type = "n", axes = F, xlab = "", ylab = "")
		abline(a = 1, b = -1)
		axis(1); axis(2)
		to.plot = rbinom(M, 1, plot.sample.rate) == 1
		to.plot.red = to.plot & abs.diff.from.line < tol/2
		to.plot.gray = to.plot & abs.diff.from.line > tol/2
		points(sims[to.plot.red,1], sims[to.plot.red,2], pch = 19, cex = .25, col = rgb(1, .5, .5, alpha = .5))
		points(sims[to.plot.gray,1], sims[to.plot.gray,2], pch = 19, cex = .25, col = rgb(.5, .5, .5, alpha = .5))
		abline(a = (x-y)/(1 - x - y), b = 1) # this is on the 0-1 scale
		exp.vec = (alpha.vec/sum(alpha.vec)) # on the 0-1 scale
		points(exp.vec[1], exp.vec[2], pch = 19, col = "red")
	}
	
	mean(abs.diff.from.line < tol/2)/tol  
	
}






### below were used in checking a bug

sim.pr.tie.between.1.and.2.between.a.and.b = function(alpha.vec, M = 10000, tol = .005, a = 0, b = 1/2){
	
	stopifnot(length(alpha.vec) == 3)
	sims = rdirichlet(n = M, alpha = alpha.vec)
	yes.vec = abs(sims[,1] - sims[,2]) < tol & sims[,1] > a & sims[,1] < b
	list(
		sims = sims,
		yes.vec = yes.vec,
		prob = mean(yes.vec)/tol
		)	
}


sim.pr.1.trails.2.by.less.than.x = function(alpha.vec, M = 10000, x = .005, a = 0, b = 1/2){
	
	stopifnot(length(alpha.vec) == 3)
	sims = rdirichlet(n = M, alpha = alpha.vec)
	
	margin.21 = sims[,2] - sims[,1]
	trailing.vec = margin.21 > 0 & margin.21 < x & sims[,1] > a & sims[,1] < b
	within.vec = abs(sims[,1] - sims[,2]) < x & sims[,1] > a & sims[,1] < b
	list(trailing.vec = trailing.vec, within.vec = within.vec, trailing.prob = mean(trailing.vec)/x, within.prob = mean(within.vec)/x, sims = sims)
}









# utilities for making a diagnostic ternary plot 

simplex.x <- function(x){
  if(sum(x) == 0){return(.5)}
  return( (x[2] + 0.5 * x[3]) / sum(x))
}
simplex.y <- function(x){
	if(sum(x) == 0){return(sqrt(.75)*(1/3))}  
	return( (sqrt(0.75) *  x[3]) / sum(x))
} 
simplex.z <- function(x){
	sqrt(.75)*x[1]
}



add.ternary.boundary = function(){
	bca.vertex = c(1, 0, 0)
	cba.vertex = c(0, 1, 0)
	bac.vertex = c(0, 0, 1)
	bac.v.x = simplex.x(bac.vertex); bac.v.y = simplex.y(bac.vertex)
	bca.v.x = simplex.x(bca.vertex); bca.v.y = simplex.y(bca.vertex)
	cba.v.x = simplex.x(cba.vertex); cba.v.y = simplex.y(cba.vertex) 

	lines(c(bac.v.x, cba.v.x), c(bac.v.y, cba.v.y))
	lines(c(bca.v.x, cba.v.x), c(bca.v.y, cba.v.y))
	lines(c(bac.v.x, bca.v.x), c(bac.v.y, bca.v.y))
}

add.ternary.point = function(point, col = "black", pch = 19, cex = 1){
	points(x = simplex.x(point), y = simplex.y(point), pch = pch, col = col, cex = cex)
}

add.ternary.lines = function(point.1, point.2, col = "black", lwd = 1, lty = 1){
	lines(x = c(simplex.x(point.1), simplex.x(point.2)), y = c(simplex.y(point.1), simplex.y(point.2)), col = col, lwd = lwd, lty = lty)
}


plot.simulation.results = function(a, b, tol, sims.out, cex = .2){
	out = sims.out
	par(mfrow = c(1,1))
	plot(c(0,1), c(0,sqrt(3/4)), type = "n", axes = F, xlab = "", ylab = "")
	add.ternary.boundary()
	to.plot.mat = out$sims[out$yes.vec,]
	for(i in 1:nrow(to.plot.mat)){
		add.ternary.point(to.plot.mat[i,], col = rgb(.3,.3,.3, alpha = .5), cex = cex)
	}
	bottom.point = c(a, a, 1-2*a)
	top.point = c(b, b, 1-2*b)
	add.ternary.lines(point.1 = bottom.point + c(-tol/2, tol/2, 0), point.2 = bottom.point + c(tol/2, -tol/2, 0))
	add.ternary.lines(point.1 = top.point + c(-tol/2, tol/2, 0), point.2 = top.point + c(tol/2, -tol/2, 0))
	add.ternary.lines(point.1 = bottom.point + c(-tol/2, tol/2, 0), point.2 = top.point + c(-tol/2, tol/2, 0))
	add.ternary.lines(point.1 = bottom.point + c(tol/2, -tol/2, 0), point.2 = top.point + c(tol/2, -tol/2, 0))	
}





