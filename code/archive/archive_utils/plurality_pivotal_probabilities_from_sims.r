# three-candidate plurality pivotal probs from simulations. 

plurality.pivotal.probs.from.sims = function(v.vec, s, M = 100000, tol = .005){
		
	# expecting shares for A, B, C	
	stopifnot(length(v.vec) == 3)
	stopifnot(length(s) == 1)	
	sims = rdirichlet(n = M, alpha = v.vec*s)
	
	list(
		"AB" = mean(sims[,1] > sims[,3] & sims[,2] > sims[,3] & abs(sims[,1] - sims[,2]) < tol/2)/tol,
		"AB" = mean(sims[,1] > sims[,2] & sims[,3] > sims[,2] & abs(sims[,1] - sims[,3]) < tol/2)/tol,
		"BC" = mean(sims[,2] > sims[,1] & sims[,3] > sims[,1] & abs(sims[,2] - sims[,3]) < tol/2)/tol
		)
		
}
