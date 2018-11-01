# code for outputting pivotal probabilities from simulations given Dirichlet belief parameters 

av.pivotal.probs.from.sims = function(v.vec = NULL, s.vec = NULL, M = 100000, tol = .005){
		
	# expecting shares for AB, AC, BA, BC, CA, CB, AX, BX, CX	
	stopifnot(length(v.vec) == 9)
	stopifnot(length(s.vec) == 4)	
	if(length(unique(s.vec)) == 1){
		sims = rdirichlet(n = M, alpha = v.vec*s.vec[1])	
		# first preference vote shares for each candidate	
		A.fp.share = apply(sims[,c(1,2,7)], 1, sum)
		B.fp.share = apply(sims[,c(3,4,8)], 1, sum)
		C.fp.share = apply(sims[,c(5,6,9)], 1, sum)		
	}	
	else{
		fp.alpha = s.vec[1]*c(sum(v.vec[c(1,2,7)]), sum(v.vec[c(3,4,8)]), sum(v.vec[c(5,6,9)]))
		fp.sims = rdirichlet(n = M, alpha = fp.alpha)
		# first preference vote shares for each candidate
		A.fp.share = fp.sims[,1]
		B.fp.share = fp.sims[,2]
		C.fp.share = fp.sims[,3]
		# second-preference shares 
		ballot.shares.A = A.fp.share*rdirichlet(n = M, alpha = s.vec[2]*v.vec[c(1,2,7)])
		ballot.shares.B = B.fp.share*rdirichlet(n = M, alpha = s.vec[3]*v.vec[c(3,4,8)])
		ballot.shares.C = C.fp.share*rdirichlet(n = M, alpha = s.vec[4]*v.vec[c(5,6,9)])
		sims = cbind(ballot.shares.A[,1:2], ballot.shares.B[,1:2], ballot.shares.C[,1:2], ballot.shares.A[,3], ballot.shares.B[,3], ballot.shares.C[,3])
	}
		
	# get pivotal probabilities from the simulations
	
	# identify cases where there is a second round -- if not, you can't be pivotal (assuming all ballot orders get some support).  
	second.round = apply(cbind(A.fp.share, B.fp.share, C.fp.share), 1, max) < .5    # used to say "- tol" here, but this seems wrong.  
			
	# FP margins
	AB.diff.1 = A.fp.share - B.fp.share
	AC.diff.1 = A.fp.share - C.fp.share
	BC.diff.1 = B.fp.share - C.fp.share
	
	# second-round matchups 
	second.round.AB = AC.diff.1 > 0 & BC.diff.1 > 0 
	second.round.BC = AB.diff.1 < 0 & AC.diff.1 < 0 
	second.round.AC = AB.diff.1 > 0 & BC.diff.1 < 0 

	# two-party preferences
	share.preferring.A.to.B = A.fp.share + sims[,5]
	share.preferring.A.to.C = A.fp.share + sims[,3]
	share.preferring.B.to.C = B.fp.share + sims[,1]
	
	# majority margins for each pairing
	AB.diff.2 =  A.fp.share + sims[,5] - (B.fp.share + sims[,6]) 
	AC.diff.2 =  A.fp.share + sims[,3] - (C.fp.share + sims[,4]) 
	BC.diff.2 =  B.fp.share + sims[,1] - (C.fp.share + sims[,2]) 
	
	# first round ties for second: I have checked this and I think it's fine
	first.round.ab.tie = AC.diff.1 < 0 & BC.diff.1 < 0 & AB.diff.1 > -tol/2 & AB.diff.1 < tol/2 # C in first; A beats B (or vice versa) by less than tol/2 in FP votes. 
	first.round.ac.tie = AB.diff.1 < 0 & BC.diff.1 > 0 & AC.diff.1 > -tol/2 & AC.diff.1 < tol/2 # B in first; A beats C (or vice versa) by less than tol/2.
	first.round.bc.tie = AB.diff.1 > 0 & AC.diff.1 > 0 & BC.diff.1 > -tol/2 & BC.diff.1 < tol/2 # A in first; B beats C (or vice versa) by less than tol/2.
		
	piv.prob.list = list(
		# first-round pivotal events
		"AB.AB" = mean(second.round & first.round.ab.tie & AC.diff.2 > 0 & BC.diff.2 > 0),
		"AB.AC" = mean(second.round & first.round.ab.tie & AC.diff.2 > 0 & BC.diff.2 < 0),
		"AB.CB" = mean(second.round & first.round.ab.tie & AC.diff.2 < 0 & BC.diff.2 > 0),

		"AC.AC" = mean(second.round & first.round.ac.tie & AB.diff.2 > 0 & BC.diff.2 < 0),
		"AC.BC" = mean(second.round & first.round.ac.tie & AB.diff.2 < 0 & BC.diff.2 < 0),		
		"AC.AB" = mean(second.round & first.round.ac.tie & AB.diff.2 > 0 & BC.diff.2 > 0),

		"BC.BC" = mean(second.round & first.round.bc.tie & AB.diff.2 < 0 & AC.diff.2 < 0),
		"BC.AC" = mean(second.round & first.round.bc.tie & AB.diff.2 > 0 & AC.diff.2 < 0),
		"BC.BA" = mean(second.round & first.round.bc.tie & AB.diff.2 < 0 & AC.diff.2 > 0),

		# second-round pivotal events: I have checked this and it is fine. 
		"AB" = mean(second.round & second.round.AB & share.preferring.A.to.B > (1 - sims[,9] - tol/2)/2 & share.preferring.A.to.B < (1 - sims[,9] + tol/2)/2), 
		"AC" = mean(second.round & second.round.AC & share.preferring.A.to.C > (1 - sims[,8] - tol/2)/2 & share.preferring.A.to.C < (1 - sims[,8] + tol/2)/2),
		"BC`" = mean(second.round & second.round.BC & share.preferring.B.to.C > (1 - sims[,7] - tol/2)/2 & share.preferring.B.to.C < (1 - sims[,7] + tol/2)/2)
								
		)
	
	as.list(unlist(piv.prob.list)/tol)  # scale by tol. This yields approximately n times the pivotal probability, so the real pivotal probability is this divided by n.
	
}




### unused diagnostics 

if(F){
####################
# diagnostic events#
####################

list(
# first-round ties for second
"ab.1" = mean(second.round & first.round.ab.tie),
"ac.1" = mean(second.round & first.round.ac.tie),
"bc.1" = mean(second.round & first.round.bc.tie),

# TPP ties -- conditional on continuing to second round
"ab.2" = mean(abs(AB.diff.2) < tol/2),
"ac.2" = mean(abs(AC.diff.2) < tol/2),
"bc.2" = mean(abs(BC.diff.2) < tol/2),

# TPP ties another way
"ab.2a" = mean(AB.diff.2 < 0 & AB.diff.2 > -tol),
"ac.2a" = mean(AC.diff.2 < 0 & AC.diff.2 > -tol),
"bc.2a" = mean(BC.diff.2 < 0 & BC.diff.2 > -tol),

# TPP ties yet another way

"ab.2b" = mean(A.fp.share + sims[,5] > .5 - tol & A.fp.share + sims[,5] < .5),
"ac.2b" = mean(A.fp.share + sims[,3] > .5 - tol & A.fp.share + sims[,3] < .5),
"bc.2b" = mean(B.fp.share + sims[,1] > .5 - tol & B.fp.share + sims[,1] < .5),

# each party in last, and continuing to second round
"a.last" = mean(second.round & second.round.BC)*tol,
"b.last" = mean(second.round & second.round.AC)*tol,
"c.last" = mean(second.round & second.round.AB)*tol
)
}
