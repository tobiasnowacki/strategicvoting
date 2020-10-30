### 

EU_given_piv_probs_and_utility = function(piv_probs, utility, piv.events = c("AB", "AC", "BC", "AB.AB", "AB.AC", "AB.CB", "AC.AC", "AC.BC", "AC.AB", "BC.BC", "BC.AC", "BC.BA"), ballots = c("AB", "AC", "BA", "BC", "CA", "CB"), rule = "AV"){
	stopifnot(length(utility) == 3)
	stopifnot(class(piv_probs) == "list")
	if(rule == "plurality"){piv.events = piv.events[1:3]}
	piv.probs = unlist(piv_probs)[piv.events]
	outcome.mat = outcome_mat_from_utility(utility, piv.events, ballots, rule)
	as.vector(piv.probs) %*% outcome.mat
}

compare_ballots_event_by_event = function(piv_probs, utility, ballot = "AB", piv.events = c("AB", "AC", "BC", "AB.AB", "AB.AC", "AB.CB", "AC.AC", "AC.BC", "AC.AB", "BC.BC", "BC.AC", "BC.BA"), ballots = c("AB", "AC", "BA", "BC", "CA", "CB")){
	stopifnot(length(utility) == 3)
	stopifnot(class(piv_probs) == "list")
	piv.probs = unlist(piv_probs)[piv.events]
	outcome.mat = outcome_mat_from_utility(utility, piv.events, ballots)
	compare.mat = outcome.mat[,-which(ballots == ballot)]
	out = compare.mat - outcome.mat[,ballot]
	out*as.vector(piv.probs)
}

plot_expected_gain_relative_to_sincerity = function(piv_probs, utility, ballot_from, ballot_to){
	eu_diff_mat = compare_ballots_event_by_event(piv_probs, utility, ballot_from)
	y.locs = seq(nrow(eu_diff_mat), 1)
	x.vals = eu_diff_mat[,ballot_to]
	plot(range(x.vals), range(y.locs) + c(-.25, .25), type = "n", xlab = "Expected gain", ylab = "Pivotal event", axes = F)
	axis(1, at = c(-1000, 1000))
	abline(v = 0, lty = 3)
	for(i in 1:nrow(eu_diff_mat)){
		lines(c(0, x.vals[i]), rep(y.locs[i], 2), lwd = 1.5)
		points(x.vals[i], y.locs[i], pch = 19, cex = .75)
		text(0, y.locs[i] + .25, cex = .6, labels = rownames(eu_diff_mat)[i])
	}	
}

plot_expected_gain_loss_relative_to_sincerity = function(piv_probs, utility, ballot_from, ballot_to, main = NULL){
	eu_diff_mat = compare_ballots_event_by_event(piv_probs, utility, ballot_from, piv.events = piv.events, ballots = ballots)
	y.locs = seq(nrow(eu_diff_mat), 1)
	x.vals = eu_diff_mat[,ballot_to]/max(abs(eu_diff_mat[,ballot_to])) # taking scale out of it
	main.label = paste0(ifelse(is.null(main), "", main), ifelse(sum(x.vals) > 0, ": Yes", ": No"))
	plot(range(abs(x.vals)) + c(-.11, .1), range(y.locs) + c(-.25, .25), type = "n", xlab = "Expected gain/loss", ylab = "Pivotal event", axes = F, main = main.label)
	axis(1, at = c(-1000, 1000))
	abline(v = 0, lty = 3)
	for(i in 1:nrow(eu_diff_mat)){
		this.col = ifelse(x.vals[i] > 0, "green", ifelse(x.vals[i] < 0, "red", "gray"))
		lines(c(0, abs(x.vals[i])), rep(y.locs[i], 2), lwd = 1.5, col = this.col)
		points(abs(x.vals[i]), y.locs[i], pch = 19, cex = .75, col = this.col)
		text(-.02, y.locs[i], cex = .6, labels = rownames(eu_diff_mat)[i], adj = 1)
	}	
}

plot_expected_gain_loss_relative_to_sincerity_before_and_after = function(piv_probs_before, piv_probs_after, utility, ballot_from, ballot_to, main = NULL, y.inc = .2, piv.events = c("AB", "AC", "BC", "AB.AB", "AB.AC", "AB.CB", "AC.AC", "AC.BC", "AC.AB", "BC.BC", "BC.AC", "BC.BA"), ballots = c("AB", "AC", "BA", "BC", "CA", "CB")){
	eu_diff_mat_before = compare_ballots_event_by_event(piv_probs_before, utility, ballot_from, piv.events = piv.events, ballots = ballots)
	eu_diff_mat_after = compare_ballots_event_by_event(piv_probs_after, utility, ballot_from, piv.events = piv.events, ballots = ballots)
	y.locs = seq(nrow(eu_diff_mat_before), 1)
	x.vals_before = eu_diff_mat_before[,ballot_to]/max(abs(eu_diff_mat_before[,ballot_to])) # taking scale out of it
	x.vals_after = eu_diff_mat_after[,ballot_to]/max(abs(eu_diff_mat_before[,ballot_to])) # but keeping on same relative scale
	before_part = ifelse(sum(x.vals_before) > 0, ": Yes", ": No")
	after_part = ifelse(sum(x.vals_after) > 0, " -> Yes", " -> No")
	main.label = paste0(ifelse(is.null(main), paste0(ballot_from, " -> ", ballot_to), main), before_part, after_part)
	plot(range(abs(c(x.vals_before, x.vals_after))) + c(-.11, .1), range(y.locs) + c(-.25, .25), type = "n", xlab = "Expected gain/loss", ylab = "Pivotal event", axes = F, main = main.label)
	axis(1, at = c(-1000, 1000))
	abline(v = 0, lty = 3)
	for(i in 1:nrow(eu_diff_mat_before)){
		this.col = ifelse(x.vals_before[i] > 0, "green", ifelse(x.vals_before[i] < 0, "red", "gray"))
		lines(c(0, abs(x.vals_before[i])), rep(y.locs[i] + y.inc, 2), lwd = 1.5, col = this.col, lty = 2)
		lines(c(0, abs(x.vals_after[i])), rep(y.locs[i] - y.inc, 2), lwd = 1.5, col = this.col)
		points(abs(x.vals_before[i]), y.locs[i] + y.inc, pch = 21, cex = .75, col = this.col)
		points(abs(x.vals_after[i]), y.locs[i] - y.inc, pch = 19, cex = .75, col = this.col)
		text(-.02, y.locs[i], cex = .6, labels = rownames(eu_diff_mat_before)[i], adj = 1)
	}	
}



outcome_mat_from_utility = function(utility, piv.events = c("AB", "AC", "BC", "AB.AB", "AB.AC", "AB.CB", "AC.AC", "AC.BC", "AC.AB", "BC.BC", "BC.AC", "BC.BA"), ballots = c("AB", "AC", "BA", "BC", "CA", "CB"), rule = "AV"){
	uA = utility[1]; uB = utility[2]; uC = utility[3]
	if(rule == "AV"){
		out = rbind(
		# what is the outcome at each pivotal event and ballot
		c(rep(uA, 2), rep(uB, 2), uA, uB),  # AB
		c(rep(uA, 2), uA, uC, rep(uC, 2)),  # AC
		c(uB, uC, rep(uB, 2), rep(uC, 2)),  # BC
		c(rep(uA, 2), rep(uB, 2), rep((uA + uB)/2, 2)),  # AB.AB
		c(rep(uA, 2), rep(uC, 2), rep((uA + uC)/2, 2)),  # AB.AC
		c(rep(uC, 2), rep(uB, 2), rep((uB + uC)/2, 2)),  # AB.CB
		c(rep(uA, 2), rep((uA + uC)/2, 2), rep(uC, 2)),  # AC.AC
		c(rep(uB, 2), rep((uB + uC)/2, 2), rep(uC, 2)),  # AC.BC
		c(rep(uA, 2), rep((uA + uB)/2, 2), rep(uB, 2)),  # AC.AB
		c(rep((uB + uC)/2, 2), rep(uB, 2), rep(uC, 2)),  # BC.BC
		c(rep((uA + uC)/2, 2), rep(uA, 2), rep(uC, 2)),  # BC.AC
		c(rep((uA + uB)/2, 2), rep(uB, 2), rep(uA, 2))   # BC.BA
		)
	}else if(rule == "plurality"){
		piv.events = piv.events[1:3]
		ballots = c("A", "B", "C")
		out = rbind(
			c(uA, uB, (uA + uB)/2),
			c(uA, (uA + uC)/2, uC),
			c((uB + uC)/2, uB, uC)
		)
	}
	rownames(out) = piv.events
	colnames(out) = ballots
	out[piv.events, ballots]
}

P_mat_at_pivotal_events = function(piv_probs, piv.events = c("AB", "AC", "BC", "AB.AB", "AB.AC", "AB.CB", "AC.AC", "AC.BC", "AC.AB", "BC.BC", "BC.AC", "BC.BA"), ballots = c("AB", "AC", "BA", "BC", "CA", "CB"), candidates = c("A", "B", "C"), rule = "AV", normalize = FALSE){
	# we want to make a K x B matrix: probability of each candidate being elected given each ballot, focusing only on pivotal events.
	# in other words, the output of this function tells you the probability of each candidate being elected given each ballot, in expectation, conditioning on being pivotal if we normalize. 

	# we start by creating matrices that give you the probability of A, B, and C winning at each pivotal event as a function of the marginal ballot -- E rows and B columns, where E is number of pivotal events and B is number of ballots
	if(rule == "plurality"){
		pA.mat = rbind(
			c(1,0,.5),
			c(1,.5,0),
			c(0,0,0))
		pB.mat = rbind(
			c(0,1,.5),
			c(0,0,0),
			c(.5,1,0))
		pC.mat = 1 - pA.mat - pB.mat
		# the piv.probs vector is E by 1.
		piv.events = piv.events[1:3]
	}else if(rule == "AV"){
		pA.mat = rbind(
		c(rep(1, 2), rep(0, 2), 1, 0),  # AB
		c(rep(1, 2), 1, 0, rep(0, 2)),  # AC
		c(rep(0, 6)), # BC
		c(rep(1, 2), rep(0, 2), rep(1/2, 2)),  # AB.AB
		c(rep(1, 2), rep(0, 2), rep(1/2, 2)),  # AB.AC
		c(rep(0, 6)),  # AB.CB
		c(rep(1, 2), rep(1/2, 2), rep(0, 2)),  # AC.AC
		c(rep(0, 6)),  # AC.BC
		c(rep(1, 2), rep(1/2, 2), rep(0, 2)),  # AC.AB
		c(rep(0, 6)),  # BC.BC
		c(rep(1/2, 2), rep(1, 2), rep(0, 2)),  # BC.AC
		c(rep(1/2, 2), rep(0, 2), rep(1, 2))   # BC.BA
		)
		pB.mat = rbind(
		c(rep(0, 2), rep(1, 2), 0, 1),  # AB
		c(rep(0, 2), 0, 0, rep(0, 2)),  # AC
		c(1, 0, rep(1, 2), rep(0, 2)),  # BC
		c(rep(0, 2), rep(1, 2), rep(1/2, 2)),  # AB.AB
		c(rep(0, 6)),  # AB.AC
		c(rep(0, 2), rep(1, 2), rep(1/2, 2)),  # AB.CB
		c(rep(0, 6)),  # AC.AC
		c(rep(1, 2), rep(1/2, 2), rep(0, 2)),  # AC.BC
		c(rep(0, 2), rep(1/2, 2), rep(1, 2)),  # AC.AB
		c(rep(1/2, 2), rep(1, 2), rep(0, 2)),  # BC.BC
		c(rep(0, 6)),  # BC.AC
		c(rep(1/2, 2), rep(1, 2), rep(0, 2))   # BC.BA
		)
		pC.mat = 1 - pA.mat - pB.mat
	}
	# then we combine this with the pivotal probabilities, which we organize in an E x B matrix
	piv.probs = matrix(unlist(piv_probs)[piv.events], nrow = length(piv.events), ncol = length(ballots))
	# the output should be K x B 
	out = rbind(apply(piv.probs*pA.mat, 2, sum), apply(piv.probs*pB.mat, 2, sum), apply(piv.probs*pC.mat, 2, sum)) # this must be dreadfully inefficient -- can't figure it out		
	# label the matrix.
	colnames(out) = ballots
	rownames(out) = candidates
	# normalize -- this makes it conditional on pivotal events.
	if(normalize == TRUE){
	  for(j in 1:ncol(out)){
	    out[,j] = out[,j]/sum(out[,j])
	  }
	} # otherwise it's just the raw pivotal probabilities.
	out
}

# this function gets the P_mat created above and combines it with a matrix of utilities. 
EU_mat_given_piv_probs_and_u_mat = function(piv_probs, u_mat, piv.events = c("AB", "AC", "BC", "AB.AB", "AB.AC", "AB.CB", "AC.AC", "AC.BC", "AC.AB", "BC.BC", "BC.AC", "BC.BA"), ballots = c("AB", "AC", "BA", "BC", "CA", "CB"), candidates = c("A", "B", "C"), rule = "AV"){
	# K x B
	P_mat = P_mat_at_pivotal_events(piv_probs, piv.events, ballots, candidates, rule = rule)  # more arguments needed
	stopifnot(ncol(u_mat) == nrow(P_mat))
	# n x B
	out = u_mat%*%P_mat
	rownames(out) = rownames(u_mat)
	colnames(out) = ballots
	out
}

# OK so now what. a survey is taken. this gives us an indication of the distribution of preference types. 
# now we want to know what each type does in response if he is (a) level-0, (b) level-1 with the belief that $p$ other voters are level-0. 

# as the basis we have this function
# we want to adjust the v.vec in response to the optimal voting strategy. 
simulation_thing = function(v.vec, s.vec, u_mat = NULL, rule = "AV", piv.events = NULL, ballots = NULL, candidates = c("A", "B", "C"), p = .2){
	if(rule == "AV"){
		if(is.null(piv.events)){piv.events = c("AB", "AC", "BC", "AB.AB", "AB.AC", "AB.CB", "AC.AC", "AC.BC", "AC.AB", "BC.BC", "BC.AC", "BC.BA")}
		if(is.null(ballots)){ballots = c("AB", "AC", "BA", "BC", "CA", "CB")}
		piv_probs.1 = av.pivotal.event.probs.general(v.vec, s.vec = s.vec)
	}else if(rule == "plurality"){
		if(is.null(piv.events)){piv.events = c("AB", "AC", "BC")}
		if(is.null(ballots)){ballots = c("A", "B", "C")}
		if(length(v.vec) == 9){v.vec = c(sum(v.vec[c(1:2, 7)]), sum(v.vec[c(3:4, 8)]), sum(v.vec[c(5:6, 9)]))} # combining
		if(length(v.vec) == 6){v.vec = c(sum(v.vec[c(1:2)]), sum(v.vec[c(3:4)]), sum(v.vec[c(5:6)]))} # combining
		piv_probs.1 = plurality.pivotal.probabilities(v.vec, s.vec[1])
	}
	if(is.null(u_mat)){
		u_mat = matrix(data = c(c(1,.5,0), c(1,0,.5), c(.5, 1, 0), c(0,1,.5), c(.5,0,1), c(0,.5,1)), byrow = T, ncol = 3, dimnames = list(c("AB", "AC", "BA", "BC", "CA", "CB"), c("A", "B", "C")))
	}
	eu_mat.1 = EU_mat_given_piv_probs_and_u_mat(piv_probs.1, u_mat, piv.events, ballots, candidates, rule) # this tells us the utility from ballot from each type. 
	new_v_vec = adjust_v_vec(v.vec, eu_mat.1, p)
	if(rule == "AV"){
		piv_probs.2 = av.pivotal.event.probs.general(new_v_vec, s.vec = s.vec)
	}else if(rule == "plurality"){
		piv_probs.2 = plurality.pivotal.probabilities(new_v_vec, s.vec[1])
	}
	eu_mat.2 = EU_mat_given_piv_probs_and_u_mat(piv_probs.2, u_mat, piv.events, ballots, candidates, rule) # this tells us the utility from ballot from each type. 
	list(original_v_vec = v.vec, updated_v_vec = new_v_vec, p = p, eu_mat.1 = eu_mat.1, eu_mat.2 = eu_mat.2)
}

# Next: better way to summarize the eu_mats and what they say about tactical voting in this setting. 


which.is.max = function(vec){
	out = which(vec == max(vec))
	stopifnot(length(out) == 1)
	out
}

adjust_v_vec = function(v.vec, eu_mat, p){
	# u_mat is a bunch of unique types in alphabetical order, eu_mat is also based on types. so we determine what the optimal vote is for each, and adjust v.vec accordingly.
	ballots = colnames(eu_mat)
	opt_vote = ballots[apply(eu_mat, 1, which.is.max)]
	fully_adjusted_v_vec = rep(0, length(v.vec))
	for(i in 1:length(v.vec)){
		this_count = table(opt_vote)[ballots[i]]
		if(!is.na(this_count)){
			fully_adjusted_v_vec[i] = this_count*v.vec[i]			
		}
	}
	fully_adjusted_v_vec = fully_adjusted_v_vec/sum(fully_adjusted_v_vec)
	return((1-p)*v.vec + p*fully_adjusted_v_vec)
}
  