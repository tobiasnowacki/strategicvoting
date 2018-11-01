##################################################
## Project: Strategic Voting in RCV
## Script purpose: Useful functions for project
## Date: 30/10/2018
## Author:
##################################################

##################################################
# (Delete later, only necessary for toy example) #
setwd("~/github/strategicvoting/utils")
source("av_pivotal_probs_analytical_general_v2.r")
##################################################

opt_vote <- function(utility_df, p_list, type = "rcv"){
# Calculate Optimal Votes
## Takes input: utility dataframe; list of lists of pivotal probabilities.
## Gives output: matrix of expected utilities of each vote

	A <- c(1, 0, 0)
	B <- c(0, 1, 0)
	C <- c(0, 0, 1)
	AB <- c(0.5, 0.5, 0)
	AC <- c(0.5, 0, 0.5)
	BC <- c(0, 0.5, 0.5)
	w_av <- as.data.frame(matrix(data = c(A, A, B, B, A, B,                    # AB
	                                      A, A, A, C, C, C,                    # AC
	                                      B, C, B, B, C, C,                    # BC
	                                      A, A, B, B, AB, AB,                  # AB.AB
	                                      C, C, B, B, BC, BC,                  # AB.CB
	                                      A, A, C, C, AC, AC,                  # AB.AC
	                                      A, A, AB, AB, B, B,                  # AC.AB
	                                      B, B, BC, BC, C, C,                  # AC.BC
	                                      A, A, AC, AC, C, C,                  # AC.AC
	                                      AC, AC, A, A, C, C,                  # BC.AC
	                                      AB, AB, B, B, A, A,                  # BC.BA
	                                      BC, BC, B, B, C, C),                 # BC.BC
	                          byrow = TRUE, ncol = 6*3, nrow = 12))

	p_list <- c(p_list$AB, p_list$AC, p_list$BC,
              p_list$AB.AB, p_list$AB.CB, p_list$AB.AC,
              p_list$AC.AB, p_list$AC.BC, p_list$AC.AC,
              p_list$BC.AC, p_list$BC.BA, p_list$BC.BC)
  
	w_test <- w_av[, rep(1:18, nrow(utility_df))]
	u_test <- utility_df[, rep(1:3, 6)]
	u_test <- as.vector(t(u_test))

	test <- matrix((mapply(`*`,w_test,u_test,SIMPLIFY=TRUE)), nrow = 12)

	cond_utils_sum <- test[, c(T, F, F)] + test[, c(F, T, F)] +
	test[, c(F, F, T)]
	cond_utils_sum_event <- cond_utils_sum * p_list
	vote_utils <- colSums(cond_utils_sum_event)
	return(vote_utils)
}

# Andy's function(s) (for reference):

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
	return(out)
}

EU_given_piv_probs_and_utility = function(piv_probs, utility, piv.events = c("AB", "AC", "BC", "AB.AB", "AB.AC", "AB.CB", "AC.AC", "AC.BC", "AC.AB", "BC.BC", "BC.AC", "BC.BA"), ballots = c("AB", "AC", "BA", "BC", "CA", "CB"), rule = "AV"){
	stopifnot(length(utility) == 3)
	stopifnot(class(piv_probs) == "list")
	if(rule == "plurality"){piv.events = piv.events[1:3]}
	piv.probs = unlist(piv_probs)[piv.events]
	outcome.mat = outcome_mat_from_utility(utility, piv.events, ballots, rule)
	as.vector(piv.probs) %*% outcome.mat
}


# Toy example to check that functions return the same output

v_vec <- c(0.02, 0.08, 0.35, 0.2, 0.175, 0.175, 0, 0, 0)
s_vec <- rep(85, 4)
p_probs <- av.pivotal.event.probs.general(v_vec, s_vec)
u_df <- data.frame(A = 1, B = 0.5, C = 0)

opt_vote(u_df, p_probs)
EU_given_piv_probs_and_utility(p_probs, c(1, 0.5, 0))
# Replication successful!

system.time({opt_vote(u_df, p_probs)})
system.time({EU_given_piv_probs_and_utility(p_probs, c(1, 0.5, 0))})
