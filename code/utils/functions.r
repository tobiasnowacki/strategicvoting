##################################################
## Project: Strategic Voting in RCV
## Script purpose: Useful functions for project
## Date: 30/10/2018
## Author:
##################################################


opt_vote <- function(utility_df, p_list, type = "rcv"){
# Calculate Optimal Votes
## Takes input: utility dataframe; list of lists of pivotal probabilities.
## Gives output: matrix of expected utilities of each vote
	if(type == "rcv"){
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
		vote_utils <- matrix(vote_utils, ncol = 6, byrow = T)
		return(vote_utils)
	}
	if(type == "plur"){
		p_list <- c(p_list$AB, p_list$AC, p_list$BC)

		w_plurality <- as.data.frame(matrix(data = c(1, 0, 0, 0, 1, 0, .5, .5, 0,
                               1, 0, 0, .5, 0, .5, 0, 0, 1,
                               0, .5, .5, 0, 1, 0, 0, 0, 1),
                      byrow = T, nrow = 3))
		w_test <- w_plurality[ , rep(1:9, nrow(utility_df))]
		u_test <- utility_df[, rep(1:3, 3)]
		u_test <- as.vector(t(u_test))

		test <- matrix((mapply(`*`,w_test, u_test, SIMPLIFY=TRUE)), nrow = 3)
		cond_utils_sum <- test[, c(T, F, F)] + test[, c(F, T, F)] +
		test[, c(F, F, T)]
		cond_utils_sum_event <- cond_utils_sum * p_list

		vote_utils <- colSums(cond_utils_sum_event)
		vote_utils <- matrix(vote_utils, ncol = 3, byrow = T)
		return(vote_utils)
	}
}

return_sv_prop <- function(v_vec, util_df, s_breaks, full_mat = FALSE){
	# Function that takes a ballot profile, a dataframe of utilities, and a list of s values
	#	Returns: a data.frame with voters by vote type for levels of s

	v_vec <- as.numeric(v_vec)
	sin_vec <- apply(util_df, 1, sin_vote_scalar)

	# RCV exercise
	p_list <- lapply(s_breaks, function(x) av.pivotal.event.probs.general(v_vec, rep(x, 4)))
	eu_list <- lapply(p_list, function(x) opt_vote(util_df, x, type = "rcv"))
	opt_vec_list <- lapply(eu_list, function(x) opt_vote_scalar(x))
	opt_vote_prop <- lapply(opt_vec_list, function(x) sum_opt_votes(sin_vec, x, type = "rcv"))
	prop_df <- do.call(rbind, opt_vote_prop)
	prop_df <- as.data.frame(prop_df)

	# Plurality exercise
	v_vec_three <- c(v_vec[1] + v_vec[2], v_vec[3] + v_vec[4], v_vec[5] + v_vec[6])
	p_list_plur <- lapply(s_breaks, function(x) plurality.pivotal.probabilities(v_vec_three, x))
	eu_list_plur <- lapply(p_list_plur, function(x) opt_vote(util_df, x, type = "plur"))

	opt_vec_list_plur <- lapply(eu_list_plur, function(x) opt_vote_scalar(x))
	opt_vote_prop_plur <- lapply(opt_vec_list_plur, function(x) sum_opt_votes(sin_vec, x, type = "plur"))

	prop_df_plur <- do.call(rbind, opt_vote_prop_plur)

	if(full_mat == TRUE){
		return(list(opt_vec_list, opt_vec_list_plur))
	}

	# Add plurality to DF
	prop_df$plurality_first <- as.numeric(prop_df_plur[, 1])
	prop_df$plurality_second <- as.numeric(prop_df_plur[, 2] )
	prop_df[, 1:3] <- t(apply(prop_df[, 1:3], 1, function(x) x / sum(x)))
	prop_df[, 4:5] <- t(apply(prop_df[, 4:5], 1, function(x) x / sum(x)))
	prop_df$s <- unlist(s_breaks)
	names(prop_df)[1:3] <- c("first", "second", "third")
	return(prop_df)

}

# Function: take ballot profile, DF of utilities, and returns vector of tactical incentives at different
return_sv_tau <- function(v_vec, util_df, s_breaks){
	v_vec <- as.numeric(v_vec)
	sin_vec <- apply(util_df, 1, sin_vote_scalar)
	
	# RCV part
	p_list <- lapply(s_breaks, function(x) av.pivotal.event.probs.general(v_vec, rep(x, 4)))
	eu_list <- lapply(p_list, function(x) opt_vote(util_df, x, type = "rcv"))
	tau_list <- lapply(eu_list, function(x) calculate_tau(x, sin_vec))

	# Plurality part
	# Do the same for plurality
	sin_vec_plur <- sin_vote_plur_transform(sin_vec)
	v_vec_three <- c(v_vec[1] + v_vec[2], v_vec[3] + v_vec[4], v_vec[5] + v_vec[6])
	p_list_plur <- lapply(s_breaks, function(x) plurality.pivotal.probabilities(v_vec_three, x))
	eu_list_plur <- lapply(p_list_plur, function(x) opt_vote(util_df, x, type = "plur"))
	tau_list_plur <- lapply(eu_list_plur, function(x) calculate_tau(x, sin_vec_plur))

	# Render optimal votes
	opt_list <- lapply(eu_list, function(x) opt_vote_scalar(x))
	opt_vec <- unlist(opt_list)

	opt_list_plur <- lapply(eu_list_plur, function(x) opt_vote_scalar(x))
	opt_vec_plur <- unlist(opt_list_plur)

	# Merge into one big data-frame
	n <- length(s_breaks)
	# return(tau_list[[1]][1:3])
	tau_vec_rcv <- unlist(tau_list)
	# return(tau_vec_rcv)
	tau_vec_plur <- unlist(tau_list_plur)
	out_df <- as.data.frame(cbind(rep(sin_vec, n), rep(sin_vec_plur, n), tau_vec_rcv, tau_vec_plur, opt_vec, opt_vec_plur, rep(unlist(s_breaks), each = nrow(util_df))))
	#out_df <- apply(out_df, 2, function(x) unlist(x))
	names(out_df) <- c("sin_rcv", "sin_plur", "tau_rcv", "tau_plur", "opt_rcv", "opt_plur", "s")
	return(out_df)

}

# Transforms optimal vote scalar from RCV to plurality (aggregates into three parties)
sin_vote_plur_transform <- function(x){
	x[x %in% c(1, 2)] <- 1
	x[x %in% c(3, 4)] <- 2
	x[x %in% c(5, 6)] <- 3
	return(x)
}


# Function: Take EU df and calculate optimal vote

opt_vote_scalar <- function(eu_df){
	# Function that takes a dataframe of expected utilities (6 or 3 times no. of respondents)
	# 	Returns: vector of optimal votes (as scalars)
	apply(eu_df, 1, function(x) which(x == max(x))) 
}

sin_vote_scalar <- function(util){
	# Input: DF of utilities
	# Output: Vector of sincere preferences (as scalars) -- under AV
	stopifnot(length(util) == 3)
	max <- which(util == max(util))
	min <- which(util == min(util))

	if (max == 1){
	if (min == 3){
	  out <- 1
	}
	if (min == 2){
	  out <- 2
	}
	}
	if (max == 2){
	if (min == 3){
	  out <- 3
	}
	if (min == 1){
	  out <- 4
	}
	}
	if (max == 3){
	if (min == 2){
	  out <- 5
	}
	if (min == 1){
	  out <- 6
	}
	}
	return(out)
}

calculate_tau <- function(eu_df, sin_vec){
	df <- cbind(eu_df, sin_vec)
	df_ncol <- ncol(df)
	tau <- apply(df, 1, function(x){
		sin <- x[df_ncol]
		eus <- x[-df_ncol]
		tau <- max(eus[-sin]) - max(eus[sin])
		return(tau)
	})
	return(tau)
}

sum_opt_votes <- function(sin_vec, opt_vec, type = "rcv"){
	# Input: vector of sincere preferences, vector of opt. votes
	# Output: Table
	if(type == "rcv"){
		sin_vec <- factor(sin_vec, levels = 1:6)
		opt_vec <- factor(opt_vec, levels = 1:6)
		tab <- table(sin_vec, opt_vec)
		out <- rbind(c(tab[1, 1], tab[1, 3], tab[1, 5]),
					 c(tab[2, 2], tab[2, 5], tab[2, 3]),
					 c(tab[3, 3], tab[3, 1], tab[3, 6]),
					 c(tab[4, 4], tab[4, 6], tab[4, 1]),
					 c(tab[5, 5], tab[5, 2], tab[5, 4]),
					 c(tab[6, 6], tab[6, 4], tab[6, 2]))
		return(colSums(out))
	}
	if(type == "plur"){
		sin_vec <- factor(sin_vec, levels = 1:6)
		opt_vec <- factor(opt_vec, levels = 1:3)
		tab <- table(sin_vec, opt_vec)
		out <- rbind(c(tab[1, 1], tab[1, 2]),
					 c(tab[2, 1], tab[2, 3]),
					 c(tab[3, 2], tab[3, 1]),
					 c(tab[4, 2], tab[4, 3]),
					 c(tab[5, 3], tab[5, 1]),
					 c(tab[6, 3], tab[6, 2]))
		return(colSums(out))
	}

}

qq_function <- function(const_bp_no_trunc, utils, s){
	# Takes dataframe of ballot proportions, dataframe of utilities, and level of uncertainty
	# Returns dataframe of coordiantes for qq-plot, grouped by constituency
	const_taus <- list()
	const_taus_qq <- list()
	for(i in 1:nrow(const_bp_no_trunc)){
			# print(const_bp_no_trunc[i, ])
			v_vec <- const_bp_no_trunc[i, 2:10]
			const_taus[[i]] <- return_sv_tau(as.numeric(v_vec), utils, s)
			const_taus_qq[[i]] <- as.data.frame(qqplot(x = unlist(const_taus[[i]]$tau_vec_plur[const_taus[[i]]$V4 == s]), 
				y = unlist(const_taus[[i]]$tau_vec_rcv[const_taus[[i]]$V4 == s]), plot.it = FALSE))
	}
	const_taus_qq_df <- do.call(rbind, const_taus_qq)
	const_taus_qq_df$const <- rep(const_bp_no_trunc$district, each = nrow(aes_utils))
	return(const_taus_qq_df)
}

# Andy's function(s) (for reference):

# outcome_mat_from_utility = function(utility, piv.events = c("AB", "AC", "BC", "AB.AB", "AB.AC", "AB.CB", "AC.AC", "AC.BC", "AC.AB", "BC.BC", "BC.AC", "BC.BA"), ballots = c("AB", "AC", "BA", "BC", "CA", "CB"), rule = "AV"){
# 	uA = utility[1]; uB = utility[2]; uC = utility[3]
# 	if(rule == "AV"){
# 		out = rbind(
# 		# what is the outcome at each pivotal event and ballot
# 		c(rep(uA, 2), rep(uB, 2), uA, uB),  # AB
# 		c(rep(uA, 2), uA, uC, rep(uC, 2)),  # AC
# 		c(uB, uC, rep(uB, 2), rep(uC, 2)),  # BC
# 		c(rep(uA, 2), rep(uB, 2), rep((uA + uB)/2, 2)),  # AB.AB
# 		c(rep(uA, 2), rep(uC, 2), rep((uA + uC)/2, 2)),  # AB.AC
# 		c(rep(uC, 2), rep(uB, 2), rep((uB + uC)/2, 2)),  # AB.CB
# 		c(rep(uA, 2), rep((uA + uC)/2, 2), rep(uC, 2)),  # AC.AC
# 		c(rep(uB, 2), rep((uB + uC)/2, 2), rep(uC, 2)),  # AC.BC
# 		c(rep(uA, 2), rep((uA + uB)/2, 2), rep(uB, 2)),  # AC.AB
# 		c(rep((uB + uC)/2, 2), rep(uB, 2), rep(uC, 2)),  # BC.BC
# 		c(rep((uA + uC)/2, 2), rep(uA, 2), rep(uC, 2)),  # BC.AC
# 		c(rep((uA + uB)/2, 2), rep(uB, 2), rep(uA, 2))   # BC.BA
# 		)
# 	}else if(rule == "plurality"){
# 		piv.events = piv.events[1:3]
# 		ballots = c("A", "B", "C")
# 		out = rbind(
# 			c(uA, uB, (uA + uB)/2),
# 			c(uA, (uA + uC)/2, uC),
# 			c((uB + uC)/2, uB, uC)
# 		)
# 	}
# 	rownames(out) = piv.events
# 	colnames(out) = ballots
# 	return(out)
# }

# EU_given_piv_probs_and_utility = function(piv_probs, utility, piv.events = c("AB", "AC", "BC", "AB.AB", "AB.AC", "AB.CB", "AC.AC", "AC.BC", "AC.AB", "BC.BC", "BC.AC", "BC.BA"), ballots = c("AB", "AC", "BA", "BC", "CA", "CB"), rule = "AV"){
# 	stopifnot(length(utility) == 3)
# 	stopifnot(class(piv_probs) == "list")
# 	if(rule == "plurality"){piv.events = piv.events[1:3]}
# 	piv.probs = unlist(piv_probs)[piv.events]
# 	outcome.mat = outcome_mat_from_utility(utility, piv.events, ballots, rule)
# 	as.vector(piv.probs) %*% outcome.mat
# }


# Toy example to check that functions return the same output

# v_vec <- c(0.02, 0.08, 0.35, 0.2, 0.175, 0.175, 0, 0, 0)
# s_vec <- rep(85, 4)
# p_probs <- av.pivotal.event.probs.general(v_vec, s_vec)
# u_df <- data.frame(A = 1, B = 0.5, C = 0)

# opt_vote(u_df, p_probs)
# EU_given_piv_probs_and_utility(p_probs, c(1, 0.5, 0))
# # Replication successful!

# system.time({opt_vote(u_df, p_probs)})
# system.time({EU_given_piv_probs_and_utility(p_probs, c(1, 0.5, 0))})
