##################################################
## Project: Strategic Voting in RCV
## Script purpose: Useful functions for project
## Date: 30/10/2018
## Author:
##################################################
library(questionr)

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
		u_df <- matrix(rep(u_test, 12), nrow = 12, byrow = T)

		test <- w_test * u_df

		cond_utils_sum <- test[, c(T, F, F)] + test[, c(F, T, F)] +
		test[, c(F, F, T)]
		#return(cond_utils_sum)
		vote_utils <- t(as.matrix(p_list)) %*% as.matrix(cond_utils_sum)
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

		u_df <- matrix(rep(u_test, 3), nrow = 3, byrow = T)
		test <- w_test * u_df

		cond_utils_sum <- test[, c(T, F, F)] + test[, c(F, T, F)] +
		test[, c(F, F, T)]
		#return(cond_utils_sum)
		vote_utils <- t(as.matrix(p_list)) %*% as.matrix(cond_utils_sum)
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

sv_prop <- function(tau_obj, weights = "fill"){
	tau_list <- split(tau_obj, tau_obj$s)
  if(weights == "fill") {
    weights <- rep(1, length(tau_list[[1]]$sin_rcv))
  }
	
	# RCV proportions
	opt_vote_prop <- lapply(tau_list, function(x) sum_opt_votes(x$sin_rcv, x$opt_rcv, type = "rcv", weights = weights))
	prop_df_rcv <- do.call(rbind, opt_vote_prop)

	opt_vote_prop_plur <- lapply(tau_list, function(x) sum_opt_votes(x$sin_rcv, x$opt_plur, type = "plur", weights = weights))
	prop_df_plur <- do.call(rbind, opt_vote_prop_plur)

	prop_df <- cbind(prop_df_rcv, prop_df_plur)
	return(prop_df)
}

# Function: take ballot profile, DF of utilities, and returns vector of tactical incentives at different
return_sv_tau <- function(v_vec, util_df, s_breaks){
	v_vec <- as.numeric(v_vec) / sum(v_vec)
	sin_vec <- apply(util_df, 1, sin_vote_scalar)
	
	# RCV part
	p_list <- lapply(s_breaks, function(x) av.pivotal.event.probs.general(v_vec, rep(x, 4)))
	eu_list <- lapply(p_list, function(x) opt_vote(util_df, x, type = "rcv"))
	tau_list <- lapply(eu_list, function(x) calculate_tau(x, sin_vec))

	# Plurality part
	# Do the same for plurality
	sin_vec_plur <- sin_vote_plur_transform(sin_vec)
	v_vec_three <- c(v_vec[1] + v_vec[2] + v_vec[7], v_vec[3] + v_vec[4] + v_vec[8], v_vec[5] + v_vec[6] + v_vec[9])
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
	scalars <- apply(eu_df, 1, function(x){ 
	  out <- which(x == max(x))
	  if(length(out) > 1){out <- NA}
	  return(out)
	    }
	  )
	if(class(scalars) == "list"){
	  
	}
	return(scalars)
}

sin_vote_scalar <- function(util){
	# Input: DF of utilities
	# Output: Vector of sincere preferences (as scalars) -- under AV
	stopifnot(length(util) == 3)
	max <- which(util == max(util))
	min <- which(util == min(util))
	if(length(max) > 1){
	  print(util)
	}
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

sum_opt_votes <- function(sin_vec, opt_vec, type = "rcv", weights = rep(1, length(opt_vec))){
	# Input: vector of sincere preferences, vector of opt. votes
  # todo: add weights through wtd.table()
	# Output: Table
  
	if(type == "rcv"){
		sin_vec <- factor(sin_vec, levels = 1:6)
		opt_vec <- factor(opt_vec, levels = 1:6)
		tab <- wtd.table(sin_vec, opt_vec, weights = weights)
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
		tab <- wtd.table(sin_vec, opt_vec, weights = weights)
		out <- rbind(c(tab[1, 1], tab[1, 2]),
					 c(tab[2, 1], tab[2, 3]),
					 c(tab[3, 2], tab[3, 1]),
					 c(tab[4, 2], tab[4, 3]),
					 c(tab[5, 3], tab[5, 1]),
					 c(tab[6, 3], tab[6, 2]))
		return(colSums(out))
	}

}

vote_matrix <- function(df, type = "rcv"){
	# Input: dataframe output from return_sv_tau
	# Output: 6x6 (RCV) or 6x3 (Plurality) matrix of votes
	df$sin_rcv <- factor(df$sin_rcv, levels = 1:6)
	if(type == "rcv"){
		df$opt_rcv <- factor(df$opt_rcv, levels = 1:6)
		tab <- tapply(df$opt_rcv, df$sin_rcv, table)
		tab <- do.call(rbind, tab)
		return(tab)
	}
	if(type == "plur"){
		df$opt_plur <- factor(df$opt_plur, levels = 1:3)
		tab <- tapply(df$opt_plur, df$sin_rcv, table)
		tab <- do.call(rbind, tab)
		return(tab)
	}
}

vote_matrix_weighted <- function(df, type = "rcv", 
                                 weights = rep(1, nrow(df))){
  # Input: dataframe output from return_sv_tau
  # Output: 6x6 (RCV) or 6x3 (Plurality) matrix of votes
  df$sin_rcv <- factor(df$sin_rcv, levels = 1:6)
  if(type == "rcv"){
    df$opt_rcv <- factor(df$opt_rcv, levels = 1:6)
    tab <- matrix(wtd.table(df$sin_rcv, df$opt_rcv, weights), ncol = 6)
    return(tab)
  }
  if(type == "plur"){
    df$opt_plur <- factor(df$opt_plur, levels = 1:3)
    tab <- matrix(wtd.table(df$sin_rcv, df$opt_plur, weights), ncol = 3)
    return(tab)
  }
}

# Functions for qq-plot

qq_function_two <- function(tau_obj, utils){
  tau_list <- split(tau_obj, tau_obj$s)
  qq_list <- lapply(tau_list, function(x) as.data.frame(qqplot(x = unlist(x$tau_plur), y = unlist(x$tau_rcv), plot.it = FALSE)))
  qq_df <- as.data.frame(do.call(rbind, qq_list))
  qq_df$s <- rep(unlist(s_list), each = nrow(utils))
  return(qq_df)
}

# Functions for interdependence

new_v_vec <- function(vote_mat, v_vec_init_weighted, lambda_list, type = "rcv"){
	v_vec_init_weighted_plur <- c(v_vec_init_weighted[1] + v_vec_init_weighted[2] + v_vec_init_weighted[7],
		v_vec_init_weighted[3] + v_vec_init_weighted[4] + v_vec_init_weighted[8], 
		v_vec_init_weighted[5] + v_vec_init_weighted[6] + v_vec_init_weighted[9])

	new_vec <- v_vec_init_weighted[1:6] %*% vote_mat
	new_vec <- new_vec / sum(new_vec)
	
	if(type == "rcv"){
	  new_vec_rcv <- c(new_vec, 0, 0, 0)
		new_vec <- lapply(lambda_list, function(lambda) c(lambda * new_vec_rcv + (1 - lambda) * v_vec_init_weighted))
		return(new_vec)
	}
	if(type == "plur"){
	  new_vec_plur <- c(rep(new_vec, each = 2) / 2, 0, 0, 0)
		new_vec <- lapply(lambda_list, function(lambda) {
			x <- lambda * new_vec_plur + (1 - lambda) * v_vec_init_weighted
			return(x)
			})
		return(new_vec)
	}
}

level_two_props <- function(v_vec, lambda_list, util, sv_df, s){
	# Create 6x6 and 6x3 matrices
	vote_mat_rcv <- vote_matrix(sv_df, type = "rcv")
	vote_mat_plur <- vote_matrix(sv_df, type = "plur")
	v_vec_init_weighted <- as.numeric(v_vec[1:9] / sum(v_vec[1:9]))

	v_vec_rcv <- new_v_vec(vote_mat_rcv, v_vec_init_weighted, lambda_list, type = "rcv")
	v_vec_plur <- new_v_vec(vote_mat_plur, v_vec_init_weighted, lambda_list, type = "plur")
	
	sv_lvl_two_rcv <- lapply(v_vec_rcv, function(x) return_sv_tau(x, util, s))
	lvl_two_summary_rcv <- lapply(sv_lvl_two_rcv, function(x) return_lvl_two_prop(sv_df, x, type = "rcv"))
	lvl_two_summary_rcv <- do.call(rbind, lvl_two_summary_rcv)

	sv_lvl_two_plur <- lapply(v_vec_plur, function(x) return_sv_tau(x, util, s))
	lvl_two_summary_plur <- lapply(sv_lvl_two_plur, function(x) return_lvl_two_prop(sv_df, x, type = "plur"))
	lvl_two_summary_plur <- do.call(rbind, lvl_two_summary_plur)

	df_out <- as.data.frame(cbind(lvl_two_summary_rcv, lvl_two_summary_plur))
	names(df_out) <- c("L1RCV", "L0RCV", "L1PLUR", "L0PLUR")
	df_out$s <- s
	df_out$lambda <- unlist(lambda_list)

	return(df_out)
}

return_lvl_two_prop <- function(sv_df, lvl_2, type = "rcv"){
	n <- nrow(sv_df)
	if(type == "rcv"){
		vs_lvl_one_rcv <- 1 - sum(lvl_2$opt_rcv == sv_df$opt_rcv) / n
		vs_sin_rcv <- 1 - sum(lvl_2$opt_rcv == lvl_2$sin_rcv) / n
		return(c(vs_lvl_one_rcv, vs_sin_rcv))
	}
	if(type == "plur"){
		vs_lvl_one_plur <- 1 - sum(lvl_2$opt_plur == sv_df$opt_plur) / n
		vs_sin_plur <- 1 - sum(lvl_2$opt_plur == lvl_2$sin_plur) / n
		return(c(vs_lvl_one_plur, vs_sin_plur))
	}
}

convert_andy_to_sv_item <- function(list_item, s, v_vec){
  # Generates sv object from Andy's function and converts it into my data structure -- much faster!
  out_rcv <- sv(U = list_item$U, weights = list_item$weights, s = s, rule = "AV", v.vec = v_vec)
  out_plur <- sv(U = list_item$U, weights = list_item$weights, s = s, v.vec = v_vec)
  sin_rcv <- apply(out_rcv$V0, 1, function(x) which(x == 1))
  sin_plur <- apply(out_plur$V0, 1, function(x) which(x == 1))
  tau_rcv <- as.numeric(out_rcv$tau)
  tau_plur <- as.numeric(out_plur$tau)
  opt_rcv <- apply(out_rcv$V.mat, 1, function(x) which(x == 1))
  opt_plur <- apply(out_plur$V.mat, 1, function(x) which(x == 1))
  s <- rep(s, nrow(list_item$U))
  df <- as.data.frame(cbind(sin_rcv, sin_plur, tau_rcv, tau_plur, opt_rcv, opt_plur, s))
  return(df)
}

# write new function in terms of level_two_props_cases
# need to re-write this
convert_andy_to_sv_item_two <- function(U, w, s, v_vec){ 
  # Generates sv object from Andy's function and converts it into my data structure -- much faster!
	n_obs <- nrow(U)
  out_rcv <- sv(U = U, weights = w, s = s, rule = "AV", v.vec = v_vec)
  v_vec_plur <- c(v_vec[1] + v_vec[2], v_vec[3] + v_vec[4], v_vec[5] + v_vec[6])
  out_plur <- sv(U = U, weights = w, s = s, v.vec = v_vec_plur)
  psum_rcv <- sum(unlist(out_rcv$piv.probs))
	rcvpp_col <- matrix(unlist(out_rcv$piv.probs), nrow = n_obs, ncol = 12, byrow = T, dimnames = list(NULL, names(unlist(out_rcv$piv.probs))))
  psum_plur <- sum(unlist(out_plur$piv.probs))
	plurpp_col <- matrix(unlist(out_plur$piv.probs), nrow = n_obs, ncol = 3, byrow = T, dimnames = list(NULL, c("ABp", "ACp", "BCp")))
  sin_rcv <- apply(out_rcv$V0, 1, function(x) which(x == 1))
  sin_plur <- apply(out_plur$V0, 1, function(x) which(x == 1))
  tau_rcv <- as.numeric(out_rcv$tau)
  tau_plur <- as.numeric(out_plur$tau)
  tau_tilde_rcv <- as.numeric(tau_rcv/ psum_rcv)
  tau_tilde_plur <- as.numeric(tau_plur/ psum_plur)
  eu_mat_rcv <- out_rcv$eu.mat
  eu_mat_plur <- out_plur$eu.mat
  colnames(eu_mat_plur) <- c("EU_A", "EU_B", "EU_C")
  sin_mat <- out_rcv$V0
  rcv_mat <- out_rcv$V.mat
  opt_rcv <- apply(rcv_mat, 1, function(x) which(x == 1)[1])
  # Throw error message
  opt_plur <- apply(out_plur$V.mat, 1, function(x) which(x == 1)[1])
  # opt_rcv <- out_rcv$opt.votes.strategic
  # opt_plur <- out_plur$opt.votes.strategic
  s <- rep(s, nrow(U))
  df <- as.data.frame(cbind(sin_rcv, sin_plur, tau_rcv, tau_plur, tau_tilde_rcv, tau_tilde_plur, opt_rcv, opt_plur, s, rcvpp_col, plurpp_col, w, eu_mat_rcv, eu_mat_plur, U))
  return(df)
}

	

retrieve_opt <- function(mat_row){
	if(sum(mat_row) == 1){
		out <- which(mat_row == 1)[1]
	}
	if(sum(mat_row > 1)){

	}
}

level_two_props_cses <- function(v_vec, lambda_list, util, sv_df, s, w = c(rep(1, nrow(util)))){
  # Create 6x6 and 6x3 matrices
  vote_mat_rcv <- vote_matrix_weighted(sv_df, type = "rcv", w)
  vote_mat_plur <- vote_matrix_weighted(sv_df, type = "plur", w)
  v_vec_init_weighted <- as.numeric(v_vec[1:9] / sum(v_vec[1:9]))

  v_vec_rcv <- new_v_vec(vote_mat_rcv, v_vec_init_weighted, lambda_list, type = "rcv")
  v_vec_plur <- new_v_vec(vote_mat_plur, v_vec_init_weighted, lambda_list, type = "plur")

  
  # sv_lvl_two_rcv <- lapply(v_vec_rcv, function(x) return_sv_tau(x, util, s))
  sv_lvl_two_rcv <- lapply(v_vec_rcv, function(x) convert_andy_to_sv_item_two(util, w, s, x[1:6]))
  lvl_two_summary_rcv <- lapply(sv_lvl_two_rcv, function(x) return_lvl_two_prop(sv_df, x, type = "rcv"))
  lvl_two_summary_rcv <- do.call(rbind, lvl_two_summary_rcv)
  
  sv_lvl_two_plur <- lapply(v_vec_plur, function(x) convert_andy_to_sv_item_two(util, w, s, x[1:6]))
  lvl_two_summary_plur <- lapply(sv_lvl_two_plur, function(x) return_lvl_two_prop(sv_df, x, type = "plur"))
  lvl_two_summary_plur <- do.call(rbind, lvl_two_summary_plur)
  
  df_out <- as.data.frame(cbind(lvl_two_summary_rcv, lvl_two_summary_plur))
  names(df_out) <- c("L1RCV", "L0RCV", "L1PLUR", "L0PLUR")
  df_out$s <- s
  df_out$lambda <- unlist(lambda_list)
  
  return(df_out)
}

non_monoton <- function(sv_obj, piv_probs, weights = rep(1, nrow(sv_obj))){
  V0_tab <- matrix(0, ncol = 6, nrow = nrow(sv_obj))
  for(i in 1:nrow(V0_tab)){
    x <- sv_obj$sin_rcv[[i]]
    V0_tab[[i, x]] <- 1
  }
  no_show <- cbind(rep(piv_probs$AB.CB, nrow(sv_obj)), piv_probs$AC.BC, piv_probs$AB.AC, piv_probs$BC.AC, piv_probs$AC.AB, piv_probs$BC.BA)
  no_show_prob <- apply(no_show * V0_tab, 1, sum, na.rm = T)
  no_show_sum <- sum(no_show_prob * weights, na.rm = T) / sum(weights[!is.na(no_show_prob)], na.rm = T)


  nonmon1 <- cbind(rep(piv_probs$BC.AC, nrow(sv_obj)), piv_probs$BC.BA, piv_probs$AC.BC, piv_probs$AC.AB, piv_probs$AB.CB, piv_probs$AB.AC)  
  nonmon1_prob <- apply(nonmon1 * V0_tab, 1, sum, na.rm = T)
  nonmon1_sum <- sum(nonmon1_prob * weights, na.rm = T) / sum(weights[!is.na(nonmon1_prob)], na.rm = T)
  
  nonmon2 <- nonmon1[, c(2, 1, 4, 3, 6, 5)]
  nonmon2_prob <- apply(nonmon2 * V0_tab, 1, sum, na.rm = T)
  nonmon2_sum <- sum(nonmon2_prob * weights, na.rm = T) / sum(weights[!is.na(nonmon2_prob)], na.rm = T)
  
  total <- sum(unlist(piv_probs))
  
  return(list(no_show = no_show_sum, nonmon1 = nonmon1_sum, nonmon2 = nonmon2_sum, total = total))
  }

wasted_vote <- function(sv_obj, piv_probs, weights = rep(1, nrow(sv_obj))){
  V0_tab <- matrix(0, ncol = 3, nrow = nrow(sv_obj))
  for(i in 1:nrow(V0_tab)){
    x <- sv_obj$sin_plur[[i]]
    V0_tab[[i, x]] <- 1
  }
  wasted <- cbind(rep(piv_probs$BC, nrow(sv_obj)), piv_probs$AC, piv_probs$AB)
  #return(list(wasted, weights))
  wasted_prob <- apply(wasted * V0_tab, 1, sum, na.rm = T)
  #return(list(wasted_prob, weights))
  wasted_sum <- sum(wasted_prob * weights, na.rm = T) / sum(weights[!is.na(wasted_prob)], na.rm = T)
  
  total <- sum(unlist(piv_probs))
  return(list(wasted = wasted_sum, total = total))
}

## CSES-specific functions

remove_nas <- function(x){
  mat <- cbind(x$U, x$weights)
  mat <- na.omit(mat)
  return(list(U = mat[, 1:3], weights = as.numeric(mat[, 4])))
}

create_v_vec <- function(x){
  x$U <- x$U + runif(nrow(x$U) * ncol(x$U), min = 0, max = 0.001)
  sin_vote <- as.numeric(apply(x$U, 1, function(x) sin_vote_scalar(x)))
  num_list <- c(1:6)
  sin_df <- (sapply(num_list, function(x) as.numeric(sin_vote == x)))
  sin_df <- sin_df * x$weights
  sin_vec <- colSums(sin_df)
  return(sin_vec / sum(sin_vec))
}

# Classification function (from Andy)

classify.vec = function(v.vec, the.floor = .6, neutral.max = .1){
  v.vec = v.vec/sum(v.vec)
  m.vec = c(v.vec[1]/sum(v.vec[1:2]), v.vec[3]/sum(v.vec[3:4]), v.vec[5]/sum(v.vec[5:6]))
  most.neutral = which(abs(m.vec - .5) == min(abs(m.vec - .5)))
  if(most.neutral == 1){
    if(m.vec[2] > the.floor & m.vec[3] > the.floor){return("SP(A)")}
    if(m.vec[2] < 1 - the.floor & m.vec[3] < 1 - the.floor){return("DM(A)")}
    if(abs(m.vec[2] - .5) < neutral.max & abs(m.vec[3] - .5) < neutral.max){return("N(A)")}
  }else if(most.neutral == 2){
    if(m.vec[1] > the.floor & m.vec[3] < 1 - the.floor){return("SP(B)")}
    if(m.vec[1] < 1 - the.floor & m.vec[3] > the.floor){return("DM(B)")}
    if(abs(m.vec[1] - .5) < neutral.max & abs(m.vec[3] - .5) < neutral.max){return("N(B)")}
  }else if(most.neutral == 3){
    if(m.vec[1] < 1 - the.floor & m.vec[2] < 1 - the.floor){return("SP(C)")}
    if(m.vec[1] > the.floor & m.vec[2] > the.floor){return("DM(C)")}
    if(abs(m.vec[1] - .5) < neutral.max & abs(m.vec[2] - .5) < neutral.max){return("N(C)")}
  }
  "O"
}

# Functions for Condorcet winner
gen_v_zero <- function(sin_rcv){
# takes RCV (6 factors) vector as input
  v_zero_mat_rcv <- matrix(0, nrow = length(sin_rcv), ncol = 6)
  v_zero_mat_plur <- matrix(0, nrow = length(sin_rcv), ncol = 3)
  sin_plur <- sin_vote_plur_transform(sin_rcv)
  for(i in 1:length(sin_rcv)){
    v_zero_mat_rcv[i, sin_rcv[i]] <- 1
    v_zero_mat_plur[i, sin_plur[i]] <- 1
  }
  return(list(v_zero_mat_rcv, v_zero_mat_plur))
}

gen_v_zero_rcv <- function(sin_rcv){
# takes RCV (6 factors) vector as input
  v_zero_mat_rcv <- matrix(0, nrow = length(sin_rcv), ncol = 6)
  sin_plur <- sin_vote_plur_transform(sin_rcv)
  for(i in 1:length(sin_rcv)){
    v_zero_mat_rcv[i, sin_rcv[i]] <- 1
    }
  return(v_zero_mat_rcv)
}

gen_v_zero_plur <- function(sin_plur){
# Takes plurality vector (3 factors) as input
	v_zero_mat_plur <- matrix(0, nrow = length(sin_plur), ncol = 3)
	for(i in 1:length(sin_plur)){
		v_zero_mat_plur[i, sin_plur[i]] <- 1
	}
	return(v_zero_mat_plur)
}

cw_prop <- function(U, v_mat, v_zero, lambda, weights, size = 300, iterations = 1000, rule = "AV"){
	main_df <- t(replicate(iterations, {cw_prop_one(U, v_mat, v_zero, lambda, weights, size)}))
	cw_mat <- main_df[, 1:3]
	v_vec_mat <- main_df[, -c(1:3)]
	winners_mat <- victory.probs.from.sims(v_vec_mat, rule = rule, return.matrix = T)
	return(sum(winners_mat * cw_mat) / iterations)
}

cw_prop_one <- function(U, v_mat, v_zero, lambda, weights, size = 300){
	samp_ind <- sample(1:nrow(U), size = size, replace = T, prob = weights)
	samp_u <- U[samp_ind, ]
	cw_one <- condorcet.winner(samp_u, weights = weights[samp_ind])
	strat_ind <- c(rep(FALSE, floor((1 - lambda) * size)), rep(TRUE, ceiling(lambda * size)))
	strat_ind <- sample(strat_ind)
	new_vote_mat <- v_zero[samp_ind,]
	new_vote_mat[strat_ind, ] <- v_mat[samp_ind, ][strat_ind, ]
	new_v_vec <- ballot.props.from.vote.mat.and.weights(new_vote_mat, weights[samp_ind])
	return(c(cw_one, new_v_vec))
}

remove_nas <- function(x){
  mat <- cbind(x$U, x$weights)
  mat <- na.omit(mat)
  return(list(U = mat[, 1:3], weights = as.numeric(mat[, 4])))
}

create_v_vec <- function(x){
  x$U <- x$U + runif(nrow(x$U) * ncol(x$U), min = 0, max = 0.001)
  sin_vote <- as.numeric(apply(x$U, 1, function(x) sin_vote_scalar(x)))
  num_list <- c(1:6)
  sin_df <- (sapply(num_list, function(x) as.numeric(sin_vote == x)))
  sin_df <- sin_df * x$weights
  sin_vec <- colSums(sin_df)
  return(sin_vec / sum(sin_vec))
}

return_cw_df <- function(cw_win_rcv, cw_win_plur, lambdas, country_weight){
cw_rcv_df <- as.data.frame(do.call(rbind, cw_win_rcv))
cw_rcv_df[, 1:2] <- apply(cw_rcv_df[, 1:2], 2, as.numeric)
cw_rcv_df$cweight <- rep(country_weight, each = 5)

cw_rcv <- as.data.frame(t(sapply(lambdas, function(x) {
  z <- (boot(cw_rcv_df[cw_rcv_df$lambdas == x, c(1, 4)], boot_wmean, 1000) %>% boot.ci(type = "perc"))
  ci <- z[[4]][4:5]
  return(c(z[2], ci))
})))
cw_rcv$type <- "IRV"
cw_rcv$lambda <- lambdas

cw_plur_df <- as.data.frame(do.call(rbind, cw_win_plur))
cw_plur_df[, 1:2] <- apply(cw_plur_df[, 1:2], 2, as.numeric)
cw_plur_df$cweight <- rep(country_weight, each = 5)

cw_plur <- as.data.frame(t(sapply(lambdas, function(x) {
  z <- (boot(cw_plur_df[cw_plur_df$lambdas == x, c(1, 4)], boot_wmean, 1000) %>% boot.ci(type = "perc"))
  ci <- z[[4]][4:5]
  return(c(z[2], ci))
})))
cw_plur$type <- "Plurality"
cw_plur$lambda <- lambdas

cw_df <- rbind(cw_rcv, cw_plur)
cw_df[, 1:3] <- apply(cw_df[, 1:3], 2, as.numeric)
names(cw_df)[1:3] <- c("mu", "lower", "upper")
return(cw_df)
}

boot_wmean <- function(x, d){
  weighted.mean(x[d, 1], x[d, 2])
}

## STUFF FOR INTERDEPENDENCE

iteration_return_summary <- function(object, v_vec, lambda, s){
  obj <- object
  tab <- convert_andy_to_sv_item_two(obj$U, obj$weights, s, v_vec)

  prev_rcv <- weighted.mean(tab$tau_rcv > 0, weights = obj$weights)
  mag_rcv <- weighted.mean(tab$tau_rcv[tab$tau_rcv > 0], weights = obj$weights[tab$tau_rcv > 0])
  eb_rcv <- prev_rcv * mag_rcv
  strat_vec_rcv <- as.numeric(table(factor(tab$opt_rcv, 1:6))/nrow(tab))
  new_vec_rcv <- lambda * strat_vec_rcv + (1 - lambda) * v_vec

  prev_plur <- weighted.mean(tab$tau_plur > 0, weights = obj$weights)
  weight_mag_plur <- sum(obj$weights[tab$tau_plur > 0])
  # print(weight_mag_plur)
  mag_plur <- weighted.mean(tab$tau_plur[tab$tau_plur > 0], weights = obj$weights[tab$tau_plur > 0])
  eb_plur <- prev_plur * mag_plur

  strat_vec_plur <- rep(as.numeric(table(factor(tab$opt_plur, 1:3))/nrow(tab)), each = 2) /2
  new_vec_plur <- lambda * strat_vec_plur + (1 - lambda) * v_vec

  return(list(rcv_sum = c(prev_rcv, mag_rcv, eb_rcv),
              rcv_vec = new_vec_rcv,
              plur_sum = c(prev_plur, mag_plur, eb_plur),
              plur_vec = new_vec_plur))
}

iteration_wrapper <- function(object, v_vec, lambda, s, k){
    rcv_sum_list <- list()
    rcv_v_vec_list <- list(v_vec)
    plur_sum_list <- list()
    plur_v_vec_list <- list(v_vec)

    # RCV loop
    for (j in 1:k) {
       out <- iteration_return_summary(object, rcv_v_vec_list[[j]], lambda, s)
       rcv_sum_list[[j]] <- out$rcv_sum
       rcv_v_vec_list[[j + 1]] <- out$rcv_vec
    }
    rcv_sum <- as.data.frame(do.call(rbind, rcv_sum_list))
    names(rcv_sum) <- c("Prevalence", "Magnitude", "ExpBenefit")
    rcv_sum$k <- 1:k
    rcv_v_vec <- as.data.frame(do.call(rbind, rcv_v_vec_list))

    # Plurality loop
    for (j in 1:k) {
       out <- iteration_return_summary(object, plur_v_vec_list[[j]], lambda, s)
       plur_sum_list[[j]] <- out$plur_sum
       plur_v_vec_list[[j + 1]] <- out$plur_vec
    }
    plur_sum <- as.data.frame(do.call(rbind, plur_sum_list))
    names(plur_sum) <- c("Prevalence", "Magnitude", "ExpBenefit")
    plur_sum$k <- 1:k
    plur_v_vec <- as.data.frame(do.call(rbind, plur_v_vec_list))

    return(list(rcv_sum, rcv_v_vec, plur_sum, plur_v_vec))
}

## NEW VERSION AS OF MAY 2019

one_iteration <- function(object, v_vec, lambda, s){
  obj <- object
  tab <- convert_andy_to_sv_item_two(obj$U, obj$weights, s, v_vec)

  strat_vec_rcv <- wtd.table(x = factor(tab$opt_rcv, 1:6),
                             weights = obj$weights) %>%
  								as.numeric
  strat_vec_rcv <- strat_vec_rcv / sum(strat_vec_rcv)
  new_vec_rcv <- lambda * strat_vec_rcv + (1 - lambda) * v_vec

  strat_plur_min <- wtd.table(x = factor(tab$opt_plur, 1:3),
                             weights = obj$weights) %>%
  								as.numeric
  strat_vec_plur <- rep(strat_plur_min / sum(strat_plur_min)
  						, each = 2) /2
  new_vec_plur <- lambda * strat_vec_plur + (1 - lambda) * v_vec

  return(list(df = tab,
              rcv_vec = new_vec_rcv,
              rcv_best_response = strat_vec_rcv,
              plur_vec = new_vec_plur,
              plur_best_response = strat_vec_plur))
}

many_iterations <- function(object, v_vec, lambda, s, k){
    rcv_df_list <- list()
    rcv_v_vec_list <- list(v_vec)
    plur_df_list <- list()
    plur_v_vec_list <- list(v_vec)

    # RCV loop
    for (j in 1:k) {
       out <- one_iteration(object, rcv_v_vec_list[[j]], lambda, s)
       rcv_df_list[[j]] <- out$df %>% mutate(iter = j)
       rcv_v_vec_list[[j + 1]] <- out$rcv_vec
    }
    rcv_sum <- as.data.frame(do.call(rbind, rcv_df_list))
    rcv_v_vec <- as.data.frame(do.call(rbind, rcv_v_vec_list))

    # Plurality loop
    for (j in 1:k) {
       out <- one_iteration(object, plur_v_vec_list[[j]], lambda, s)
       plur_df_list[[j]] <- out$df %>% mutate(iter = j)
       plur_v_vec_list[[j + 1]] <- out$plur_vec
    }
    plur_sum <- as.data.frame(do.call(rbind, plur_df_list))
    plur_v_vec <- as.data.frame(do.call(rbind, plur_v_vec_list))

    return(list(rcv_sum, rcv_v_vec, plur_sum, plur_v_vec))
}

many_iterations_until_convergence <- function(object, v_vec, lambda, s, thresh, max_iter){
    rcv_df_list <- list()
    rcv_v_vec_list <- list(v_vec)
    rcv_br_v_vec <- list(v_vec)
    plur_df_list <- list()
    plur_v_vec_list <- list(v_vec)
    plur_br_v_vec <- list(v_vec)

    w <- object$weights

    # RCV loop
    k_rcv <- 0
    epsilon_rcv <- 1
    thresh_ind <- 0
    while (k_rcv < max_iter) {
       k_rcv <- k_rcv + 1	
       out <- one_iteration(object, rcv_v_vec_list[[k_rcv]], lambda, s)
       rcv_df_list[[k_rcv]] <- out$df %>% mutate(iter = k_rcv, 
                                                 converged = thresh_ind)
       rcv_v_vec_list[[k_rcv + 1]] <- out$rcv_vec
       rcv_br_v_vec[[k_rcv + 1]] <- out$rcv_best_response
       epsilon_rcv <- sqrt(sum((rcv_v_vec_list[[k_rcv + 1]] - rcv_br_v_vec[[k_rcv + 1]])^2))
       if (epsilon_rcv < thresh){
       	thresh_ind <- 1
       }

    }
    rcv_sum <- as.data.frame(do.call(rbind, rcv_df_list))
    rcv_v_vec <- as.data.frame(do.call(rbind, rcv_v_vec_list))
    rcv_convg <- thresh_ind

    k_plur <- 0
    epsilon_plur <- 1
    thresh_ind <- 0

    # Plurality loop
    while (k_plur < max_iter) {
       k_plur <- k_plur + 1	
       out <- one_iteration(object, plur_v_vec_list[[k_plur]], lambda, s)
       plur_df_list[[k_plur]] <- out$df %>% mutate(iter = k_plur,
                                                   converged = thresh_ind)
       plur_v_vec_list[[k_plur + 1]] <- out$plur_vec
       plur_br_v_vec[[k_plur + 1]] <- out$plur_best_response
       epsilon_plur <- sqrt(sum((plur_v_vec_list[[k_plur + 1]] - plur_br_v_vec[[k_plur + 1]])^2))
       if(epsilon_plur < thresh){
       	thresh_ind <- 1
       }
    }

    plur_sum <- as.data.frame(do.call(rbind, plur_df_list))
    plur_v_vec <- as.data.frame(do.call(rbind, plur_v_vec_list))
    plur_convg <- thresh_ind

    return(list(rcv_df = rcv_sum, 
                rcv_v_vec = rcv_v_vec, 
                plur_df = plur_sum, 
                plur_v_vec = plur_v_vec, 
                rcv_br = rcv_br_v_vec, 
                plur_br = plur_br_v_vec,
                rcv_convg = rcv_convg,
                plur_convg = plur_convg,
                weights = w))
}

# This is a slimmed down version that only computes the v_vec path for RCV iterations.
iteration_convergence_slim <- function(object, v_vec, lambda, s, max_iter){
	rcv_v_vec_list <- list(v_vec)
	rcv_br_list <- list(v_vec)
	# iteration loop
	for(i in 1:max_iter){
		out <- one_iteration(object, rcv_v_vec_list[[i]], lambda, s)
		rcv_v_vec_list[[i + 1]] <- out$rcv_vec
		rcv_br_list[[i + 1]] <- out$rcv_br
	}

	return(rcv_v_vec_list)
}

# Slimmed down functions for sensitivity analysis

one_iteration_rcv_only <- function(object, v_vec, lambda, s){
  obj <- object
  tab <- convert_andy_to_sv_item_two(obj$U, obj$weights, s, v_vec)

  strat_vec_rcv <- as.numeric(table(factor(tab$opt_rcv, 1:6))/nrow(tab))
  new_vec_rcv <- lambda * strat_vec_rcv + (1 - lambda) * v_vec

  return(list(df = tab,
              rcv_vec = new_vec_rcv))
}

many_iterations_rcv_only <- function(object, v_vec, lambda, s, k){
    rcv_df_list <- list()
    rcv_v_vec_list <- list(v_vec)

    # RCV loop
    for (j in 1:k) {
       out <- one_iteration_rcv_only(object, 
                                     rcv_v_vec_list[[j]], 
                                     lambda, 
                                     s)
       rcv_df_list[[j]] <- out$df %>% mutate(iter = j)
       rcv_v_vec_list[[j + 1]] <- out$rcv_vec
    }
    rcv_sum <- as.data.frame(do.call(rbind, rcv_df_list))
    rcv_v_vec <- as.data.frame(do.call(rbind, rcv_v_vec_list))

    return(list(rcv_sum, rcv_v_vec))
}

many_iterations_rcv_only_light <- function(object, v_vec, lambda, s, k){
    rcv_df_list <- list()
    rcv_v_vec_list <- list(v_vec)

    # RCV loop
    for (j in 1:k) {
       out <- one_iteration_rcv_only(object, 
                                     rcv_v_vec_list[[j]], 
                                     lambda, 
                                     s)
       rcv_v_vec_list[[j + 1]] <- out$rcv_vec
    }
    rcv_v_vec <- as.data.frame(do.call(rbind, rcv_v_vec_list))

    return(rcv_v_vec)
}

piv_ratio <- function(x){
	pprobs <- x[1, ] %>% select("AB":"BCp")
	prat <- sum(pprobs[13:15]) / sum(pprobs[1:12])
	return(prat)
}


# Write a function to draw 100 v.vecs with distance smaller than original

draw_restricted_vvecs <- function(orig_vvec, eqm_vvec, n){
	orig_dist <- sqrt(sum((eqm_vvec - orig_vvec)^2))
	first_draw <- rdirichlet(10 * n, as.numeric(orig_vvec) * 85)
	smaller_dist <- apply(first_draw, 1, function(x) sqrt(sum((x - orig_vvec)^2)) < orig_dist)
	new_df <- first_draw[smaller_dist, ]
	i <- 1
	while(nrow(new_df) < n){
		print(nrow(new_df))
		print(i)
		i <- i + 1
		another_draw <- rdirichlet(10 * n,  as.numeric(orig_vvec) * 85)
		another_dist <- apply(another_draw, 1, function(x) sqrt(sum((x - orig_vvec)^2)) < orig_dist)
		new_df <- rbind(new_df, another_draw[another_dist, ])
	}
	if(nrow(new_df) > n){
		new_df <- new_df[1:n, ]
	}
	return(new_df)
}

# Euclidean distances

euclid <- function(x) {
  # Euclidean distance between one iteration and next
  if (isdf <- is.data.frame(x)) {
    x <- data.matrix(x)
  }
  dij <- c(NA, sqrt(rowSums((tail(x, -1) - head(x, -1))^2)))
  x <- cbind(x, diff = dij)
  if (isdf) {
    x <- as.data.frame(x)
  }
  x
}

euclid_first <- function(df){
  # Euclidean distance between one iteration and first
  df_follow <- df[2:nrow(df), ]
  dist <- apply(df_follow, 1, function(x) sqrt(sum((x - df[1, ])^2)))
  return(dist)
}

euclid_together <- function(df){
  a <- euclid(df)
  b <- c(NA, euclid_first(df))
  return(cbind(a, b))
}


### FUNCTIONS FOR BIG LOOPS

get_vote_mat_df <- function(case_iter){
  opt <- optimal.vote.from.V.mat(case_iter$V.mat) %>% as.vector
  sin <- optimal.vote.from.V.mat(case_iter$sincere.vote.mat) %>% as.vector

  tab <- table(sin, factor(opt, c("abc", "acb", "bac", "bca", "cab", "cba"))) %>% as.data.frame 
  return(tab)
}

iter_wrapper <- function(case){
  list_wrap <- lapply(case, function(x) get_vote_mat_df(x))
  for(i in 1:length(list_wrap)){
    list_wrap[[i]]$iter <- i
    names(list_wrap[[i]])[2] <- "opt"
  }
  df_wrap <- do.call(rbind, list_wrap) %>% mutate(vote = interaction(sin, opt))
  return(df_wrap)
}

calc_lag_dist <- function(df, n_lag, avg_span){
	# calculates
	sqsums <- rowSums(df^2)
	max_dist <- nrow(df) - n_lag
	d <- rep(NA, n_lag)
	for(i in 1:max_dist){
		avg_d <- sum(colMeans(df[i:(i + avg_span), ])^2) 
		avg_d <- sqsums[i + n_lag - 1]
		d[i + n_lag] <- sqrt(abs(avg_d - sqsums[i + n_lag]))
	}
	return(d)
}

plot_v_vec_distance <- function(obj, filepath, 
                                custom_red = 0.0001, 
                                n_lag = NULL,
                                avg_span){

################################################
# Takes the full object and creates the following
# six distance plots:
# 	"/dist_irv" -- v-vec(t - 1) to v-vec(t) under IRV
# 	"/dist_plur" -- v-vec(t - 1) to v-vec(t) under Plur
# 	"/dist_irv_lagged" -- v-vec(lag) to v-vec(t) under IRV
#   "/dist_plur_lagged" -- v-vec(lag) to v-vec(t) under Plur
# 	"/dist_br" -- v-vec(t) to BR(t) under both
#	"/dist_br_lagged" -- v-vec(t) to BR(t) under both
################################################

	cat("V.VEC TO V.VEC DISTANCE \n")
	cat("Calculating distances. \n")

	conv_dist <- list()

	for (j in 1:length(obj)){
		# print(j)
		df <- obj[[j]][[2]] 
		df_dist <- euclid_together(df) %>%
		mutate(system = 'IRV',
		       case = names(obj)[j],
		       iter = 1:nrow(.),
		       converged = obj[[j]][[7]] %>% as.factor)
		if(is.numeric(n_lag) == TRUE){
			df_dist$avg_dist <- calc_lag_dist(df, n_lag, avg_span)
		}
		conv_dist[[(2 * j) - 1]] <- df_dist

		df <- obj[[j]][[4]] 
		df_dist <- euclid_together(df) %>%
		mutate(system = 'Plurality',
		       case = names(obj)[j],
		       iter = 1:nrow(.),
		       converged = obj[[j]][[8]] %>% as.factor)
		if(is.numeric(n_lag) == TRUE){
			df_dist$avg_dist <- calc_lag_dist(df, n_lag, avg_span)
		}
		conv_dist[[2 * j]] <- df_dist
	}

	conv_dist_df <- do.call(rbind, conv_dist)

	cat("Distances done. Plotting now. \n")

	# Plot IRV convergence
	p1 <- ggplot(conv_dist_df %>% filter(system == "IRV"), aes(iter, log(diff))) +
		geom_line(aes(group = interaction(case, system)), 
		          alpha = 0.1,
		          color = "blue") +
		theme_sv() +
		labs(x = "Iteration (Strategic Responsiveness)",
		     y = "ln(Distance)")
	ggsave(here(paste0(filepath, "/dist_irv.pdf")), 
	       p1, 
	       device = cairo_pdf,
	       width = 4, 
	       height = 4)

	# Plot distance to lagged average
	p1_lag <- ggplot(conv_dist_df %>% filter(system == "IRV"), 
	                 aes(iter, log(avg_dist))) +
		geom_line(aes(group = interaction(case, system)), 
		          alpha = 0.1,
		          color = "blue") +
		theme_sv() +
		labs(x = "Iteration (Strategic Responsiveness)",
		     y = "ln(Distance to lagged average)")
	ggsave(here(paste0(filepath, "/dist_irv_lagged.pdf")), 
	       p1_lag, 
	       device = cairo_pdf,
	       width = 4,
	       height = 4)

	# Plot plurality convergence
	p2 <- ggplot(conv_dist_df %>% filter(system == "Plurality"), aes(iter, log(diff))) +
		geom_line(aes(group = interaction(case, system)), 
		          alpha = 0.1, 
		          color = "blue") +
		theme_sv() + 
		labs(x = "Iteration (Strategic Responsiveness)",
		     y = "ln(Distance)")
	ggsave(here(paste0(filepath, "/dist_plur.pdf")), 
	       p2, 
	       device = cairo_pdf,
	       width = 4,
	       height = 4)

	################################################

	cat("V.VEC TO BEST RESPONSE DISTANCE \n")
	cat("Distances done. Plotting now. \n")

	# compare best response vector to polling v_vec
	br_dist_list <- list()
	for(c in 1:length(cases_converge)){
		case_name <- names(big_list_na_omit)[c]
		br_dist_list[[c]] <- br_distance(cases_converge[[c]],
		                                 n_lag = 20,
		                                 avg_span = 10,
		                                 case = case_name)
	}
	br_dist <- do.call(rbind, br_dist_list)
	names(br_dist)[1] <- 'd'
	
	# plot results

	# t to t-1
	plot_br_1 <- ggplot(br_dist %>% filter(iter > 0), aes(iter, log(d))) +
		geom_line(aes(group = case), alpha = 0.1) +
		facet_wrap(. ~ system) +
		labs(x = "Iteration", y = "ln(Distance)") +
		theme_sv()
	ggsave(here(paste0(filepath, "/dist_br.pdf")),
	            plot_br_1,
	            height = 4,
	            width = 8,
	            device = cairo_pdf)

	# t to avg(t-10 : t-20)
	plot_br_2 <- ggplot(br_dist %>% filter(iter > 0), aes(iter, log(l_avg))) +
		geom_line(aes(group = case), alpha = 0.1) +
		facet_wrap(. ~ system) +
		labs(x = "Iteration", y = "ln(Lagged Distance)") +
		theme_sv()
	ggsave(here(paste0(filepath, "/dist_br_lag.pdf")),
	            plot_br_2,
	            height = 4,
	            width = 8,
	            device = cairo_pdf)

	# linear, rather than logged
	plot_br_3 <- ggplot(br_dist %>% filter(iter > 0), aes(iter, d)) +
		geom_line(aes(group = case), alpha = 0.1) +
		facet_wrap(. ~ system) +
		labs(x = "Iteration", y = "Distance") +
		theme_sv()
	ggsave(here(paste0(filepath, "/dist_br_lin.pdf")),
	            plot_br_3,
	            height = 4,
	            width = 8,
	            device = cairo_pdf)

	plot_br_4 <- ggplot(br_dist %>% filter(iter > 0), aes(iter, l_avg)) +
		geom_line(aes(group = case), alpha = 0.1) +
		facet_wrap(. ~ system) +
		labs(x = "Iteration", y = "Lagged Distance") +
		theme_sv()
	ggsave(here(paste0(filepath, "/dist_br_lag_lin.pdf")),
	            plot_br_4,
	            height = 4,
	            width = 8,
	            device = cairo_pdf)

	save(br_dist, file = here(paste0(filepath, "v_vecs_br.Rdata")))

	# if(is.numeric(n_lag) == TRUE ){
	# 	return(conv_dist_df)
	# }

}

# Combined v_vecs (IRV and plurality)
joint_v_vec_plot <- function(obj, filepath){
	rcv_df <- list()
	plur_df <- list()
	# Pull cases together into one DF
	for (i in 1:length(obj)){
		rcv_df[[i]] <- obj[[i]][[2]] %>% mutate(iter = 1:nrow(.),
		                                   case = names(obj)[i],
		                                   system = "IRV",
		                                   state = c("first", rep("middle", nrow(.) - 2), "last"))
		plur_df[[i]] <- obj[[i]][[4]] %>% mutate(iter = 1:nrow(.),
		                                   case = names(obj)[i],
		                                   system = "Plurality",
		                                   state = c("first", rep("middle", nrow(.) - 2), "last"))
	}
	rcv_df <- do.call(rbind, rcv_df)
	plur_df <- do.call(rbind, plur_df)
	v_vec_df <- rbind(rcv_df, plur_df)

	# Create plot
	p1 <- ggtern(v_vec_df, aes(V1 + V2, V3 + V4, V5 + V6)) +
	  geom_line(aes(group = interaction(system, case)),
	            alpha = 0.2) +
	  geom_point(data = v_vec_df %>% filter(state %in% c("first", "last")),
	             aes(colour = state),
	             size = 0.5,
	             alpha = 0.6) +
	  facet_wrap(~ system) +
	  theme_sv() +
	  theme(legend.position = "bottom") +
	  labs(x = "A", y = "B", z = "C",
	       colour = "Iteration")
	 ggsave(here(paste0(filepath, "/v_vec_path.pdf")), 
	       p1, 
	       device = cairo_pdf,
	       height = 4,
	       width = 6)
}

non_conv_v_vec_plot <- function(obj, filepath, max_iter){
	which_cases <- which(sapply(obj, function(x) x[[7]] == 0))
	rcv_df <- list()

	if(length(which_cases) < 1){
		return(cat("No non-convergent cases"))
	}

	for(i in which_cases){
		rcv_df[[i]] <- obj[[i]][[2]] %>% mutate(iter = 1:nrow(.),
		                                   case = names(obj)[i],
		                                   system = "IRV",
		                                   state = c("first", rep("middle", nrow(.) - 2), "last"))
	}
	rcv_df <- do.call(rbind, rcv_df)
	ggtern(rcv_df, aes(V1 + V2, V3 + V4, V5 + V6)) +
		geom_line(alpha = .2) +
		geom_point(data = rcv_df %>% filter(state %in% c("first", "last")), aes(colour = state)) +
		facet_wrap(~ case) +
		theme_sv() +
		labs(x = "A", y = "B", z = "C",
		     phase = "Iteration")

}

non_conv_strat_votes <- function(obj, filepath, max_iter, all_df = FALSE){

	vote_tab <- lapply(obj, function(x){
		x[[1]] %>% 
			group_by(sin_rcv, opt_rcv, iter) %>% 
			summarise(freq = sum(opt_rcv)) %>% 
			filter(sin_rcv != opt_rcv)
	})
	for(i in 1:length(obj)){
		vote_tab[[i]] <- vote_tab[[i]] %>% 
			mutate(case = names(obj)[i])
	}
	vote_tab <- do.call(rbind, vote_tab)

	if(all_df == FALSE){
		which_cases <- which(sapply(obj, function(x) x[[7]] == 0))
			if(length(which_cases) < 1){
		return(cat("No non-convergent cases"))
	}
		vote_tab <- vote_tab %>% filter(case %in% names(obj)[which_cases])
	}
	return(vote_tab)
}

# TODO: 1. create list of cases and weights;
# 2. add weighted means
get_sum_stats <- function(obj){
	list <- lapply(obj, function(x){
		w <- x[[9]]
		out <- list()
		out[["irv"]] <- x[[1]] %>% 
			group_by(iter) %>% 
			summarise(Prevalence = wtd.mean(tau_rcv > 0, w), 
			          Magnitude = wtd.mean(tau_rcv[tau_rcv > 0], w[tau_rcv > 0]), 
			          ExpBenefit = Prevalence * Magnitude,
			          System = "IRV")
		out[["plur"]] <- x[[3]] %>%	
			group_by(iter) %>% 
			summarise(Prevalence = wtd.mean(tau_plur > 0, w), 
			          Magnitude = wtd.mean(tau_plur[tau_plur > 0], w[tau_plur > 0]), 
			          ExpBenefit = Prevalence * Magnitude,
			          System = "Plurality")
		out <- do.call(rbind, out)
	})
	for(i in 1:length(obj)){
		list[[i]] <- list[[i]] %>% mutate(case = names(obj)[i])
	}
	list <- do.call(rbind, list)
	return(list)
}

# can I do this somehow more efficiently?
# necessary for conjectures (adjusted to cases_converge)

# TODO: fix utilities (since now appended to main 
rearrange <- function(x){
  voter <- as.numeric(x[1])
  order_utils <- unlist(c(cmat_u[voter, ]))
  utils <- x[35:37]
  utils <- utils[order_utils]
  order_pprobs <- unlist(c(cmat_plur[voter, ]))
  order_pprobs_rcv <- unlist(c(cmat_rcv[voter, ]))
  pprobs_plur <- x[22:24]
  pprobs_plur <- unlist(c(pprobs_plur[order_pprobs]))
  pprobs_rcv <- x[10:21]
  pprobs_rcv <- unlist(c(pprobs_rcv[order_pprobs_rcv]))
  out <- unlist(c(utils, pprobs_plur, pprobs_rcv))
  names(out) <- c("uA", "uB", "uC", "ABpo", "ACpo", "BCpo", "ABo", "ACo", "BCo", "AB_ABo", "AB_ACo", "AB_CBo", "AC_ACo", "AC_ABo", "AC_BCo", "BC_BCo", "BC_BAo", "BC_ACo")
  return(out)
}

# produce a vector of distance elements for each row of a given v_vec DF
get_distance_by_row <- function(df){
	out <- df^2 %>% rowSums %>% sqrt
	return(out)
}

# produce DF with averaged values
get_lagged_avg <- function(df, avg_span, n_lag){
	new_mat <- matrix(NA, nrow = nrow(df), ncol = ncol(df))
	for(i in 1:(nrow(df) - n_lag)){
		new_mat[i + n_lag, ] <- colMeans(df[i:(i + avg_span), ])
	}
	return(new_mat)
}

# Neat function to calculate distances: vvec and best response
# could modify to accommodate vvec to vvec too...
br_distance <- function(obj, n_lag, avg_span, case = "NA"){
	br_rcv <- obj[[5]] %>% do.call(rbind, .)
	d_rcv <- (obj[[2]] - br_rcv) %>%
		get_distance_by_row(.) %>%
		as.data.frame %>%
		mutate(iter = 0:(length(.) - 1),
		       system = "IRV",
		       case = case)
	avg_rcv <- get_lagged_avg(obj[[2]], 
	                          avg_span = avg_span, 
	                          n_lag = n_lag)
	d_rcv_avg <- (avg_rcv - br_rcv) %>% 
		get_distance_by_row(.)
	d_rcv$l_avg <- d_rcv_avg

	br_plur <- obj[[6]] %>% do.call(rbind, .)
	d_plur <- (obj[[4]] - br_plur) %>%
		get_distance_by_row(.) %>%
		as.data.frame %>%
		mutate(iter = 0:(length(.) - 1),
		       system = "Plurality",
		       case = case)	  	
	avg_plur <- get_lagged_avg(obj[[4]],
	                           avg_span = avg_span,
	                           n_lag = n_lag)
	d_plur$l_avg <- (avg_plur - br_plur) %>%
		get_distance_by_row(.)
	return(rbind(d_rcv, d_plur))
}