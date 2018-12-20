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
convert_andy_to_sv_item_two <- function(U, w, s, v_vec){
  # Generates sv object from Andy's function and converts it into my data structure -- much faster!
  out_rcv <- sv(U = U, weights = w, s = s, rule = "AV", v.vec = v_vec)
  v_vec_plur <- c(v_vec[1] + v_vec[2], v_vec[3] + v_vec[4], v_vec[5] + v_vec[6])
  out_plur <- sv(U = U, weights = w, s = s, v.vec = v_vec_plur)
  sin_rcv <- apply(out_rcv$V0, 1, function(x) which(x == 1))
  sin_plur <- apply(out_plur$V0, 1, function(x) which(x == 1))
  tau_rcv <- as.numeric(out_rcv$tau)
  tau_plur <- as.numeric(out_plur$tau)
  opt_rcv <- apply(out_rcv$V.mat, 1, function(x) which(x == 1))
  opt_plur <- apply(out_plur$V.mat, 1, function(x) which(x == 1))
  s <- rep(s, nrow(U))
  df <- as.data.frame(cbind(sin_rcv, sin_plur, tau_rcv, tau_plur, opt_rcv, opt_plur, s))
  return(df)
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