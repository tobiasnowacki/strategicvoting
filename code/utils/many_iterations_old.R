convert_andy_to_sv_item_two_old <- function(U, w, s, v_vec){ 
  # Generates sv object from Andy's function and converts it into my data structure -- much faster!
  n_obs <- nrow(U)
  out_rcv <- sv_old(U = U, weights = w, s = s, rule = "AV", v.vec = v_vec)
  v_vec_plur <- c(v_vec[1] + v_vec[2], v_vec[3] + v_vec[4], v_vec[5] + v_vec[6])
  out_plur <- sv_old(U = U, weights = w, s = s, v.vec = v_vec_plur)
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
  return(list(df, out_plur$p.mat, out_plur$piv.probs))
}


one_iteration_old <- function(object, v_vec, lambda, s, ae){
  obj <- object
  tab <- convert_andy_to_sv_item_two_old(obj$U, obj$weights, s, v_vec)
  tabmat = tab[[2]]
  tabpiv = tab[[3]]
  tab = tab[[1]]
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
              plur_best_response = strat_vec_plur,
              p_mat = tabmat,
              p_pr = tabpiv))
}


many_iterations_until_convergence_old <- function(object, v_vec, lambda, s, thresh, max_iter){
    rcv_df_list <- list()
    rcv_v_vec_list <- list(v_vec)
    rcv_br_v_vec <- list(v_vec)
    plur_df_list <- list()
    plur_v_vec_list <- list(v_vec)
    plur_br_v_vec <- list(v_vec)
    plur_out_pmat_list <- list()
    plur_out_prob_list <- list()


    w <- object$weights

    # RCV loop
    k_rcv <- 0
    epsilon_rcv <- 1
    thresh_ind <- 0
    while (k_rcv < max_iter) {
       k_rcv <- k_rcv + 1   
       out <- one_iteration_old(object, rcv_v_vec_list[[k_rcv]], lambda, s, ae)
       rcv_df_list[[k_rcv]] <- out$df %>% mutate(iter = k_rcv, 
                                                 converged = thresh_ind)
       rcv_v_vec_list[[k_rcv + 1]] <- out$rcv_vec
       rcv_br_v_vec[[k_rcv + 1]] <- out$rcv_best_response
       epsilon_rcv <- sqrt(sum((rcv_v_vec_list[[k_rcv + 1]] - rcv_br_v_vec[[k_rcv + 1]])^2))
       # cat(c(k_rcv, epsilon_rcv, thresh, "\n"))
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
       out <- one_iteration_old(object, plur_v_vec_list[[k_plur]], lambda, s, ae)
       plur_df_list[[k_plur]] <- out$df %>% mutate(iter = k_plur,
                                                   converged = thresh_ind)
       plur_v_vec_list[[k_plur + 1]] <- out$plur_vec
       plur_br_v_vec[[k_plur + 1]] <- out$plur_best_response
       plur_out_pmat_list[[k_plur]] <- out$p_mat
       plur_out_prob_list[[k_plur]] <- out$p_pr
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
                weights = w,
                p_mat_list = plur_out_pmat_list,
                p_pr_list = plur_out_prob_list))
}