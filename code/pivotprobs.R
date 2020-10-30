library(tidyverse)
library(devtools)
devtools::install_github("aeggers/pivotprobs")
library(pivotprobs)

alpha3 <- c(.4, .35, .25)*85
electorate_size <- 15000

plurality_election(k = 3, n = electorate_size) %>% 
  election_event_probs(method = "ev", alpha = alpha3) %>% 
  combine_P_matrices()

U_test = big_list_na_omit[[7]]$U
v_vec_test = big_list_na_omit[[7]]$v_vec

names(big_list_na_omit[[1]])
head(U_test)

sv = function(U, weights = NULL, v.vec = NULL, s, rule = "plurality", V0 = NULL, ae_pack = TRUE){
  
  stopifnot(!is.null(colnames(U)))
  candidates = sort(colnames(U))
  U = U[, candidates]
  
  stopifnot(rule %in% c("plurality", "AV"))
  stopifnot(ncol(U) == 3) # we can't do more than three with this technique -- expanding would require generalizing the P_mat_at_pivotal_events step below 
  
  if(is.null(weights)){weights = rep(1, nrow(U))}
  
  # get sincere vote mat (V0)
  
  if(is.null(V0)){
    V0 = sincere.vote.mat.from.U(U, rule = rule, candidates = candidates)
  }
  
  if(is.null(v.vec)){
    v.vec = ballot.props.from.vote.mat.and.weights(V0, weights = weights)
  }
  
  if(rule %in% c("AV")){
    ballots = apply(permutations(n = length(candidates), r = length(candidates), v = candidates, repeats.allowed = F), 1, paste, collapse = "")
  }else{
    ballots = candidates
  }
  
  # get pivotal probabilities -- depends on the rule
  if(rule == "AV"){
    if(ae_pack == FALSE){
        pps = av.pivotal.event.probs.general(v.vec = c(v.vec, 0,0,0), s.vec = rep(s, 4))  # adding truncated ballots
        P.mat = t(P_mat_at_pivotal_events(pps, rule = rule, ballots = ballots))
        #  return(P.mat)
    }
    if(ae_pack == TRUE){
        pps = irv_election(n = 1000000) %>%
            election_event_probs(method = "en", alpha = (v.vec * 85))
        P.mat <- pps %>%
            combine_P_matrices() %>% t()
        pps = map(pps, "integral")
        pps = pps[c("a_b", "a_c", "b_c", "a_b|ab", "a_b|ac", "a_b|cb", "a_c|ac", "a_c|ab", "a_c|bc", "b_c|bc", "b_c|ba", "b_c|ac")]
        names(pps) = c("AB", "AC", "BC", "AB.AB", "AB.AC", "AB.CB", "AC.AC", "AC.AB", "AC.BC", "BC.BC", "BC.BA", "BC.AC")
        P.mat = P.mat[-7, ]
        #  return(P.mat)
    }

  }else if(rule == "plurality"){
    if(ae_pack == FALSE){
        v_vec_test_3 = c(v.vec[1] + v.vec[2], v.vec[3] + v.vec[4], v.vec[5] + v.vec[6])
        pps = plurality.pivotal.probabilities(v.vec = v_vec_test_3, s = s) # a touch of noise to prevent any alpha component from being 0 
        P.mat = t(P_mat_at_pivotal_events(pps, rule = rule, ballots = ballots))
        # return(P.mat)

    # return(pps)
    } else if(ae_pack == TRUE){
        v_vec_test_3 = c(v.vec[1] + v.vec[2], v.vec[3] + v.vec[4], v.vec[5] + v.vec[6])
        pps = plurality_election(k = 3, n = 1000000) %>%
            election_event_probs(method = "ev", alpha = (v_vec_test_3 * 85))
        P.mat = pps %>% combine_P_matrices() %>% t()
        pps = map(pps, "integral")
        pps = pps[c("a_b", "a_c", "b_c")]
        names(pps) = c("AB", "AC", "BC")
        P.mat = P.mat[-4, ]
        # return(P.mat)
    }
  }
  # v_vecplur = c(v.vec[1] + v.vec[2], v.vec[3] + v.vec[4], v.vec[5] + v.vec[6])
  # P.mat = plurality_election(k = 3, n = 1000000) %>%
  #   election_event_probs(method = "ev", alpha = v_vecplur * s) %>%
  #   combine_P_matrices()
  # return(P.mat)
  # P.mat = t(P_mat_at_pivotal_events(pps, rule = rule, ballots = ballots, normalize = FALSE))
  # return(P.mat)
    # return(P.mat)
  rownames(P.mat) = ballots
  colnames(P.mat) = candidates

  ballot.prop.mat = V1.ballot.prop.mat = V0.ballot.prop.mat = NULL
  
  eu.by.ballot = t(P.mat%*%t(U))
  colnames(eu.by.ballot) = ballots
  V.mat = ballot.mat.from.eu.mat(eu.by.ballot,
                                 break.ties.with.sincerity = TRUE,
                                 sincere.vote.mat = V0)

  # optimal vote
  opt.votes.strategic = optimal.vote.from.V.mat(V.mat)
  opt.votes.sincere = optimal.vote.from.V.mat(V0)

  # calculate tau 
  tau = get.tau.from.eu.by.ballot.and.V0(eu.by.ballot, V0)
  
  # now output 
  list(opt.votes.strategic = opt.votes.strategic, opt.votes.sincere = opt.votes.sincere, piv.probs = pps, tau = tau, weights = weights, V.mat = V.mat, V0 = V0, eu.mat = eu.by.ballot)
  
}


out = sv(U_test, s = 85, v.vec = v_vec_test, rule = "plurality", ae_pack = TRUE)
out$piv.probs 
out$p.mat / 2
rowSums(out$p.mat)
colSums(out$V.mat)

out1 = sv(U_test, s = 85, v.vec = v_vec_test, rule = "plurality", ae_pack = FALSE)
out1$piv.probs
out1$p.mat
colSums(out1$V.mat)

out = sv(U_test, s = 85, v.vec = v_vec_test_3, rule = "plurality", ae_pack = TRUE)
out * 1000000
out$piv.probs
colSums(out$V.mat)

out1 = sv(U_test, s = 85, v.vec = v_vec_test_3, rule = "plurality", ae_pack = FALSE)
out1

out / out[1, 1]
out1 / out1[1, 1]




head(out1)
out1$piv.probs
colSums(out1$V.mat)

v_vec_test_3 = c(v_vec_test[1] + v_vec_test[2], v_vec_test[3] + v_vec_test[4], v_vec_test[5] + v_vec_test[6])
plurality3 <- plurality_election(k = 3, n = 400) %>%
    election_event_probs(method = "ev", alpha = (v_vec_test_3 * 85)) %>% combine_P_matrices()

plurality3
test = plurality3 %>% map("integral")
test = test[c(1, 2, 4)] %>% unlist + test[c(3, 5, 6)] %>% unlist
test
names(test) = c("AB", "AC", "BC")


plurality3 <- plurality_election(k = 3, n = 15000) 
alpha3 <- c(.4, .35, .25)*85

# compute pivot event probabilities using the Eggers-Vivyan method 
plurality3 %>% 
  election_event_probs(method = "ev", alpha = alpha3) -> pps
pps %>% combine_P_matrices


P_mat_at_pivotal_events(test, rule = "plurality", ballots = c("A", "B", "C"))




head(out)
out$piv.probs
out1 = sv(U_test, s = 85, v.vec = v_vec_test)
head(out1)


# Redefine function
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