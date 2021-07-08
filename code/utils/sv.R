#### sv

# function to get strategic voting info given utility matrix, weights, and ingredients of (Dirichlet) belief
# a slimmed-down version of iteration_simulation() from general_iteration_simulation_approach.R


# Updated as of October 2020 to include pivotprobs package from AE.
sv <- function(
  U, 
  weights = NULL, 
  v.vec = NULL, 
  s, 
  rule = "plurality", 
  V0 = NULL, 
  ae_pack = TRUE
) {
  stopifnot(!is.null(colnames(U)))
  candidates <- sort(colnames(U))
  U <- U[, candidates]

  stopifnot(rule %in% c("plurality", "AV"))
  stopifnot(ncol(U) == 3) # we can't do more than three with this technique -- expanding would require generalizing the P_mat_at_pivotal_events step below

  if (is.null(weights)) {
    weights <- rep(1, nrow(U))
  }

  # get sincere vote mat (V0)

  if (is.null(V0)) {
    V0 <- sincere.vote.mat.from.U(U, rule = rule, candidates = candidates)
  }

  if (is.null(v.vec)) {
    v.vec <- ballot.props.from.vote.mat.and.weights(V0, weights = weights)
  }

  if (rule %in% c("AV")) {
    ballots <- apply(
      permutations(
        n = length(candidates),
        r = length(candidates),
        v = candidates, 
        repeats.allowed = FALSE),
      1, 
      paste,
      collapse = ""
    )
  } else {
    ballots <- candidates
  }

  # get pivotal probabilities -- depends on the rule
  if (rule == "AV") {
    v.vec.tri <- c(v.vec[1] + v.vec[2], v.vec[3] + v.vec[4], v.vec[5] + v.vec[6])
    if (sum(v.vec.tri > 0.005) == 3) {
      pps <- irv_election(n = 1) %>%
        election_event_probs(method = "en", alpha = (v.vec * s))
    }
    if (sum(v.vec.tri > 0.005) < 3) {
      # if close to edge of ternary, run MC simulation method instead
      pps <- irv_election(n = 10000) %>%
        election_event_probs(
          method = "mc",
          num_sims = 400000,
          alpha = (v.vec * s)
        )
      # multiply integral by n
      pps <- lapply(pps, function(x) {
        x$integral <- x$integral * 10000
        return(x)
      })
    }
    P.mat <- pps %>%
      combine_P_matrices()
    P.mat <- P.mat[, -7]
    pps <- map(pps, "integral")
    pps <- pps[c("a_b", "a_c", "b_c", "a_b|ab", "a_b|ac", "a_b|cb", "a_c|ac", "a_c|ab", "a_c|bc", "b_c|bc", "b_c|ba", "b_c|ac")]
    names(pps) <- c("AB", "AC", "BC", "AB.AB", "AB.AC", "AB.CB", "AC.AC", "AC.AB", "AC.BC", "BC.BC", "BC.BA", "BC.AC")
  } else if (rule == "plurality") {
    pps <- plurality_election(k = 3, n = 1) %>%
      election_event_probs(method = "ev", alpha = (v.vec * s))
    P.mat <- pps %>% combine_P_matrices()
    P.mat <- P.mat[, -4]
    pps <- map(pps, "integral")
  }
  probability.pivotal <- sum(P.mat[, 1])
  normalized.P.mat <- P.mat / probability.pivotal

  ballot.prop.mat <- V1.ballot.prop.mat <- V0.ballot.prop.mat <- NULL

  eu.by.ballot <- as.matrix(U) %*% normalized.P.mat
  colnames(eu.by.ballot) <- ballots

  rownames(P.mat) <- candidates
  colnames(P.mat) <- ballots
  colnames(eu.by.ballot) <- ballots
  sincere_pref_mat <- sincere_pref_mat_from_U(U, rule = rule)
  V.mat <- ballot_mat_from_eu_mat(eu.by.ballot,
    break_ties_with_sincerity = TRUE,
    sincere_mat = sincere_pref_mat
  )
  # V.mat = ballot.mat.from.eu.mat(eu.by.ballot,
  # break.ties.with.sincerity = TRUE,
  # sincere.vote.mat = V0)

  # optimal vote
  opt.votes.strategic <- optimal.vote.from.V.mat(V.mat)
  opt.votes.sincere <- optimal.vote.from.V.mat(V0)

  # calculate tau
  tau <- get.tau.from.eu.by.ballot.and.V0(eu.by.ballot, V0, sincere_mat = sincere_pref_mat)

  # now output
  list(opt.votes.strategic = opt.votes.strategic, opt.votes.sincere = opt.votes.sincere, piv.probs = pps, tau = tau, weights = weights, V.mat = V.mat, V0 = V0, eu.mat = eu.by.ballot, p.mat = P.mat)
}


# one question is whether we want to assess the success of the *ex ante* CW, or the ex post.
# depends on how we think about the uncertainty.
## one view is that we know the true type proportions (properties of the DGP), and results differ from these due to e.g. turnout variation. But this isn't the right way to do it. We have the results of a poll, and we think the true type proportions might be different from what's in the poll -- if the poll was wrong you don't want the system that recovers the CW in the poll.

evaluate_success_of_CW_given_U_and_V.mat <- function(U, V.mat, V0 = NULL, lambdas = c(.1, .2, .5), weights = NULL, rule = "plurality", m = 300, M = 1000) {
  if (is.null(weights)) {
    weights <- rep(1, nrow(U))
  }

  stopifnot(!is.null(colnames(U)))
  if (is.null(V0)) {
    V0 <- sincere.vote.mat.from.U(U, rule = rule, candidates = colnames(U))
  }

  stopifnot(nrow(U) == nrow(V.mat))
  if (rule == "AV") {
    stopifnot(ncol(V.mat) == 6)
  }
  if (rule == "plurality") {
    stopifnot(ncol(V.mat) == 3)
  }

  cw.mat <- matrix(NA, nrow = 0, ncol = 3)

  lambda.v.vec.mats <- list()
  for (lambda in lambdas) {
    lambda.v.vec.mats[[paste0("lambda_", lambda)]] <- matrix(NA, nrow = 0, ncol = ifelse(rule == "AV", 6, 3))
  }

  for (i in 1:M) {
    indexes <- sample(1:nrow(U), size = m, replace = T, prob = weights)
    this.U <- U[indexes, ]
    cw.vec <- condorcet.winner(this.U, weights = weights[indexes])
    cw.mat <- rbind(cw.mat, as.integer(cw.vec))
    for (lambda in lambdas) {
      votes.strategically <- rep(F, m)
      if (lambda > 0) {
        votes.strategically[1:floor(lambda * m)] <- T
      }
      votes.strategically <- sample(votes.strategically)
      strategic.vote.mat <- V0[indexes, ] # start with sincere votes
      strategic.vote.mat[votes.strategically, ] <- V.mat[indexes, ][votes.strategically, ] # stick in strategic votes
      strategic.v.vec <- ballot.props.from.vote.mat.and.weights(strategic.vote.mat, weights[indexes])
      # add a bit of noise to break ties
      strategic.v.vec <- strategic.v.vec + runif(length(strategic.v.vec), min = 0, max = .00001)
      lambda.v.vec.mats[[paste0("lambda_", lambda)]] <- rbind(lambda.v.vec.mats[[paste0("lambda_", lambda)]], strategic.v.vec)
    }
  }

  out <- list()
  for (lam in 1:length(lambdas)) {
    strategic.vps.mat <- victory.probs.from.sims(lambda.v.vec.mats[[lam]], rule = rule, return.matrix = T)
    cw.winner <- c(as.numeric(apply(strategic.vps.mat * cw.mat, 1, function(x) 1 %in% x)))
    mu <- mean(cw.winner)
    out[[lam]] <- c(mu)
  }
  out <- do.call(rbind, out)
  # out[, 4] <- lambdas
  return(cbind(out, lambdas))
}

no_show_non_mon_from_sv_object <- function(sv.obj) {
  ns.mat <- cbind(rep(sv.obj$piv.probs$AB.CB, nrow(sv.obj$V0)), sv.obj$piv.probs$AC.BC, sv.obj$piv.probs$AB.AC, sv.obj$piv.probs$BC.AC, sv.obj$piv.probs$AC.AB, sv.obj$piv.probs$BC.BA)
  ns.prob.vec <- apply(ns.mat * sv.obj$V0, 1, sum, na.rm = T)
  avg.ns.prob <- sum(ns.prob.vec * sv.obj$weights, na.rm = T) / sum(sv.obj$weights[!is.na(ns.prob.vec)], na.rm = T)

  nm1.mat <- cbind(rep(sv.obj$piv.probs$BC.AC, nrow(sv.obj$V0)), sv.obj$piv.probs$BC.BA, sv.obj$piv.probs$AC.BC, sv.obj$piv.probs$AC.AB, sv.obj$piv.probs$AB.CB, sv.obj$piv.probs$AB.AC)
  nm1.prob.vec <- apply(nm1.mat * sv.obj$V0, 1, sum, na.rm = T)
  avg.nm1.prob <- sum(nm1.prob.vec * sv.obj$weights, na.rm = T) / sum(sv.obj$weights[!is.na(nm1.prob.vec)], na.rm = T)

  nm2.mat <- nm1.mat[, c(2, 1, 4, 3, 6, 5)]
  nm2.prob.vec <- apply(nm2.mat * sv.obj$V0, 1, sum, na.rm = T)
  avg.nm2.prob <- sum(nm2.prob.vec * sv.obj$weights, na.rm = T) / sum(sv.obj$weights[!is.na(nm2.prob.vec)], na.rm = T)

  list(avg.ns.prob = avg.ns.prob, avg.nm1.prob = avg.nm1.prob, avg.nm2.prob = avg.nm2.prob, prob.sum = sum(unlist(sv.obj$piv.probs)), fr.prob.sum = sum(unlist(sv.obj$piv.probs)) - sum(unlist(sv.obj$piv.probs))[c("AB", "AC", "BC")])
}

plurality_wasted_vote_from_sv_object <- function(sv.obj) {
  wv.mat <- cbind(rep(sv.obj$piv.probs$BC, nrow(sv.obj$V0)), sv.obj$piv.probs$AC, sv.obj$piv.probs$AB)
  wv.prob.vec <- apply(wv.mat * sv.obj$V0, 1, sum, na.rm = T)
  avg.wv.prob <- sum(wv.prob.vec * sv.obj$weights, na.rm = T) / sum(sv.obj$weights[!is.na(wv.prob.vec)], na.rm = T)

  list(avg.wv.prob = avg.wv.prob, prob.sum = sum(unlist(sv.obj$piv.probs)))
}