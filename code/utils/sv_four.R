#### sv
## Adapted for four candidates
# function to get strategic voting info given
# utility matrix, weights, and ingredients of (Dirichlet) belief
library(pivotprobs)

source("code/utils/new_sv_helpers.R")

sv_four <- function(U,
                    weights = NULL,
                    v.vec = NULL,
                    s,
                    rule = "plurality",
                    V0 = NULL,
                    ae_pack = TRUE) {
    stopifnot(!is.null(colnames(U)))
    candidates <- sort(colnames(U))
    U <- U[, candidates]

    stopifnot(rule %in% c("plurality", "AV"))
    stopifnot(ncol(U) == 4)

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
                repeats.allowed = FALSE
            ),
            1,
            paste,
            collapse = ""
        )
    } else {
        ballots <- candidates
    }

    # IRV pivprobs
    if (rule == "AV") {
        # Step 1: Simulate results
        mc_sims <- simulate_ordinal_results_from_dirichlet(
            k = 4,
            n = 200000,
            alpha = v.vec * s
        )
        # Step 2: Calculate pprobs
        pps <- mc_sims %>%
<<<<<<< HEAD
<<<<<<< HEAD
            irv_pivot_probs_four_cands(n = 1, reporting = 0)
=======
=======
>>>>>>> 95425d540f4efbe579af50da3c43f21f50cc9383
            irv_pivot_probs_four_cands(n = 1000, reporting = 0)

        pps <- lapply(pps, function(x) {
            x$integral <- x$integral * 1000
            return(x)
        })
        
<<<<<<< HEAD
>>>>>>> 95425d540f4efbe579af50da3c43f21f50cc9383
=======
>>>>>>> 95425d540f4efbe579af50da3c43f21f50cc9383
        # Step 3: Get pmat
        P.mat <- pps %>%
            combine_P_matrices()
        P.mat <- P.mat[, -25] # delete abstention
    } else if (rule == "plurality") {
        pps <- plurality_election(k = 4, n = 1000) %>%
            election_event_probs(
                method = "sc",
                alpha = (v.vec * s),
                drop_dimension = TRUE,
                merge_adjacent_pivot_events = TRUE
            )
        
        # Divide  to get pivot for each voter
        pps <- lapply(pps, function(x) {
            x$integral <- x$integral * 1000
            return(x)
        })
      
        P.mat <- pps %>% combine_P_matrices()
        P.mat <- P.mat[, -5]
        pps <- map(pps, "integral")
    }

    # Normalize (this is old code)
    probability.pivotal <- sum(P.mat[, 1])
    normalized.P.mat <- P.mat

    # Set up ballot matrices
    ballot.prop.mat <- V1.ballot.prop.mat <- V0.ballot.prop.mat <- NULL
    eu.by.ballot <- as.matrix(U) %*% normalized.P.mat
    colnames(eu.by.ballot) <- ballots
    rownames(P.mat) <- candidates
    colnames(P.mat) <- ballots
    colnames(eu.by.ballot) <- ballots
    sincere_pref_mat <- sincere_pref_mat_from_U(
        U,
        rule = rule, 
        candidates = candidates
    )
    V.mat <- ballot_mat_from_eu_mat(
        eu.by.ballot,
        break_ties_with_sincerity = TRUE,
        sincere_mat = sincere_pref_mat
    )

    # determine optimal vote
    opt.votes.strategic <- optimal.vote.from.V.mat(V.mat)
    opt.votes.sincere <- optimal.vote.from.V.mat(V0)

    # calculate tau
    tau <- get.tau.from.eu.by.ballot.and.V0(
        eu.by.ballot,
        V0,
        sincere_mat = sincere_pref_mat
    )

      # get best response v.vec
    br.v.vec <- ballot.props.from.vote.mat.and.weights(
      V.mat,
      weights = weights
    )
      
    # now output
    return(
        list(
           opt.votes.strategic = opt.votes.strategic,
           opt.votes.sincere = opt.votes.sincere, 
           piv.probs = pps, 
           tau = tau, 
           weights = weights, 
           V.mat = V.mat, 
           V0 = V0, 
           eu.mat = eu.by.ballot, 
           best.response.v.vec = br.v.vec,
           v.vec.before = v.vec,
<<<<<<< HEAD
<<<<<<< HEAD
           p.mat <- P.mat
=======
           p.mat = P.mat
>>>>>>> 95425d540f4efbe579af50da3c43f21f50cc9383
=======
           p.mat = P.mat
>>>>>>> 95425d540f4efbe579af50da3c43f21f50cc9383
        )
    )
}


