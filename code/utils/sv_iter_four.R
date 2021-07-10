# Function to run iterations

source("code/utils/sv_four.R")
source("code/utils/new_sv_helpers.R")

sv_iter <- function(U,
                    s,
                    starting.v.vec = NULL,
                    weights = NULL,
                    rule = "plurality",
                    lambda = .1,
                    epsilon = .001,
                    max.iterations = 200,
                    until.convergence = TRUE,
                    sincere.proportion = 0,
                    candidates = c("a", "b", "c", "d"),
                    ballots = c("abc", "acb", "bac", "bca", "cab", "cba"),
                    the.floor = .01,
                    noisy = FALSE) {
    
    # default handling
    if (rule == "plurality") {
        ballots <- candidates
    }
    if (is.null(weights)) {
        weights <- rep(1, nrow(U))
    }

    sincere.vote.mat <- sincere.vote.mat.from.U(
        U,
        rule = rule, 
        candidates = candidates
    )
    sincere.pref.mat <- sincere_pref_mat_from_U(
        U,
        rule = rule, 
        candidates = candidates
    )
    
    # get starting vvec
    if (is.null(starting.v.vec)) {
        sincere.v.vec <- ballot.props.from.vote.mat.and.weights(
            sincere.vote.mat, 
            weights
        )
    } else {
        sincere.v.vec <- starting.v.vec
    }
    
    # Run first iteration
    out <- list()
    out[[1]] <- sv_four(
        U = U,
        weights = weights,
        v.vec = sincere.v.vec,
        s = s,
        rule = rule,
        V0 = sincere.vote.mat,
    )

    # Run subsequent iterations
    for (i in 2:max.iterations) {
        if (noisy) {
            cat("Iteration ", i, "\n", sep = "")
        }
        # response vectors
        strategic.v.vec.i <- lambda * out[[i - 1]]$best.response.v.vec +
            (1 - lambda) * out[[i - 1]]$v.vec.before
        overall.v.vec.i <- sincere.proportion * sincere.v.vec +
            (1 - sincere.proportion) * strategic.v.vec.i

        # Replace this below with sv function
        out[[i]] <- sv_four(
            U = U,
            weights = weights,
            v.vec = overall.v.vec.i,
            s = s,
            rule = rule,
            V0 = sincere.vote.mat,
        )
    }
    return(out)
}

