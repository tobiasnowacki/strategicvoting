### July 2021
### Script to get utility data frame for four parties
### 

# Path to CSES data
cses.location <- "data/cses/data/"  

# Load in original CSES data
c4 <- read.csv(paste0(cses.location, "cses4.csv"))
c3 <- read.csv(paste0(cses.location, "cses3.csv"))
c2 <- read.csv(paste0(cses.location, "cses2.csv"))
c1 <- read.csv(paste0(cses.location, "cses1.csv"))

# Set up variables
prefix <- c("A", "B", "C", "D") # prefixes for variables in the CSES waves
varnames <- c("3020", "3037", "3009", "3011") # the varname for like-dislike in each survey
ubound <- .5 # utility noise bound
big_list <- list()
bigger_list <- list()

##
## Define helper functions (usually in sv_helpers.R) 
##

# Function to get 1/0 mat from utilities
#   U: utility matrix
#   rule: "AV" or "plurality"
#   candidates: vector of cand names
sincere.vote.mat.from.U = function(
    U, 
    rule = "plurality", 
    candidates = c("A", "B", "C", "D")) {

    # do so under AV first
    if (rule %in% c("AV")) {
        # Create all possible ballot permutations
        ballots <- apply(
            permutations(
                n = length(candidates),
                r = length(candidates),
                v = candidates,
                repeats.allowed = F
            ),
            1,
            paste,
            collapse = ""
        )
        
        # Create sincere P.mat
        sincere.P <- matrix(
            NA,
            nrow = length(ballots),
            ncol = length(candidates)
        )
        
        # Cycle through ballot-candidate options
        for (i in 1:length(ballots)) {
            for (j in 1:length(candidates)) {
                sincere.P[i, j] = length(candidates) - which(
                    grepl(candidates[j], str_split(ballots[i], "")[[1]])
                )
            }
        }
    # Deal with plurality case instead
    } else { 
        ballots <- candidates
        sincere.P <- diag(length(candidates))
    }

    # Create EU matrix
    sincere.eu.by.ballot <- t(sincere.P %*% t(U))
    colnames(sincere.eu.by.ballot) <- ballots # rename cols

    # Recreate ballot matrix from EU mat (and return)
    ballot.mat.from.eu.mat(
        sincere.eu.by.ballot,
        break.ties.with.sincerity = FALSE
    )
}

# Legacy function to create EU mat -> Ballot mat
ballot.mat.from.eu.mat = function(
    eu.mat, 
    break.ties.with.sincerity = TRUE, 
    sincere.vote.mat = NULL) {

    # If numerical ties w/ sincerity ought to be resolved
    if (break.ties.with.sincerity) {
        if (is.null(sincere.vote.mat)) {
            cat("you must pass a sincere.vote.mat!\n")
        }
        # add a little bit of extra weight to sincere option
        eu.mat <- eu.mat + sincere.vote.mat * (10^(-10))
    }

    # Maximum EU for each voter
    max.eus = apply(eu.mat, 1, max, na.rm = TRUE)
    # Create ballot matrix
    ballot.mat = matrix(NA, ncol = ncol(eu.mat), nrow = nrow(eu.mat))
    # Fill in ballot matrix (candidate by candidate)
    for (j in 1:ncol(eu.mat)) {
        ballot.mat[, j] = as.integer(eu.mat[, j] == max.eus)
    }

    # Additional trouble handling 
    # not relevant in the CSES utils case since all we want is the sincere ballot mat.
    if (break.ties.with.sincerity) {
        if (is.null(sincere.vote.mat)) {
            cat("you must pass a sincere.vote.mat!\n")
        }
        # get rows where ballot.mat has multiple entries
        multi_ballot <- apply(ballot.mat, 1, sum) > 1
        # replace with sincere vote
        ballot.mat[multi_ballot, ] <- sincere.vote.mat[multi_ballot, ]
        # problem here: no guarantee that ties include sincere option...
    }

    # Give columns appropriate names
    colnames(ballot.mat) = colnames(eu.mat)
    return(ballot.mat)
}


# Function to get new v_vec from vote matrix and weights
ballot.props.from.vote.mat.and.weights <- function(V, weights) {
    
    # Create weight matrix
    weight.mat <- matrix(weights,
        nrow = nrow(V),
        ncol = ncol(V),
        byrow = FALSE
    )
    
    # Fill in with 0 where NA
    weight.mat[is.na(V)] <- 0
    out <- apply(V * weight.mat, 2, sum, na.rm = TRUE) / apply(weight.mat, 2, sum, na.rm = TRUE)
    out / sum(out) # not sure why I have to normalize, but there it is
}

## 
## Run through data
##

# For-loop
for(j in 4:1){

    # Set up data and grab cases
    d = get(paste0("c", j))
    p = prefix[j]
    case_name = paste0(p, "1004")
    varname = varnames[j]
    cases = unique(as.character(d[[case_name]]))

    # For the cases within partial data, run loop
    for(i in 1:length(cases)){
        # Get utility data
        X = d[
                d[[case_name]] == cases[i],
                paste0(p, varname, "_", c("A", "B", "C", "D"))
        ]
        X[X > 10] = NA # replacing 99 with NA
        # Get weights
        weights = d[d[[case_name]] == cases[i], paste0(p, "1010_3")]
        # add some noise
        X = X + runif(n = nrow(X)*ncol(X), min = -ubound, max = ubound)
        # determine order of finish
        V0 = sincere.vote.mat.from.U(
            X,
            rule = "plurality", 
            candidates = c("A", "B", "C", "D")
        )
        # get v.vec
        v.vec = ballot.props.from.vote.mat.and.weights(V0, weights = weights) + runif(n = 4, min = 0, max = .0001) # tie always possible
        v.vec = v.vec/sum(v.vec)
        # reorganize so columns are in order of finish
        X = X[, order(v.vec, decreasing = T)]
        # rename so A is the expected winner
        colnames(X) = c("A", "B", "C", "D")
        if(length(unique(X[, 4])) == 1){next}
        # print status
        cat(" == ", cases[i], " == ")
        big_list[[cases[i]]] = list("U" = X, "weights" = weights)
    }
    cat("\n")
}

# Finishing touches and export
cases <- sort(names(big_list))
save(big_list, cases, file = "data/cses/data/big_list_4_parties.RData")

