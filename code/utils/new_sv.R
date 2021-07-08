source("code/utils/new_sv_helpers.R")

# Updated as of October 2020 to include pivotprobs package from AE.
sv = function(U, 
              weights = NULL, 
              v.vec = NULL, 
              s, 
              rule = "plurality", 
              V0 = NULL,
              sin_pref_mat = NULL, 
              ae_pack = NULL){
  
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
    # Combine v.vec for checking what method to run
    v.vec.tri = c(v.vec[1] + v.vec[2], v.vec[3] + v.vec[4], v.vec[5] + v.vec[6])

    # Default: Eggers-Nowacki method
    if(sum(v.vec.tri > 0.005) == 3){
     pps = irv_election(n = 1) %>%
             election_event_probs(method = "en",  alpha = (v.vec * s))
    }
    # if close to edge of ternary, run MC simulation method instead     
    if(sum(v.vec.tri > 0.005) < 3){
      
      pps = irv_election(n = 10000) %>%
           election_event_probs(method = "mc",  
                                num_sims = 400000,
                                alpha = (v.vec * s))
      
      # multiply integral by n
      pps = lapply(pps, function(x) {
        x$integral = x$integral * 10000
        return(x)
        } 
      )
    }
    
    # Combine P matrix
    P.mat <- pps %>%
        combine_P_matrices()
    
    # Get rid of abstentions
    P.mat = P.mat[, -7]

    # Save pivot probabilities
    pps = map(pps, "integral")
    pps = pps[c("a_b", "a_c", "b_c", "a_b|ab", "a_b|ac", "a_b|cb", "a_c|ac", "a_c|ab", "a_c|bc", "b_c|bc", "b_c|ba", "b_c|ac")]
    names(pps) = c("AB", "AC", "BC", "AB.AB", "AB.AC", "AB.CB", "AC.AC", "AC.AB", "AC.BC", "BC.BC", "BC.BA", "BC.AC")

  }else if(rule == "plurality"){
    pps = plurality_election(k = 3, n = 1) %>%
        election_event_probs(method = "ev", 
                             alpha = (v.vec * s))
    P.mat = pps %>% combine_P_matrices()
    P.mat = P.mat[, -4]
    pps = map(pps, "integral")
  }

  # Normalise
  # probability.pivotal = sum(P.mat[,1])
  normalized.P.mat = P.mat # /probability.pivotal

  ballot.prop.mat = V1.ballot.prop.mat = V0.ballot.prop.mat = NULL
  
  # Get EU matrix
  eu.by.ballot = as.matrix(U) %*% normalized.P.mat
  colnames(eu.by.ballot) = ballots

  rownames(P.mat) = candidates
  colnames(P.mat) = ballots
  colnames(eu.by.ballot) = ballots

  # Get optimal vote matrix (w/ sincerity adjustment)
  if(is.null(sin_pref_mat)){
    sin_pref_mat = sincere_pref_mat_from_U(U, rule = rule)
  }
  V.mat = ballot_mat_from_eu_mat(eu.by.ballot,  break_ties_with_sincerity = TRUE, 
    sincere_mat = sin_pref_mat)

  # optimal vote
  opt.votes.strategic = optimal.vote.from.V.mat(V.mat)
  opt.votes.sincere = optimal.vote.from.V.mat(V0)

  # calculate tau 
  tau = get.tau.from.eu.by.ballot.and.V0(eu.by.ballot, V0, sincere_mat = sin_pref_mat)

  # get best response v.vec
  br.v.vec = ballot.props.from.vote.mat.and.weights(V.mat, weights = weights)
  
  # now output 
  return(
    list(opt.votes.strategic = opt.votes.strategic, 
       opt.votes.sincere = opt.votes.sincere, 
       piv.probs = pps, 
       tau = tau, 
       weights = weights, 
       V.mat = V.mat, 
       V0 = V0, 
       eu.mat = eu.by.ballot, 
       best.response.v.vec = br.v.vec,
       v.vec.before = v.vec,
       p.mat = P.mat)
  )
}
