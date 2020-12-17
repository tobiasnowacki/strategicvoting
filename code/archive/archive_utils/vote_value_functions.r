how.many.at.max = function(vec){
	length(which(vec == max(vec, na.rm = T)))
}

which.is.highest = function(vec){
    out = which(vec == max(vec, na.rm = T))
    if(length(out) == 1){
        out
    }else if(length(out) == 0){
        NA
    }else{
        sample(out)[1]
    }
}

convert_pp_vec_to_pp_mat = function(pp.vec, parties = c("lib", "con", "ndp", "bq", "grn")){
  out = matrix(0, nrow = length(parties), ncol = length(parties))
  colnames(out) = rownames(out) = parties
  for(i in 1:length(parties)){
    for(j in 1:length(parties)){
      if(i == j){next}
      this.name = paste0("p.", parties[i], ".", parties[j])
      if(this.name %in% names(pp.vec)){
        this = as.numeric(pp.vec[this.name])      
      }else{
        this.name = paste0("p.", parties[j], ".", parties[i])
        this = as.numeric(pp.vec[this.name])
      }
      out[i,j] = this 
    }
  }
  out
}

get_vote_values_and_tau_and_vote_code_for_voter_given_utilities_and_pivotal_prob_mat_and_vote_indicator_vec = function(u, ppm, viv, parties = c("lib", "con", "ndp", "bq", "grn")){
  stopifnot(rownames(ppm) == parties & colnames(ppm) == parties)
  stopifnot(names(u) == rownames(ppm))
  u = as.numeric(u)
  u.diff.mat = matrix(NA, nrow = length(u), ncol = 0)
  for(j in 1:length(u)){u.diff.mat = cbind(u.diff.mat, u - u[j])}
  vv = apply(u.diff.mat*ppm, 1, sum, na.rm = T)
  tau = max(vv[!is.na(u) & u != max(u, na.rm = T)], na.rm = T) - vv[!is.na(u) & u == max(u, na.rm = T)]
  # 1: voted for fave
  # 2: voted for max-v option among non-faves
  # 3: other 
  vote_code = as.integer(ifelse(u[viv == 1] == max(u, na.rm = T), 1, ifelse(vv[viv == 1] == max(vv[u != max(u, na.rm = T)], na.rm = T), 2, 3)))
  list(vv = vv, tau = tau, vote_code = vote_code)
}


calculate.vote.values.and.tau.and.vote.codes = function(parties = c("con", "lab", "ld", "snp", "pc", "ukip", "grn"), 	
	pivotality.suffix = "", feel.suffix = "post", output.suffix = ".post", df.name = "D5", old.tau.too = F){
	DF = get(df.name)
	# We create a new dataframe to fill up and then (outside of function) merge with D
	DD = as.data.frame(matrix(NA, nrow = nrow(DF), ncol = 0))
	# Calculate vote values for each party
	for(party in parties){
		vname = paste0("V.", party, output.suffix) 
		DD[[vname]] = 0
		for(party2 in parties){
			### get the pivotal probability
			if(party == party2){next}
			this.name = paste0("p.", party, ".", party2, pivotality.suffix)
			if(!this.name %in% names(DF)){this.name = paste0("p.", party2, ".", party, pivotality.suffix)}  # reverse if not in there -- e.g. if con.lab is not in there, try lab.con
			if(!this.name %in% names(DF)){cat("No pivotality measure for", party2, "and", party, "\n"); 1/0}  # if still not in there, raise an error
			### calculate p_{jk}(u_j - u_k)
			this.increment = DF[[this.name]]*(DF[[paste0(party, "feel", feel.suffix)]] - DF[[paste0(party2, "feel", feel.suffix)]])
			this.increment[is.na(this.increment)] = 0  # if one of the feeling scores is NA, we add or subtract nothing. 
			DD[[vname]] = DD[[vname]] + this.increment
		}
		DD[[vname]][is.na(DF[[paste0(party, "feel", feel.suffix)]])] = NA
	}
	# Now having calculated the values we calculate the maximum values
	# we'll use two matrices: the matrix of values
	V.mat = DD[,paste0("V.", parties, output.suffix)]  # this is just DD at this point.
	# and the matrix of indicators for which party is the respondent's fave, and which party the respondent voted for
	party1.mat = vote.mat = max.Vnot1.mat = matrix(NA, nrow = nrow(DF), ncol = 0)
	for(party in parties){
		# party1feel is the highest feeling score given. 
		# party1.mat is a logical matrix with n rows and k (number of parties) columns -- indicating which party is party1 for each individual 
		party1.mat = cbind(party1.mat, DF[[paste0("party1", feel.suffix)]] == party)
		vote.mat = cbind(vote.mat, DF$reg.vote.post == party)
	}
	DD[[paste0("highest.value.party", output.suffix)]] = parties[apply(V.mat, 1, which.is.highest)]
	DD[[paste0("how.many.at.highest.value", output.suffix)]] = apply(V.mat, 1, how.many.at.max)
	### V1 is the value of the respondent's favorite party
	V.mat2 = V.mat
	V.mat2[!party1.mat] = NA
	# this is NA for party's other than the respondent's fave. 
	DD[[paste0("V1", output.suffix)]] = apply(V.mat2, 1, sum, na.rm = T)
	### max.Vnot1 is the value of the highest-value party other than the respondent's fave. 
	vname = paste0("max.Vnot1", output.suffix)
	V.mat3 = V.mat
	V.mat3[party1.mat] = NA # set to NA for fave
	DD[[vname]] = apply(V.mat3, 1, max, na.rm = T) # the value of your favorite
	DD[[vname]][is.infinite(DD[[vname]])] = NA   # I think this is what happens when we have no non-NA values?
	DD[[paste0("highest.not.first.party", output.suffix)]] = parties[apply(V.mat3, 1, which.is.highest)] # the identiy of your favorite
	# Just count up how many are tied for first among non-fave parties
	DD[[paste0("number.highest.not.first", output.suffix)]] = apply(V.mat3, 1, how.many.at.max)
	for(party in parties){
		max.Vnot1.mat = cbind(max.Vnot1.mat, DD[[paste0("V.", party, output.suffix)]] == DD[[vname]])  # a matrix with T where the value for this party is equal to DD$max.Vnot1post
	}
	
	# I think these things are not needed anymore
	if(old.tau.too){
		# Now just get the elements to calculate tau
		max.Vnot1 = DD[[paste0("max.Vnot1", output.suffix)]]
		v1 = DD[[paste0("V1", output.suffix)]]
		tau.vname = paste0("tau", output.suffix)
		DD[[tau.vname]] = (max.Vnot1 - v1)/(max.Vnot1 + v1)  # calculate tau.
		DD[[tau.vname]][!is.na(DD[[tau.vname]]) & DD[[tau.vname]] < -1] = -1 # The argument is that you can always actually vote for someone with no chance -- spoil a ballot, whatever. Could equivalently set max.Vnot1 to 0 when it is negative. 
		DD[[tau.vname]][is.infinite(DD[[tau.vname]])] = NA		
	}
	
	## OK now let's deal with vote codes 
	DD[[paste0("vote.first", output.suffix)]] = apply(party1.mat*vote.mat, 1, sum, na.rm = T) # did you vote for your fave?  # this one will actually be the same across value measures
	DD[[paste0("vote.highest.not.first", output.suffix)]] = apply(max.Vnot1.mat*vote.mat, 1, sum, na.rm = T) # did you vote for a party that was tied for first among those NOT your fave? 
	DD[[paste0("vote.first", output.suffix)]][is.na(DF$reg.vote.post)] = DD[[paste0("vote.highest.not.first", output.suffix)]][is.na(DF$reg.vote.post)] = NA
	# turn the vote into a vote code
	DD[[paste0("vote.code", output.suffix)]] = NA
	DD[[paste0("vote.code", output.suffix)]][DD[[paste0("vote.first", output.suffix)]] == 1] = 1
	DD[[paste0("vote.code", output.suffix)]][DD[[paste0("vote.highest.not.first", output.suffix)]] == 1] = 2
	DD[[paste0("vote.code", output.suffix)]][DD[[paste0("vote.first", output.suffix)]] == 0 & DD[[paste0("vote.highest.not.first", output.suffix)]] == 0 & !is.na(DF$reg.vote.post)] = 3
	# outputs this thing, which we then merge with the DF.
	DD
}

test = F
if(test){
	
	cat("testing the calculation of vote values")
	# get the dataset 
	piv.df = "empty"
	analysis.type = "analytical"
	s = 85
	type = "actual"
	
	# get the dataset -- just part for testing
	dv5 = read.csv("./../data/processed/D_from_strategic_heterogeneity_v5.csv", as.is = T)
	dv5 = dv5[!is.na(dv5$confeelpost), ][1:100,]
	
	# get the pivs together
    for(year in c(2005, 2010, 2015)){
   
        these.pivs = read.csv(paste0("./../data/processed/pivotality_measures/", ifelse(analysis.type == "simulation", n, ""), "s_", s, "_", type, "_", year, ifelse(analysis.type == "analytical", "_analytical", ""), ".csv"), as.is = T)
		these.pivs$year = year
        if(class(piv.df) == "character"){piv.df = these.pivs}else{piv.df = rbind(piv.df, these.pivs)}
   
    }

	# merge them
    dddd = merge(dv5, piv.df, by = c("year", "refno"), all.x = T, suffixes = c("", paste0(".", s)))  
    
	# calculate the vote values
    dv5.a = cbind(dddd, calculate.vote.values.and.tau.and.vote.codes(pivotality.suffix = "", output.suffix = ".post", df.name = "dddd"))

	# now some tests
	stopifnot(sum(dv5.a$highest.not.first.party == dv5.a$party1post, na.rm = T) == 0)

}
