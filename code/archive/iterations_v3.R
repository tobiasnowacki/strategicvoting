library(here)
# getting code
source(here("code", "utils", "general_iteration_simulation_approach.r")) 
source(here("code", "utils", "sv.r"))
source(here("code/utils/functions.r"))


# country VAP and number of cases 
vap = read.table(here("/data/case_vap.csv"), header = T)
vap$case_weight = normalize.to.1(vap$VAP)/vap$Freq # this is the weight that should be assigned to each election for this country

# get the CSES data 
load(here("output/big_list_2.RData"))
date = "20190612"

sorted.cases = sort(names(big_list))
sorted.cases = sorted.cases[!sorted.cases %in% c("BLR_2008", "LTU_1997")]

# number of countries
length(unique(gsub("_\\d+", "", sorted.cases)))

# describe sample sizes excluding people with missing data   
samps = rep(NA, length(sorted.cases))
for(i in 1:length(sorted.cases)){
  samps[i] = nrow(na.omit(big_list[[sorted.cases[i]]]$U))
}

case_vap_weight = rep(NA, length(sorted.cases))
for(i in 1:length(sorted.cases)){
  cc = paste(strsplit(sorted.cases[i], "")[[1]][1:3], collapse = "")
  case_vap_weight[i] = vap$case_weight[vap$cntry == cc]
}

# set parameters
ss = seq(10, 85, 15)
rules = c("plurality", "AV")

lambda = .05
m = 60
the.floor = 0

# vap stuff
vap_weight_mat = matrix(case_vap_weight, nrow = m, ncol = length(sorted.cases), byrow = T)
vap_weight_mat_uk_only = vap_weight_mat
for(j in 1:ncol(vap_weight_mat)){
  if(grepl("GBR", sorted.cases[j])){
    vap_weight_mat_uk_only[,j] = 1/3
  }else{
    vap_weight_mat_uk_only[,j] = 0
  }
}

# run loop
big_out = list()

for(case in sorted.cases){
  cat("=== ", case, " ===\n", sep = "")
  this.U = as.matrix(big_list[[case]]$U)
  colnames(this.U) = c("a", "b", "c")
  this.weight = big_list[[case]]$weights
  big_out[[case]] = list()
  for(rule in rules){
    cat(rule, ": ", sep = "")
    big_out[[case]][[rule]] = list()
    for(s in ss){
      cat(s, " ")
      big_out[[case]][[rule]][[paste0("s=", s)]] = iterated_best_response_sequence(U = this.U, s = s, weights = this.weight, rule = rule, lambda = lambda, until.convergence = F, max.iterations = m, sincere.proportion = 0, candidates = c("a", "b", "c"), ballots = c("abc", "acb", "bac", "bca", "cab", "cba"), the.floor = the.floor, noisy = F)
    }
    cat("\n")
  }
  cat("\n")
}


