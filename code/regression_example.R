##################################################
## Project: Strategic Voting in RCV
## Script purpose: Describing cases
## Date: 16/12/2018
## Author:
##################################################

### 
### Dependencies
###

library(here)
library(gtools)
library(ggtern)
library(stargazer)

source(here("utils/functions.r"))

### 
### Create ballot profiles
###
n <- 1000 # No. of cases
s_list <- list(85) # Level of info

hyper_v_vec <- c(0.30, 0.13, 0.185, 0.185, 0.03, 0.17)
v_vec_df <- as.data.frame(rdirichlet(n, 40 * hyper_v_vec))

ggtern(v_vec_df, aes(V1 + V2, V3 + V4, V5 + V6)) +
  geom_point(alpha = 0.05)

v_vec_df$mAB <- v_vec_df[, 1] / (v_vec_df[, 1] + v_vec_df[, 2])
v_vec_df$mBA <- v_vec_df[, 3] / (v_vec_df[, 3] + v_vec_df[, 4])
v_vec_df$mCB <- v_vec_df[, 6] / (v_vec_df[, 5] + v_vec_df[, 6])

### Check strategic voting incentives


voter_vec <- c("ABC", "ACB", "BAC", "BCA", "CAB", "CBA")
beta <- 0.5

voter_df <- data.frame(A = c(1, 1, beta, 0, beta, 0), 
                       B = c(beta, 0, 1, 1, 0, beta), 
                       C = c(0, beta, 0, beta, 1, 1))

case_list <- list()
for(i in 1:n){
  print(i)
  obj <- return_sv_tau(as.numeric(c(v_vec_df[i, 1:6], 0, 0, 0)), voter_df, s_list)
  obj$voter <- voter_vec
  case_list[[i]] <- obj
}

case_df <- as.data.frame(do.call(rbind, case_list))

### Build dataframe
v_vec_df[, 10:15] <- sapply(voter_vec, function(x) case_df$tau_rcv[case_df$voter == x] > 0)
names(v_vec_df) <- c("vABC", "vACB", "vBAC", "vBCA", "vCAB", "vCBA", "mAB", "mBA", "mCB", "svABC", "svACB", "svBAC", "svBCA", "svCAB", "svCBA")

### Run regression

mod_ABC <- lm("svABC ~ I(vABC + vACB) + I(vBAC + vBCA) + I((vABC + vACB) * (vBAC + vBCA))  + mAB + mCB", v_vec_df)
summary(mod_ABC)

mod_ACB <- lm("svACB ~ I(vABC + vACB) + I(vBAC + vBCA) + I((vABC + vACB) * (vBAC + vBCA)) + mAB + mCB", v_vec_df)
summary(mod_ACB)

ggtern(v_vec_df, aes(vABC + vACB, vBAC + vBCA, vCAB + vCBA)) +
  geom_point(aes(colour = svABC), alpha = 0.1)
ggsave("../output/figures/regression/abc.pdf")

ggtern(v_vec_df, aes(vABC + vACB, vBAC + vBCA, vCAB + vCBA)) +
  geom_point(aes(colour = svACB), alpha = 0.1)
ggsave("../output/figures/regression/acb.pdf")