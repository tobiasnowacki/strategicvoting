## Script to generate conjecture tests in results section.
## 11 September 2019.
## Toby Nowacki

### -------------------------
### LOAD DEPENDENCIES -------

# install.packages("here")
library(here) 						# to get dir
source(here("code/full_header.R")) 	# fn's and data
source(here("code/prep_cses.R")) 	# data prep


### ------------------
### SET VALUES -------

# Precision parameter (s)
s_list <- as.list(c(10, 25, 40, 55, 70, 85))

# Strategicness parameter (λ)
lambda_list <- as.list(c(0.05, 0.1, 0.01))

# Convergence distance criterion (ε)
epsilon_thresh <- 0.001

# Maximum number of iterations
k_max <- 300

### ---
### RUN ITERATIONS

# for(lambda_val in 1:3){
for(lambda_val in 1){
	lambda <- lambda_list[[lambda_val]]

	cat(paste0("=== Beginning loop for lambda == ", lambda, ". === \n"))

	# for(s_loop_val in c(1)){ 
	for(s_loop_val in c(6)){ 

	  s_val <- s_list[[s_loop_val]]

	  cat(paste0("=== Beginning loop for s == ", s_val, ". === \n"))

	  which_cases <- 1:160   # default is 1:160
	  # max_iter_val <- 250	
	  max_iter_val <- 60
	  # which_cases <- 1:3        # for test purposes only!

	  cat("Running parallel iterations. \n")

	  cl <- makeCluster(4)
	  registerDoParallel(cl)
	  # parallelise
	  cases_converge <- foreach(case = which_cases, 
	              .packages = c("gtools", "stringr", "tidyverse",
	                            "questionr")
	              ) %dopar% {
	              out <- many_iterations_until_convergence(
	                    big_list_na_omit[[case]], 
	                    big_list_na_omit[[case]]$v_vec, 
	                    lambda, 
	                    s_val, 
	                    0.001, 
	                    max_iter_val)
	              # sink("log_convergence.txt", append = TRUE)
	              # cat(paste(case, " done.", "\n"))
	              out
	  }

	  stopCluster(cl)

	  names(cases_converge) <- names(big_list_na_omit)[which_cases]


	  cat("Done. \n")

	}

}


# Conjecture tests
big_rcv_sum <- list()
big_plur_sum <- list()
for(i in 1:160){
	print(i)
	mini_list1 <- cases_converge[[i]][[1]]
	mini_list1$case <- names(big_list_na_omit)[i]
	mini_list1$system <- "IRV"
	big_rcv_sum[[i]] <- mini_list1
	mini_list2 <- cases_converge[[i]][[3]]
	mini_list2$case <- names(big_list_na_omit)[i]
	mini_list2$system <- "Plurality"
	big_plur_sum[[i]] <- mini_list2
}

big_rcv_sum <- do.call(rbind, big_rcv_sum)
big_plur_sum <- do.call(rbind, big_plur_sum)


## CONJECTURE 1a
# Summed pivotal probabilities
rcv_summary <- big_rcv_sum %>% group_by(case, s, iter) %>% summarise_at(vars(AB:BCp), first) 
rcv_summary$rcv_sum <- rowSums(rcv_summary[, 4:16])
rcv_summary$plur_sum <- rowSums(rcv_summary[, 7:19])
rcv_summary <- rcv_summary %>% mutate(diff = plur_sum - rcv_sum)

# Do the same thing for plurality
plur_summary <- big_plur_sum %>% group_by(case, s, iter) %>% summarise_at(vars(AB:BCp), first) 
plur_summary$rcv_sum <- rowSums(plur_summary[, 4:16])
plur_summary$plur_sum <- rowSums(plur_summary[, 7:19])
plur_summary <- plur_summary %>% mutate(diff = plur_sum - rcv_sum)

cmat_plur <- matrix(c(1, 2, 3, 2, 1, 3, 1, 3, 2, 3, 1, 2, 2, 3, 1, 3, 2, 1), byrow = T, ncol = 3)

cmat_rcv <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 2, 1, 3, 7, 8, 9, 4, 5, 6, 10, 12, 11, 1, 3, 2, 4, 6, 5, 10, 11, 12, 7, 8, 9, 3, 1, 2, 10, 11, 12, 4, 6, 5, 7, 9, 8, 2, 3, 1, 7, 9, 8, 10, 12, 11, 4, 5, 6, 3, 2, 1, 10, 12, 11, 7, 9, 8, 4, 6, 5), byrow = T, ncol = 12)

cmat_u <- matrix(c(1, 2, 3, 1, 3, 2, 2, 1, 3, 2, 3, 1, 3, 1, 2, 3, 2, 1), byrow = T, ncol = 3)

# re-order RCV DF
rcvdf85 <- big_rcv_sum %>% filter(s == 85)
cat("About to rearrange.")
ordered_u_probs <- t(apply(rcvdf85, 1, rearrange))
cat("Done rearranging.")
ordered_u_probs <- apply(ordered_u_probs, 2, as.numeric)

rcvdf852 <- as.data.frame(ordered_u_probs)
rcvdf852 <- rcvdf852 %>% mutate(
  ben_rcv2 = AB_CBo * (uB - uC) + BC_BCo * ((uB - uC)/2) + BC_ACo * ((uA - uC)/2), 
  cost_rcv2 = ABo * (uA - uB) + AB_ABo * (uA - uB) + AB_ACo * (uA - uC) + AC_ACo * ((uA - uC) / 2) + AC_ABo * ((uA - uB) / 2) + AC_BCo * ((uB - uC)/2) + BC_BAo * ((uA - uB)/2),
  ben_rcv3 = AB_CBo * ((uB - uC) / 2) + BC_BAo * ((uA - uB) / 2),
  cost_rcv3 = ACo * (uA - uC) + BCo * (uB - uC) + AB_ABo * ((uA - uB) / 2) + AB_ACo * (uA - uC)/2 + AC_ACo * (uA - uC)/2 + AC_ABo * (uA - uB) + BC_BCo * (uC - uB)/2 + BC_ACo * (uA - uC)/2,
  pprob2 = AB_CBo + BC_BCo + BC_ACo, 
  pprob3 = AB_CBo + BC_BAo)

# compute correlations
rcvdf853 <- cbind(rcvdf85, rcvdf852)
rcvdf853 <- rcvdf853 %>% group_by(case, iter) %>% summarise(corr_rcv2 = cor(ben_rcv2, cost_rcv2), corr_rcv3 = cor(ben_rcv3, cost_rcv3),pprob2 = weighted.mean(pprob2, weight = w), pprob3 = weighted.mean(pprob3, weight = w), path = "IRV")

# Tidy DF for plotting
rcvdf854 <- rcvdf853 %>% gather(key = "quantity", value = "value", corr_rcv2:pprob3) 
rcvdf854 <- rcvdf854 %>% mutate(type = ifelse(quantity %in% c("corr_rcv2", "pprob2"), "IRV_second", "IRV_third"), quantity = ifelse(quantity %in% c("corr_rcv2", "corr_rcv3"), "correlation", "pprob"))
rcvdf854 <- spread(rcvdf854, key = quantity, value = value)

# Re-order plurality DF
plurdf85 <- big_plur_sum %>% filter(s == 85)
cat("About to re-arrange.")
ordered_u_probs_plur <- t(apply(plurdf85, 1, rearrange))
ordered_u_probs_plur <- apply(ordered_u_probs_plur, 2, as.numeric)

# Tidy DF for plotting
plurdf852 <- as.data.frame(ordered_u_probs_plur)
plurdf852 <- plurdf852 %>% mutate(ben_p = BCpo * ((uB - uC) / 2), cost_p = ABpo * (uA - uB) + ACpo * ((uB - uC)/2), pprob_plur = BCpo
)

# compute correlations
plurdf853 <- cbind(plurdf85, plurdf852)
plurdf853 <- plurdf853 %>% group_by(case, iter) %>% summarise(correlation = cor(ben_p, cost_p), pprob = weighted.mean(pprob_plur, weight = w), type = "Plurality", path = "Plurality")

# get country weights
# sv_list <- list()
n <- length(big_list_na_omit)
country_weight <- matrix(nrow = n, ncol = 2)
for(i in 1:n){
  progress(i)
  if(i == n) cat("Done! \n")
  this_list <- big_list_na_omit[[i]]
  # df_list <- lapply(s_list, function(x) convert_andy_to_sv_item_two(this_list$U, this_list$weights, x, this_list$v_vec))
  # df <- as.data.frame(do.call(rbind, df_list))
  dfcase <- names(big_list_na_omit)[[i]]
  # df$weight <- big_list_na_omit[[i]]$weights
  dfcountry <- substr(dfcase, 1, 3)
  dfweight_sum <- sum(big_list_na_omit[[i]]$weights)
  dfVAP <- vap$VAP[vap$cntry == dfcountry]
  dfm <- vap$Freq[vap$cntry == dfcountry]
  # df$weight_rep <- df$weight * (df$VAP / (df$weight_sum * df$m))
  #df <- apply(df, 2, as.numeric)
  country_weight[i, 1] <- names(big_list_na_omit)[[i]]
  country_weight[i, 2] <- dfVAP / dfm
  # sv_list[[i]] <- df
}



# Bind DF together and tidy
weightdf <- as.data.frame(country_weight) 
names(weightdf) <- c("case", "ctryweight")
weightdf$ctryweight <- as.numeric(as.character(weightdf$ctryweight))

conjdf <- rbind(plurdf853, rcvdf854) %>% left_join(weightdf)
conjdf_quant <- conjdf %>% group_by(iter, type) %>% summarise(
  corr_q025 = wtd.quantile(correlation, q = 0.025, weight = ctryweight),
  corr_q25 = wtd.quantile(correlation, q = 0.25, weight = ctryweight),
  corr_q50 = wtd.quantile(correlation, q = 0.5, weight = ctryweight),
  corr_q75 = wtd.quantile(correlation, q = 0.75, weight = ctryweight),
  corr_q975 = wtd.quantile(correlation, q = 0.975, weight = ctryweight),
  corr_mean = wtd.mean(correlation, weight = ctryweight),
  pprob_q025 = wtd.quantile(pprob, q = 0.025, weight = ctryweight),
  pprob_q25 = wtd.quantile(pprob, q = 0.25, weight = ctryweight),
  pprob_q50 = wtd.quantile(pprob, q = 0.5, weight = ctryweight),
  pprob_q75 = wtd.quantile(pprob, q = 0.75, weight = ctryweight),
  pprob_q975 = wtd.quantile(pprob, q = 0.975, weight = ctryweight),
  pprob_mean = wtd.mean(pprob, weight = ctryweight))

# save conjdf
save(conjdf, file = here("output/conjdf.Rdata"))

# Pprobs plot (raw)
ggplot(conjdf, aes(x = iter)) +
geom_line(aes(y = pprob, group = interaction(case, type), colour = type), alpha = 0.1) +
geom_line(data = conjdf_quant, aes(x = iter, y = pprob_mean, group = type, colour = type), lwd = 2) + 
geom_hline(yintercept = 0, lty = "dashed") +
labs(x = "Iteration", y = "Probability vote is beneficial * electorate size") +
theme_sv()  +
ylim(c(0, 0.5))
ggsave(here("output/figures/conj1.pdf"), device = cairo_pdf,
       width = 4,
       height = 4)

# Pprobs plot (quantiles)
ggplot(conjdf_quant, aes(x = iter)) +
geom_line(aes(y = pprob_q50, group = type, colour = type), alpha = 1) +
geom_ribbon(aes(ymin = pprob_q25, ymax = pprob_q75, group = type, fill = type), alpha = 0.25) + 
geom_hline(yintercept = 0, lty = "dashed") +
theme_sv() 

# Test

ggplot(plurdf853, aes(x = iter)) +
geom_line(aes(y = correlation, group = interaction(case, type), colour = type), alpha = 0.05) +
geom_hline(yintercept = 0, lty = "dashed") +
theme_sv()

# Correlations plot
ggplot(conjdf, aes(x = iter)) +
geom_line(aes(y = correlation, group = interaction(case, type), colour = type), alpha = 0.05) +
geom_line(data = conjdf_quant, aes(x = iter, y = corr_mean, group = type, colour = type), lwd = 2) + 
geom_hline(yintercept = 0, lty = "dashed") +
theme_sv()
ggsave(here("output/figures/conj2.pdf"), device = cairo_pdf,
       height = 4,
       width = 4)

# Correlations plot (quantiles)
ggplot(conjdf_quant, aes(x = iter)) +
geom_line(aes(y = corr_q50, group = type, colour = type), alpha = 1) +
geom_ribbon(aes(ymin = corr_q25, ymax = corr_q75, group = type, fill = type), alpha = 0.25) + 
geom_hline(yintercept = 0, lty = "dashed") +
theme_sv() 
