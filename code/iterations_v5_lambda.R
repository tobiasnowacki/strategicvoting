## Full script to generate results for paper
## latest version: 11/Sep/2019
## Toby Nowacki

### -------------------------
### TO-DO -------
### 	* adjust number of iterations
### 	* parameter values

### -------------------------
### LOAD DEPENDENCIES -------

# install.packages("here")
library(here) 						# to get dir
source(here("code/full_header.R")) 	# fn's and data
source(here("code/prep_cses.R")) 	# data prep

### ------------------
### SET VALUES -------

s_list <- as.list(c(10, 25, 40, 55, 70, 85)) # precision (s)
lambda_list <- as.list(c(0.05, 0.1, 0.01))	 # responsiveness ()
epsilon_thresh <- 0.001						 # threshold (epsilon)
max_iter_val <- 250							 # no of iterations
which_cases <- 1:160						 # which cases?

### ---
### RUN ITERATIONS

# For loop over values of s
save_out <- list()

# for(lambda_val in 1:3){
for(lambda_val in 3){
	lambda <- lambda_list[[lambda_val]]

	cat(paste0("=== Beginning loop for lambda == ", lambda, ". === \n"))

	save_out[[lambda_val]] <- list()

	# for(s_loop_val in c(1)){ 
	for(s_loop_val in c(1, 4, 6)){ 

	  s_val <- s_list[[s_loop_val]]

	  cat(paste0("=== Beginning loop for s == ", s_val, ". === \n"))

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
	                    epsilon_thresh, 
	                    max_iter_val)
	              out
	  }

	  stopCluster(cl)

	  cat("Done. \n")

	  path <- paste0("output/figs_v2/", lambda_val, "/", s_val)
	  path_files <- paste0("output/files/", lambda_val, "/", s_val)

	  # name list obj
	  names(cases_converge) <- names(big_list_na_omit)[which_cases]

	  # save v_vec paths somewhere...
	  out <- lapply(cases_converge, function(x) x[c(2, 4)])
	  save(out, file = here(paste0(path_files, "v_vecs_", lambda_val, "_", s_val, ".Rdata")))
	  cat("Saved v_vecs!")

	  # produce iteration plot(s)
	  plot_v_vec_distance(cases_converge, path, 
	                                  n_lag = 20, 
	                                  avg_span = 10)

	  cat("v_vec distance plotted: complete. \n")


	  ### PLAYGROUND
	  # cl <- makeCluster(6)
	  # registerDoParallel(cl)
	  # # parallelise
	  # test <- foreach(case = which_cases, 
	  #             .packages = c("gtools", "stringr", "tidyverse",
	  #                           "questionr")
	  #             ) %dopar% {
	  #             out <- iteration_convergence_slim(
	  #                   big_list_na_omit[[case]], 
	  #                   big_list_na_omit[[case]]$v_vec, 
	  #                   lambda, 
	  #                   s_val, 
	  #                   max_iter_val)
	  #             out <- do.call(rbind, out)
	  #             # sink("log_convergence.txt", append = TRUE)
	  #             # cat(paste(case, " done.", "\n"))
	  #             out
	  # }

	  # stopCluster(cl)

	  # test_item <- iteration_convergence_slim(big_list_na_omit[[1]], big_list_na_omit[[1]]$v_vec, 0.05, 85, 250)
	  # do.call(rbind, test_item)
	  ###

	  # produce v_vec path plot(s)
	  # all together
	  # todo: correct plots' size
	  # todo: given I save them previously I could do this later, no?
	  joint_v_vec_plot(cases_converge, path)

	  cat("v_vec paths plotted. \n")

	  # # rcv separate ones by non_convergence
	  # non_conv_v_vec_plot(cases_converge, path, max_iter_val)
	  # vote_tally <- non_conv_strat_votes(cases_converge, path, max_iter_val)

	  # # rcv plot non_convergent
	  # unique_nc_cases <- unique(vote_tally$case)

	  # # re-label as votes
	  # vote_tally$sin_rcv <- recode(vote_tally$sin_rcv, "1" = "ABC", 
	  # 	       			"2" = "ACB", 
	  # 	       			"3" = "BAC", 
	  # 	       			"4" = "BCA", 
	  # 	       			"5" = "CAB", 
	  # 	       			"6" = "CBA")

	  # vote_tally$opt_rcv <-	recode(vote_tally$opt_rcv, "1" = "ABC", 
	  # 	       			"2" = "ACB", 
	  # 	       			"3" = "BAC", 
	  # 	       			"4" = "BCA", 
	  # 	       			"5" = "CAB", 
	  # 	       			"6" = "CBA")

	  # # plot non-convergent cases in detail
	  # for(i in unique_nc_cases){
	  #   ggplot(vote_tally %>% filter(case == i), 
	  #          aes(iter, freq)) +
	  #     geom_line(aes(group = opt_rcv, colour = opt_rcv %>% as.factor)) +
	  #     facet_grid(~ sin_rcv) +
	  #     theme_sv() +
	  #     labs(x = "Iteration",
	  #          y = "Frequency", 
	  #          colour = "Optimal Vote",
	  #          title = i) +
	  #     theme(legend.position = "bottom")
	  #     ggsave(here(paste0(path, "/nc_", i, ".pdf")),
	  #            device = cairo_pdf,
	  #            width = 5,
	  #            height = 5)
	  # }

	  # cat("non-convergent cases plotted. \n")

	  # # group expected benefit and the like
	  # summary_stats <- get_sum_stats(cases_converge)
	  # # todo here: (a) correct weights (just weighted means)
	  # #        (b) produce averages across cases.
	  # summary_stats_wide <- summary_stats %>% 
	  #   gather(., key = "Statistic", 
	  #          value = "Value", "Prevalence":"ExpBenefit")

	  # summary_agg <- summary_stats_wide %>% group_by(iter, Statistic, System) %>%
	  #   summarise(Value = mean(Value))

	  # # Summary statistics
	  # ggplot(summary_stats_wide, aes(iter, Value)) +
	  #   geom_line(aes(group = interaction(System, case, Statistic),
	  #                 colour = System), alpha = 0.3) +
	  #   geom_line(data = summary_agg %>% filter(System == "Plurality"),
	  #             aes(group = interaction(System, Statistic)), alpha = 1,
	  #             color = "#CC6600", lwd = 1.1) +
	  #   geom_line(data = summary_agg %>% filter(System == "IRV"),
	  #             aes(group = interaction(System, Statistic)), alpha = 1,
	  #             color = "#004C99", lwd = 1.1) +
	  #   scale_color_manual(values = cbbPalette[c(3, 2)]) +
	  #   facet_wrap(. ~ Statistic, scales = "free_y") +
	  #   theme_sv() +
	  #   guides(colour = guide_legend(override.aes = list(alpha = 1))) +
	  #                 theme(legend.position = "bottom", legend.direction = "horizontal") +
	  #   labs(x = "Degree of Strategicness (Iterations)")
	  # ggsave(here(paste0(path, "/main_results.pdf")), 
	  #        device = cairo_pdf)

	  # cat("Summary statistics plotted. \n")

	}
}



# ## MESSY STUFF: CONJECTURE TESTS



# # Conjecture tests
# big_rcv_sum <- list()
# big_plur_sum <- list()
# for(i in 1:160){
# 	print(i)
# 	mini_list1 <- cases_converge[[i]][[1]]
# 	mini_list1$case <- names(big_list_na_omit)[i]
# 	mini_list1$system <- "IRV"
# 	big_rcv_sum[[i]] <- mini_list1
# 	mini_list2 <- cases_converge[[i]][[3]]
# 	mini_list2$case <- names(big_list_na_omit)[i]
# 	mini_list2$system <- "Plurality"
# 	big_plur_sum[[i]] <- mini_list2
# }

# big_rcv_sum <- do.call(rbind, big_rcv_sum)
# big_plur_sum <- do.call(rbind, big_plur_sum)


# ## CONJECTURE 1a
# # Summed pivotal probabilities
# rcv_summary <- big_rcv_sum %>% group_by(case, s, iter) %>% summarise_at(vars(AB:BCp), first) 
# rcv_summary$rcv_sum <- rowSums(rcv_summary[, 4:16])
# rcv_summary$plur_sum <- rowSums(rcv_summary[, 7:19])
# rcv_summary <- rcv_summary %>% mutate(diff = plur_sum - rcv_sum)

# # Do the same thing for plurality
# plur_summary <- big_plur_sum %>% group_by(case, s, iter) %>% summarise_at(vars(AB:BCp), first) 
# plur_summary$rcv_sum <- rowSums(plur_summary[, 4:16])
# plur_summary$plur_sum <- rowSums(plur_summary[, 7:19])
# plur_summary <- plur_summary %>% mutate(diff = plur_sum - rcv_sum)

# cmat_plur <- matrix(c(1, 2, 3, 2, 1, 3, 1, 3, 2, 3, 1, 2, 2, 3, 1, 3, 2, 1), byrow = T, ncol = 3)

# cmat_rcv <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 2, 1, 3, 7, 8, 9, 4, 5, 6, 10, 12, 11, 1, 3, 2, 4, 6, 5, 10, 11, 12, 7, 8, 9, 3, 1, 2, 10, 11, 12, 4, 6, 5, 7, 9, 8, 2, 3, 1, 7, 9, 8, 10, 12, 11, 4, 5, 6, 3, 2, 1, 10, 12, 11, 7, 9, 8, 4, 6, 5), byrow = T, ncol = 12)

# cmat_u <- matrix(c(1, 2, 3, 1, 3, 2, 2, 1, 3, 2, 3, 1, 3, 1, 2, 3, 2, 1), byrow = T, ncol = 3)

# # re-order RCV DF
# rcvdf85 <- big_rcv_sum %>% filter(s == 85)
# ordered_u_probs <- t(apply(rcvdf85, 1, rearrange))
# ordered_u_probs <- apply(ordered_u_probs, 2, as.numeric)

# rcvdf852 <- as.data.frame(ordered_u_probs)
# rcvdf852 <- rcvdf852 %>% mutate(
#   ben_rcv2 = AB_CBo * (uB - uC) + BC_BCo * ((uB - uC)/2) + BC_ACo * ((uA - uC)/2), 
#   cost_rcv2 = ABo * (uA - uB) + AB_ABo * (uA - uB) + AB_ACo * (uA - uC) + AC_ACo * ((uA - uC) / 2) + AC_ABo * ((uA - uB) / 2) + AC_BCo * ((uB - uC)/2) + BC_BAo * ((uA - uB)/2),
#   ben_rcv3 = AB_CBo * ((uB - uC) / 2) + BC_BAo * ((uA - uB) / 2),
#   cost_rcv3 = ACo * (uA - uC) + BCo * (uB - uC) + AB_ABo * ((uA - uB) / 2) + AB_ACo * (uA - uC)/2 + AC_ACo * (uA - uC)/2 + AC_ABo * (uA - uB) + BC_BCo * (uC - uB)/2 + BC_ACo * (uA - uC)/2,
#   pprob2 = AB_CBo + BC_BCo + BC_ACo, 
#   pprob3 = AB_CBo + BC_BAo)

# # compute correlations
# rcvdf853 <- cbind(rcvdf85, rcvdf852)
# rcvdf853 <- rcvdf853 %>% group_by(case, iter) %>% summarise(corr_rcv2 = cor(ben_rcv2, cost_rcv2), corr_rcv3 = cor(ben_rcv3, cost_rcv3),pprob2 = weighted.mean(pprob2, weight = w), pprob3 = weighted.mean(pprob3, weight = w), path = "IRV")

# # Tidy DF for plotting
# rcvdf854 <- rcvdf853 %>% gather(key = "quantity", value = "value", corr_rcv2:pprob3) 
# rcvdf854 <- rcvdf854 %>% mutate(type = ifelse(quantity %in% c("corr_rcv2", "pprob2"), "IRV_second", "IRV_third"), quantity = ifelse(quantity %in% c("corr_rcv2", "corr_rcv3"), "correlation", "pprob"))
# rcvdf854 <- spread(rcvdf854, key = quantity, value = value)

# # Re-order plurality DF
# plurdf85 <- big_plur_sum %>% filter(s == 85)
# ordered_u_probs_plur <- t(apply(plurdf85, 1, rearrange))
# ordered_u_probs_plur <- apply(ordered_u_probs_plur, 2, as.numeric)

# # Tidy DF for plotting
# plurdf852 <- as.data.frame(ordered_u_probs_plur)
# plurdf852 <- plurdf852 %>% mutate(ben_p = BCpo * ((uB - uC) / 2), cost_p = ABpo * (uA - uB) + ACpo * ((uB - uC)/2), pprob_plur = BCpo
# )

# # compute correlations
# plurdf853 <- cbind(plurdf85, plurdf852)
# plurdf853 <- plurdf853 %>% group_by(case, iter) %>% summarise(correlation = cor(ben_p, cost_p), pprob = weighted.mean(pprob_plur, weight = w), type = "Plurality", path = "Plurality")

# # get country weights
# # sv_list <- list()
# n <- length(big_list_na_omit)
# country_weight <- matrix(nrow = n, ncol = 2)
# for(i in 1:n){
#   progress(i)
#   if(i == n) cat("Done! \n")
#   this_list <- big_list_na_omit[[i]]
#   # df_list <- lapply(s_list, function(x) convert_andy_to_sv_item_two(this_list$U, this_list$weights, x, this_list$v_vec))
#   # df <- as.data.frame(do.call(rbind, df_list))
#   dfcase <- names(big_list_na_omit)[[i]]
#   # df$weight <- big_list_na_omit[[i]]$weights
#   dfcountry <- substr(dfcase, 1, 3)
#   dfweight_sum <- sum(big_list_na_omit[[i]]$weights)
#   dfVAP <- vap$VAP[vap$cntry == dfcountry]
#   dfm <- vap$Freq[vap$cntry == dfcountry]
#   # df$weight_rep <- df$weight * (df$VAP / (df$weight_sum * df$m))
#   #df <- apply(df, 2, as.numeric)
#   country_weight[i, 1] <- names(big_list_na_omit)[[i]]
#   country_weight[i, 2] <- dfVAP / dfm
#   # sv_list[[i]] <- df
# }



# # Bind DF together and tidy
# weightdf <- as.data.frame(country_weight) 
# names(weightdf) <- c("case", "ctryweight")
# weightdf$ctryweight <- as.numeric(as.character(weightdf$ctryweight))

# conjdf <- rbind(plurdf853, rcvdf854) %>% left_join(weightdf)
# conjdf_quant <- conjdf %>% group_by(iter, type) %>% summarise(
#   corr_q025 = wtd.quantile(correlation, q = 0.025, weight = ctryweight),
#   corr_q25 = wtd.quantile(correlation, q = 0.25, weight = ctryweight),
#   corr_q50 = wtd.quantile(correlation, q = 0.5, weight = ctryweight),
#   corr_q75 = wtd.quantile(correlation, q = 0.75, weight = ctryweight),
#   corr_q975 = wtd.quantile(correlation, q = 0.975, weight = ctryweight),
#   corr_mean = wtd.mean(correlation, weight = ctryweight),
#   pprob_q025 = wtd.quantile(pprob, q = 0.025, weight = ctryweight),
#   pprob_q25 = wtd.quantile(pprob, q = 0.25, weight = ctryweight),
#   pprob_q50 = wtd.quantile(pprob, q = 0.5, weight = ctryweight),
#   pprob_q75 = wtd.quantile(pprob, q = 0.75, weight = ctryweight),
#   pprob_q975 = wtd.quantile(pprob, q = 0.975, weight = ctryweight),
#   pprob_mean = wtd.mean(pprob, weight = ctryweight))

# # save conjdf
# save(conjdf, file = here(paste0(path_files, "/conjdf.Rdata")))

# # TODO: fix these.

# # Pprobs plot (raw)
# ggplot(conjdf, aes(x = iter)) +
# geom_line(aes(y = pprob, group = interaction(case, type), colour = type), alpha = 0.1) +
# geom_line(data = conjdf_quant, aes(x = iter, y = pprob_mean, group = type, colour = type), lwd = 2) + 
# geom_hline(yintercept = 0, lty = "dashed") +
# labs(x = "Iteration", y = "Probability vote is beneficial * electorate size") +
# theme_sv()  +
# ylim(c(0, 0.5))
# ggsave(here("output/figures/conj1.pdf"), device = cairo_pdf,
#        width = 4,
#        height = 4)

# # Pprobs plot (quantiles)
# ggplot(conjdf_quant, aes(x = iter)) +
# geom_line(aes(y = pprob_q50, group = type, colour = type), alpha = 1) +
# geom_ribbon(aes(ymin = pprob_q25, ymax = pprob_q75, group = type, fill = type), alpha = 0.25) + 
# geom_hline(yintercept = 0, lty = "dashed") +
# theme_sv() 

# # Test

# ggplot(plurdf853, aes(x = iter)) +
# geom_line(aes(y = correlation, group = interaction(case, type), colour = type), alpha = 0.05) +
# geom_hline(yintercept = 0, lty = "dashed") +
# theme_sv()

# # Correlations plot
# ggplot(conjdf, aes(x = iter)) +
# geom_line(aes(y = correlation, group = interaction(case, type), colour = type), alpha = 0.05) +
# geom_line(data = conjdf_quant, aes(x = iter, y = corr_mean, group = type, colour = type), lwd = 2) + 
# geom_hline(yintercept = 0, lty = "dashed") +
# theme_sv() +
# labs(x = "Iteration", y = "Correlation") +
# theme(legend.position = "bottom",
#       legend.direction = "horizontal")
# ggsave(here("output/figures/conj2a.pdf"), device = cairo_pdf,
#        width = 4,
#        height= 4)

# # Correlations plot (quantiles)
# ggplot(conjdf_quant, aes(x = iter)) +
# geom_line(aes(y = corr_q50, group = type, colour = type), alpha = 1) +
# geom_ribbon(aes(ymin = corr_q25, ymax = corr_q75, group = type, fill = type), alpha = 0.25) + 
# geom_hline(yintercept = 0, lty = "dashed") +
# theme_sv() 


# TODO: conjecture tests (for lambda and s)


# TODO still: 'diagnostic plots' (case-by-case)
# Another question is: what about random starting points?
# Should I include this into the script?

### ROBUSTNESS CHECKS
### ---------------
# I should check if I get the same with Andy's function.
# tie_test_list <- list()
# for(i in 1:160){
#   print(i)
#   tie_test_list[[i]] <- iterated_best_response_sequence(U = big_list_na_omit[[i]]$U %>% as.matrix, 
#                                 s = s_val, 
#                                 weights = big_list_na_omit[[i]]$weights, 
#                                 rule = "AV", 
#                                 lambda = lambda,
#                                 until.convergence = F, 
#                                 max.iterations = 250, 
#                                 sincere.proportion = 0, 
#                                 candidates = c("a", "b", "c"), 
#                                 ballots = c("abc", "acb", "bac", "bca", "cab", "cba"), 
#                                 the.floor = 0, 
#                                 noisy = F)
# }

# Create distance vector
# dist_list <- list()
# for(i in 1:160){
#   vec_dist <- do.call(rbind, lapply(tie_test_list[[i]], function(x) x$distance.from.last))
#   dist_df <- data.frame(dist = vec_dist, iter = 1:length(vec_dist), case = names(big_list_na_omit)[i])
#   dist_list[[i]] <- dist_df
# }

# big_dist_df <- do.call(rbind, dist_list)

# # Plot vvec distances
# ggplot(big_dist_df, aes(iter, log(dist))) +
#   geom_line(aes(group = case), alpha = 0.1) +
#   geom_line(data = big_dist_df %>% filter(case == "SWE_2014"), aes(group = case), colour = "red") +
#   theme_sv()


# # Get expected benefit
# test <- expected.benefit.mat.from.ibrs(tie_test_list[[2]])
# test %>% colMeans()

# summary_stats_wide %>% 
# 	filter(case == "AUT_2013" & Statistic == "ExpBenefit") # Comparison

# # Get prevalence
# ebm.av.int <- test
# ebm.av.int[test > 0] <- 1
# ebm.av.int %>% colMeans()

# summary_stats_wide %>% filter(case == "AUT_2013" & Statistic == "Prevalence")

# Get magnitude.


## I'm getting the same results!
## But the distance plots are different.
## That's probably because they weight different things.
## Weighting.

### -------------------------
