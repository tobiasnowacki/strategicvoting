### -------------------------
### LOAD DEPENDENCIES -------

# install.packages("here")
library(here) 						# to get dir
source(here("code/full_header.R")) 	# fn's and data
source(here("code/prep_cses.R")) 	# data prep
library(devtools)
source_url("https://raw.githubusercontent.com/tobiasnowacki/RTemplates/master/plottheme.R")

### ------------------
### SET VALUES -------

cmd_line_args <- commandArgs(trailingOnly = TRUE)
cat(cmd_line_args, sep = "n")

start_time <- proc.time()

s_list <- as.list(c(10, 25, 40, 55, 70, 85)) # precision (s)
lambda_list <- as.list(c(0.05, 0.1, 0.01))	 # responsiveness ()
epsilon_thresh <- 0.001						 # threshold (epsilon)
max_iter_val <- 60				 # no of iterations
which_cases <- 1:2				 # which cases?

s_choice = 6
l_choice = 1

# ifelse(length(cmd_line_args >= 1),
#        s_choice <- as.numeric(cmd_line_args[1]),
#        s_choice <- 6)
# ifelse(length(cmd_line_args >= 2),
#        l_choice <- as.numeric(cmd_line_args[2]),
#        l_choice <- 1)

### ---
### RUN ITERATIONS

# test_vec = big_list_na_omit[[1]]$v_vec
# test_3 = c(test_vec[1] + test_vec[2], test_vec[3] + test_vec[4], test_vec[5] + test_vec[6])

# test1 = sv(big_list_na_omit[[1]]$U, v.vec = test_vec, s = 85, rule = "AV")
# test1

# test2 = sv(big_list_na_omit[[1]]$U, v.vec = test_3, s = 85, ae_pack = FALSE)
# str(test1)
# str(test2)

# colSums(test1$V.mat)
# colSums(test2$V.mat)

# plot(density(test2$tau))

# test_old = one_iteration(big_list_na_omit[[1]], 
#                      big_list_na_omit[[1]]$v_vec, 0.05, 85)

# head(test[[1]])

# plot(density(test[[1]]$tau_plur))
# plot(density(test[[1]]$tau_rcv))
# plot(density(test_old[[1]]$tau_rcv))
# plot(density(test_old[[1]]$tau_plur))

# plot(density(test[[1]]$tau_rcv / test[[1]]$tau_plur))
# plot(density(test_old[[1]]$tau_rcv / test_old[[1]]$tau_plur))



# test[[4]]
# test_old[[4]]
# head(test)

# names(test$plur_best_response)

# for(lambda_val in 1:3){
for(lambda_val in l_choice){
	lambda <- lambda_list[[lambda_val]]

	method_setting = TRUE

	cat(paste0("=== Beginning loop for lambda == ", lambda, ". === \n"))

	for(s_loop_val in s_choice){ 

	  s_val <- s_list[[s_loop_val]]

	  cat(paste0("=== Beginning loop for s == ", s_val, ". === \n"))

	  cat("Running parallel iterations. \n")

	  clno <- detectCores()
	  cl <- makeCluster(clno)
	  registerDoParallel(cl)
	  # parallelise
	  # cases_converge_1 <- foreach(case = which_cases, 
	  #             .packages = c("gtools", "stringr", "tidyverse", "pivotprobs",
	  #                           "questionr")
	  #             ) %dopar% {
	  #             out <- many_iterations_until_convergence(
	  #                   big_list_na_omit[[case]], 
	  #                   big_list_na_omit[[case]]$v_vec, 
	  #                   lambda, 
	  #                   s_val, 
	  #                   epsilon_thresh, 
	  #                   max_iter_val, ae = TRUE)
	  #             out
	  # }


	  # cat("Second code base. \n")

	  # parallelise
	  cases_converge <- foreach(case = which_cases, 
	              .packages = c("gtools", "stringr", "tidyverse", "pivotprobs",
	                            "questionr")
	              ) %dopar% {
	              out <- many_iterations_until_convergence(
	                    big_list_na_omit[[case]], 
	                    big_list_na_omit[[case]]$v_vec, 
	                    lambda, 
	                    s_val, 
	                    epsilon_thresh, 
	                    max_iter_val, ae = method_setting)
	              out
	  }

	  stopCluster(cl)

	  cat("Done. \n")

	  add_str = "ae"
	  if(method_setting == FALSE){
	  	add_str = "old"
	  }

	  path <- paste0("output/figures/", lambda_val, "/", s_val)
	  path_files <- paste0("output/files/", lambda_val, "/", s_val)

	  # # name list obj
	  names(cases_converge) <- names(big_list_na_omit)[which_cases]
	  # names(cases_converge_1) <- names(big_list_na_omit)[which_cases]

	  # Investigate
	  # head(cases_converge[[1]])[[1]][1:10, ]
	  # table(cases_converge[[1]][[1]]$sin_rcv, cases_converge[[1]][[1]]$opt_rcv)

	  # # save v_vec paths somewhere...
	  out <- lapply(cases_converge, function(x) x[c(2, 4)])
	  save(out, file = here(paste0(path_files, "v_vecs_", lambda_val, "_", s_val, add_str, ".Rdata")))
	  # cat("Saved v_vecs!")

	  # # produce iteration plot(s)
	  # plot_v_vec_distance(cases_converge, path, 
	  #                                 n_lag = 20, 
	  #                                 avg_span = 10)

	  # cat("v_vec distance plotted: complete. \n")

	  # # produce v_vec path plot(s)
	  v_mut = joint_v_vec_plot(cases_converge, path)

	  # cat("v_vec paths plotted. \n")

	  # # group expected benefit and the like
	  summary_stats <- get_sum_stats(cases_converge)

	  # Identify cases where IRV prevalence is high
	  # suspect_cases = summary_stats %>% filter(iter == 50, 
	  #                          System == "IRV",
	  #                          Prevalence > 0.25) %>%
	  # 	distinct(case) %>% unlist

	  # Append v_mut to signal these cases
	  v_mut = v_mut %>%
	  	mutate(irv_strange = case == "SWE_2014")

	  swe_check = cases_converge$SWE_2014$rcv_df
	  swe_check %>% filter(iter == 12) %>% head
	  iter10 = swe_check %>% filter(iter == 12)
	  table(iter10$sin_rcv, iter10$opt_rcv)


	  # Plot cases
	  ggplot(v_mut %>% filter(system == "IRV", irv_strange == TRUE), aes(x = C, y = B)) +
	      geom_line(aes(
	  				colour = case),
	                alpha = 0.5) +
	      # geom_point(data = v_mut %>% filter(state %in% c("first", "last")),
	      #            aes(colour = state),
	      #            size = 0.5,
	      #            alpha = 0.6) +
	      coord_fixed() +
	      geom_ternary_boundary() +
	      theme_tn()
	   # ggsave(paste0("output/figures/irv_trouble_v_vec", add_str, ".pdf"), device = cairo_pdf)

	  # Get minimum vote share for all cases and iters (in IRV)
	  min_share = out %>% map(~ apply(.x$rcv_v_vec[, 5:6], 1, sum))
	  # min_share = out %>% map(~ .x$rcv_v_vec[, 5])

	  min_share_df = min_share %>% map_dfr(~ .x %>% as.data.frame %>% 
	                        mutate(iter = 1:61), .id = "case") %>%
	  	mutate(irv_strange = case == "SWE_2014")
	  names(min_share_df)[2] = "value"

	  # Identify strange case
	  min_share_df %>% filter(irv_strange == TRUE, iter == 61) %>%
	  	arrange(value)

	  ggplot(min_share_df, aes(iter, value)) +
	  	geom_line(aes(group = case, colour = irv_strange)) +
	  	theme_tn()
	  # ggsave(paste0("output/figures/irv_trouble_v_vec_min", add_str, ".pdf"), device = cairo_pdf)
	  # summary_stats_old = summary_stats
	  # todo here: (a) correct weights (just weighted means)
	  #        (b) produce averages across cases.
  	  # save(summary_stats, file = here(paste0(path_files, "summary_", lambda_val, "_", s_val, ".Rdata")))
	  summary_stats_wide <- summary_stats %>% 
	    gather(., key = "Statistic", 
	           value = "Value", "Prevalence":"ExpBenefit")

	  summary_stats %>%
	  	filter(iter > 20 & System == "IRV" & Prevalence > 0.3)

	  summary_stats %>%
	  	filter(iter %in% 10:20 & System == "IRV" & case == "SWE_2014")

	  # weight by case weight
	  summary_agg <- summary_stats_wide %>% group_by(iter, Statistic, System) %>%
	    summarise(Value = wtd.mean(Value, case_weight_tbl$case_weight))

	  # Summary statistics
	  ggplot(summary_stats_wide, aes(iter, Value)) +
	    geom_line(aes(group = interaction(System, case, Statistic),
	                  colour = System), alpha = 0.3) +
	    geom_line(data = summary_agg %>% filter(System == "Plurality"),
	              aes(group = interaction(System, Statistic)), alpha = 1,
	              color = "#CC6600", lwd = 1.1) +
	    geom_line(data = summary_agg %>% filter(System == "IRV"),
	              aes(group = interaction(System, Statistic)), alpha = 1,
	              color = "#004C99", lwd = 1.1) +
	    scale_color_manual(values = cbbPalette[c(3, 2)]) +
	    facet_wrap(. ~ Statistic, scales = "free_y") +
	    theme_sv() +
	    guides(colour = guide_legend(override.aes = list(alpha = 1))) +
	                  theme(legend.position = "bottom", legend.direction = "horizontal") +
	    labs(x = "Degree of Strategicness (Iterations)")
	  # ggsave(here(paste0(path, "/main_results.pdf")), 
	         # device = cairo_pdf)

	  # cat("Summary statistics plotted. \n")

	  # if(s_choice == 6 & l_choice == 1){

	  # 	# rcv separate ones by non_convergence
	  # 	non_conv_v_vec_plot(cases_converge, path, max_iter_val)
	  # 	vote_tally <- non_conv_strat_votes(cases_converge, path, max_iter_val)
	  # 	save(vote_tally, file = here(paste0(path_files, "votetally_", lambda_val, "_", s_val, ".Rdata")))


	  # 	# rcv plot non_convergent
	  # 	unique_nc_cases <- unique(vote_tally$case)

	  # 	if(!is.null(vote_tally)){

		 #  	# re-label as votes
		 #  	vote_tally$sin_rcv <- recode(vote_tally$sin_rcv, "1" = "ABC", 
		 #  		       			"2" = "ACB", 
		 #  		       			"3" = "BAC", 
		 #  		       			"4" = "BCA", 
		 #  		       			"5" = "CAB", 
		 #  		       			"6" = "CBA")

		 #  	vote_tally$opt_rcv <-	recode(vote_tally$opt_rcv, "1" = "ABC", 
		 #  		       			"2" = "ACB", 
		 #  		       			"3" = "BAC", 
		 #  		       			"4" = "BCA", 
		 #  		       			"5" = "CAB", 
		 #  		       			"6" = "CBA")

		 #  	# plot non-convergent cases in detail
		 #  	for(i in unique_nc_cases){
		 #  	  ggplot(vote_tally %>% filter(case == i), 
		 #  	         aes(iter, freq)) +
		 #  	    geom_line(aes(group = opt_rcv, colour = opt_rcv %>% as.factor)) +
		 #  	    facet_grid(~ sin_rcv) +
		 #  	    theme_sv() +
		 #  	    labs(x = "Iteration",
		 #  	         y = "Frequency", 
		 #  	         colour = "Optimal Vote",
		 #  	         title = i) +
		 #  	    theme(legend.position = "bottom")
		 #  	    ggsave(here(paste0(path, "/nc_", i, ".pdf")),
		 #  	           device = cairo_pdf,
		 #  	           width = 5,
		 #  	           height = 5)
		 #  	}

	  # 	}

	  # 	cat("non-convergent cases plotted. \n")

	  # 	# Conjecture tests

	  # 	cat("starting conjecture tests. \n")

	  # 	# Prepare data
	  # 	big_rcv_sum <- list()
	  # 	big_plur_sum <- list()
	  # 	for(i in which_cases){
	  # 		mini_list1 <- cases_converge[[i]][[1]] %>% filter(iter < 61)
	  # 		mini_list1$case <- names(big_list_na_omit)[i]
	  # 		mini_list1$system <- "IRV"
	  # 		big_rcv_sum[[i]] <- mini_list1
	  # 		mini_list2 <- cases_converge[[i]][[3]] %>% filter(iter < 61)
	  # 		mini_list2$case <- names(big_list_na_omit)[i]
	  # 		mini_list2$system <- "Plurality"
	  # 		big_plur_sum[[i]] <- mini_list2
	  # 	}

	  # 	big_rcv_sum <- do.call(rbind, big_rcv_sum)
	  # 	big_plur_sum <- do.call(rbind, big_plur_sum)

	  # 	## CONJECTURE 1a
	  # 	# Summed pivotal probabilities
	  # 	rcv_summary <- big_rcv_sum %>% group_by(case, s, iter) %>% summarise_at(vars(AB:BCp), first) 
	  # 	rcv_summary$rcv_sum <- rowSums(rcv_summary[, 4:16])
	  # 	rcv_summary$plur_sum <- rowSums(rcv_summary[, 7:19])
	  # 	rcv_summary <- rcv_summary %>% mutate(diff = plur_sum - rcv_sum)

	  # 	# Do the same thing for plurality
	  # 	plur_summary <- big_plur_sum %>% group_by(case, s, iter) %>% summarise_at(vars(AB:BCp), first) 
	  # 	plur_summary$rcv_sum <- rowSums(plur_summary[, 4:16])
	  # 	plur_summary$plur_sum <- rowSums(plur_summary[, 7:19])
	  # 	plur_summary <- plur_summary %>% mutate(diff = plur_sum - rcv_sum)

	  # 	cmat_plur <- matrix(c(1, 2, 3, 2, 1, 3, 1, 3, 2, 3, 1, 2, 2, 3, 1, 3, 2, 1), byrow = T, ncol = 3)

	  # 	cmat_rcv <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 2, 1, 3, 7, 8, 9, 4, 5, 6, 10, 12, 11, 1, 3, 2, 4, 6, 5, 10, 11, 12, 7, 8, 9, 3, 1, 2, 10, 11, 12, 4, 6, 5, 7, 9, 8, 2, 3, 1, 7, 9, 8, 10, 12, 11, 4, 5, 6, 3, 2, 1, 10, 12, 11, 7, 9, 8, 4, 6, 5), byrow = T, ncol = 12)

	  # 	cmat_u <- matrix(c(1, 2, 3, 1, 3, 2, 2, 1, 3, 2, 3, 1, 3, 1, 2, 3, 2, 1), byrow = T, ncol = 3)

	  # 	# re-order RCV DF
	  # 	rcvdf85 <- big_rcv_sum %>% filter(s == 85)
	  # 	ordered_u_probs <- t(apply(rcvdf85, 1, rearrange))
	  # 	ordered_u_probs <- apply(ordered_u_probs, 2, as.numeric)

	  # 	rcvdf852 <- as.data.frame(ordered_u_probs)
	  # 	rcvdf852 <- rcvdf852 %>% mutate(
	  # 	  ben_rcv2 = AB_CBo * (uB - uC) + BC_BCo * ((uB - uC)/2) + BC_ACo * ((uA - uC)/2), 
	  # 	  cost_rcv2 = ABo * (uA - uB) + AB_ABo * (uA - uB) + AB_ACo * (uA - uC) + AC_ACo * ((uA - uC) / 2) + AC_ABo * ((uA - uB) / 2) + AC_BCo * ((uB - uC)/2) + BC_BAo * ((uA - uB)/2),
	  # 	  ben_rcv3 = AB_CBo * ((uB - uC) / 2) + BC_BAo * ((uA - uB) / 2),
	  # 	  cost_rcv3 = ACo * (uA - uC) + BCo * (uB - uC) + AB_ABo * ((uA - uB) / 2) + AB_ACo * (uA - uC)/2 + AC_ACo * (uA - uC)/2 + AC_ABo * (uA - uB) + BC_BCo * (uC - uB)/2 + BC_ACo * (uA - uC)/2,
	  # 	  pprob2 = AB_CBo + BC_BCo + BC_ACo, 
	  # 	  pprob3 = AB_CBo + BC_BAo)

	  # 	# compute correlations
	  # 	rcvdf853 <- cbind(rcvdf85, rcvdf852)
	  # 	rcvdf853 <- rcvdf853 %>% group_by(case, iter) %>% summarise(corr_rcv2 = cor(ben_rcv2, cost_rcv2), corr_rcv3 = cor(ben_rcv3, cost_rcv3),pprob2 = weighted.mean(pprob2, weight = w), pprob3 = weighted.mean(pprob3, weight = w), path = "IRV")

	  # 	# Tidy DF for plotting
	  # 	rcvdf854 <- rcvdf853 %>% gather(key = "quantity", value = "value", corr_rcv2:pprob3) 
	  # 	rcvdf854 <- rcvdf854 %>% mutate(type = ifelse(quantity %in% c("corr_rcv2", "pprob2"), "IRV_second", "IRV_third"), quantity = ifelse(quantity %in% c("corr_rcv2", "corr_rcv3"), "correlation", "pprob"))
	  # 	rcvdf854 <- spread(rcvdf854, key = quantity, value = value)

	  # 	# Re-order plurality DF
	  # 	plurdf85 <- big_plur_sum %>% filter(s == 85)
	  # 	ordered_u_probs_plur <- t(apply(plurdf85, 1, rearrange))
	  # 	ordered_u_probs_plur <- apply(ordered_u_probs_plur, 2, as.numeric)

	  # 	# Tidy DF for plotting
	  # 	plurdf852 <- as.data.frame(ordered_u_probs_plur)
	  # 	plurdf852 <- plurdf852 %>% mutate(ben_p = BCpo * ((uB - uC) / 2), cost_p = ABpo * (uA - uB) + ACpo * ((uB - uC)/2), pprob_plur = BCpo
	  # 	)

	  # 	# compute correlations
	  # 	plurdf853 <- cbind(plurdf85, plurdf852)
	  # 	plurdf853 <- plurdf853 %>% group_by(case, iter) %>% summarise(correlation = cor(ben_p, cost_p), pprob = weighted.mean(pprob_plur, weight = w), type = "Plurality", path = "Plurality")

	  # 	# get country weights
	  # 	# sv_list <- list()
	  # 	n <- length(big_list_na_omit)
	  # 	country_weight <- matrix(nrow = n, ncol = 2)
	  # 	for(i in 1:n){
	  # 	  progress(i)
	  # 	  if(i == n) cat("Done! \n")
	  # 	  this_list <- big_list_na_omit[[i]]
	  # 	  dfcase <- names(big_list_na_omit)[[i]]
	  # 	  dfcountry <- substr(dfcase, 1, 3)
	  # 	  dfweight_sum <- sum(big_list_na_omit[[i]]$weights)
	  # 	  dfVAP <- vap$VAP[vap$cntry == dfcountry]
	  # 	  dfm <- vap$Freq[vap$cntry == dfcountry]
	  # 	  country_weight[i, 1] <- names(big_list_na_omit)[[i]]
	  # 	  country_weight[i, 2] <- dfVAP / dfm
	  # 	}

	  # 	# Bind DF together and tidy
	  # 	weightdf <- as.data.frame(country_weight) 
	  # 	names(weightdf) <- c("case", "ctryweight")
	  # 	weightdf$ctryweight <- as.numeric(as.character(weightdf$ctryweight))

	  # 	conjdf <- rbind(plurdf853, rcvdf854) %>% left_join(weightdf)
	  # 	conjdf_quant <- conjdf %>% group_by(iter, type) %>% summarise(
	  # 	  corr_q025 = wtd.quantile(correlation, q = 0.025, weight = ctryweight),
	  # 	  corr_q25 = wtd.quantile(correlation, q = 0.25, weight = ctryweight),
	  # 	  corr_q50 = wtd.quantile(correlation, q = 0.5, weight = ctryweight),
	  # 	  corr_q75 = wtd.quantile(correlation, q = 0.75, weight = ctryweight),
	  # 	  corr_q975 = wtd.quantile(correlation, q = 0.975, weight = ctryweight),
	  # 	  corr_mean = wtd.mean(correlation, weight = ctryweight),
	  # 	  pprob_q025 = wtd.quantile(pprob, q = 0.025, weight = ctryweight),
	  # 	  pprob_q25 = wtd.quantile(pprob, q = 0.25, weight = ctryweight),
	  # 	  pprob_q50 = wtd.quantile(pprob, q = 0.5, weight = ctryweight),
	  # 	  pprob_q75 = wtd.quantile(pprob, q = 0.75, weight = ctryweight),
	  # 	  pprob_q975 = wtd.quantile(pprob, q = 0.975, weight = ctryweight),
	  # 	  pprob_mean = wtd.mean(pprob, weight = ctryweight))

	  # 	# save conjdf
	  # 	save(conjdf, file = here(paste0(path_files, "conjdf.Rdata")))

	  # 	cat("plotting conjecture figures. \n")

	  # 	# Pprobs plot (raw)
	  # 	ggplot(conjdf, aes(x = iter)) +
	  # 	geom_line(aes(y = pprob, group = interaction(case, type), colour = type), alpha = 0.1) +
	  # 	geom_line(data = conjdf_quant, aes(x = iter, y = pprob_mean, group = type, colour = type), lwd = 2) + 
	  # 	geom_hline(yintercept = 0, lty = "dashed") +
	  # 	labs(x = "Iteration", y = "Probability vote is beneficial * electorate size") +
	  # 	theme_sv()  +
	  # 	ylim(c(0, 0.5))
	  # 	ggsave(here(paste0(path, "/conj1.pdf")), device = cairo_pdf,
	  # 	       width = 4,
	  # 	       height = 4)

	  # 	# Correlations plot
	  # 	ggplot(conjdf, aes(x = iter)) +
	  # 	geom_line(aes(y = correlation, group = interaction(case, type), colour = type), alpha = 0.05) +
	  # 	geom_line(data = conjdf_quant, aes(x = iter, y = corr_mean, group = type, colour = type), lwd = 2) + 
	  # 	geom_hline(yintercept = 0, lty = "dashed") +
	  # 	theme_sv() +
	  # 	labs(x = "Iteration", y = "Correlation") +
	  # 	theme(legend.position = "bottom",
	  # 	      legend.direction = "horizontal")
	  # 	ggsave(here(paste0(path, "/conj2.pdf")), device = cairo_pdf,
	  # 	       width = 4,
	  # 	       height= 4)
	  # }


	}
}

# cases_converge[[1]]$plur_v_vec[2, ]
# cases_converge_1[[1]]$plur_v_vec[2, ]

# cases_converge[[1]]$p_pr_list[[2]]
# cases_converge_1[[1]]$p_pr_list[[2]]

# cases_converge[[1]]$p_mat_list[2]
# cases_converge_1[[1]]$p_mat_list[2]

# cases_converge[[1]]$plur_v_vec[60, ]
# cases_converge_1[[1]]$plur_v_vec[60, ]

# cases_converge[[1]]$p_pr_list[[60]]
# cases_converge[[1]]$p_pr_list[[100]]
# cases_converge_1[[1]]$p_pr_list[[60]]
# cases_converge_1[[1]]$p_pr_list[[100]]

# piv_prob_1 = map(cases_converge[[1]]$p_pr_list, ~ .x[["AB"]]) %>% unlist
# piv_prob_2 = map(cases_converge_1[[1]]$p_pr_list, ~ .x[["AB"]]) %>%
# 	unlist

# AB_prob_tibble = tibble(iter = 1:150, 
#                         old = piv_prob_1, 
#                         new = piv_prob_2) %>%
# 	pivot_longer(old:new)

# ggplot(AB_prob_tibble, aes(iter, value)) +
# 	geom_line(aes(colour = name))
# ggsave("output/figures/piv_prob_analytics.pdf")


# (cases_converge[[1]]$p_mat_list[100])[[1]] 
# (cases_converge_1[[1]]$p_mat_list[60])[[1]] / 2

# cases_converge_old = cases_converge
# out_old = out
# summary_stats_old = summary_stats

# head(summary_stats_old)
# head(summary_stats)

# head(summary_stats_old %>% filter(System == "Plurality"))
# head(summary_stats %>% filter(System == "Plurality"))

# cases_converge[[1]][[1]] %>% head
# cases_converge_old[[1]][[1]] %>% head

# out[[1]][[2]][39:49, ]
# out_old[[1]][[2]][39:49, ]

# str(summary_stats)


