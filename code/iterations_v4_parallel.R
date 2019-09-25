## Full script to generate results for paper
## 7 August 2019.
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
lambda_list <- as.list(c(0.05, 0.1, 0.5))

# Convergence distance criterion (ε)
epsilon_thresh <- 0.001

# Maximum number of iterations
k_max <- 300

### ---
### RUN ITERATIONS

# For loop over values of s

# writeLines(c(""), "log_convergence.txt")

save_out <- list()

# temporary (later replaced by loop)
# lambda_counter <- 1
# for(lambda_val in c(1, 2)){
for(lambda_val in c(1)){
	lambda <- lambda_list[[lambda_val]]

	cat(paste0("=== Beginning loop for lambda == ", lambda, ". === \n"))

	save_out[[lambda_val]] <- list()

	# for(s_loop_val in c(1, 4, 6)){ 
	for(s_loop_val in c(6)){ 

	  s_val <- s_list[[s_loop_val]]

	  cat(paste0("=== Beginning loop for s == ", s_val, ". === \n"))

	  which_cases <- 1:160   # default is 1:160
	  # max_iter_val <- 250	
	  max_iter_val <- 250	
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

	  save_out[[lambda_val]][[paste0(s_val)]] <- cases_converge

	  cat("Done. \n")

	  path <- paste0("output/figs_v2/", lambda_val, "/", s_val)

	  # name list obj
	  names(cases_converge) <- names(big_list_na_omit)[which_cases]

	  # produce iteration plot(s)
	  # will need to adjust euclid because now the distance is compared to best response v_vec!
	  test_out <- plot_v_vec_distance(cases_converge, path, 
	                                  n_lag = 20, 
	                                  avg_span = 10)

	  cat("v_vec distance plotted. \n")

	  # produce v_vec path plot(s)
	  # all together
	  joint_v_vec_plot(cases_converge, path)

	  # distance to t-50 avg
	  # early_late_v_vec_plot()

	  cat("v_vec paths plotted. \n")

	  # rcv separate ones by non_convergence
	  non_conv_v_vec_plot(cases_converge, path, max_iter_val)
	  vote_tally <- non_conv_strat_votes(cases_converge, path, max_iter_val)

	  # rcv plot non_convergent
	  unique_nc_cases <- unique(vote_tally$case)

	  for(i in unique_nc_cases){
	    ggplot(vote_tally %>% filter(case == i), 
	           aes(iter, freq)) +
	      geom_line(aes(group = opt_rcv, colour = opt_rcv %>% as.factor)) +
	      facet_grid(~ sin_rcv) +
	      theme_sv()
	      ggsave(here(paste0(path, "/nc_", i, ".pdf")),
	             device = cairo_pdf)
	  }

	  cat("non-convergent cases plotted. \n")

	  # group expected benefit and the like
	  summary_stats <- get_sum_stats(cases_converge)
	  # todo here: (a) correct weights (just weighted means)
	  #        (b) produce averages across cases.
	  summary_stats_wide <- summary_stats %>% 
	    gather(., key = "Statistic", 
	           value = "Value", "Prevalence":"ExpBenefit")

	  summary_agg <- summary_stats_wide %>% group_by(iter, Statistic, System) %>%
	    summarise(Value = mean(Value))

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
	  ggsave(here(paste0(path, "/main_results.pdf")), 
	         device = cairo_pdf)

	  cat("Summary statistics plotted. \n")

	}
}

# for each s in s_list
## run loop over cases => cases_converge
## get iteration plot
## get v_vec paths
## get non-convergent ones (SV type)
## get expected benefit, magnitude, etc
## Done!.

# TODO still: 'diagnostic plots' (case-by-case)
# Another question is: what about random starting points?
# Should I include this into the script?

### ROBUSTNESS CHECKS
### ---------------
# # I should check if I get the same with Andy's function.
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

# # Create distance vector
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

# # Get magnitude.


# ## I'm getting the same results!
# ## But the distance plots are different.
# ## That's probably because they weight different things.
# ## Weighting.

### -------------------------
