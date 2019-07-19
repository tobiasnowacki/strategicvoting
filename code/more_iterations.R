## Note that right now the code runs with iterations 1-50 in sequential computing, and 51-100 in parallel computing. 
# I need to fix this when running the code again (make all of it parallel).

#===================================================
# Dependencies
#===================================================

# Load packages
requiredPackages <- c("here", "ggplot2", "ggtern", "dplyr", "purrr", "tidyr", "lmtest", "sandwich", "plm", "extrafont", "RColorBrewer", "boot", "svMisc", "ggtern", "reldist", "gridExtra", "ggpubr")
ipak <- function(pkg){
        new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
        if (length(new.pkg))
                install.packages(new.pkg, dependencies = TRUE)
        sapply(pkg, require, character.only = TRUE)
}
ipak(requiredPackages)

# Load functions
source(here("code/utils/functions.r"))
source(here("code/utils/av_pivotal_probs_analytical_general_v2.r"))
source(here("code/utils/plurality_pivotal_probabilities_analytical.r"))
source(here("code/utils/general_iteration_simulation_approach.r"))
source(here("code/utils/sv.r"))

# Load existing data
load(here("output/manyiterations.RData"))

# Load ggplot theme
theme_sv <- function(){
  theme_bw(base_size=11, base_family="Roboto Light") %+replace%
  theme(
    panel.grid.major =  element_line(
      colour = "grey50",
      size = 0.2,
      linetype = "dotted"),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "grey97"),
    plot.margin = unit(c(0.2, 1, 0.2, 1), "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.title = element_text(size = 10, family = "Roboto Medium", face = "bold"),
    strip.background = element_rect(fill= NULL, colour = "white", linetype = NULL),
    strip.text = element_text(colour = 'grey50', size = 9, vjust = 0.5, family = "Roboto Medium")
  )
}

# Packages for parallelisation.
library(foreach)
library(doParallel)


#===================================================
# Check convergence beyond 60 iterations.
#===================================================

big_rcv_sum <- list()
big_plur_sum <- list()
big_rcv_vec <- list()
big_plur_vec <- list()
piv_ratio_rcv <- list()
piv_ratio_plur <- list()

# Choice parameters
lambda <- 0.05
k <- 60
writeLines(c(""), "log_convergence.txt")

cl <- makeCluster(4)
registerDoParallel(cl)

# Loop over cases
prec <- 6
s_val <- s_list[[prec]]
rcv_sum <- list()
plur_sum <- list()
rcv_vec <- list()
plur_vec <- list()
rcv_piv <- list()
plur_piv <- list()
cases_converge <- foreach(case = 1:160, 
  	        .packages = c("gtools", "stringr", "tidyverse")
  	        ) %dopar% {
  		# sink("log_convergence.txt", append = TRUE)
	    # cat(paste(case, " . ", "\n"))
	    # need to adjust function (k is undefined)
	    out <- many_iterations_until_convergence(big_list_na_omit[[case]], big_list_na_omit[[case]]$v_vec, lambda, s_val, 0.0001, 150)
	    # rcv_sum <- out[[1]] %>% mutate(case = names(big_list_na_omit)[[case]])
	    # plur_sum <- out[[3]] %>% mutate(case = names(big_list_na_omit)[[case]])
	    # rcv_vec <- out[[2]]
	    # plur_vec <- out[[4]]
	    list(out)
}

stopCluster(cl)


# Compute Euclidean distances

conv_dist <- list()

for (j in 1:160){
	print(j)
	df <- cases_converge[[j]][[1]][[2]] 
	df_dist <- euclid_together(df) %>%
	mutate(system = 'IRV',
	       case = names(big_list_na_omit)[j],
	       iter = 1:nrow(.))
	conv_dist[[(2 * j) - 1]] <- df_dist

	df <- cases_converge[[j]][[1]][[4]] 
	df_dist <- euclid_together(df) %>%
	mutate(system = 'Plurality',
	       case = names(big_list_na_omit)[j],
	       iter = 1:nrow(.))
	conv_dist[[2 * j]] <- df_dist
}

conv_dist_df <- do.call(rbind, conv_dist)

# IRV convergence
ggplot(conv_dist_df %>% filter(system == "IRV"), aes(iter, log(diff))) +
	geom_line(aes(group = interaction(case, system)), 
	          alpha = 0.5) +
	geom_hline(yintercept = log(0.0001), colour = "red") +
	theme_sv()
ggsave(here("output/figures/convergence_irv.pdf"), device = cairo_pdf)

# Plurality convergence
ggplot(conv_dist_df %>% filter(system == "Plurality"), aes(iter, log(diff))) +
	geom_line(aes(group = interaction(case, system)), 
	          alpha = 0.5) +
	geom_hline(yintercept = log(0.0001), colour = "red") +
	theme_sv()
ggsave(here("output/figures/convergence_plur.pdf"), device = cairo_pdf)

# Pull out cases that did not converge
non_conv_cases <- conv_dist_df %>% filter(iter == 151) %>% select(case) %>% unlist
non_conv_ids <- which(names(big_list_na_omit) %in% non_conv_cases)

# Run 300+ iterations for these cases
cl <- makeCluster(4)
registerDoParallel(cl)

prec <- 6
s_val <- s_list[[prec]]

non_convergent <- foreach(case = non_conv_ids, 
  	        .packages = c("gtools", "stringr", "tidyverse")
  	        ) %dopar% {
  		# sink("log_convergence.txt", append = TRUE)
	    # cat(paste(case, " . ", "\n"))
	    # need to adjust function (k is undefined)
	    out <- many_iterations_until_convergence(big_list_na_omit[[case]], big_list_na_omit[[case]]$v_vec, lambda, s_val, 0.0001, 300)
	    return(out)
}

# add distance & iteration column to DF
non_conv_dist <- lapply(non_convergent, function(x){
	vvecs <- x[[2]]
	dist_df <- euclid_together(vvecs) %>% mutate(iter = 1:nrow(.))
	return(dist_df)
})

still_not <- which(lapply(non_conv_dist, nrow) == 301)
non_conv_dist <- non_conv_dist[still_not]

# add case number
for(i in 1:length(non_conv_dist)){
	non_conv_dist[[i]] <- non_conv_dist[[i]] %>% mutate(case = i)
}

non_conv_dist <- do.call(rbind, non_conv_dist)

# plot path
ggtern(non_conv_dist, aes(V1 + V2, V3 + V4, V5 + V6)) +
	geom_line(aes(colour = iter, alpha = iter)) +
	facet_wrap(~ case) +
	theme_sv()
ggsave(here("output/figures/non_conv_path.pdf"), device = cairo_pdf)

# plot distance
ggplot(non_conv_dist, aes(iter, log(diff))) +
	geom_line(aes(group = case,  colour = case %>% as.factor)) +
	theme_sv()
ggsave(here("output/figures/non_conv_dist.pdf"), device = cairo_pdf)



#===================================================
# Other starting points analysis.
#===================================================

uniform_ternary <- rdirichlet(100, rep(1, 6))

# Set up lists.
big_rcv_sum_sense <- list()
big_rcv_vec_sense <- list()

# Choice parameters (already in loaded .RData)
# lambda <- 0.05
# k <- 60

# This is marginally faster after I cut the functions down.

# For each randomly generated v_vec, use this as a starting point and follow the IRV learning path.
# m <- 100 # number of random draws

# for(rand_iter in 11:50){
#   prec <- 85
#   s_val <- 85
#   cat(paste0("\n === starting point = ", rand_iter, " =============== \n"))
#   # rcv_sum <- list()
#   rcv_vec <- list()
#   # rcv_piv <- list()
#   rand_v_vec <- uniform_ternary[rand_iter, ] %>% as.numeric
#   for (case in 1:160) {
#     cat(paste0(case, ": ", names(big_list_na_omit)[case], "   "))
#     out <- many_iterations_rcv_only(big_list_na_omit[[case]], rand_v_vec, lambda, s_val, k)
#     # rcv_sum[[case]] <- cbind(out[[1]], names(big_list_na_omit)[[case]], rand_iter)
#     rcv_vec[[case]] <- cbind(out, names(big_list_na_omit)[[case]],
#                        "IRV", 
#                        1:61)
#     # rcv_piv[[case]] <- piv_ratio(out[[1]])
#   }
#   # rcv_sum <- do.call(rbind, rcv_sum)
#   # big_rcv_sum_sense[[rand_iter]] <- rcv_sum
#   big_rcv_vec_sense[[rand_iter]] <- rcv_vec
# }

# # convert list into DF
# iter_df <- do.call(rbind, lapply(big_rcv_vec_sense, function(x) do.call(rbind, x))) %>% as.data.frame
# iter_df$rand_iter <- rep(1:length(big_rcv_vec_sense), each = 9760)

# # convert shares into numeric format
# iter_df[, c(1:6, 9)] <- apply(iter_df[, c(1:6, 9)], 2, as.numeric)
# names(iter_df)[c(7, 9)] <- c("case", "iter")


cl <- makeCluster(7)
registerDoParallel(cl)


# For loop -- every iteration is a random vvec
out <- foreach(rand_iter = 1:100,
               .packages = c("gtools", "stringr")
               ) %dopar% {
	prec <- 85
	s_val <- 85
	cat(paste0("\n === starting point = ", rand_iter, " =============== \n"))
	rcv_vec <- list()
	rand_v_vec <- uniform_ternary[rand_iter, ] %>% as.numeric
	for (case in 1:160) {
	    cat(paste0(case, ": ", names(big_list_na_omit)[case], "   "))
	    out <- many_iterations_rcv_only(big_list_na_omit[[case]], rand_v_vec, lambda, s_val, k)
	    rcv_vec[[case]] <- cbind(out, names(big_list_na_omit)[[case]],
	                       "IRV", 
	                       1:61, rand_iter)
	}
	rcv_vec
}
stopCluster(cl)

# Unpack
out2 <- lapply(out, function(x) do.call(rbind, x)) %>% do.call(rbind, .)

out3 <- rbind(iter_df, out2)
names(out3)[7:10] <- c("case", "system", "iter", "rand_iter") 
out3[, c(1:6, 9)] <- apply(out3[, c(1:6, 9)], 2, as.numeric)


# original paths (code from paper_code_eqm.r)
# requires pre-loaded data (intermediate3.RData)
rcv_vec_df <- big_rcv_vec[[6]] %>% do.call(rbind, .) %>% mutate(case = rep(names(big_list_na_omit), each = 61), iter = rep(1:61, 160), A = V1 + V2, B = V3 + V4, C = V5 + V6)

# =====================================================
# Plot iteration paths
# =====================================================

ggtern(out3, aes(V1 + V2, V3 + V4, V5 + V6)) +
	# Plot starting points
	geom_line(aes(group = interaction(case, rand_iter)), 
	          colour = "grey", 
	          alpha = 0.3) +
	geom_point(data = out3 %>% filter(iter == 1),
				colour = "red", 
				alpha = 0.3,
				size = 1) +
	# Plot end points
	geom_point(data = out3 %>% filter(iter == 61), 
	           colour = "blue", 
	           alpha = 0.3, 
	           size = 1) +
	#Plot starting points (true)
  	geom_point(data = rcv_vec_df[rcv_vec_df$iter == 1, ], 
             alpha = .7,
             colour = "orange", 
             size = 1.3) + 
  	# Plot end points (true)
  	geom_point(data = rcv_vec_df[rcv_vec_df$iter == 61, ], 
             size = 1.3, 
             colour = "#00BCD6", 
             alpha = .7) + 
  	# Mark "true" path with red line
  	geom_line(data = rcv_vec_df, 
  	          aes(group = case), colour = "red") +
	facet_wrap(~ case) +
	theme_sv() +
	labs(x = "A", y = "B", z = "C")

# Export
ggsave(here("output/figures/random_starting_iteration.pdf"), width = 20, height = 20, device = cairo_pdf)

save.image(here("output/manyiterations.Rdata"))


# =====================================================
# Random vvecs with smaller distance
# =====================================================

cl <- makeCluster(7)
registerDoParallel(cl)

case_list <- names(big_list_na_omit)
writeLines(c(""), "log.txt")

restricted_iters <- foreach(i = 1:160,
        .packages = c("gtools", "stringr", "tidyverse")) %dopar% {
	sink("log.txt", append = TRUE)
	orig_vvec <- rcv_vec_df %>% filter(case == case_list[i] & iter == 1) %>% select(1:6)
	eqm_vvec <- rcv_vec_df %>% filter(case == case_list[i] & iter == 61) %>% select(1:6)
	draw_vecs <- draw_restricted_vvecs(orig_vvec, eqm_vvec, 50)
	case_iters <- list()
	for(j in 1:nrow(draw_vecs)){
		cat(paste(case_list[i], " . ", j,"\n"))
		out <- many_iterations_rcv_only(big_list_na_omit[[i]], draw_vecs[j, ] %>% as.numeric, lambda, s_val, k)
	    # rcv_sum[[case]] <- cbind(out[[1]], names(big_list_na_omit)[[case]], rand_iter)
	    case_iters[[j]] <- cbind(out, case_list[i],
	                       "IRV", 
	                       1:61, j)
	}
	case_iters
}

stopCluster(cl)

#===================================================
# Sensitivity to choice of lambda
#===================================================

# TODO still
