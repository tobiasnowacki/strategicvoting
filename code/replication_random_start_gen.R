# paste code from more_iterations.R
# run iterations
# check distances (numerically)
# 	a. mean distance between 60th iteration
#	b. mean distance to starting point
#	c. mean distance to 'true' 60th iteration

cmd_line_args <- commandArgs(trailingOnly = TRUE)

# load dependencies
library(here)						# to get dir
source(here("code/full_header.R")) 	# fn's and data
source(here("code/prep_cses.R")) 	# data prep

set.seed(cmd_line_args[1])

names_vec <- names(big_list_na_omit)
v_vec_list <- list()

uniform_ternary <- rdirichlet(10, rep(1, 6))

# Set up lists.
big_rcv_sum_sense <- list()
big_rcv_vec_sense <- list()

# Set up parameters.
lambda <- 0.05
s_val <- 85

cl <- makeCluster(7)
registerDoParallel(cl)

# For loop -- every iteration is a random vvec
out <- foreach(rand_iter = 1:10,
               .packages = c("gtools", "stringr", "tidyverse")
               ) %dopar% {
	prec <- 85
	s_val <- 85
	cat(paste0("\n === starting point = ", rand_iter, " =============== \n"))
	rcv_vec <- list()
	rand_v_vec <- uniform_ternary[rand_iter, ] %>% as.numeric
	for (case in 1:160) {
	    cat(paste0(case, ": ", names(big_list_na_omit)[case], "   "))
	    out <- many_iterations_rcv_only_light(big_list_na_omit[[case]], rand_v_vec, lambda, s_val, 60)
	    rcv_vec[[case]] <- cbind(out, names(big_list_na_omit)[[case]],
	                       "IRV", 
	                       1:61, rand_iter)
	}
	rcv_vec
}
stopCluster(cl)

cat("Done.")

save.image(here(paste0("output/files/random_", cmd_line_args[1], ".Rdata")))


# ---
