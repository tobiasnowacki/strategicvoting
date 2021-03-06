# paste code from more_iterations.R
# run iterations
# check distances (numerically)
# 	a. mean distance between 60th iteration
#	b. mean distance to starting point
#	c. mean distance to 'true' 60th iteration

# load dependencies
library(here)						# to get dir
source(here("code/full_header.R")) 	# fn's and data
source(here("code/prep_cses.R")) 	# data prep

names_vec <- names(big_list_na_omit)
v_vec_list <- list()

uniform_ternary <- rdirichlet(100, rep(1, 6))

# Set up lists.
big_rcv_sum_sense <- list()
big_rcv_vec_sense <- list()

# Set up parameters.
lambda <- 0.05
s_val <- 85

cl <- makeCluster(7)
registerDoParallel(cl)

# For loop -- every iteration is a random vvec
out <- foreach(rand_iter = 1:2,
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

save.image(here("output/files/test.Rdata"))
# ---

# put together into one DF
file_names <- c("31", "36", "37", "38", "40", "51", "60")
merge_list <- list()
for(name in file_names){
	load(here(paste0("output/files/test_", name, ".Rdata")))
	for(i in 1:length(out)){
		out[[i]] <- do.call(rbind, out[[i]])
		names(out[[i]])[7:10] <- c("case", "system", "iter", "rand_iter")
	}
	full_df <- do.call(rbind, out) %>% filter(iter == 61) %>% mutate(partition = name)
	merge_list[[name]] <- full_df
}
full_df <- do.call(rbind, merge_list)

# load baseline case and put into one table
load(here("output/files/1/85v_vecs_1_85.Rdata"))

v_vec_list <- list()
for(n in 1:160){
	item1 <- out[[n]][[1]] %>% mutate(iter = 1:251,
							 case = names_vec[n],
							 system = "IRV",
							 s = 85,
							 lambda = 0.05)
	v_vec_list[[names_vec[n]]] <- item1
}
baseline_df <- do.call(rbind, v_vec_list)

# compare the two
dist_list <- list()
for(j in names_vec){
	rands <- full_df %>% filter(case == j)
	base <- baseline_df %>% filter(case == j & iter == 61) %>% as.numeric
	d_obj <- rbind(base[1:6], rands[, 1:6])
	d_mat <- dist(d_obj) %>% as.matrix
	rands <- rands %>% mutate(base_d = d_mat[-1, 1])
	dist_list[[j]] <- rands
}

dist_df <- do.call(rbind, dist_list)

# plot boxplots
ggplot(dist_df, aes(case, base_d)) +
	geom_boxplot(aes(group = case),
		fill = "grey80",
		outlier.size = 0.5) +
	coord_flip() +
	theme_sv() +
	theme(axis.text.y = element_text(size = rel(0.5))) +
	labs(x = "Distances to baseline case", y = "Case")
ggsave(here("output/figs_v2/random_dist.pdf"),
	device = cairo_pdf,
	height = 14,
	width = 6)