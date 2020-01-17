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

file_names <- 1:10 %>% as.character
merge_list <- list()
for(name in file_names){
	load(here(paste0("output/files/random_", name, ".Rdata")))
	for(i in 1:length(out)){
		out[[i]] <- do.call(rbind, out[[i]])
		names(out[[i]])[7:10] <- c("case", "system", "iter", "rand_iter")
	}
	full_df <- do.call(rbind, out) %>% mutate(partition = name)
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

# for each case, get the distance...

# compare the two
dist_list <- list()
for(j in names_vec){
	print(j)
	rands <- full_df %>% filter(case == j)
	base <- baseline_df %>% filter(case == j & iter == 251) %>% as.numeric
	d_obj <- rbind(base[1:6], rands[, 1:6])
	d_mat <- dist(d_obj) %>% as.matrix
	rands <- rands %>% mutate(base_d = d_mat[-1, 1])
	dist_list[[j]] <- rands
}

dist_df <- do.call(rbind, dist_list)

# print case by case
for(j in names_vec){
	print(j)
	ggplot(dist_df %>% filter(case == j), aes(iter, base_d)) +
	geom_boxplot(aes(group = iter), 
		fill = "grey80",
		outlier.size = 0.5) +
	theme_sv() +
	labs(x = "Iteration", 
		y = "Distance to voteshare at 250th (baseline)")
	ggsave(here(paste0("output/figures/algconv/", j, ".pdf")),
	device = cairo_pdf)
}

# print altogether
ggplot(dist_df, aes(iter, base_d)) +
geom_boxplot(aes(group = iter), 
	fill = "grey80",
	outlier.size = 0.5) +
theme_sv() +
labs(x = "Iteration", 
	y = "Distance to voteshare at 250th (baseline)") +
facet_wrap(~ case)
ggsave(here("output/figures/algconv/joint.pdf"),
device = cairo_pdf)

# summarise quantiles and plot
quant_df <- dist_df %>% group_by(case, iter) %>%
	summarise(mean_d = mean(base_d), Median = median(base_d),
		"90th Percentile" = quantile(base_d, 0.9),
		q95 = quantile(base_d, 0.95),
		"99th Percentile" = quantile(base_d, 0.99)) %>%
	pivot_longer(mean_d:'99th Percentile')

ggplot(quant_df, aes(x = iter)) +
	geom_point(data = quant_df %>% filter(name %in% c("Median", "90th Percentile", "99th Percentile")),
		aes(group = case, y = value, colour = name), alpha = 0.1,
		lty = "dashed") +
	theme_sv() +
	labs(x = "Iteration", y = "Distance to baseline case at 250th iteration",
		colour = "Summary statistic") +
	theme(legend.position = "bottom") +
	guides(colour = guide_legend(override.aes = list(alpha = 1)))

ggsave(here("output/figures/algconv/quantile_summary.pdf"),
	device = cairo_pdf)


# plot boxplots
ggplot(dist_df, aes(case, base_d)) +
	geom_boxplot(aes(group = case),
		fill = "grey80",
		outlier.size = 0.5) +
	coord_flip() +
	theme_sv() +
	theme(axis.text.y = element_text(size = rel(0.5))) +
	labs(x = "Distances to baseline case", y = "Case")
ggsave(here("output/figures/random_dist.pdf"),
	device = cairo_pdf,
	height = 14,
	width = 6)


# print selected cases
select_cases <- tibble(case = c("AUS_2013", "GBR_2015", "DEU_2005", "FRA_2012"),
	full_case = c("Australia (2013)", "United Kingdom (2015)", "Germany (2005)", "France (2012)"))

select_df <- full_df %>% 
	right_join(select_cases) 

ggtern(select_df, aes(V1 + V2, V3 + V4, V5 + V6)) +
	geom_line(aes(group = interaction(rand_iter, partition)),
		alpha = 0.1) +
	geom_point(data = select_df %>% filter(iter == 61),
		colour = "blue", alpha = 0.3, size = 1) +
	geom_point(data = select_df %>% filter(iter == 1),
		colour = "red", alpha = 0.3, size = 1) +
	theme_sv() +
	labs(x = "A", y = "B", z = "C") +
	facet_wrap(. ~ full_case)
ggsave(here("output/figures/random_select.pdf"),
	device = cairo_pdf,
	height = 8,
	width = 8)