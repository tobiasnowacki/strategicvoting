### PRELIMINARIES
# install.packages("here")
library(here) 						# to get dir
source(here("code/full_header.R")) 	# fn's and data
source(here("code/prep_cses.R")) 	# data prep

s_param <- c("10", "55", "85")
lambda_param <- c("1", "2", "3")
names_vec <- names(big_list_na_omit)


### V_VEC DISTANCES

v_vec_list <- list()
for(i in s_param){
	for(j in lambda_param){
		# load file
		datapath <- paste0("output/files/", j, "/", i, "v_vecs_", j, "_", i, ".Rdata")
		load(here(datapath))
		# append cases w/ relevant data
		for(n in 1:160){
			item1 <- out[[n]][[1]] %>% mutate(iter = 1:251,
									 case = names_vec[n],
									 system = "IRV",
									 s = i,
									 lambda = j)
			item2 <- out[[n]][[2]] %>% mutate(iter = 1:251,
									 case = names_vec[n],
									 system = "Plurality",
									 s = i,
									 lambda = j)
			v_vec_list[[paste0(i, "_", j, "_", names_vec[n])]] <- rbind(item1, item2)
		}
	}
}


### EXPECTED BENEFIT SUMMARY ---

summary_list <- list()
for(i in s_param){
	for(j in lambda_param){
		# load file
		datapath <- paste0("output/files/", j, "/", i, "summary_", j, "_", i, ".Rdata")
		load(here(datapath))
		# append cases w/ relevant data
		summary_list[[paste0(j, "_", i)]] <- summary_stats %>% 
			mutate(s = i, 
				lambda = j)
		}
	}	

# For each element in summary_list, perform X




# first plot param configs

# v_vec_distances
plot_v_vec_distance(cases_converge, path, 
                                n_lag = 20, 
                                avg_span = 10)

# exp benefit summaries
summary_stats_wide <- summary_stats %>% 
  gather(., key = "Statistic", 
         value = "Value", "Prevalence":"ExpBenefit")

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
ggsave(here(paste0(path, "/main_results.pdf")), 
       device = cairo_pdf)

cat("Summary statistics plotted. \n")