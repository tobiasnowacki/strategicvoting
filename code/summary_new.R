### LOAD DEPENDENCIES -------

# install.packages("here")
library(here) 						# to get dir
source(here("code/full_header.R")) 	# fn's and data
source(here("code/prep_cses.R")) 	# data prep

# 
s_param <- c("10", "55", "85")
lambda_param <- c("1", "2", "3")


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

hidden_scale_adj <- data.frame(iter = 3, value = 0.69, name = "Prevalence", s = "s = 10")

summary_df <- do.call(rbind, summary_list) %>% 
	mutate(cntry = substr(case, 1, 3)) %>% 
	inner_join(case_weight_tbl) %>% 
	pivot_longer(Prevalence:ExpBenefit) %>% 
	mutate(s = recode(s, `10` = "s = 10", `55` = "s = 55", `85` = "s = 85"))

agg_df <- summary_df %>% 
	group_by(iter, s, lambda, name, System) %>% 
	summarise(value = wtd.mean(value, case_weight))	

main_plot <- ggplot(summary_df %>% filter(lambda == 1 & iter < 61), aes(iter, value)) +
	geom_line(aes(group = interaction(case, System), colour = System), alpha = 0.1) +
	scale_color_manual(values = cbbPalette[c(3, 2)]) +
	geom_point(data = hidden_scale_adj, alpha = 0) +
	geom_line(data = agg_df %>% filter(System == "IRV" & lambda == 1 & iter < 61), 
	          color = "#004C99", lwd = 1.1) + 
	geom_line(data = agg_df %>% filter(System == "Plurality" & lambda == 1 & iter < 61), 
	          color = "#CC6600", lwd = 1.1) + 
	facet_wrap(s ~ name, scales = "free_y") +
	theme_sv() +
	guides(colour = guide_legend(override.aes = list(alpha = 1))) +
	              theme(legend.position = "bottom", legend.direction = "horizontal") +
	labs(x = "Degree of Strategicness (Iterations)", y = "Value")
ggsave(here("output/figs_v2/iterated_complete.pdf"), main_plot, device = cairo_pdf, width = 6, height = 6.5)
