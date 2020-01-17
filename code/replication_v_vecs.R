# New file to get v_vec data together and create new plots / distance tables

# load dependencies
library(here)						# to get dir
source(here("code/full_header.R")) 	# fn's and data
source(here("code/prep_cses.R")) 	# data prep

s_param <- c("10", "55", "85")
lambda_param <- c("1", "2", "3")
names_vec <- names(big_list_na_omit)
v_vec_list <- list()

# load v_vec data and put into one big table
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

v_vec_df <- do.call(rbind, v_vec_list) %>% mutate(params = paste0(s, "_", lambda))
table(v_vec_df$params)

# first: plot v_vec paths
for(val in unique(v_vec_df$params)){
	print(val)
	sub_df <- v_vec_df %>% filter(params == val & iter %in% c(1, 251)) %>% 
		mutate(Iteration = ifelse(iter == 1, "First", "Last"))
	ll <- substr(val, 4, 4)
	ss <- substr(val, 1, 2)
	figpath <- paste0("output/figures/", ll, "/", ss, "/v_vec_path_v2.pdf")
	ggtern(sub_df, aes(V1 + V2, V3 + V4, V5 + V6)) +
		geom_line(aes(group = interaction(case)),
			alpha = 0.1) +
		geom_point(data = sub_df,
			aes(colour = Iteration),
			size = 1.2,
			alpha = .5) +
		facet_wrap(~ system) +
		#theme_sv() +
		labs(x = "A", y = "B", z = "C") +
		scale_colour_manual(values = c("red", "blue")) +
		theme(panel.spacing = unit(0, "cm"),
			plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
			legend.position = "bottom",
			legend.direction = "horizontal")
	ggsave(here(figpath), width = 8, height = 3)

}

# second: compare distances between jth iteration and baseline
irv_d_list <- list()
plur_d_list <- list()
comp_iter <- c(61, 101, 251)
avg_lag <- 0

for(case in 1:160){
	casename <- names_vec[case]
	for(j in comp_iter){
		for(s_val in c(10, 55, 85)){
			list_case <- paste0(case, "_", j)
		# get baseline vec -- always 250th iteration...
		first_vec <- v_vec_df %>% 
			filter(iter == 251 & 
					params == "85_1" &
					case == casename & 
					system == "IRV")
		# compute distance for IRV iterations
		irv_vecs <- v_vec_df %>% 
			filter(iter %in% ((j-avg_lag):j),
					case == casename &
					system == "IRV") %>%
			group_by(s, lambda) %>%
			summarise(V1 = mean(V1),
				V2 = mean(V2),
				V3 = mean(V3),
				V4 = mean(V4),
				V5 = mean(V5),
				V6 = mean(V6))
		joint_mat <- rbind(first_vec[1:6], irv_vecs[, 3:8])
		d <- dist(joint_mat) %>% as.matrix
		irv_vecs <- irv_vecs %>% as.data.frame %>% mutate(d_start = d[2:10, 1],
			comp_iter = j)
		irv_d_list[[list_case]] <- irv_vecs
		}
		

	}	
}

# SNIPPET
# compute distance for plurality iterations
# plur_vecs <- v_vec_df %>% 
# 	filter(iter == j,
# 			case == casename &
# 			system == "Plurality")
# joint_mat <- rbind(first_vec[1:6], plur_vecs[, 1:6])
# d <- dist(joint_mat) %>% as.matrix
# plur_vecs <- plur_vecs %>% mutate(d_start = d[2:10, 1],
# 	comp_iter = j)
# plur_d_list[[list_case]] <- plur_vecs

irv_d <- do.call(rbind, irv_d_list) %>% mutate(s = recode(s, "10" = "s = 10", "55" = "s = 55", "85" = "s = 85"), lambda = recode(lambda, "1" = "lambda = 0.05", "2" = "lambda = 0.1", "3" = "lambda = 0.01"))
plur_d <- do.call(rbind, plur_d_list)

# check that baseline case is unit mass at zero
# no longer applies if we're taking distance to 250th iteration.
# irv_d %>% filter(s == "s = 85", lambda == "lambda = 0.05") %>% summarise(sd(d_start))

# facet plot
ggplot(irv_d %>% filter(!(s == "s = 85" & lambda == "lambda = 0.05")),
						aes(d_start)) +
	geom_density(aes(y = ..scaled..,
		fill = as.factor(comp_iter)), 
	alpha = 0.5) +
	facet_grid(s ~ lambda) +
	theme_sv() +
	labs(x = "Distance to baseline case",
		y = "(Normalised) density",
		fill = "After ... iterations") +
	theme(legend.position = "bottom",
		legend.direction = "horizontal")
ggsave(here("output/figures/distance_to_baseline.pdf"))

# third: compare distances between 250th iterations for each case, conditional on s.
# effectively just modified code from above (not great, need to fix in the future)
v_vec_df %>% filter(iter %in% c(1, 251)) %>% head()

irv_d_list <- list()
plur_d_list <- list()
comp_iter <- c(61, 101, 251)
avg_lag <- 0

for(case in 1:160){
	casename <- names_vec[case]
	for(j in comp_iter){
		for(s_val in c(10, 55, 85)){
			list_case <- paste0(case, "_", s_val, "_", j)
			# get baseline vec -- always 250th iteration...
			first_vec <- v_vec_df %>% 
				filter(iter == 251 & 
						s == s_val &
						lambda == 1 &
						case == casename & 
						system == "IRV")
			# compute distance for IRV iterations
			irv_vecs <- v_vec_df %>% 
				filter(iter %in% ((j-avg_lag):j),
						s == s_val &
						case == casename &
						system == "IRV") %>%
				group_by(s, lambda) %>%
				summarise(V1 = mean(V1),
					V2 = mean(V2),
					V3 = mean(V3),
					V4 = mean(V4),
					V5 = mean(V5),
					V6 = mean(V6))
			joint_mat <- rbind(first_vec[1:6], irv_vecs[, 3:8])
			d <- dist(joint_mat) %>% as.matrix
			irv_vecs <- irv_vecs %>% as.data.frame %>% mutate(d_start = d[2:4, 1],
				comp_iter = j, case = casename)
			irv_d_list[[list_case]] <- irv_vecs
		}
	}	
}

irv_d <- do.call(rbind, irv_d_list) %>% mutate(s = recode(s, "10" = "s = 10", "55" = "s = 55", "85" = "s = 85"), 
	lambda = recode(lambda, "1" = "lambda = 0.05", "2" = "lambda = 0.1", "3" = "lambda = 0.01"))

# quick check that cases with lambda = 0.05 have mass unit at zero.
irv_d %>% filter(lambda == "lambda = 0.05" & comp_iter == "251") %>% head

# facet plot
ggplot(irv_d %>% filter(!lambda == "lambda = 0.05"),
						aes(d_start)) +
	geom_density(aes(y = ..scaled..,
		fill = as.factor(comp_iter)), 
	alpha = 0.5) +
	facet_grid(s ~ lambda) +
	theme_sv() +
	labs(x = "Distance to baseline with same precision (s)",
		y = "(Normalised) density",
		fill = "After ... iterations") +
	theme(legend.position = "bottom",
		legend.direction = "horizontal")
ggsave(here("output/figures/distance_to_baseline_by_s.pdf"))

