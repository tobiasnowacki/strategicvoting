# New file to get v_vec data together and create new plots / distance tables

source("code/prep_cses.R") 	# data prep
source("code/utils/sv_theme_template.R")

s_param <- c("10", "55", "85")
lambda_param <- c("1", "2", "3")
names_vec <- names(big_list_na_omit)
v_vec_list <- list()

# load v_vec data and put into one big table
# I don't think I need the plurality results here...?
for(i in s_param){
  for(j in lambda_param){
    # load file
    datapath <- paste0("output/files/", j, "/", i, "_vvec.Rdata")
    load(datapath)
    # append cases w/ relevant data
    for(n in 1:160){
      item1 <- vvecdf[[n]][[1]] %>% mutate(iter = 1:250,
                                        case = names_vec[n],
                                        system = "IRV",
                                        s = i,
                                        lambda = j)
      # item2 <- vvecdf[[n]][[2]] %>% mutate(iter = 1:250,
      #                                   case = names_vec[n],
      #                                   system = "Plurality",
      #                                   s = i,
      #                                   lambda = j)
      v_vec_list[[paste0(i, "_", j, "_", names_vec[n])]] <- item1
    }
  }
}



v_vec_df <- do.call(rbind, v_vec_list) %>% 
  mutate(params = paste0(s, "_", lambda))
table(v_vec_df$params)

# first: plot v_vec paths not necessary since done in other script!
# for(val in unique(v_vec_df$params)){
#   print(val)
#   sub_df <- v_vec_df %>% filter(params == val & iter %in% c(1, 251)) %>% 
#     mutate(Iteration = ifelse(iter == 1, "First", "Last"))
#   ll <- substr(val, 4, 4)
#   ss <- substr(val, 1, 2)
#   figpath <- paste0("output/figures/", ll, "/", ss, "/v_vec_path_v2.pdf")
#   ggtern(sub_df, aes(abc + acb, bac + bca, cab + cba)) +
#     geom_line(aes(group = interaction(case)),
#               alpha = 0.1) +
#     geom_point(data = sub_df,
#                aes(colour = Iteration),
#                size = 1.2,
#                alpha = .5) +
#     facet_wrap(~ system) +
#     #theme_sv() +
#     labs(x = "A", y = "B", z = "C") +
#     scale_colour_manual(values = c("red", "blue")) +
#     theme(panel.spacing = unit(0, "cm"),
#           plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
#           legend.position = "bottom",
#           legend.direction = "horizontal") +
#     theme_sv()
#   ggsave(here(figpath), width = 8, height = 3)
#   
# }

# second: compare distances between jth iteration and baseline
irv_d_list <- list()
plur_d_list <- list()
comp_iter <- c(61, 101, 250)
avg_lag <- 0

for(case in 1:160){
  casename <- names_vec[case]
  for(j in comp_iter){
    for(s_val in c(10, 55, 85)){
      list_case <- paste0(case, "_", j)
      # get baseline vec -- always 250th iteration...
      first_vec <- v_vec_df %>% 
        filter(iter == 250 & 
                 params == "85_1" &
                 case == casename & 
                 system == "IRV")
      # compute distance for IRV iterations
      irv_vecs <- v_vec_df %>% 
        filter(iter %in% ((j-avg_lag):j),
               case == casename &
                 system == "IRV") %>%
        group_by(s, lambda) %>%
        summarise(V1 = mean(abc),
                  V2 = mean(acb),
                  V3 = mean(bac),
                  V4 = mean(bca),
                  V5 = mean(cab),
                  V6 = mean(cba))
      names(irv_vecs)[3:8] = names(first_vec)[1:6]
      joint_mat <- rbind(first_vec[1:6], irv_vecs[, 3:8])
      d <- dist(joint_mat) %>% 
        as.matrix
      irv_vecs <- irv_vecs %>% 
        as.data.frame %>% 
        mutate(d_start = d[2:10, 1],
                                                        comp_iter = j)
      irv_d_list[[list_case]] <- irv_vecs
    }
    
    
  }	
}

irv_d <- do.call(rbind, irv_d_list) %>% 
  mutate(s = recode(s, "10" = "s = 10", "55" = "s = 55", "85" = "s = 85"), 
         lambda = recode(lambda, "1" = "lambda = 0.05", "2" = "lambda = 0.1", "3" = "lambda = 0.01"))
plur_d <- do.call(rbind, plur_d_list)


# facet plot
ggplot(irv_d %>% filter(!(s == "s = 85" & lambda == "lambda = 0.05")),
       aes(d_start)) +
  geom_density(aes(y = ..scaled..,
                   fill = as.factor(comp_iter)), 
               alpha = 0.3) +
  facet_grid(s ~ lambda) +
  theme_tn() +
  labs(x = "Distance to baseline case",
       y = "(Normalized) density",
       fill = "After ... iterations") +
  theme(legend.position = "bottom",
        legend.direction = "horizontal")
ggsave("output/figures/distance_to_baseline.pdf")

# third: compare distances between 250th iterations for each case, conditional on s.
# effectively just modified code from above (not great, need to fix in the future)
# v_vec_df %>% filter(iter %in% c(1, 250)) %>% head()

irv_d_list <- list()
plur_d_list <- list()
comp_iter <- c(61, 101, 250)
avg_lag <- 0

for(case in 1:160){
  casename <- names_vec[case]
  for(j in comp_iter){
    for(s_val in c(10, 55, 85)){
      list_case <- paste0(case, "_", s_val, "_", j)
      # get baseline vec -- always 250th iteration...
      first_vec <- v_vec_df %>% 
        filter(iter == 250 & 
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
        summarise(abc = mean(abc),
                  acb = mean(acb),
                  bac = mean(bac),
                  bca = mean(bca),
                  cab = mean(cab),
                  cba = mean(cba))
      joint_mat <- rbind(first_vec[1:6], irv_vecs[, 3:8])
      d <- dist(joint_mat) %>% as.matrix
      irv_vecs <- irv_vecs %>% 
        as.data.frame %>% 
        mutate(d_start = d[2:4, 1],
               comp_iter = j, case = casename)
      irv_d_list[[list_case]] <- irv_vecs
    }
  }	
}

irv_d <- do.call(rbind, irv_d_list) %>% mutate(s = recode(s, "10" = "s = 10", "55" = "s = 55", "85" = "s = 85"), 
                                               lambda = recode(lambda, "1" = "lambda = 0.05", "2" = "lambda = 0.1", "3" = "lambda = 0.01"))

# quick check that cases with lambda = 0.05 have mass unit at zero.
irv_d %>% filter(lambda == "lambda = 0.05" & comp_iter == "250") %>% head

# facet plot
ggplot(irv_d %>% filter(!lambda == "lambda = 0.05"),
       aes(d_start)) +
  geom_density(aes(y = ..scaled..,
                   fill = as.factor(comp_iter)), 
               alpha = 0.5) +
  facet_grid(s ~ lambda) +
  theme_tn() +
  labs(x = "Distance to baseline with same precision (s)",
       y = "(Normalized) density",
       fill = "After ... iterations") +
  theme(legend.position = "bottom",
        legend.direction = "horizontal")
ggsave("output/figures/distance_to_baseline_by_s.pdf")



