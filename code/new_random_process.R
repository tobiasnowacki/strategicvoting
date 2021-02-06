# process iterations with random starting points
# produce figures:
#   Distance to baseline convergence
#   Quantile summary
#   Selected cases (main paper figure)

# LOAD DEPENDENCIES ---------------------------------------------
# Load functions
source("code/utils/new_sv_iter.R")
source("code/utils/sv_theme_template.R")
source("code/utils/ternary_functions.R")

# Load data
load("output/big_list_2.RData")
vap <- read.csv("data/case_vap.csv", sep = "") # voting age pop.
cat("Data imported. \n")
source("code/prep_cses.R")  # prep CSES

# Load random start iterations
file_names <- 1:10 %>% as.character
merge_list <- list()
for(name in file_names){
  load(paste0("output/files/random_", name, ".Rdata"))
  for(i in 1:length(out)){
    out[[i]] <- do.call(rbind, out[[i]])
    names(out[[i]])[7:10] <- c("case", "system", "iter", "rand_iter")
  }
  full_df <- do.call(rbind, out) %>% mutate(partition = name)
  merge_list[[name]] <- full_df
}
full_df <- do.call(rbind, merge_list)

# load baseline case and put into one table
load("output/files/1/85_vvec.Rdata")
v_vec_list <- list()
for(n in 1:160){
  item1 <- vvecdf[[n]][[1]] %>% 
    mutate(iter = 1:250,
          case = names_vec[n],
          system = "IRV",
          s = 85,
          lambda = 0.05)
  v_vec_list[[names_vec[n]]] <- item1
}
baseline_df <- do.call(rbind, v_vec_list)

# COMPARE RANDOM STARTING ITERS TO BASELINE -----------------
dist_list <- list()
for(j in names_vec){
  print(j)
  rands <- full_df %>% filter(case == j)
  base <- baseline_df %>% filter(case == j & iter == 250) %>% as.numeric
  d_obj <- rbind(base[1:6], rands[, 1:6])
  d_mat <- dist(d_obj) %>% as.matrix
  rands <- rands %>% mutate(base_d = d_mat[-1, 1])
  dist_list[[j]] <- rands
}
dist_df <- do.call(rbind, dist_list)

# ANALYSIS AND PLOTTING ------------------------------------
# Plot distance to baseline case-by-case
for(j in names_vec){
  print(j)
  ggplot(dist_df %>% filter(case == j), aes(iter, base_d)) +
    geom_boxplot(aes(group = iter), 
                 fill = "grey80",
                 outlier.size = 0.5) +
    theme_tn() +
    labs(x = "Iteration", 
         y = "Distance to voteshare at 250th (baseline)")
  ggsave(paste0("output/figures/algconv/", j, ".pdf"),
         device = cairo_pdf)
}

# # print all together
# ggplot(dist_df, aes(iter, base_d)) +
#   geom_boxplot(aes(group = iter), 
#                fill = "grey80",
#                outlier.size = 0.5) +
#   theme_tn() +
#   labs(x = "Iteration", 
#        y = "Distance to voteshare at 250th (baseline)") +
#   facet_wrap(~ case)
# ggsave("output/figures/algconv/joint.pdf",
#        device = cairo_pdf)

# summarise quantiles and plot
quant_df <- dist_df %>% group_by(case, iter) %>%
  summarise(mean_d = mean(base_d), Median = median(base_d),
            "90th Percentile" = quantile(base_d, 0.9),
            q95 = quantile(base_d, 0.95),
            "99th Percentile" = quantile(base_d, 0.99)) %>%
  pivot_longer(mean_d:'99th Percentile')

p_quant = ggplot(quant_df, aes(x = iter)) +
  geom_point(data = quant_df %>% filter(name %in% c("Median", "90th Percentile", "99th Percentile")),
             aes(group = case, y = value, colour = name), alpha = 0.05,
             lty = "dashed") +
  theme_tn() +
  labs(x = "Iteration", y = "Distance to baseline case at 250th iteration",
       colour = "Summary statistic") +
  theme(legend.position = "bottom") +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))

p_quant = ggplot2:::print.ggplot(p_quant)
ggsave(p_quant, "output/figures/algconv/quantile_summary.pdf")


# print selected cases
select_cases <- tibble(case = c("AUS_2013", "GBR_2015", "DEU_2005", "FRA_2012"),
                       full_case = c("Australia (2013)", "United Kingdom (2015)", "Germany (2005)", "France (2012)"))
select_df <- full_df %>% 
  right_join(select_cases) 

library(ggtern) # can use here because no more ggplots...
ggtern(select_df, aes(abc + acb, bac + bca, cab + cba)) +
  geom_line(aes(group = interaction(rand_iter, partition)),
            lwd = 0.5,
            alpha = 0.1) +
  geom_point(data = select_df %>% filter(iter == 60),
             colour = "blue", alpha = 0.3, size = 1) +
  geom_point(data = select_df %>% filter(iter == 1),
             colour = "red", alpha = 0.3, size = 1) +
  theme_tn() +
  labs(x = "A", y = "B", z = "C") +
  facet_wrap(. ~ full_case)
ggsave("output/figures/random_select.pdf",
       height = 8,
       width = 8)
