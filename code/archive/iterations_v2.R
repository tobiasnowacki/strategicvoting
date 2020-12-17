# Next step(s): 
#   1. wrap everything in bigger functions;
#   2. run for different levels of lambda;
#   3. include random starting points;
#   4. clean up code / refactor

# 
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
load(here("output/big_list_2.RData"))
load(here("output/manyiterations2.RData"))


source(here("code/utils/refactored_functions.r"))



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

# Get na_omit and v_vecs
big_list_na_omit <- lapply(big_list, function(x) remove_nas(x))
big_list_na_omit[[39]] <- NULL
big_list_na_omit[[144]] <- NULL

# get v_vecs
sin_vote_list <- lapply(big_list_na_omit, function(x) sincere.vote.mat.from.U(x$U, rule = "AV"))
v_vec_list <- list()
# Append to big list
for(i in 1:length(sin_vote_list)){
  weights <- big_list_na_omit[[i]]$weights
  v_vec <- ballot.props.from.vote.mat.and.weights(sin_vote_list[[i]], weights)
  big_list_na_omit[[i]]$v_vec <- as.numeric(v_vec)
}


#
# ITERATIONS UNTIL CONVERGENCE
#

# set values
lambda <- 0.5

# Run Andy's function (with adapted fix) over all cases
tie_test_list <- list()
for(i in 1:160){
  print(i)
  tie_test_list[[i]] <- iterated_best_response_sequence(U = big_list_na_omit[[i]]$U %>% as.matrix, 
                                s = 85, 
                                weights = big_list_na_omit[[i]]$weights, 
                                rule = "AV", 
                                lambda = lambda,
                                until.convergence = T, 
                                max.iterations = 500, 
                                sincere.proportion = 0, 
                                candidates = c("a", "b", "c"), 
                                ballots = c("abc", "acb", "bac", "bca", "cab", "cba"), 
                                the.floor = 0, 
                                noisy = F)
}

# Create distance vector
dist_list <- list()
for(i in 1:160){
  vec_dist <- do.call(rbind, lapply(tie_test_list[[i]], function(x) x$distance.from.last))
  dist_df <- data.frame(dist = vec_dist, iter = 1:length(vec_dist), case = names(big_list_na_omit)[i])
  dist_list[[i]] <- dist_df
}

big_dist_df <- do.call(rbind, dist_list)

# Plot vvec distances
ggplot(big_dist_df %>% filter(iter < 150), aes(iter, log(dist))) +
  geom_line(aes(group = case), alpha = 0.1) +
  geom_line(data = big_dist_df %>% filter(case == "SWE_2014"), aes(group = case), colour = "red") +
  theme_sv()
# save plot
ggsave(here("output/figures_new/conv_cleaned.pdf"), device = cairo_pdf)

# get cases with no convergence
non_conv_ids <- which(sapply(dist_list, function(x) nrow(x)) == 499)
names(big_list_na_omit)[non_conv_ids]

# generate v_vec paths for cases with no convergence
vvec_list <- list()
for(i in non_conv_ids){
  vvecs <- do.call(rbind, lapply(tie_test_list[[i]], function(x) x$v.vec)) %>% as.data.frame
  vvec_list[[i]] <- vvecs %>% mutate(iter = 1:nrow(vvecs), case = names(big_list_na_omit)[i])
}

vvec_df <- do.call(rbind, vvec_list)

# Plot v_vec paths for cases with no convergence
ggtern(vvec_df, aes(abc + acb, bac + bca, cab + cba)) +
  geom_line(alpha = .2) +
  geom_point(data = vvec_df %>% filter(iter == 1), colour = "blue", size = 0.3) +
  facet_wrap(~ case) +
  theme_sv() +
  labs(A = "A", B = "B", C = "C")
ggsave(here("output/figures_new/non_conv_path.pdf"), device = cairo_pdf)

# Tabulate sincere and strategic votes for non_convergent cases
for(i in non_conv_ids){
  case_optimal <- sapply(tie_test_list[[i]], function(x) optimal.vote.from.V.mat(x$V.mat)) %>% as.vector
  case_sincere <- sapply(tie_test_list[[i]], function(x) optimal.vote.from.V.mat(x$sincere.vote.mat)) %>% as.vector

  print(table(case_sincere, case_optimal))
}

# Write function to get table of sincere and optimal votes...
get_vote_mat_df <- function(case_iter){
  opt <- optimal.vote.from.V.mat(case_iter$V.mat) %>% as.vector
  sin <- optimal.vote.from.V.mat(case_iter$sincere.vote.mat) %>% as.vector

  tab <- table(sin, factor(opt, c("abc", "acb", "bac", "bca", "cab", "cba"))) %>% as.data.frame 
  return(tab)
}

iter_wrapper <- function(case){
  list_wrap <- lapply(case, function(x) get_vote_mat_df(x))
  for(i in 1:length(list_wrap)){
    list_wrap[[i]]$iter <- i
    names(list_wrap[[i]])[2] <- "opt"
  }
  df_wrap <- do.call(rbind, list_wrap) %>% mutate(vote = interaction(sin, opt))
  return(df_wrap)
}

# Apply getting table to all cases
vote_matrix <- list()
for(i in 1:160){
  vote_matrix[[i]] <- iter_wrapper(tie_test_list[[i]])
}

for(i in non_conv_ids){
  print(i)
  ex_case <- i
  cleaned <- unique(vote_matrix[[ex_case]]$vote[vote_matrix[[ex_case]]$Freq > 0])
  cleaned_mat <- vote_matrix[[ex_case]] %>% filter(vote %in% cleaned & sin != opt)

  # Plot vote changes over iterations
  ggplot(cleaned_mat %>% filter(iter < 300), aes(iter, Freq)) +
    geom_line(aes(group = vote, colour = vote)) +
    facet_grid(~ sin) +
    theme_sv()
  # get name
  case_name <- names(big_list_na_omit)[i]
  filename <- paste0("output/figures_new/", lambda, "_", case_name, "_nonsincere.pdf")  
  ggsave(here(filename), 
         device = cairo_pdf)
}


# Examine SWE_2014 case
lapply(tie_test_list[[32]][1:10], function(x) optimal.vote.from.V.mat(x$V.mat) %>% table)

lapply(tie_test_list[[32]][1:10], function(x) optimal.vote.from.V.mat(x$V.mat) %>% table)

# strange case w/ ABC --> ACB (if left unchecked)
abc <- tie_test_list[[32]][[10]]$eu.by.ballot[126, ]

#
# RUN OWN FUNCTION TO SEE IF CODE MATCHES
#
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
      out <- many_iterations_until_convergence(big_list_na_omit[[case]], big_list_na_omit[[case]]$v_vec, lambda, s_val, 0.001, 150)
      # rcv_sum <- out[[1]] %>% mutate(case = names(big_list_na_omit)[[case]])
      # plur_sum <- out[[3]] %>% mutate(case = names(big_list_na_omit)[[case]])
      # rcv_vec <- out[[2]]
      # plur_vec <- out[[4]]
      list(out)
}

stopCluster(cl)

### COMPARISON -----
# SWE_2014 case
out <- many_iterations_until_convergence(big_list_na_omit[[32]], big_list_na_omit[[32]]$v_vec, lambda, s_val, 0.001, 150)

# check they have the same no. of iterations
tie_test_list[[32]] %>% length
nrow(out[[2]])

# check the v_vecs are the same
tie_test_list[[32]][[71]]$v.vec
out[[2]][70:72, ]
### COMPARISON ENDS ---

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

# Plot IRV convergence -- looks the same!
ggplot(conv_dist_df %>% filter(system == "IRV"), aes(iter, log(diff))) +
  geom_line(aes(group = interaction(case, system)), 
            alpha = 0.5) +
  geom_hline(yintercept = log(0.0001), colour = "red") +
  theme_sv()
ggsave(here("output/figures/convergence_irv_toby.pdf"), device = cairo_pdf)



