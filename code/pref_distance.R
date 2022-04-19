#-----
# Script to analyse similarity of preferences between cases.
# (April 2022)
#-----

# Dependencies ----------------

library(tidyverse)
library(gtools)
library(fields)
library(ggsci)
library(ggtern)
library(devtools)
source_url('https://raw.githubusercontent.com/tobiasnowacki/RTemplates/master/plottheme.R')

# Load functions
source("code/utils/new_sv_iter.R")
# Load data
load("output/big_list_2.RData")
vap <- read.csv("data/case_vap.csv", sep = "") # voting age pop.
cat("Data imported. \n")
source("code/prep_cses.R")  # data prep
systems <- read.csv("data/systems.csv")

# Functions -------------------

#' Function to return distance dataframe
#' data should include 'case' and 'type' along with utils
return_distdf <- function(
  data, type, name_vec, case_supplied = FALSE){
  abc_type <- data[data[["type"]] == type, ]

  if (case_supplied){
    name_vec <- abc_type[["case"]]
  }

  abc_dist <- rdist(
    abc_type[, c("A", "B", "C")],
    abc_type[, c("A", "B", "C")])

  rownames(abc_dist) <- name_vec
  colnames(abc_dist) <- name_vec

  ind <- which(upper.tri(abc_dist, diag = TRUE), arr.ind = TRUE)
  nn <- dimnames(abc_dist)
  dist_mat <- data.frame(row = nn[[1]][ind[, 1]],
            col = nn[[2]][ind[, 2]],
            val = abc_dist[ind])

  dist_mat_ext <- dist_mat |>
    left_join(systems, by = c("row" = "case")) |>
    mutate(system_left = system) |>
    dplyr::select(row, col, val, system_left) |>
    left_join(systems, by = c("col" = "case")) |>
    mutate(system_right = system) |>
    dplyr::select(row, col, val, system_left, system_right) |>
    mutate(
      system_comp = paste0(system_left, ":", system_right),
      type = type
    ) |>
    filter(row != col)

  return(dist_mat_ext)
}

# Distance between random draws ------------

# Draw 500-voter sample from every case
util_list <- list()
for (i in 1:160) {
  n <- nrow(big_list_na_omit[[i]]$U)
  indices <- sample(1:n, 500, replace = TRUE)
  util_list[[i]] <- big_list_na_omit[[i]]$U[indices, ]
}

# similarity between kth case and all others
k <- 49
dist_mats <- list()
for (i in 1:160){
  dist_mats[[i]] <- rdist(util_list[[k]], util_list[[i]])
}
names(dist_mats) <- names(big_list_na_omit)[1:160]

plot_df <- map_dfr(
  dist_mats,
  ~ as.vector(.x),
  .id = "case"
) %>%
  pivot_longer(cols = everything()) %>%
  left_join(systems, by = c("name" = "case")) %>%
  mutate(is_self = name == names(big_list_na_omit[k]))

ggplot(plot_df, aes(x = value)) +
  geom_line(aes(group = name, colour = interaction(system, is_self)), stat = "density", alpha = 0.2) +
  scale_colour_nejm()

# Distance between voters of the same type --------------

# Put all utils into one DF
all_utils <- map_dfr(
  big_list_na_omit,
  ~ .x$U,
  .id = "case"
)

# Identify voter types and label utils data
prefmat <- sincere_pref_mat_from_U(all_utils[, 2:4])
pref_string <- apply(prefmat, 1, function(x) paste0(x[1], x[2], x[3]))
pref_string[pref_string == "123"] <- "CBA"
pref_string[pref_string == "132"] <- "BCA"
pref_string[pref_string == "213"] <- "CAB"
pref_string[pref_string == "312"] <- "ACB"
pref_string[pref_string == "321"] <- "ABC"
pref_string[pref_string == "231"] <- "BAC"
all_utils$type <- pref_string

all_utils$type <- fct_relevel(all_utils$type, c(
  "BCA", "BAC", "ABC", "ACB", "CAB", "CBA"))

# Compute the average profile by case and type
summary_df <- all_utils |>
  group_by(case, type) |>
  summarise(A = mean(A), B = mean(B), C = mean(C)) |>
  left_join(systems) |>
  filter(system == "fptp" | system == "rcv") |>
  arrange(case, type)

# Plot the average preference profiles (by system)
ggtern(summary_df, aes(B, A, C)) +
  geom_path(aes(colour = system, group = case), alpha = 0.3) +
  geom_point(aes(colour = system), alpha = 0.3)

# Calculate distance between cases (by same voter type)
out_map <- map_dfr(
    c("BCA", "BAC", "ABC", "ACB", "CAB", "CBA"),
    ~ return_distdf(summary_df, type = .x, case_supplied = TRUE)
  )

# The distance between average voter profiles isn't that different between FPTP and PR!
ggplot(out_map, aes(x = val)) +
  geom_density(aes(colour = system_comp)) +
  facet_wrap(~ type) +
  theme_tn()
