# Building new dataset for analysis
library(tidyverse)
library(rio)
library(gtools)
library(questionr)
library(gridExtra)
# library(devtools)

# Uncomment the following when running on server

    # # Load data -------------------------------
    # load("output/files/1/85_winners_tbl.Rdata") # add win proportions
    # load("output/files/1/85_winners.Rdata") # add expected utilities
    # load("output/files/1/85_sum.Rdata") # add SV analyses
    # load("output/files/1/85_vvec.Rdata") # add vvecs
    # load("output/big_list_2.RData") ## add original data
    # case_df <- import("data/systems.csv") # add original electoral systems
    # vap <- read.csv("data/case_vap.csv", sep = "") # voting age pop.

    # # Load scripts and templates ------------
    # source("code/utils/sv_theme_template.R")
    # source("code/utils/new_sv_iter.R")
    # source("code/prep_cses.R") # data prep

    # # Now that everything is loaded, save intermediate file
    # save.image(file = "output/files/1/85_workbench.Rdata")

# On client, load the data as follows:
load("output/files/1/85_workbench.Rdata")

# For SV statistics, convert iteration to numeric
sum_df <- sum_df %>% mutate(iter = as.numeric(iter))

# Function to find condorcet winner from 6-item vvec
return_condorcet <- function(vvec){
    a_b <- sum(vvec[c(1, 2, 5)]) - sum(vvec[c(3, 4, 6)])
    a_c <- sum(vvec[c(1, 2, 3)]) - sum(vvec[c(4, 5, 6)])
    b_c <- sum(vvec[c(1, 3, 4)]) - sum(vvec[c(2, 5, 6)])

    if (a_b > 0 & a_c > 0) {
       return("A")
    } else if (a_c < 0 & b_c < 0) {
       return("C")
    } else if (a_b < 0 & b_c > 0) {
       return("B")
    } else {
       return("NA")
    }
}

# Apply to vvecs
c_winner <- map_chr(
    big_list_na_omit, ~
    return_condorcet(.x$v_vec)
    ) %>% 
    enframe() %>%
    rename(case = name, cwinner = value)

# Merge winners, EU and summary stats together
ddf <- winners_df %>%
    left_join(win_df) %>%
    left_join(sum_df) %>%
    left_join(case_df) %>%
    left_join(c_winner) %>%
    rename(
        winpr_a = V1,
        winpr_b = V2,
        winpr_c = V3,
        avg_exp_eu = value
    ) %>%
    mutate(
        system_group = case_when(
            system %in% c("pr", "mixed") ~ "pr",
            system %in% c("fptp", "tworound", "parallel") ~ "plur",
            system %in% c("rcv") ~ "rcv", 
        TRUE ~ "other"
        ),
        cwin_prob = case_when(
            cwinner == "A" ~ winpr_a,
            cwinner == "B" ~ winpr_b,
            cwinner == "C" ~ winpr_c
        )
    )


# Compare differences in AEU
diff_df_extra <- ddf %>%
    group_by(case) %>%
    summarise(
        diff_p1_rcv1 = avg_exp_eu[system == "plurality" & iter == 1] - avg_exp_eu[system == "rcv" & iter == 1],
        diff_p60_rcv1 = avg_exp_eu[system == "plurality" & iter == 60] - avg_exp_eu[system == "rcv" & iter == 1],
        diff_p60_rcv60 = avg_exp_eu[system == "plurality" & iter == 60] - avg_exp_eu[system == "rcv" & iter == 60],
        diff_p60_p1 = avg_exp_eu[system == "plurality" & iter == 60] - avg_exp_eu[system == "plurality" & iter == 1],
        diff_rcv60_rcv1 = avg_exp_eu[system == "rcv" & iter == 60] - avg_exp_eu[system == "rcv" & iter == 1]
    ) %>%
    left_join(
        case_weight_tbl # Merge in case weights
    )

# Compute weighted means of AEU
avg_p1_rcv1 <- wtd.mean(diff_df_extra$diff_p1_rcv1, diff_df_extra$case_weight)
avg_p60_rcv60 <- wtd.mean(diff_df_extra$diff_p60_rcv60, diff_df_extra$case_weight)
avg_p60_p1 <- wtd.mean(diff_df_extra$diff_p60_p1, diff_df_extra$case_weight)
avg_rcv60_rcv1 <- wtd.mean(diff_df_extra$diff_rcv60_rcv1, diff_df_extra$case_weight)


# Difference in AEU Plot --------------------------------

# Function to standardise the two plot panels
# Return as list to append to existing ggplot(s)
joint_plot_elements <- function() {
    list(
        geom_vline(xintercept = 0, lty = "dotted", colour = "#333333"),
        geom_hline(yintercept = 0, lty = "dotted", colour = "#333333"),
        geom_point(colour = "#727272", alpha = 0.3),
        lims(x = c(-0.65, 0.65), y = c(-0.65, 0.35)),
        theme_tn()
    )
}

# Left panel
px1 <- ggplot(diff_df_extra, aes(diff_p1_rcv1, diff_p60_rcv60)) +
    geom_hline(yintercept = avg_p60_rcv60, colour = "blue") +
    geom_vline(xintercept = avg_p1_rcv1, colour = "blue") +
    labs(
        x = expression(paste(Delta, " EU (P 1st, RCV 1st)")),
        y = expression(paste(Delta, " EU (P 60th, RCV 60st)"))
    ) +
    joint_plot_elements() +
    ggtitle("Across Systems") 

# Right panel
px2 <- ggplot(diff_df_extra, aes(diff_p60_p1, diff_rcv60_rcv1)) +
    geom_hline(yintercept = avg_rcv60_rcv1, colour = "blue") +
    geom_vline(xintercept = avg_p60_p1, colour = "blue") +
    labs(
        x = expression(paste(Delta, " EU (P 60st, P 1st)")),
        y = expression(paste(Delta, " EU (RCV 60th, RCV 1th)"))
    ) +
    ggtitle("Across Iterations") +
    joint_plot_elements()

# Combine and save
px_comb <- gridExtra::grid.arrange(px1, px2, ncol = 2)
ggsave(
    px_comb, 
    file = "output/figures/1/85_aeu.pdf", 
    width = 6, height = 3)

# Additional Analyses ----------------------------------------

# Only retain 1st and 60th iteration
ddf_pruned <- ddf %>% filter(iter == 1 | iter == 60)

# Density of Condorcet Winner Plot
cwin_density <- ggplot(ddf_pruned, aes(x = cwin_prob)) +
    geom_density(aes(fill = as.factor(iter)), alpha = 0.5) +
    facet_wrap( ~ system) +
    labs(x = "Pr(Condorcet Winner Win)", y = "Density", fill = "Iteration") +
    theme_tn()
ggsave(cwin_density, file = "output/figures/1/85_cwin_density.pdf", 
    width = 6, height = 3)

# Boxplot of C candidate winning
last_win_box <- ggplot(ddf_pruned, aes(x = winpr_c, y = as.factor(iter))) +
    geom_boxplot(aes(fill = as.factor(iter)), 
        alpha = 0.5, width = 0.5) +
    facet_wrap( ~ system) +
    labs(x = "Pr(C Candidate Wins)", y = "Iteration", fill = "Iteration") +
    theme_tn() +
    scale_y_discrete(labels = c("1st", "60th"))
ggsave(last_win_box, file = "output/figures/1/85_lastwin_box.pdf", 
    width = 6, height = 3)
