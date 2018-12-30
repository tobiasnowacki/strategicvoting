##################################################
## Project: Strategic Voting in RCV
## Script purpose: Checking theoretical predictions
## Date: 19/12/2018
## Author:
##################################################

### 
### Dependencies
###

library(here)
library(gtools)
library(ggtern)
library(stargazer)
library(tidyr)
library(rworldmap)

source(here("utils/functions.r"))
source(here("utils", "av_pivotal_probs_analytical_general_v2.r"))
sim_appr2 <- here("utils", "general_iteration_simulation_approach.r")
source(sim_appr2)
sv_file <-  here("utils/sv.r")
source(sv_file)
source(here("utils", "plot_results_on_ternary_mirror.r"))

load(here("..", "output", "cses_big_list_2.RData"))

### 
### Analysis
###

# Create list with v_vecs from CSES utility dfs:
big_list_na_omit <- lapply(big_list, function(x) remove_nas(x))
sin_vote_list <- lapply(big_list_na_omit, function(x) sincere.vote.mat.from.U(x$U, rule = "AV"))
v_vec_list <- list()
for(i in 1:length(sin_vote_list)){
  weights <- big_list_na_omit[[i]]$weights
  v_vec <- ballot.props.from.vote.mat.and.weights(sin_vote_list[[i]], weights)
  big_list_na_omit[[i]]$v_vec <- as.numeric(v_vec)
}

# Drop NA cases
names(big_list)[[39]]
names(big_list)[[145]]
big_list_na_omit[[39]] <- NULL
big_list_na_omit[[144]] <- NULL

# Set level of uncertainty
s <- 85

# For each case:
# 1. get v.vec 
# 1. get "transition matrix"
# 2. get pivotal probabilities

results <- list()
for(i in 1:length(big_list_na_omit)){
	print(i)
	this_list <- big_list_na_omit[[i]]
	pprobs <- av.pivotal.event.probs.general(
		c(this_list$v_vec, 0, 0, 0),
		rep(s, 4)
		)
	class <- classify.vec(this_list$v_vec)
	sv_obj <- convert_andy_to_sv_item_two(
		this_list$U, 
		this_list$weights, 
		s, 
		this_list$v_vec)
	trans_matrix <- vote_matrix_weighted(
		sv_obj,
		type = "rcv", 
		weights = big_list_na_omit[[i]]$weights
		)
	results[[i]] <- list(this_list$v_vec, pprobs, sv_obj, trans_matrix, class)
}

# Mapping cases
cntry <- sapply(names(big_list_na_omit), function(x) substr(x, 1, 3))
cntry_list <- as.data.frame(table(cntry))
map_df <- joinCountryData2Map(cntry_list, joinCode = "ISO3", nameJoinColumn = "cntry")
pdf(here("../output/figures/case_map.pdf"))
mapCountryData(map_df, nameColumnToPlot = "Freq", catMethod = "categorical", mapTitle = "Number of cases by country")
dev.off()

# Classification

# Create matrix of v_vecs
pprobs_df <- as.data.frame(t(sapply(results, `[[`, 2)))
pprobs_df$type <- sapply(results, `[[`, 5)
pprobs_df$case <- c(names(big_list_na_omit))

prob_df_long <- gather(pprobs_df, key = "event", value = "prob", 1:12)
prob_df_long$event <- as.factor(prob_df_long$event)
prob_df_long$prob <- as.numeric(prob_df_long$prob)

v_vec_df <- t(sapply(results, `[[`, 1))

get_strat_vote_vec <- function(vote_mat){
	vote_mat <- t(apply(vote_mat, 1, function(x) x / sum(x)))
	return(c(vote_mat[1, 3], vote_mat[1, 5], vote_mat[2, 3], vote_mat[2, 5], 
		vote_mat[3, 1], vote_mat[3, 6], vote_mat[4, 1], vote_mat[4, 6], 
		vote_mat[5, 2], vote_mat[5, 4], vote_mat[6, 2], vote_mat[6, 4]))
}

v_vec_df <- as.data.frame(cbind(v_vec_df, t(sapply(results, function(x) get_strat_vote_vec(x[[4]])))))
names(v_vec_df) <- c("ABC", "ACB", "BAC", "BCA", "CAB", "CBA", 
	"ABCBAC", "ABCCAB", "ACBBAC", "ACBCAB", 
	"BACABC", "BACCBA", "BCAABC", "BCACBA",
	"CABACB", "CABBCA", "CBAACB", "CBABCA")

v_vec_df$type <- sapply(results, `[[`, 5)
v_vec_df$case <- c(names(big_list_na_omit))

df_long <- gather(v_vec_df[, 7:20], key = "svtype", value = "prop", ABCBAC:CBABCA)

# Summary statistics

mean(unlist(lapply(big_list_na_omit, function(x) nrow(x$U))))
sd(unlist(lapply(big_list_na_omit, function(x) nrow(x$U))))
min(unlist(lapply(big_list_na_omit, function(x) nrow(x$U))))
max(unlist(lapply(big_list_na_omit, function(x) nrow(x$U))))
table(v_vec_df$type)

# First preference distribution (handled better in replication_cses.r)
v_vec_df$def_party <- NA
v_vec_df$def_party[v_vec_df$type %in% c("SP(A)", "DM(A)")] <- "A"
v_vec_df$def_party[v_vec_df$type %in% c("SP(B)", "DM(B)")] <- "B"
v_vec_df$def_party[v_vec_df$type %in% c("SP(C)", "DM(C)")] <- "C"
v_vec_df$def_class <- "Other"
v_vec_df$def_class[v_vec_df$type %in% c("SP(A)", "SP(B)", "SP(C)")] <- "single-peaked"
v_vec_df$def_class[v_vec_df$type %in% c("DM(A)", "DM(B)", "DM(C)")] <- "divided-majority"
v_vec_df$def_class[v_vec_df$type %in% c("N")] <- "neutral"

ggtern(v_vec_df, aes(ABC + ACB, BAC + BCA, CAB + CBA)) +
	geom_point(aes(colour = def_party), alpha = 0.25) +
	facet_wrap(~ def_class) +
	theme_bw() +
	labs(x = "A", y = "B", z = "C", colour = "Defining Party") +
	theme(legend.position = "bottom")
ggsave(here("../output/figures/cses_fp.pdf"), width = 6, height = 3)


# Check "dominance" of predicted event.
pprobs_df$dominant <- NA
pprobs_df$dominant[v_vec_df$type == "SP(A)"] <- pprobs_df$"AC.AB"[v_vec_df$type == "SP(A)"]
pprobs_df$dominant[v_vec_df$type == "DM(A)"] <- pprobs_df$"BC.BC"[v_vec_df$type == "DM(A)"]
pprobs_df$dominant[v_vec_df$type %in% c("SP(C)", "DM(B)")] <- pprobs_df$"BC.AC"[v_vec_df$type %in% c("SP(C)", "DM(B)")]
pprobs_df$dominant[v_vec_df$type %in% c("SP(B)", "DM(C)")] <- pprobs_df$"BC.BA"[v_vec_df$type %in% c("SP(B)", "DM(C)")]

pprobs_df$ratio <- apply(pprobs_df[, c(4:12, 15)], 1, function(x) {
		vec <- as.numeric(x)
		max(vec[vec != x$dominant]) / x$dominant
	})

ggplot(pprobs_df) +
	geom_histogram(aes(x = ratio), bins = 300)
ggsave(here("../output/figures/prediction/pivotal_ratios.pdf"))

ggplot(pprobs_df) +
	geom_histogram(aes(x = ratio), bins = 50) +
	xlim(0, 1.5)
ggsave(here("../output/figures/prediction/pivotal_ratios_truncated.pdf"))

# Checking massive outlier (ratio > 300)
which.max(pprobs_df$ratio)
pprobs_df$type[119]
plot.av.result(c(unlist(v_vec_df[119, 1:6])), add.fp.result = T)
pprobs_df[119, ]

# Plotting strategic votes -- SINGLE PEAKED CASES

# Plotting pivotal probabilities
ggplot(prob_df_long[prob_df_long$type == "SP(A)" & !(prob_df_long$event %in% c("AB", "AC", "BC")), ]) +
	geom_boxplot(aes(x = event, y = prob)) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	labs(x = "Event", y = "Pr(Event)", caption = "A+")
ggsave(here("../output/figures/prediction/pprob_sp_a.pdf"), height = 4, width = 4)

ggplot(prob_df_long[prob_df_long$type == "SP(B)" & !(prob_df_long$event %in% c("AB", "AC", "BC")), ]) +
	geom_boxplot(aes(x = event, y = prob)) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	labs(x = "Event", y = "Pr(Event)", caption = "B+")
ggsave(here("../output/figures/prediction/pprob_sp_b.pdf"), height = 4, width = 4)


ggplot(prob_df_long[prob_df_long$type == "SP(C)" & !(prob_df_long$event %in% c("AB", "AC", "BC")), ]) +
	geom_boxplot(aes(x = event, y = prob)) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	labs(x = "Event", y = "Pr(Event)", caption = "C+")
ggsave(here("../output/figures/prediction/pprob_sp_c.pdf"), height = 4, width = 4)


# Plotting strategic vote incidences
ggplot(df_long[df_long$type == "SP(A)", ]) +
	geom_boxplot(aes(x = svtype, y = prop, fill = svtype)) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	labs(x = "Incentive Type", y = "Proportion of all voters",
		caption = "A+") +
	scale_fill_manual(values = c("white", "white", "white", "grey",
								 "white", "light grey", "white", "light grey",
								 "light grey", "white", "white", "white"),
						guide = FALSE)
ggsave(here("../output/figures/prediction/svinc_sp_a.pdf"), height = 4, width = 4)

ggplot(df_long[df_long$type == "SP(B)", ]) +
	geom_boxplot(aes(x = svtype, y = prop, fill = svtype)) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	labs(x = "Incentive Type", y = "Proportion of all voters",
		caption = "B+") +
	scale_fill_manual(values = c("white", "dark grey", "white", "dark grey",
								 "white", "white", "white", "white",
								 "white", "white", "white", "orange"),
						guide = FALSE)
ggsave(here("../output/figures/prediction/svinc_sp_b.pdf"), height = 4, width = 4)

ggplot(df_long[df_long$type == "SP(C)", ]) +
	geom_boxplot(aes(x = svtype, y = prop, fill = svtype)) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	labs(x = "Incentive Type", y = "Proportion of all voters",
		caption = "C+") +
	scale_fill_manual(values = c("dark grey", "white", "dark grey", "white",
								 "white", "white", "white", "orange",
								 "white", "white", "white", "white"),
						guide = FALSE)
ggsave(here("../output/figures/prediction/svinc_sp_c.pdf"), , height = 4, width = 4)

# Plotting strategic votes -- DIVIDED MAJORITY CASES
# Pivotal probabilities

ggplot(prob_df_long[prob_df_long$type == "DM(A)" & !(prob_df_long$event %in% c("AB", "AC", "BC")), ]) +
	geom_boxplot(aes(x = event, y = prob)) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	labs(x = "Event", y = "Pr(Event)", caption = "A-")
ggsave(here("../output/figures/prediction/pprob_dm_a.pdf"), height = 4, width = 4)

ggplot(prob_df_long[prob_df_long$type == "DM(B)" & !(prob_df_long$event %in% c("AB", "AC", "BC")), ]) +
	geom_boxplot(aes(x = event, y = prob)) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	labs(x = "Event", y = "Pr(Event)", caption = "B-")
ggsave(here("../output/figures/prediction/pprob_dm_B.pdf"), height = 4, width = 4)

ggplot(prob_df_long[prob_df_long$type == "DM(C)" & !(prob_df_long$event %in% c("AB", "AC", "BC")), ]) +
	geom_boxplot(aes(x = event, y = prob)) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	labs(x = "Event", y = "Pr(Event)", caption = "C-")
ggsave(here("../output/figures/prediction/pprob_dm_c.pdf"), height = 4, width = 4)

# Strategic incentives

ggplot(df_long[df_long$type == "DM(A)", ]) +
	geom_boxplot(aes(x = svtype, y = prop, fill = svtype)) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	labs(x = "Incentive Type", y = "Proportion of all voters",
		caption = "A-") +
	scale_fill_manual(values = c("dark grey", "white", "white", "dark grey",
								 "white", "white", "white", "white",
								 "white", "white", "white", "white"),
						guide = FALSE)
ggsave(here("../output/figures/prediction/svinc_dm_a.pdf"), height = 4, width = 4)

ggplot(df_long[df_long$type == "DM(B)", ]) +
	geom_boxplot(aes(x = svtype, y = prop, fill = svtype)) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	labs(x = "Incentive Type", y = "Proportion of all voters",
		caption = "B-") +
	scale_fill_manual(values = c("grey", "white", "grey", "white",
								 "white", "white", "white", "grey",
								 "white", "white", "white", "white"),
						guide = FALSE)
ggsave(here("../output/figures/prediction/svinc_dm_b.pdf"), height = 4, width = 4)

ggplot(df_long[df_long$type == "DM(C)", ]) +
	geom_boxplot(aes(x = svtype, y = prop, fill = svtype)) +
	theme_bw() +
	theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
	labs(x = "Incentive Type", y = "Proportion of all voters",
		caption = "C-") +
	scale_fill_manual(values = c("white", "grey", "white", "grey",
								 "white", "white", "white", "white",
								 "white", "white", "white", "grey"),
						guide = FALSE)
ggsave(here("../output/figures/prediction/svinc_dm_c.pdf"), height = 4, width = 4)







# Other

# Incidence of main strategic vote in C+, given relative frequency of "conflicting" event (AC.BC)
ggplot(v_vec_df[v_vec_df$type == "SP(C)", ]) +
	geom_point(aes(x = as.numeric(pprobs_df[v_vec_df$type == "SP(C)", "AC.BC"]) / as.numeric(pprobs_df[v_vec_df$type == "SP(C)", "BC.AC"]), 
		y = ACBCAB)) +
	theme_bw() +
	labs(x = "Pr(AC.BC) / Pr(BC.AC)", y = "Pr(ACBCAB)")
ggsave("../output/figures/prediction/sp_c_odd.pdf", height = 4, width = 4)

# Distribution of sincere FP in single-peaked cases.
ggtern(v_vec_df[(v_vec_df$type %in% c("SP(A)", "SP(B)", "SP(C)")), c(1:6, 19)], aes(ABC + ACB, BAC + BCA, CAB + CBA)) +
	geom_point(aes(colour = type)) +
	labs(x = "A", y = "B", z = "C") +
	theme_bw()
ggsave("../output/figures/prediction/sp_cases_tern.pdf", height = 4, width = 4)


ggplot(v_vec_df[v_vec_df$type == "SP(B)", ]) +
	geom_density(aes(CBABCA)) + 
	geom_density(aes(ABCCAB), lty = "dotted") + 
	geom_density(aes(ACBCAB), lty = "dashed")


# The higher the probability of these pivotal events, the lower should the incidence of these strategic votes be.
plot(pprobs_df[v_vec_df$type == "SP(B)", "BC.AC"], 
	v_vec_df$ABCCAB[v_vec_df$type == "SP(B)"])
plot(pprobs_df[v_vec_df$type == "SP(B)", "BC.AC"], 
	v_vec_df$ACBCAB[v_vec_df$type == "SP(B)"])

ggplot(v_vec_df[v_vec_df$type == "SP(C)", ]) +
	geom_density(aes(BCACBA)) + 
	geom_density(aes(ABCBAC), lty = "dotted") + 
	geom_density(aes(ACBBAC), lty = "dashed")	

ggplot(v_vec_df[v_vec_df$type == "DM(B)", ]) +
	geom_density(aes(BCACBA)) + 
	geom_density(aes(ABCBAC), lty = "dotted") + 
	geom_density(aes(ACBBAC), lty = "dashed")

ggplot(v_vec_df[v_vec_df$type == "SP(A)", ]) +
	geom_density(aes(BACCBA)) + 
	geom_density(aes(BCACBA), lty = "dotted") + 
	geom_density(aes(CABACB), lty = "dashed")

ggplot(v_vec_df[v_vec_df$type == "DM(A)", ]) +
	geom_density(aes(ABCBAC)) + 
	geom_density(aes(ACBCAB), lty = "dotted") + 
	geom_density(aes(CABACB), lty = "dashed")

hist(v_vec_df[v_vec_df$type == "SP(B)", ]$BAC / (v_vec_df[v_vec_df$type == "SP(B)", ]$BAC + v_vec_df[v_vec_df$type == "SP(B)", ]$BCA))