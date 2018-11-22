### PRELIMINARIES ###

setwd("~/Dropbox/strategic_voting")
library(here)


library(haven)
df <- read_sav("data/australia/AES/aes_2013_01259.sav")

# NSW (without National preferences) still has 736 "full preference" respondents
df <- df[df$stateab == "NSW", ]


### FUNCTIONS ###


### ANALYSIS ###

short_df <- data.frame(id = df$uniqueid,
                       alp = df$b17alp,
                       lib = df$b17lib,
                       grn = df$b17grn,
                       nat = df$b17nat)

# How many have set ALP and GRN equal; how many LIB and NAT?
sum(short_df$alp == short_df$grn, na.rm = TRUE)
sum(short_df$lib == short_df$nat, na.rm = TRUE)
#sum(is.na(short_df$lib))
#sum(is.na(short_df$nat))

# combine LIB and NAT vote
short_df$lib[is.na(short_df$lib)] <- short_df$nat[is.na(short_df$lib)]
short_df$nat[is.na(short_df$nat)] <- short_df$lib[is.na(short_df$nat)]
short_df$coal <- (short_df$lib + short_df$nat)  / 2
#sum(is.na(short_df$coal))

short_df$lib <- NULL
short_df$nat <- NULL

# Check how many cases are missing and get rid of missing data
table(complete.cases(short_df)) # Answer: Not that many: 266/3955 or 0.0672 (if taking full AES sample).
short_df <- short_df[complete.cases(short_df), ]

# Find first preference
short_df$first <- apply(short_df[, c(2:4)], 1, return_max_name)

# Proportion of voters with more than one clear first preference
sum(short_df$first == 99)
# OK, quite a significant proportion of voters who seem to have equal like-dislike scores...
# Subset on those who have a clear preference.

# Get their second and third preference. (this could be done with a more general function)
short_df$second <- apply(short_df[, c(2:4)], 1, return_sec_name)
short_df$third <- apply(short_df[, c(2:4)], 1, return_third_name)

# What respondents have assigned 3 different scores?
short_df$full <-
  !(apply(short_df[, c("first", "second", "third")], 1, function(x) 99 %in% x))
sum(short_df$full)

# Calculate differences between first, second, and third preference (VNM utilities)
short_df$fs_diff <- apply(short_df[, c(2:4)], 1, first_second_difference)
short_df$st_diff <- apply(short_df[, c(2:4)], 1, function(x) first_second_difference(x[-which(x == max(x))]))

short_df$beta <- short_df$st_diff / (short_df$fs_diff + short_df$st_diff)

# pass on a matrix of utility vectors.

order_df <- short_df[short_df$full == TRUE, c(5:7, 11)]
utility_df <- matrix(NA, nrow = nrow(order_df), ncol = 3)
for (i in 1:nrow(order_df)){
     pos1 <- order_df[i, 1]
     utility_df[i, pos1] <- 1
     pos2 <- order_df[i, 2]
     utility_df[i, pos2] <- order_df[i, 4]
     pos3 <- order_df[i, 3]
     utility_df[i, pos3] <- 0
}

# Data.frame and append ID
utility_df <- as.data.frame(utility_df)
names(utility_df) <- c("LAB", "GRN", "COAL")
utility_df$ID <- short_df$id[short_df$full == TRUE]
utility_df$beta <- short_df$beta[short_df$full == TRUE]

# Get voter type
utility_df <- utility_df[, c("GRN", "COAL", "LAB", "ID", "beta")]

utility_df$type[utility_df$GRN > utility_df$COAL & utility_df$COAL > utility_df$LAB] <- "GRN.COAL.LAB"
utility_df$type[utility_df$GRN > utility_df$LAB & utility_df$LAB > utility_df$COAL] <- "GRN.LAB.COAL"
utility_df$type[utility_df$COAL > utility_df$GRN & utility_df$GRN > utility_df$LAB] <- "COAL.GRN.LAB"
utility_df$type[utility_df$COAL > utility_df$LAB & utility_df$LAB > utility_df$GRN] <- "COAL.LAB.GRN"
utility_df$type[utility_df$LAB > utility_df$GRN & utility_df$GRN > utility_df$COAL] <- "LAB.GRN.COAL"
utility_df$type[utility_df$LAB > utility_df$COAL & utility_df$COAL > utility_df$GRN] <- "LAB.COAL.GRN"
utility_df$type2 <- factor(utility_df$type, levels = c("GRN.COAL.LAB", "GRN.LAB.COAL", "COAL.GRN.LAB", "COAL.LAB.GRN",
                                                       "LAB.GRN.COAL", "LAB.COAL.GRN"))
utility_df$type2 <- factor(utility_df$type2, levels(utility_df$type2)[c(1, 3, 5, 2, 4, 6)])

# Plot frequency of beta by voter type

library(plyr)
util_count <- ddply(.data=utility_df, 
                 .(type2), 
                 summarize, 
                 n=paste("n =", length(GRN)))

ggplot(data = utility_df) +
     geom_histogram(aes(x = beta), bins = 10) +
     facet_wrap(~ type2) +
     geom_text(data = util_count, aes(x = 0.8, y = 65, label = n), inherit.aes = F) +
     theme_bw() +
     scale_y_continuous(limits = c(0, 70), expand = c(0, 0)) +
     theme(panel.grid = element_blank()) +
     labs(y = "Respondents", 
          x = expression(beta))
ggsave("code/figs/aes_respondents_beta_distribution.pdf", width = 16, height = 10, units = "cm")

# Export for simulations
write.csv(utility_df, "data/AES_utility.csv")
