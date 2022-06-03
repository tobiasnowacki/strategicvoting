# calculate average distance between top and bottom preference(s)

# Load packages and data
library(tidyverse)
library(gtools)
library(fields)
library(ggsci)
library(ggtern)
library(devtools)
library(modelsummary)
library(kableExtra)
source_url("https://raw.githubusercontent.com/tobiasnowacki/RTemplates/master/plottheme.R")

# Load functions
source("code/utils/new_sv_iter.R")
# Load data
load("output/big_list_2.RData")
vap <- read.csv("data/case_vap.csv", sep = "") # voting age pop.
cat("Data imported. \n")
source("code/prep_cses.R") # data prep
systems <- read.csv("data/systems.csv")

return_dat <- function(X) {
  diff <- mean(apply(X$U, 1, function(y) max(y) - min(y)))
  top <- mean(apply(X$U, 1, function(y) max(y)))
  bottom <- mean(apply(X$U, 1, function(y) min(y)))

  prefmat <- sincere_pref_mat_from_U(X$U)
  pref_string <- apply(prefmat, 1, function(x) paste0(x[1], x[2], x[3]))
  polarization <- min(
    c(
      mean(pref_string %in% c("123", "132")),
      mean(pref_string %in% c("231", "321")),
      mean(pref_string %in% c("213", "312"))
    )
  )
  return(
    data.frame(
      Difference = diff,
      Top = top,
      Bottom = bottom,
      Unidimensionality = 1 - polarization
    )
  )
}

easy_summary <- map_dfr(big_list_na_omit, ~ return_dat(.x), .id = "case")
easy_summary_mg <- easy_summary |>
  pivot_longer(
    !case,
    names_to = "stat",
    values_to = "value"
  ) |>
  left_join(systems, by = c("case")) |>
  filter(system %in% c("pr", "fptp"))

easy_summary_mg$stat <- fct_relevel(
  easy_summary_mg$stat, c("Top", "Bottom", "Difference", "Unidimensionality")
)

# Run KS between PR and FPTP
measures <- c("Difference", "Unidimensionality")

for (m in measures) {

  prdist <- easy_summary_mg$value[
    easy_summary_mg$system == "pr" & easy_summary_mg$stat == m
  ]

  fptpdist <- easy_summary_mg$value[
    easy_summary_mg$system == "fptp" & easy_summary_mg$stat == m
  ]

  print(paste0("Test for Equality of Distributions of ", m))

  print(ks.test(prdist, fptpdist))

}


# RCV is a little bit of an outlier - but this could just because of the much smaller sample...
ggplot(easy_summary_mg, aes(x = value)) +
  geom_density(aes(fill = system), alpha = 0.5) +
  facet_wrap(~ stat, scales = "free") +
  scale_fill_nejm() +
  labs(x = "Value", y = "Density") +
  theme_tn()

ggsave(
  "output/figures/profile_comp.pdf",
  height = 3,
  width = 8,
  device = cairo_pdf
)

# Let's try a regression
regdf <- easy_summary |>
  left_join(systems, by = "case") |>
  filter(system %in% c("pr", "fptp", "rcv"))

reglist <- list()

for (v in c("Top", "Bottom", "Difference", "Unidimensionality")){
  fml <- as.formula(paste0(v, " ~ system"))
  reglist[[v]] <- lm(fml, data = regdf)
}

cmap <- c(
   '(Intercept)' ='Plurality baseline',
   "systempr" = "PR",
   "systemrcv" = "RCV"
)

modelsummary(
  reglist,
  output = "latex",
  booktabs = TRUE,
  coef_map = cmap,
  gof_omit = "",
  caption = "\\label{tab:prefreg} Statistics describing distribution of preferences.",
  align = "ldddd"
  ) |>
  kable_styling(
    latex_options = "hold_position"
  ) |>
  footnote(
    threeparttable = TRUE,
    fixed_small_size = FALSE,
    footnote_as_chunk = TRUE,
    general_title = "",
    general = c("Each unit of observation is a CSES case. 'Top' describes the average utility of the top-ranked candidate; 'Bottom' describes the average utility of the lowest-ranked candidate. 'Difference' measures the average difference in utility between the highest- and lowest-ranked candidate. Polarization measures the smallest share of voters who put a candidate last")
  ) |>
  save_kable(file = "output/tables/pref_distance.tex")
