# List of packages required for this script
requiredPackages <- c("here", 
                      "ggplot2",
                      "ggtern",
                      "dplyr", 
                      "purrr", 
                      "tidyr", 
                      "lmtest", 
                      "sandwich", 
                      "plm", 
                      "extrafont", 
                      "RColorBrewer", 
                      "boot", 
                      "svMisc", 
                      "ggtern", 
                      "reldist", 
                      "gridExtra", 
                      "ggpubr", 
                      "doParallel", 
                      "foreach",
                      "questionr")

# Function to either install or call packages
ipak <- function(pkg){
        new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
        if (length(new.pkg))
                install.packages(new.pkg, 
                                 dependencies = TRUE,
                                 repos='http://cran.us.r-project.org')
        suppressMessages(sapply(pkg, require, 
                                character.only = TRUE))
}

# Load packages.
ipak(requiredPackages)

cat("Packages loaded. \n")

# Load functions
source(here("code/utils/functions.r"))
source(here("code/utils/av_pivotal_probs_analytical_general_v2.r"))
source(here("code/utils/plurality_pivotal_probabilities_analytical.r"))
source(here("code/utils/general_iteration_simulation_approach.r"))
source(here("code/utils/sv.r"))
source(here("code/utils/refactored_functions.r"))

cat("Functions loaded. \n")


# Load data
load(here("output/big_list_2.RData"))
vap <- read.csv(here("data/case_vap.csv"), sep = "") # voting age pop.

cat("Data imported. \n")

# Load fonts
# font_import()

# Load ggplot theme
source(here("code/utils/sv_theme_template.R"))

# Colourblind palette for plots etc.
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

cat("Ready for analysis. \n")
