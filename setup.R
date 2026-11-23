
# load libraries ####
suppressMessages(require(BayesianTools))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(ggplot2))
# suppressPackageStartupMessages(require(tidyverse))
suppressPackageStartupMessages(require(data.table))
suppressMessages(require(patchwork))
suppressMessages(require(parallel))
suppressPackageStartupMessages(require(readr))
suppressMessages(require(readODS))
suppressPackageStartupMessages(require(readxl))
suppressPackageStartupMessages(require(stringr))
suppressPackageStartupMessages(require(tidyr))
suppressMessages(require(scales))
suppressMessages(require(viridis))
options(dplyr.summarise.inform = FALSE)

library(RcppArmadillo)


# define directories ####
Dir_Local <- "C:/Users/MartinWohlfender/Documents/documents_msw/projects/winter_pressures/local/influenza_model_age_imd"

SubDir_Code <- "scripts"
Subdir_Data <- "data"
Subdir_Output <- "output"



# source R functions ####
list_functions_R <- dir(path = file.path(SubDir_Code, "functions"), recursive = TRUE)
invisible(lapply(X = list_functions_R, FUN = function(x) {source(file.path(SubDir_Code, "functions", x), echo = FALSE)}))
remove(list_functions_R)



# define global variables ####

## number of years of data
Years <- 2023:2025 # 2023-24 to 2025-26

## define parameters ####
Age_Groups_Breaks <- c(0, 5, 12, 18, 30, 50, 65, 75, 85)
Age_Groups_Names <- paste0(Age_Groups_Breaks, "-", lead(Age_Groups_Breaks) - 1)
Age_Groups_Names[length(Age_Groups_Names)] <- "85+"
Age_Groups_Names <- factor(x = Age_Groups_Names, levels = Age_Groups_Names)

# numer of age groups
N_AGE <- length(Age_Groups_Names)

# number of IMD groups (quintiles)
N_IMD <- 5

Table_Age_Groups_IMD_Quintiles <- data.table(age_grp = rep(x = Age_Groups_Names, times = N_IMD),
                                             imd_quintile = rep(x = 1:N_IMD, each = length(Age_Groups_Names)))

Broad_Ages_Names <- factor(x = c("children", "adults", "older_adults"), levels = c("children", "adults", "older_adults"))

Table_Broad_Ages_Age_Groups <- data.table(age_grp = Age_Groups_Names, 
                                          broad_age = fcn_assign_ages("children", "adults", "older_adults", Age_Groups_Names))

# variation of risk by broad age and IMD
Table_Broad_Ages_Risk_IMD <- data.table(broad_age = unique(Table_Broad_Ages_Age_Groups$broad_age),
                                        var_risk_imd = c(0.15, 0.25, 0.08))




## define LAD for which we model an Influenza season (code based on 2021 LADs)
Lad_Code <- "E08000032"





# source R scripts ####
source(file.path("scripts", "setup", "colors.R"))




# source(file.path("scripts", "setup", "load_pop_data.R"))
# source(file.path("scripts", "setup", "load_contact_data.R"))
# source(file.path("scripts", "setup", "subtype_setup.R"))
# source(file.path("scripts", "dummy_data", "produce_known_parameters.R"))
# source(file.path("scripts", "dummy_data", "produce_unknown_parameters.R"))


# source(file.path("scripts", "seir_model.R"))
# source(file.path("scripts", "dummy_mcmc", "mcmc_functions.R"))


