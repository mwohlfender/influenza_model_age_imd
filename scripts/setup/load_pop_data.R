## load population data ##

#### SETUP ####
.args <- if (interactive()) c(
  file.path(Dir_Local, Subdir_Data, "population", "imd_2025.xlsx"),
  file.path(Dir_Local, Subdir_Data, "population", "lsoa_to_region.csv"),
  file.path(Dir_Local, Subdir_Data, "inputs", "imd_age_pop.rds") ,
  file.path(Dir_Local, Subdir_Data, "inputs", "imd_age_pop_lad.rds"),
  file.path(Dir_Local, Subdir_Data, "inputs", "imd_age_pop_gp.rds")
) else commandArgs(trailingOnly = TRUE)


## load data
lsoa_imd <- data.table(read_xlsx(.args[1], sheet = 2) %>% 
                         select(starts_with("LSOA"),starts_with("Overall"),starts_with("Index")))
colnames(lsoa_imd) <- c("lsoa21cd","lsoa21nm","imd_rank","imd_decile")
lsoa_imd[, imd_quintile := ceiling(imd_decile/2)]

## load regions of England
# https://geoportal.statistics.gov.uk/datasets/ons::lsoa-2021-to-bua-to-lad-to-region-december-2022-best-fit-lookup-in-ew-v2/about
lsoa_to_reg <- data.table(read_csv(.args[2], show_col_types = F))[, c(1,11)]
colnames(lsoa_to_reg) <- c("lsoa21cd","p_engreg")

## load age-specific population in each LSOA
# https://www.ons.gov.uk/peoplepopulationandcommunity/populationandmigration/populationestimates/datasets/lowersuperoutputareamidyearpopulationestimates
# using 2024 population, but can go back to 2019 with these datasets (vary `years` and `sheet_input`)
years <- c("20192022","20222024")[2]
sheet_input <- ifelse(years == "20192022", 8, 7)
lsoa_pop <- data.table(read_xlsx(file.path(Dir_Local, Subdir_Data, "population",paste0("sapelsoasyoa", years, ".xlsx")), 
                                 sheet = sheet_input, skip = 3))
colnames(lsoa_pop) <- c("lad21cd","lad21nm","lsoa21cd","lsoa21nm","total",paste0("F_", 0:90), paste0("M_", 0:90))

## merge:
lsoa_dat <- lsoa_pop %>% 
  left_join(lsoa_imd %>% select(lsoa21cd, imd_quintile), by = "lsoa21cd") %>% 
  left_join(lsoa_to_reg, by = "lsoa21cd") %>%
  filter(substr(lsoa21cd, 1, 1) == "E")

## aggregate:
imd_dat <- data.table(lsoa_dat %>% select(!starts_with("lsoa")))
imd_dat <- imd_dat[, lapply(.SD, sum), by = c("lad21cd", "lad21nm", "p_engreg", "imd_quintile")]

## make
imd_dat_long <- imd_dat %>% 
  select(!total) %>% 
  pivot_longer(cols = !c("lad21cd", "lad21nm", "p_engreg", "imd_quintile")) %>% 
  mutate(age = as.numeric(substr(name, 3, 4))) %>% 
  group_by(lad21cd, lad21nm, p_engreg, imd_quintile, age) %>% 
  summarise(pop = sum(value))

imd_age <- imd_dat_long %>% 
  mutate(age_grp = cut(age, c(Age_Groups_Breaks, Inf), right = F, labels = Age_Groups_Names)) %>% 
  group_by(lad21cd, lad21nm, p_engreg, imd_quintile, age_grp) %>% 
  summarise(pop = sum(pop))

## save 
write_rds(imd_age, .args[3])



imd_age_pop_lad <- imd_age %>%
  filter(lad21cd == Lad_Code) %>%
  ungroup()

imd_age_pop_lad_name <- unique(imd_age_pop_lad$lad21nm)
imd_age_pop_lad_region <- unique(imd_age_pop_lad$p_engreg)

# make sure there is data about all age groups and imd quintiles in 'imd_age_pop_lad'
# fill missing population values with 0
imd_age_pop_lad <- Table_Age_Groups_IMD_Quintiles  %>%
  left_join(imd_age_pop_lad, by = c("age_grp", "imd_quintile")) %>%
  tidyr::replace_na(list(lad21cd = Lad_Code,
                         lad21nm = imd_age_pop_lad_name,
                         p_engreg = imd_age_pop_lad_region,
                         pop = 0)) %>%
  arrange(imd_quintile, age_grp) %>%
  relocate(lad21cd, lad21nm, p_engreg, age_grp, imd_quintile, pop)

remove(imd_age_pop_lad_name, imd_age_pop_lad_region)

## save 
write_rds(imd_age_pop_lad, .args[4])




# dummy population data
n_pop <- 50000
imd_age_pop_gp <- data.table(age_grp = sample(x = Age_Groups_Names, size = n_pop, prob = c(0.02, 0.03, 0.03, 0.2, 0.35, 0.22, 0.1, 0.04, 0.01), replace = TRUE),
       imd_quintile = sample(x = 1:N_IMD, size = n_pop, prob = c(0.3, 0.1, 0.15, 0.2, 0.25), replace = TRUE)) %>%
  group_by(age_grp, imd_quintile) %>%
  summarize(pop = n())

## risk dummy data
risk_group_kids <- 0.12
risk_group_adults <- 0.08
risk_group_elderly <- 0.4

risk_group_age_vector <- fcn_assign_ages(risk_group_kids, 
                                         risk_group_adults,
                                         risk_group_elderly,
                                         Age_Groups_Names)

table_age_groups_risk <- data.table(age_grp = Age_Groups_Names,
                                    risk_age_grp = risk_group_age_vector)
table_age_groups_imd_risk <- Table_Age_Groups_IMD_Quintiles %>%
  left_join(table_age_groups_risk, by = c("age_grp")) %>%
  left_join(Table_Broad_Ages_Age_Groups, by = c("age_grp")) %>%
  left_join(Table_Broad_Ages_Risk_IMD, by = c("broad_age")) %>%
  mutate(rel_risk = 1 - (imd_quintile - 3) * var_risk_imd) %>%
  mutate(risk_proportion = risk_age_grp * rel_risk) %>%
  dplyr::select(-c("rel_risk", "risk_age_grp", "var_risk_imd")) %>%
  arrange(imd_quintile, age_grp)

imd_age_pop_gp <- imd_age_pop_gp %>%
  left_join(table_age_groups_imd_risk, by = c("age_grp", "imd_quintile")) %>%
  mutate(risk_pop = round(risk_proportion * pop))

imd_age_pop_risk_gp <- imd_age_pop_gp %>%
  mutate(pop = pop - risk_pop,
         risk_level = "low") %>%
  bind_rows(imd_age_pop_gp %>%
              mutate(pop = risk_pop,
                     risk_level = "high")) %>%
  dplyr::select(-c("risk_proportion", "risk_pop"))


## dummy vaccination data
table_broad_age_groups_risk_vacc <- data.table(broad_age = c("children", "children", "adults", "adults", "older_adults", "older_adults"),
                                               risk_level = c("low", "high", "low", "high", "low", "high"),
                                               vacc_coverage = c(0.2, 0.35, 0.6, 0.6, 0.8, 0.9))

imd_age_pop_risk_gp <- imd_age_pop_risk_gp %>%
  left_join(table_broad_age_groups_risk_vacc, by = c("broad_age", "risk_level")) %>%
  mutate(pop_vacc = round(vacc_coverage * pop))


imd_age_pop_risk_vacc_gp <- imd_age_pop_risk_gp %>%
  mutate(pop = pop - pop_vacc,
         vaccinated = FALSE) %>%
  bind_rows(imd_age_pop_risk_gp %>%
              mutate(pop = pop_vacc,
                     vaccinated = TRUE)) %>%
  dplyr::select(-c("vacc_coverage", "pop_vacc")) %>%
  relocate(c("age_grp", "broad_age", "imd_quintile", "risk_level", "vaccinated", "pop"))



## save 
write_rds(imd_age_pop_risk_vacc_gp, .args[5])


