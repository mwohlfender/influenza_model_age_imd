## KNOWN PARAMETERS ##

#### SETUP ####
.args <- if (interactive()) c(
  file.path(Dir_Local, Subdir_Data, "inputs", "imd_age_pop_lad.rds"),
  file.path(Dir_Local, Subdir_Data, "inputs", "imd_age_pop_gp.rds"),
  file.path(Dir_Local, Subdir_Data, "inputs", "subtype_years.rds"),
  file.path(Dir_Local, Subdir_Data, "population", "risk_group_population_data.rds"),
  file.path(Dir_Local, Subdir_Data, "dummy_data", "known_parameters.rds")
) else commandArgs(trailingOnly = TRUE)

if(!dir.exists(file.path(Subdir_Data, "dummy_data"))){dir.create(file.path("data", "dummy_data"))}


## read in population data
imd_age_pop_lad <- read_rds(.args[1])
pop_gp_age_imd_risk_vacc <- read_rds(.args[2])


#### MAKE PARAMETERS ####

subtype_years <- read_rds(.args[3]) %>% select(season, subtype) %>% unique()
subtype_vec <- unique(subtype_years$subtype)

#### PROPORTIONS IN RISK GROUPS ####

## risk group proportions

USING_TRUE_RISK_DAT <- T

if (USING_TRUE_RISK_DAT) {
  
  ## loading in true risk group data
  risk_group_dat <- readRDS(.args[4]) %>% 
    mutate(.by = interval,
           proportion = count/sum(count)) %>% 
    filter(risk_group %like% 'high') %>%
    arrange(lower)
  
  ## matching up age groups via midpoint
  age_labels_midpoints_list <- (str_split(gsub('[+]','-100',Age_Groups_Names), pattern = '-'))
  age_labels_midpoints <- c(); corresponding_intervals <- c(); risk_group_age_vector <- c()
  for (interval in 1:N_AGE) {
    age_labels_midpoints[interval] <- mean(as.numeric(age_labels_midpoints_list[[interval]]))
    corresponding_intervals[interval] <- sum(age_labels_midpoints[interval] > risk_group_dat$lower)
    risk_group_age_vector[interval] <- risk_group_dat$proportion[corresponding_intervals[interval]]
  }
  
} else {
  
  ## dummy data
  risk_group_kids <- 0.12
  risk_group_adults <- 0.08
  risk_group_elderly <- 0.4
  
  risk_group_age_vector <- fcn_assign_ages(risk_group_kids, 
                                           risk_group_adults,
                                           risk_group_elderly,
                                           Age_Groups_Names)
  
}

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

risk_group_pop <- table_age_groups_imd_risk %>%
  left_join(imd_age_pop_lad, by = c("age_grp", "imd_quintile")) %>%
  mutate(risk_population = round(risk_proportion * pop)) %>%
  relocate("age_grp", "broad_age", "imd_quintile",
           "lad21cd", "lad21nm", "p_engreg",
           "risk_proportion", "risk_population", "pop") %>%
  arrange(imd_quintile, age_grp)



## PLOT
risk_group_pop %>% 
  ggplot() +
  geom_line(aes(x = age_grp, y = risk_population, colour = as.factor(imd_quintile), group = imd_quintile),
            lwd = 0.8) + 
  scale_color_manual(values = imd_quintile_colors) +
  labs(x = "", y = "Population in risk group", col = "IMD") +
  theme_lg()


risk_group_pop %>% 
  ggplot() +
  geom_line(aes(x = age_grp, y = risk_proportion, colour = as.factor(imd_quintile), group = imd_quintile),
            lwd = 0.8) + 
  scale_color_manual(values = imd_quintile_colors) +
  labs(x = "", y = "Proportion in risk group", col = "IMD") +
  theme_lg()


#### VACCINATION COVERAGE #### 

# https://www.gov.uk/government/statistics/seasonal-influenza-vaccine-uptake-in-children-of-school-age-monthly-data-2025-to-2026
vaccination_coverage_kids <- 0.55
# https://www.gov.uk/government/statistics/seasonal-influenza-vaccine-uptake-in-gp-patients-monthly-data-2025-to-2026
vaccination_coverage_elderly <- 0.75
vaccination_coverage_risk <- 0.4

vaccination_coverage_age_vector <- fcn_assign_ages(
  vaccination_coverage_kids, 
  vaccination_coverage_risk,
  vaccination_coverage_elderly,
  Age_Groups_Names)

## add some variation by IMD
vaccination_imd_variation <- seq(0.9, 1.1, by = 0.05)

vaccinated_pop_1 <- data.table(vcav = rep(vaccination_coverage_age_vector, N_IMD), 
                               vciv = rep(vaccination_imd_variation, each = N_AGE)) %>% 
  mutate(age_grp = rep(Age_Groups_Names, N_IMD),
         imd_quintile = rep(1:N_IMD, each = N_AGE),
         vaccinated_proportion = vcav*vciv) %>% 
  left_join(risk_group_pop, by = c('age_grp','imd_quintile')) 

# split by risk groups
vaccinated_pop <- rbind(vaccinated_pop_1 %>% mutate(risk_level = 'low', pop = pop - risk_population),
                        vaccinated_pop_1 %>% mutate(risk_level = 'high', pop = risk_population)) %>% 
  select(!c(risk_population, vcav, vciv)) %>% 
  mutate(nrow = n()) %>% 
  ## duplicate across years for now
  slice(rep(row_number(), 3)) %>% 
  mutate(start_of_season = rep(Years, each = nrow[1])) %>% select(!nrow) %>% 
  mutate(vaccinated_population = case_when(age_grp %in% c('18-25', '26-34', '35-49', '50-69') & risk_level == 'low' ~ 0,
                                           T ~ round(vaccinated_proportion*pop))) %>% 
  ## ENSURE ALL POP ORDERS ARE IMD THEN AGE
  arrange(start_of_season, desc(risk_level), imd_quintile, age_grp)


# create plot
vaccinated_pop %>%
  group_by(age_grp, imd_quintile, start_of_season) %>% 
  summarise(vaccinated_population = sum(vaccinated_population)) %>% 
  ggplot() +
  geom_line(aes(x = age_grp, y = vaccinated_population, colour = imd_quintile, group = imd_quintile)) +
  facet_wrap(~start_of_season) +
  theme_lg()

# create plot
vaccinated_pop %>% 
  ggplot() +
  geom_line(aes(age_grp, vaccinated_population, colour = imd_quintile, lty = risk_level, group = interaction(imd_quintile, risk_level))) +
  facet_wrap(~start_of_season) +
  theme_lg()

# create plot
vacc_plot <- vaccinated_pop %>% 
  ggplot() +
  geom_line(aes(x = age_grp, group = interaction(as.factor(imd_quintile), risk_level), 
                col = as.factor(imd_quintile), y = vaccinated_population/pop), lwd = 1) +
  geom_point(aes(x = age_grp, group = interaction(as.factor(imd_quintile), risk_level), 
                 y = vaccinated_population/pop, shape = risk_level), 
             col='white', size = 3) +
  geom_point(aes(x = age_grp, group = interaction(as.factor(imd_quintile), risk_level), 
                 col = as.factor(imd_quintile), y = vaccinated_population/pop, 
                 shape = risk_level), stroke=1.5, size = 3) +
  scale_color_manual(values = imd_quintile_colors) + ylim(c(0,NA)) +
  scale_shape_manual(values = c(1, 2)) +
  theme_lg() + labs(x = 'Age group', col = 'IMD quintile', 
                    y = 'Simulated vaccination coverage',
                    shape = 'Risk level') +
  theme(text = element_text(size = 14)); vacc_plot

# create plot
vaccinated_pop %>% 
  mutate(key_group = case_when(
    age_grp %in% c('0-4','5-11','12-17') ~ 'Children',
    !(age_grp %in% c('0-4','5-11','12-17','70-79','80+')) & risk_level == 'high' ~ 'Risk group (18-69)',
    age_grp %in% c('70-79','80+') ~ 'Older adults',
    T~ NA)) %>% 
  group_by(key_group, imd_quintile) %>% 
  summarise(vaccinated_population = sum(vaccinated_population),
            pop = sum(pop)) %>% filter(!is.na(key_group)) %>% 
  ggplot() + 
  geom_bar(aes(x = key_group, group = imd_quintile, 
               fill = as.factor(imd_quintile), y = 100*vaccinated_population/pop),
           stat = 'identity', position = 'dodge') +
  scale_fill_manual(values = imd_quintile_colors) + 
  theme_bw() + labs(x = '', fill = 'IMD quintile', 
                    y = 'Simulated vaccine uptake (%)') +
  theme(text = element_text(size = 14))

# create plot
risk_plot <- vaccinated_pop %>% 
  ggplot() + 
  geom_line(aes(x = age_grp, group = as.factor(imd_quintile), 
                col = as.factor(imd_quintile), y = risk_proportion), lwd = 1) +
  geom_point(aes(x = age_grp, group = as.factor(imd_quintile), y = risk_proportion), 
             col='white', size = 3) +
  geom_point(aes(x = age_grp, group = as.factor(imd_quintile), 
                 col = as.factor(imd_quintile), y = risk_proportion), 
             shape = 1, stroke=2, size = 3) +
  scale_color_manual(values = imd_quintile_colors) + ylim(c(0,NA)) +
  theme_lg() + labs(x = 'Age group', col = 'IMD quintile', 
                    y = 'Simulated percentage in clinical risk group') +
  theme(text = element_text(size = 14),
        legend.position = 'none')

vacc_plot + risk_plot + plot_layout(nrow = 2, guides = 'collect')

ggsave(file.path('output','figures','dummy_infections','dummy_vacc_risk_props.png'),
       width = 10, height = 10)

#### VACCINE EFFICACY ####
## (age-dependent, annual, subtype-specific)

#### AGAINST INFECTION ####

## FOR NOW USING MADE UP DATA
## TODO UPDATE WHEN DATA AVAILABLE

vaccination_efficacy_infection_AH1N1 <- cross_join(data.table(age_grp = Age_Groups_Names,  
                                                              VE_INF = fcn_assign_ages(0.45, 
                                                                                       0.35,
                                                                                       0.3,
                                                                                       Age_Groups_Names)),
                                                   data.table(start_of_season = Years, subtype = "AH1N1")) 

vaccination_efficacy_infection_AH3N2 <- cross_join(data.table(age_grp = Age_Groups_Names,  
                                                              VE_INF = fcn_assign_ages(0.45, 
                                                                                       0.35,
                                                                                       0.3,
                                                                                       Age_Groups_Names)),
                                                   data.table(start_of_season = Years, subtype = "AH3N2")) 

vaccination_efficacy_infection_B <- cross_join(data.table(age_grp = Age_Groups_Names,  
                                                          VE_INF = fcn_assign_ages(0.45, 
                                                                                   0.35,
                                                                                   0.3,
                                                                                   Age_Groups_Names)),
                                               data.table(start_of_season = Years, subtype = "B")) 

vaccination_efficacy_infection <- rbind(vaccination_efficacy_infection_AH1N1,
                                        vaccination_efficacy_infection_AH3N2,
                                        vaccination_efficacy_infection_B)

#### AGAINST HOSPITALISATION ####
# 2023/24: https://onlinelibrary.wiley.com/doi/epdf/10.1111/irv.70194
# 2024/25: https://www.gov.uk/government/statistics/influenza-in-the-uk-annual-epidemiological-report-winter-2024-to-2025/influenza-in-the-uk-annual-epidemiological-report-winter-2024-to-2025#secondary-care-surveillance

# from: https://doi.org/10.1111/irv.70194

VE_2_17_AH1N1_2022_2023 <- 0.61
VE_2_17_AH1N1_2023_2024 <- 0.47
VE_2_17_AH3N2_2022_2023 <- 0.63
VE_2_17_AH3N2_2023_2024 <- 0.67
VE_2_17_B_2022_2023 <- 0.87
VE_2_17_B_2023_2024 <- 0.84
VE_18_64_AH1N1_2022_2023 <- 0.34
VE_18_64_AH1N1_2023_2024 <- 0.48
VE_18_64_AH3N2_2022_2023 <- 0.27
VE_18_64_AH3N2_2023_2024 <- 0.36
VE_18_64_B_2022_2023 <- 0.53
VE_18_64_B_2023_2024 <- 0.70
VE_65_AH1N1_2022_2023 <- 0.24
VE_65_AH1N1_2023_2024 <- 0.37
VE_65_AH3N2_2022_2023 <- 0.25
VE_65_AH3N2_2023_2024 <- 0.20
VE_65_B_2022_2023 <- 0.28
VE_65_B_2023_2024 <- 0.39

# from: https://www.gov.uk/government/statistics/influenza-in-the-uk-annual-epidemiological-report-winter-2024-to-2025/influenza-in-the-uk-annual-epidemiological-report-winter-2024-to-2025#vaccination
# (Figure 56)

ukhsa_dat_2024_25 <- read_ods(file.path(Dir_Local, Subdir_Data, "ukhsa","annual_influenza_2024_2025.ods"),
                              sheet = 59, skip = 3)
colnames(ukhsa_dat_2024_25) <- c('age_group','flu_subtype','NA1','NA2','NA3','NA4','VE','VE_CIL','VE_CIU')
ukhsa_dat_2024_25 <- ukhsa_dat_2024_25 %>% mutate(VE = VE/100)

VE_2_17_AH1N1_2024_2025 <- (ukhsa_dat_2024_25 %>% filter(age_group %like% '2 to 17', flu_subtype %like% 'H1N1'))$VE
VE_2_17_AH3N2_2024_2025 <- (ukhsa_dat_2024_25 %>% filter(age_group %like% '2 to 17', flu_subtype %like% 'H3N2'))$VE
VE_2_17_B_2024_2025 <- (ukhsa_dat_2024_25 %>% filter(age_group %like% '2 to 17', flu_subtype %like% 'Influenza B'))$VE
VE_18_64_AH1N1_2024_2025 <- (ukhsa_dat_2024_25 %>% filter(age_group %like% '18 to 64', flu_subtype %like% 'H1N1'))$VE
VE_18_64_AH3N2_2024_2025 <- (ukhsa_dat_2024_25 %>% filter(age_group %like% '18 to 64', flu_subtype %like% 'H3N2'))$VE
VE_18_64_B_2024_2025 <- (ukhsa_dat_2024_25 %>% filter(age_group %like% '18 to 64', flu_subtype %like% 'Influenza B'))$VE
VE_65_AH1N1_2024_2025 <- (ukhsa_dat_2024_25 %>% filter(age_group %like% '65', flu_subtype %like% 'H1N1'))$VE
VE_65_AH3N2_2024_2025 <- (ukhsa_dat_2024_25 %>% filter(age_group %like% '65', flu_subtype %like% 'H3N2'))$VE
VE_65_B_2024_2025 <- (ukhsa_dat_2024_25 %>% filter(age_group %like% '65', flu_subtype %like% 'Influenza B'))$VE

# from: https://www.gov.uk/government/statistics/influenza-in-the-uk-annual-epidemiological-report-winter-2025-to-2026/influenza-in-the-uk-annual-epidemiological-report-winter-2025-to-2026#vaccination
# (Figure 56)

ukhsa_dat_2025_26 <- read_ods(file.path(Dir_Local, Subdir_Data, "ukhsa","annual_influenza_2025_2026.ods"),
                              sheet = 59, skip = 3)
colnames(ukhsa_dat_2025_26) <- c('age_group','flu_subtype','NA1','NA2','NA3','NA4','VE')
ukhsa_dat_2025_26 <- ukhsa_dat_2025_26 %>% 
  mutate(VE = as.numeric(sub("%.*", "", VE))/100,
         flu_subtype = gsub(substr(ukhsa_dat_2025_26$flu_subtype[2], 4, 4), '', flu_subtype))

VE_2_17_AH1N1_2025_2026 <- (ukhsa_dat_2025_26 %>% filter(age_group %like% '2 to 17', flu_subtype %like% 'H1N1'))$VE
VE_2_17_AH3N2_2025_2026 <- (ukhsa_dat_2025_26 %>% filter(age_group %like% '2 to 17', flu_subtype %like% 'H3N2'))$VE
VE_2_17_B_2025_2026 <- (ukhsa_dat_2025_26 %>% filter(age_group %like% '2 to 17', flu_subtype %like% 'InfluenzaB'))$VE
VE_18_64_AH1N1_2025_2026 <- (ukhsa_dat_2025_26 %>% filter(age_group %like% '18 to 64', flu_subtype %like% 'H1N1'))$VE
VE_18_64_AH3N2_2025_2026 <- (ukhsa_dat_2025_26 %>% filter(age_group %like% '18 to 64', flu_subtype %like% 'H3N2'))$VE
VE_18_64_B_2025_2026 <- (ukhsa_dat_2025_26 %>% filter(age_group %like% '18 to 64', flu_subtype %like% 'InfluenzaB'))$VE
VE_65_AH1N1_2025_2026 <- (ukhsa_dat_2025_26 %>% filter(age_group %like% '65', flu_subtype %like% 'H1N1'))$VE
VE_65_AH3N2_2025_2026 <- (ukhsa_dat_2025_26 %>% filter(age_group %like% '65', flu_subtype %like% 'H3N2'))$VE
VE_65_B_2025_2026 <- (ukhsa_dat_2025_26 %>% filter(age_group %like% '65', flu_subtype %like% 'InfluenzaB'))$VE

## merge data
vaccination_efficacy_hospitalisation <- CJ(subtype = subtype_vec,
                                           start_of_season = Years,
                                           age_grp = Age_Groups_Names,
                                           VE_HOSP = 0)

for (subtype_i in subtype_vec) {
  for (season_i in Years) {
    vaccination_efficacy_hospitalisation[
      subtype == subtype_i & start_of_season == season_i, 
      VE_HOSP := fcn_assign_ages(
        get(paste0('VE_2_17_', subtype_i, '_', season_i, '_', season_i + 1)), 
        get(paste0('VE_18_64_', subtype_i, '_', season_i, '_', season_i + 1)),
        get(paste0('VE_65_', subtype_i, '_', season_i, '_', season_i + 1)),
        Age_Groups_Names)]
  }
}

## plot
vaccination_efficacy_hospitalisation %>% 
  ggplot() + 
  geom_line(aes(x = age_grp, y = VE_HOSP, group = interaction(start_of_season, subtype),
                col = as.factor(start_of_season)), lwd = 1) +
  geom_point(aes(x = age_grp, y = VE_HOSP, group = interaction(start_of_season, subtype)), 
             col='white', size = 3) +
  geom_point(aes(x = age_grp, y = VE_HOSP, group = interaction(start_of_season, subtype),
                 col = as.factor(start_of_season), shape = subtype), 
             stroke=1.5, size = 3) +
  scale_shape_manual(values = c(1, 2, 4)) +
  scale_color_manual(values = season_colors) +
  theme_lg() + labs(x = 'Age group', col = 'Season start',
                    y = 'VE against hospitalisation') +
  facet_grid(subtype ~.) + 
  theme(text = element_text(size = 14))

vaccination_efficacy_hospitalisation %>% 
  ggplot() + 
  geom_line(aes(x = age_grp, y = VE_HOSP, group = interaction(start_of_season, subtype),
                col = subtype), lwd = 1) +
  geom_point(aes(x = age_grp, y = VE_HOSP, group = interaction(start_of_season, subtype)), 
             col='white', size = 3) +
  geom_point(aes(x = age_grp, y = VE_HOSP, group = interaction(start_of_season, subtype),
                 col = subtype, shape = as.factor(start_of_season)), 
             stroke=1.5, size = 3) +
  ylim(c(NA, 1)) + 
  scale_shape_manual(values = c(1, 2, 4)) +
  scale_color_manual(values = subtype_colors) +
  theme_lg() + labs(x = 'Age group', col = 'Subtype', shape = 'Season',
                    y = 'VE against hospitalisation') +
  facet_grid(start_of_season ~.) +
  theme(text = element_text(size = 14))

ggsave(file.path('output','figures','dummy_infections','ukhsa_ve_hospitalisation.png'),
       width = 10, height = 8)

## join together with vaccination coverage data

vaccinated_data <- rbind(vaccinated_pop %>% mutate(subtype = 'AH1N1'),
                         vaccinated_pop %>% mutate(subtype = 'AH3N2'),
                         vaccinated_pop %>% mutate(subtype = 'B')) %>% 
  left_join(vaccination_efficacy_infection, by = c('age_grp','start_of_season', 'subtype')) %>% 
  left_join(vaccination_efficacy_hospitalisation, by = c('age_grp','start_of_season', 'subtype')) %>% 
  mutate(effectively_vaccinated_population = round(VE_INF*vaccinated_population))

#### EPI PERIODS ####

# from Baguelin et al., 2013
epid_periods <- c(0.8, 1.8) # latent and infectious periods

#### DELAYS ####
## (in weeks)
primary_care_delay <- 1
secondary_care_delay <- 1

#### OPENSAFELY COVERAGE ####
# proportion of LAD population living in GP catchment area
# stratified by age group, IMD quintile, risk group, and vaccination status

proportion_observed <- CJ(age_grp = Age_Groups_Names,
                          imd_quintile = 1:N_IMD,
                          risk_level = c('low', 'high'),
                          vaccinated = c(FALSE, TRUE),
                          year = Years)

pop_gp_age_imd_risk_lad_unvacc <- vaccinated_data %>%
  group_by(age_grp, imd_quintile, risk_level, start_of_season) %>%
  summarize(pop = unique(pop-vaccinated_population)) %>%
  ungroup() %>%
  mutate(vaccinated = FALSE)

pop_gp_age_imd_risk_lad_vacc  <- vaccinated_data %>%
  group_by(age_grp, imd_quintile, risk_level, start_of_season) %>%
  summarize(pop = unique(vaccinated_population)) %>%
  ungroup() %>%
  mutate(vaccinated = TRUE)

pop_gp_age_imd_risk_vacc_lad <- pop_gp_age_imd_risk_lad_unvacc %>%
  bind_rows(pop_gp_age_imd_risk_lad_vacc) %>%
  dplyr::rename(c("pop_lad" = "pop")) %>%
  arrange(start_of_season, imd_quintile, age_grp) %>%
  dplyr::relocate("age_grp", "imd_quintile", "risk_level", "start_of_season", "vaccinated", "pop_lad")

proportion_observed <- proportion_observed %>%
  left_join(pop_gp_age_imd_risk_vacc_lad,
              by = c("age_grp" = "age_grp", "imd_quintile" = "imd_quintile",
                     "risk_level" = "risk_level", "vaccinated" = "vaccinated", "year" = "start_of_season")) %>%
  left_join(pop_gp_age_imd_risk_vacc %>%
              dplyr::select(c("age_grp", "imd_quintile", "risk_level", "vaccinated", "pop_gp" = "pop")),
            by = c("age_grp" = "age_grp", "imd_quintile" = "imd_quintile",
                   "risk_level" = "risk_level", "vaccinated" = "vaccinated")) %>%
  mutate(prop = pop_gp / pop_lad)



#### MAKE INTO LIST ####
known_pars <- list(years = Years,
                   pop_gp_age_imd_risk_vacc_lad = pop_gp_age_imd_risk_vacc_lad,
                   vaccinated_data = vaccinated_data,
                   epid_periods = epid_periods,
                   primary_care_delay = primary_care_delay,
                   secondary_care_delay = secondary_care_delay,
                   proportion_observed = proportion_observed)

#### SAVE KNOWN PARAMETERS ####

saveRDS(known_pars, .args[5])


