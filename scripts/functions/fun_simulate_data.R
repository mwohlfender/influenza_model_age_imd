# ## CREATE DUMMY DATA (INFECTIONS) ##
# 
# #### SETUP ####
# .args <- if (interactive()) c(
#   file.path(Dir_Local, Subdir_Data, "inputs", "imd_age_pop_lad.rds"),
#   file.path(Dir_Local, Subdir_Data, "inputs", "contact_matrix.rds"),
#   file.path(Dir_Local, Subdir_Data, "dummy_data", "known_parameters.rds"),
#   file.path(Dir_Local, Subdir_Data, "dummy_data", "unknown_parameters.rds"),
#   file.path(Dir_Local, Subdir_Data, "dummy_data", "dummy_infections.rds")
# ) else commandArgs(trailingOnly = TRUE)
# 
# 
# # input_subtype <- "AH1N1"
# input_subtype <- "AH3N2"
# input_start_of_season <- 2023
# input_susceptibility_broad_age_groups <- c(0.564, 0.416, 0.396)
# 
# 
# set.seed(60)
# 
# # LOAD DATA ####
# 
# ## read in population data #### 
# ## ENSURE ALL POP ORDERS ARE IMD THEN AGE
# imd_age_pop_lad <- readRDS(.args[1]) %>%
#   arrange(imd_quintile, age_grp)
# 
# ## read in contact matrix  ####
# contact_matrix <- readRDS(.args[2])
# 
# ## check the right number of rows (no extra defining variables)
# if (nrow(contact_matrix) != nrow(contact_matrix %>% select(p_imd_q, c_imd_q, p_age_group, c_age_group) %>% unique())) {
#   warning('Number of rows in contact matrix wrong')
# }
# 
# contact_matrix$p_age_group <- factor(contact_matrix$p_age_group, levels = Age_Groups_Names)
# contact_matrix$c_age_group <- factor(contact_matrix$c_age_group, levels = Age_Groups_Names)
# 
# ## ENSURE ALL POP ORDERS ARE IMD THEN AGE
# contact_matrix <- contact_matrix %>% 
#   arrange(p_imd_q,
#           c_imd_q,
#           p_age_group,
#           c_age_group)
# 
# ## into matrix
# cm <- contact_matrix %>%
#   ungroup() %>% 
#   mutate(p_var = paste0(p_imd_q, '_', p_age_group),
#          c_var = paste0(c_imd_q, '_', c_age_group)) %>% 
#   select(p_var, c_var, n) %>% 
#   pivot_wider(names_from = c_var, values_from = n) %>% 
#   select(!p_var) %>%
#   as.matrix()
# 
# ## EXPAND CONTACT MATRIX 
# 
# ## make matrix per capita
# per_cap_matrix_45 <- t(t(cm)/imd_age_pop_lad$pop)
# 
# ## scale up for both risk groups 
# ## (assuming per capita contacts are independent of risk level)
# pc_cm <- expand_contact_matrix(per_cap_matrix_45)
# 
# ng <- nrow(per_cap_matrix_45)
# ndim <- nrow(pc_cm)
# if (!all.equal(2*ng, ndim)) {warning('dimensions not adding up')}
# 
# 
# 
# 
# ## read in known parameters ####
# known_pars <- readRDS(.args[3])
# 
# vaccinated_pop_seasonal <- known_pars$vaccinated_data %>% 
#   filter(start_of_season == input_start_of_season,
#          subtype == input_subtype) %>% 
#   arrange(start_of_season, desc(risk_level), imd_quintile, age_grp)
# 
# ## should be ordered by IMD then age
# if (vaccinated_pop_seasonal$imd_quintile[2] != 1) {warning('vaccinated_data in wrong order')}
# if (vaccinated_pop_seasonal$age_grp[2] != Age_Groups_Names[2]) {warning('vaccinated_data in wrong order')}
# 
# 
# 
# ## UNKNOWN PARAMETERS
# unknown_pars <- readRDS(.args[4])
# 
# pars <- unknown_pars$epid_parameters %>% 
#   filter(year == input_start_of_season,
#          subtype == input_subtype)
# 
# susceptibility_vec <- fcn_assign_ages(pars$children_susceptibility,
#                                       pars$adults_susceptibility,
#                                       pars$older_adults_susceptibility,
#                                       Age_Groups_Names)
# 
# 
# 
# 
# time_series <- run_model(pop = vaccinated_pop_seasonal$pop,
#                          I0 = pars$init_infected,
#                          vacc_cov = vaccinated_pop_seasonal$vaccinated_population,
#                          ve_inf = vaccinated_pop_seasonal$VE_INF,
#                          cm = pc_cm,
#                          # trans = pars$transmissibility,
#                          # trans = 0.188665,
#                          trans = 0.1889241,
#                          susc = susceptibility_vec,
#                          lat_per = known_pars$epid_periods[1],
#                          inf_per = known_pars$epid_periods[2])
# 
# # check not rising at end of time
# n_tot_subpops <- nrow(time_series[t == max(t), ])
# n_subpops <- sum(time_series[t == max(t) - 1, ]$infections < time_series[t == max(t), ]$infections)
# if (n_subpops > 0) {
#   warning(paste0('Epidemic still rising in ', n_subpops, '/', n_tot_subpops, ' subpopulations'))
# }
# 
# # add population data
# time_series <- time_series[known_pars$pop_gp_age_imd_risk_vacc_lad %>%
#                              filter(start_of_season == input_start_of_season) %>%
#                              dplyr::select(-c("start_of_season"), "pop" = "pop_lad"),
#                            on = c('age_grp', 'imd_quintile', 'risk_level', 'vaccinated')]
# 
# 
# n_subpops_min <- sum(time_series[t == max(t), ]$infections > 5)
# 
# if (n_subpops_min > 0) {
#   warning(paste0('Epidemic still above 5 in ', n_subpops_min, '/', n_tot_subpops, ' subpopulations'))
# }
# 
# time_series[, start_date := pars$start_date]
# time_series[, subtype := input_subtype]
# 
# time_series <- time_series[t %in% 0:365] ## take data from the start of each day
# 
# 
# 
# 
# 
# 
# 
# 
# dummy_infections <- readRDS(.args[5])
# dummy_infections_2023 <- dummy_infections[[1]] %>% bind_rows(dummy_infections[[2]])
# 
# # new_time_series_1 <- time_series
# # new_time_series_2 <- time_series
# # 
# # new_time_series <- new_time_series_1 %>% bind_rows(new_time_series_2)
# 
# 
# dummy_infections_2023_summarized <- dummy_infections_2023 %>%
#   filter(t > 0) %>%
#   mutate(age_grp = case_when(grepl('18|30|50', age_grp) ~ '18-64', T ~ age_grp)) %>%
#   mutate(age_grp = factor(age_grp, levels = c('0-4','5-11','12-17','18-64','65-74','75-84','85+'))) %>%
#   group_by(t, age_grp, imd_quintile) %>% 
#   summarise(inf = sum(infections), pop = sum(pop)) 
# 
# new_time_series_summarized <- new_time_series %>%
#   filter(t > 0) %>%
#   mutate(age_grp = case_when(grepl('18|30|50', age_grp) ~ '18-64', T ~ age_grp)) %>%
#   mutate(age_grp = factor(age_grp, levels = c('0-4','5-11','12-17','18-64','65-74','75-84','85+'))) %>%
#   group_by(t, age_grp, imd_quintile) %>% 
#   summarise(inf = sum(infections), pop = sum(pop)) 
# 
# 
# ggplot() +
#   geom_line(data = dummy_infections_2023_summarized, aes(t, 100000*inf/pop, col = as.character(imd_quintile)), lwd = 0.8) +
#   geom_line(data = new_time_series_summarized, aes(t, 100000*inf/pop, col = as.character(imd_quintile)), lwd = 0.8) +
#   scale_color_manual(values = imd_quintile_colors) +
#   facet_wrap(age_grp ~ ., scales = 'free') +
#   theme_bw() +
#   labs(x = 'Day of epidemic', y = 'Infections per 100,000', col = 'IMD')
# 
# 
# ggplot() +
#   geom_line(data = dummy_infections_2023_summarized, aes(t, 100000*inf/pop, col = as.character(imd_quintile)), lwd = 0.8) +
#   # geom_line(data = new_time_series_summarized, aes(t, 100000*inf/pop, col = as.character(imd_quintile)), lwd = 0.8) +
#   scale_color_manual(values = imd_quintile_colors) +
#   facet_wrap(age_grp ~ ., scales = 'free') +
#   theme_bw() +
#   labs(x = 'Day of epidemic', y = 'Infections per 100,000', col = 'IMD')
# 
# 
# ggplot() +
#   # geom_line(data = dummy_infections_2023_summarized, aes(t, 100000*inf/pop, col = as.character(imd_quintile)), lwd = 0.8) +
#   geom_line(data = new_time_series_summarized, aes(t, 100000*inf/pop, col = as.character(imd_quintile)), lwd = 0.8) +
#   scale_color_manual(values = imd_quintile_colors) +
#   facet_wrap(age_grp ~ ., scales = 'free') +
#   theme_bw() +
#   labs(x = 'Day of epidemic', y = 'Infections per 100,000', col = 'IMD')
# 
# 
# 
