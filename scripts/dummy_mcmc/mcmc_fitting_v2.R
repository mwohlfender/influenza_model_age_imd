## MCMC FITTING - UNKNOWN REPORTING RATES ##

#### SETUP ####
.args <- #if (interactive()) c(
  c(file.path(Dir_Local, Subdir_Data, "inputs", "imd_age_pop_lad.rds"),
    file.path(Dir_Local, Subdir_Data, "inputs", "imd_age_pop_gp.rds"),
    file.path(Dir_Local, Subdir_Data, "inputs", "contact_matrix.rds"),
    file.path(Dir_Local, Subdir_Data, "inputs", "subtype_years.rds"),
    file.path(Dir_Local, Subdir_Data, "dummy_data", "dummy_surveillance.rds"),
    file.path(Dir_Local, Subdir_Data, "dummy_data", "known_parameters.rds"),
    file.path(Dir_Local, Subdir_Data, "dummy_data", "unknown_parameters.rds")
  ) #else commandArgs(trailingOnly = TRUE)



set.seed(60)

#### LOAD DATA ####

## read in population data
## ENSURE ALL POP ORDERS ARE IMD THEN AGE
imd_age_pop_lad <- readRDS(.args[1]) %>%
  arrange(imd_quintile, age_grp)

imd_age_pop_gp <- readRDS(.args[2]) %>%
  arrange(imd_quintile, age_grp)

## read in contact matrix
contact_matrix_1000 <- readRDS(.args[3])

## aggregate
contact_matrix <- contact_matrix_1000 %>% 
  group_by(p_imd_q, c_imd_q, p_age_group, c_age_group) %>% 
  summarise(n = mean(n))

contact_matrix$p_age_group <- factor(contact_matrix$p_age_group,
                                     levels = Age_Groups_Names)
contact_matrix$c_age_group <- factor(contact_matrix$c_age_group,
                                     levels = Age_Groups_Names)

## ENSURE ALL POP ORDERS ARE IMD THEN AGE
contact_matrix <- contact_matrix %>% 
  arrange(p_imd_q, c_imd_q,
          p_age_group, c_age_group,)

## into matrix
cm <- contact_matrix %>% ungroup() %>% 
  mutate(p_var = paste0(p_imd_q, '_', p_age_group),
         c_var = paste0(c_imd_q, '_', c_age_group)) %>% 
  select(p_var, c_var,n) %>% 
  pivot_wider(names_from = c_var, values_from = n) %>% 
  select(!p_var) %>% as.matrix()

## make matrix per capita
per_cap_matrix_45 <- t(t(cm)/imd_age_pop_lad$pop)

## scale up for both risk groups 
## (assuming per capita contacts are independent of risk level)
pc_cm <- expand_contact_matrix(per_cap_matrix_45)

ng <- nrow(per_cap_matrix_45)
ndim <- nrow(pc_cm)
if (!all.equal(2*ng, ndim)){warning('dimensions not adding up')}
if (!all.equal(2*N_IMD*N_AGE, ndim)){warning('dimensions not adding up')}

## SUBTYPE-SEASONS
subtype_seasons <- readRDS(.args[4])

## SURVEILLANCE DATA
surveillance_data <- readRDS(.args[5])

## KNOWN PARAMETERS
known_pars <- readRDS(.args[6])

years <- known_pars$years

delays <- c(known_pars$primary_care_delay, known_pars$secondary_care_delay)
names(delays) <- c('primary','secondary')

## UNKNOWN PARAMETERS
unknown_pars <- readRDS(.args[7])

#### RUNNING MCMC ####

## MCMC pars
burn_in <- 50
thinning_value <- 1
n_samples <- 150
n_chains <- 2

mcmc_parallel <- function(i) {
  
  txt_output <<- i
  year_i <- years[i]
  season_i <- paste0(year_i, '/', substr(year_i + 1, 3, 4))
  list_subtypes <- subtype_seasons %>% filter(season == season_i) %>% pull(subtype) %>% unique() %>% sort()
  n_subtypes <- length(list_subtypes)
  
  children_primary_rate_low_risk_1 <- unknown_pars$base_care_rates %>% filter(year == year_i, subtype == list_subtypes[1], broad_age == "children", risk_level == "low") %>% pull(gp_rate)
  children_primary_rate_high_risk_1 <- unknown_pars$base_care_rates %>% filter(year == year_i, subtype == list_subtypes[1], broad_age == "children", risk_level == "high") %>% pull(gp_rate)
  adults_primary_rate_low_risk_1 <- unknown_pars$base_care_rates %>% filter(year == year_i, subtype == list_subtypes[1], broad_age == "adults", risk_level == "low") %>% pull(gp_rate)
  adults_primary_rate_high_risk_1 <- unknown_pars$base_care_rates %>% filter(year == year_i, subtype == list_subtypes[1], broad_age == "adults", risk_level == "high") %>% pull(gp_rate)
  older_adults_primary_rate_low_risk_1 <- unknown_pars$base_care_rates %>% filter(year == year_i, subtype == list_subtypes[1], broad_age == "older_adults", risk_level == "low") %>% pull(gp_rate)
  older_adults_primary_rate_high_risk_1 <- unknown_pars$base_care_rates %>% filter(year == year_i, subtype == list_subtypes[1], broad_age == "older_adults", risk_level == "high") %>% pull(gp_rate)
  
  children_secondary_rate_low_risk_1 <- unknown_pars$base_care_rates %>% filter(year == year_i, subtype == list_subtypes[1], broad_age == "children", risk_level == "low") %>% pull(hosp_rate)
  children_secondary_rate_high_risk_1 <- unknown_pars$base_care_rates %>% filter(year == year_i, subtype == list_subtypes[1], broad_age == "children", risk_level == "high") %>% pull(hosp_rate)
  adults_secondary_rate_low_risk_1 <- unknown_pars$base_care_rates %>% filter(year == year_i, subtype == list_subtypes[1], broad_age == "adults", risk_level == "low") %>% pull(hosp_rate)
  adults_secondary_rate_high_risk_1 <- unknown_pars$base_care_rates %>% filter(year == year_i, subtype == list_subtypes[1], broad_age == "adults", risk_level == "high") %>% pull(hosp_rate)
  older_adults_secondary_rate_low_risk_1 <- unknown_pars$base_care_rates %>% filter(year == year_i, subtype == list_subtypes[1], broad_age == "older_adults", risk_level == "low") %>% pull(hosp_rate)
  older_adults_secondary_rate_high_risk_1 <- unknown_pars$base_care_rates %>% filter(year == year_i, subtype == list_subtypes[1], broad_age == "older_adults", risk_level == "high") %>% pull(hosp_rate)
  
  initial_parameters_1 <- c(0.2, rep(x = 0.5, times = 3), 2,
                            children_primary_rate_low_risk_1, adults_primary_rate_low_risk_1, older_adults_primary_rate_low_risk_1, 
                            children_primary_rate_high_risk_1, adults_primary_rate_high_risk_1, older_adults_primary_rate_high_risk_1,
                            children_secondary_rate_low_risk_1, adults_secondary_rate_low_risk_1, older_adults_secondary_rate_low_risk_1, 
                            children_secondary_rate_high_risk_1, adults_secondary_rate_high_risk_1, older_adults_secondary_rate_high_risk_1)
  
  # c(transmissibility, 3x absolute susceptibility, log of initial infected, 
  #   reporting rates for primary care, reporting rates for secondary care)
  initial_parameters <- initial_parameters_1
  
  if (n_subtypes == 2) {
    
    children_primary_rate_low_risk_2 <- unknown_pars$base_care_rates %>% filter(year == year_i, subtype == list_subtypes[2], broad_age == "children", risk_level == "low") %>% pull(gp_rate)
    children_primary_rate_high_risk_2 <- unknown_pars$base_care_rates %>% filter(year == year_i, subtype == list_subtypes[2], broad_age == "children", risk_level == "high") %>% pull(gp_rate)
    adults_primary_rate_low_risk_2 <- unknown_pars$base_care_rates %>% filter(year == year_i, subtype == list_subtypes[2], broad_age == "adults", risk_level == "low") %>% pull(gp_rate)
    adults_primary_rate_high_risk_2 <- unknown_pars$base_care_rates %>% filter(year == year_i, subtype == list_subtypes[2], broad_age == "adults", risk_level == "high") %>% pull(gp_rate)
    older_adults_primary_rate_low_risk_2 <- unknown_pars$base_care_rates %>% filter(year == year_i, subtype == list_subtypes[2], broad_age == "older_adults", risk_level == "low") %>% pull(gp_rate)
    older_adults_primary_rate_high_risk_2 <- unknown_pars$base_care_rates %>% filter(year == year_i, subtype == list_subtypes[2], broad_age == "older_adults", risk_level == "high") %>% pull(gp_rate)
    
    children_secondary_rate_low_risk_2 <- unknown_pars$base_care_rates %>% filter(year == year_i, subtype == list_subtypes[2], broad_age == "children", risk_level == "low") %>% pull(hosp_rate)
    children_secondary_rate_high_risk_2 <- unknown_pars$base_care_rates %>% filter(year == year_i, subtype == list_subtypes[2], broad_age == "children", risk_level == "high") %>% pull(hosp_rate)
    adults_secondary_rate_low_risk_2 <- unknown_pars$base_care_rates %>% filter(year == year_i, subtype == list_subtypes[2], broad_age == "adults", risk_level == "low") %>% pull(hosp_rate)
    adults_secondary_rate_high_risk_2 <- unknown_pars$base_care_rates %>% filter(year == year_i, subtype == list_subtypes[2], broad_age == "adults", risk_level == "high") %>% pull(hosp_rate)
    older_adults_secondary_rate_low_risk_2 <- unknown_pars$base_care_rates %>% filter(year == year_i, subtype == list_subtypes[2], broad_age == "older_adults", risk_level == "low") %>% pull(hosp_rate)
    older_adults_secondary_rate_high_risk_2 <- unknown_pars$base_care_rates %>% filter(year == year_i, subtype == list_subtypes[2], broad_age == "older_adults", risk_level == "high") %>% pull(hosp_rate)
    
    initial_parameters_2 <- c(0.2, rep(x = 0.5, times = 3), 2,
                              children_primary_rate_low_risk_2, adults_primary_rate_low_risk_2, older_adults_primary_rate_low_risk_2, 
                              children_primary_rate_high_risk_2, adults_primary_rate_high_risk_2, older_adults_primary_rate_high_risk_2,
                              children_secondary_rate_low_risk_2, adults_secondary_rate_low_risk_2, older_adults_secondary_rate_low_risk_2, 
                              children_secondary_rate_high_risk_2, adults_secondary_rate_high_risk_2, older_adults_secondary_rate_high_risk_2)
    
    # c(transmissibility, 3x absolute susceptibility, log of initial infected, 
    #   reporting rates for primary care, reporting rates for secondary care)
    initial_parameters <- c(initial_parameters, initial_parameters_2)
    
  }
  
  imd_spline_pars_primary <- unknown_pars$imd_spline_pars %>% filter(year == year_i, subtype == list_subtypes[1]) %>% pull(primary)
  imd_spline_pars_secondary <- unknown_pars$imd_spline_pars %>% filter(year == year_i, subtype == list_subtypes[1]) %>% pull(secondary)
  
  initial_parameters <- c(initial_parameters, imd_spline_pars_primary, imd_spline_pars_secondary)
  
  ## population data
  vaccinated_data_seasonal <- known_pars$vaccinated_data
  
  vaccinated_data_seasonal <- vaccinated_data_seasonal %>% 
    filter(start_of_season == year_i) %>% 
    mutate(age_grp = factor(age_grp, levels = Age_Groups_Names)) %>%
    arrange(subtype, desc(risk_level), imd_quintile, age_grp)
  
  vaccinated_data_seasonal_no_rep <- vaccinated_data_seasonal %>% 
    filter(subtype == vaccinated_data_seasonal$subtype[1])
  
  ## should be ordered by IMD then age
  if(vaccinated_data_seasonal$imd_quintile[2] != 1){warning('vaccinated_data_seasonal in wrong order')}
  if(vaccinated_data_seasonal$age_grp[2] != Age_Groups_Names[2]){warning('vaccinated_data_seasonal in wrong order')}
  
  demography <- vaccinated_data_seasonal_no_rep %>% 
    mutate(population = pop) %>% 
    select(age_grp, imd_quintile, risk_level, population, risk_proportion) %>% 
    arrange(desc(risk_level), imd_quintile, age_grp)
  
  ## check population sum is correct
  tot_pop <- sum(imd_age_pop_lad$pop)
  if(!all.equal(sum(demography$population), tot_pop)){warning('pop not adding up')}
  
  # demography_input <- demography
  # vaccinated_input <- vaccinated_data_seasonal
  # subtype_season <- subtype_seasons %>% filter(season == season_i)
  # cm_input <- pc_cm
  # epidemic_to_fit <- surveillance_data %>% filter(index == i)
  # epid_periods_input <- known_pars$epid_periods
  # coverage_rates <- known_pars$proportion_observed %>% filter(year == year_i)
  # care_delays <- delays
  # initial_parameters <- c(rep(x = subtype_init_pars, times = n_subtypes), rep(0, 4))
  # #   subtype-specific parameters x2, IMD spline parameters x4
  # save_all <- T
  # n_samples <- n_samples
  # nburn <- burn_in
  # thinning <- thinning_value
  # n_chains <- n_chains
  # n_cores <- 1
  # txt_output <- paste0(i)
  
  output_i <- run_mcmc_inference_v2(demography_input = demography,
                                    vaccinated_input = vaccinated_data_seasonal,
                                    subtype_season = subtype_seasons %>% filter(season == season_i),
                                    cm_input = pc_cm,
                                    epidemic_to_fit = surveillance_data %>% filter(index == i),
                                    epid_periods_input = known_pars$epid_periods,
                                    coverage_rates = known_pars$proportion_observed %>% filter(year == year_i),
                                    care_delays = delays,
                                    initial_parameters = initial_parameters,
                                    # subtype-specific parameters x2, IMD spline parameters x4
                                    save_all = T,
                                    n_samples = n_samples,
                                    nburn = burn_in,
                                    thinning = thinning_value,
                                    n_chains = n_chains,
                                    txt_output = paste0(i))
  
  return(output_i)
  
}


# mcmc_results <- mclapply(1:length(years), mcmc_parallel, mc.cores = length(years))
# mcmc_results <- mclapply(1:length(years), mcmc_parallel, mc.cores = 1)
mcmc_results <- mclapply(3:3, mcmc_parallel, mc.cores = 1)



#### SAVE RESULTS ####
time_stamp <- Sys.time()
time_stamp_char <- as.character(1000*round(as.numeric(time_stamp), 3))

# save most recently run settings as a dummy save
saveRDS(object = data.table(burn_in = burn_in,
                            thinning = thinning_value,
                            n_samples = n_samples,
                            n_chains = n_chains,
                            HPC = F,
                            date = Sys.Date()),
        file = paste0(Dir_Local, "/", Subdir_Output, "/", Subdir_Data, "/", "parameters_", time_stamp_char, ".rds"))


# save actual data
saveRDS(object = mcmc_results,
        file = paste0(Dir_Local, "/", Subdir_Output, "/", Subdir_Data, "/", "output_", time_stamp_char, ".rds"))



summary(mcmc_results[[1]])
tracePlot(mcmc_results[[1]], start = 50)
gelmanDiagnostics(mcmc_results[[1]], start = 50)

samples <- as_tibble(getSample(mcmc_results[[1]], start = 50))
mean(samples %>% pull(1))
mean(samples %>% pull(2))
getCredibleIntervals(samples, quantiles = c(0.025, 0.975))



test_output <- readRDS(file = paste0(Dir_Local, "/output/data/", "output_1794960528861.rds"))
test_parameters <- readRDS(file = paste0(Dir_Local, "/output/data/", "parameters_1794960528861.rds"))
test_burn_in <- test_parameters$burn_in

summary(test_output[[1]])
tracePlot(test_output[[1]], start = test_burn_in)
gelmanDiagnostics(test_output[[1]])

test_samples <- as_tibble(getSample(test_output, start = test_burn_in))
test_mean_01 <- mean(test_samples %>% pull(1))
test_mean_02 <- mean(test_samples %>% pull(2))
test_mean_03 <- mean(test_samples %>% pull(3))
test_mean_04 <- mean(test_samples %>% pull(4))
test_mean_05 <- 10^mean(test_samples %>% pull(5))
test_credible_intervals <- getCredibleIntervals(test_samples, quantiles = c(0.025, 0.975))
test_credible_intervals[,5] <- 10^test_credible_intervals[,5]

test_true_values_trans <- unknown_pars$epid_parameters %>%
  filter(year == 2025) %>%
  pull(transmissibility)

test_true_values_children_susc <- unknown_pars$epid_parameters %>%
  filter(year == 2025) %>%
  pull(children_susceptibility)

test_true_values_adults_susc <- unknown_pars$epid_parameters %>%
  filter(year == 2025) %>%
  pull(adults_susceptibility)

test_true_values_older_adults_susc <- unknown_pars$epid_parameters %>%
  filter(year == 2025) %>%
  pull(older_adults_susceptibility)

test_true_values_init_infected <- unknown_pars$epid_parameters %>%
  filter(year == 2025) %>%
  pull(init_infected)

test_coverage <- tibble(low = test_credible_intervals[1,],
                        high = test_credible_intervals[2,],
                        estimate = c(test_mean_01, test_mean_02, test_mean_03, test_mean_04, test_mean_05),
                        true_value = c(test_true_values_trans, test_true_values_children_susc,
                                       test_true_values_adults_susc, test_true_values_older_adults_susc,
                                       test_true_values_init_infected)) %>%
  mutate(coverage = low <= true_value & true_value <= high)


