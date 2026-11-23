#### MCMC FUNCTIONS ####

TEXT_SAVE_DIR <- file.path('mcmc_output')
if(!dir.exists(TEXT_SAVE_DIR)){dir.create(TEXT_SAVE_DIR)}

## MCMC FITTING ##
run_mcmc_inference <- function(demography_input,
                               demography_gp,
                               vaccinated_input,
                               subtype_season,
                               cm_input, 
                               epidemic_to_fit, 
                               epid_periods_input,
                               coverage_rates,
                               care_delays,
                               initial_parameters,
                               save_all = F,
                               n_samples, 
                               nburn, 
                               thinning,
                               n_chains,
                               n_cores = 1,
                               txt_output = NULL) {
  
  ## SET UP DATA FRAMES ETC. ##
  
  epidemic_dt <- as.data.table(epidemic_to_fit)
  coverage_rates$imd_quintile <- factor(coverage_rates$imd_quintile)
  n_subtypes <- n_distinct(subtype_season$subtype)
  n_pars_per_subtype <- (length(initial_parameters) - 4)/n_subtypes
  
  # DEzs sampler will require 3 chains
  n_pop_internal <- 3
  
  ll_call_count <- 0
  ll_total_calls <- nburn + n_samples*thinning
  mod_val <- if (ll_total_calls < 100) 10 else if (ll_total_calls < 2000) 100 else 1000
  
  txt_out <- file.path('mcmc_output', paste0('index_',txt_output,'.txt'))
  
  LLcache <- new.env(hash = TRUE, parent = emptyenv())
  
  make_key <- function(par) paste(sprintf("%.15g", par), collapse = "_")
  
  # Define the log likelihood function
  single_ll <- function(pars) {
    
    #### PARAMETERS: EPIDEMIC 1 ####
    epidemic_1 <- unique(subtype_season$subtype)[1]
    
    transmissibility_1 <- pars[1]
    
    susceptibility_vec_1 <- fcn_assign_ages(pars[2],
                                            pars[3],
                                            pars[4],
                                            Age_Groups_Names)
    
    init_infected_1 <- 10^(pars[5])
    care_rate_df_1 <- data.frame(broad_age = rep(x = Broad_Ages_Names, times = 2),
                                 risk_level = rep(c('low', 'high'), each = 3),
                                 primary_care = pars[6:11],
                                 secondary_care = pars[12:17])
    
    care_rate_age_df_1 <- data.table(cross_join(care_rate_df_1, Table_Broad_Ages_Age_Groups))[broad_age.x == broad_age.y, ]
    care_rate_age_df_1[, c('broad_age.x', 'broad_age.y') := NULL]
    
    imd_spline_pars <- data.table(primary = pars[1:2 + n_subtypes*n_pars_per_subtype],
                                  secondary = pars[3:4 + n_subtypes*n_pars_per_subtype])
    rel_imd_rep_rates <- data.frame(imd_quintile = 1:5,
                                    rel_primary_rates = imd_spline(imd_spline_pars$primary),
                                    rel_secondary_rates = imd_spline(imd_spline_pars$secondary))
    
    reporting_rates_1 <- cross_join(care_rate_age_df_1, rel_imd_rep_rates) %>% 
      mutate(primary_care = primary_care*rel_primary_rates,
             secondary_care = secondary_care*rel_secondary_rates) %>% 
      select(!c(rel_primary_rates, rel_secondary_rates))
    
    reporting_dt_1 <- as.data.table(reporting_rates_1)
    reporting_long_1 <- melt(reporting_dt_1,
                             measure.vars = c('primary_care', 'secondary_care'),
                             variable.name = 'setting',
                             value.name = 'rate')
    
    # any out-of-bounds proposals slipping past the prior
    if (any(is.na(pars)) || 
        transmissibility_1 < min_trans || 
        sum(susceptibility_vec_1 < min_susc) > 0  || sum(susceptibility_vec_1 > max_susc) > 0 ||
        init_infected_1 < 0 || init_infected_1 > min(demography_input$population)) {
      
      return(-Inf)
      
    }
    
    if (n_subtypes == 2) {
      
      #### PARAMETERS: EPIDEMIC 2 ####
      epidemic_2 <- unique(subtype_season$subtype)[2]
      
      transmissibility_2 <- pars[1 + n_pars_per_subtype]
      
      susceptibility_vec_2 <- fcn_assign_ages(pars[2 + n_pars_per_subtype],
                                              pars[3 + n_pars_per_subtype],
                                              pars[4 + n_pars_per_subtype],
                                              Age_Groups_Names)
      
      init_infected_2 <- 10^(pars[5 + n_pars_per_subtype])
      care_rate_df_2 <- data.frame(broad_age = rep(x = Broad_Ages_Names, times = 2),
                                   risk_level = rep(c('low','high'), each = 3),
                                   primary_care = pars[6:11 + n_pars_per_subtype],
                                   secondary_care = pars[12:17 + n_pars_per_subtype])
      care_rate_age_df_2 <- data.table(cross_join(care_rate_df_2, Table_Broad_Ages_Age_Groups))[broad_age.x == broad_age.y,]
      care_rate_age_df_2[, c('broad_age.x','broad_age.y') := NULL]
      
      reporting_rates_2 <- cross_join(care_rate_age_df_2, rel_imd_rep_rates) %>% 
        mutate(primary_care = primary_care*rel_primary_rates,
               secondary_care = secondary_care*rel_secondary_rates) %>% 
        select(!c(rel_primary_rates, rel_secondary_rates))
      
      reporting_dt_2 <- as.data.table(reporting_rates_2)
      reporting_long_2 <- melt(reporting_dt_2,
                               measure.vars = c('primary_care', 'secondary_care'),
                               variable.name = 'setting',
                               value.name = 'rate')
      
      # any out-of-bounds proposals slipping past the prior
      if (any(is.na(pars)) || 
          transmissibility_2 < min_trans || 
          sum(susceptibility_vec_2 < min_susc) > 0  || sum(susceptibility_vec_2 > max_susc) > 0 ||
          init_infected_2 < 0 || init_infected_2 > min(demography_input$population) ) {
        return(-Inf)
      }
      
    }
    
    
    if (n_subtypes == 2) {
      
      #### RUNNING MULTIPLE EPIDEMICS ####
      epi_data_1 <- produce_epi_data(demography_gp = demography_gp, 
                                     vaccinated_input = vaccinated_input,
                                     subtype_input = epidemic_1,
                                     cm_input = cm_input, 
                                     init_infected = init_infected_1,
                                     transmissibility = transmissibility_1,
                                     susceptibility_vec = susceptibility_vec_1,
                                     start_of_epidemic = as.Date(paste0('01-09-', year(epidemic_dt$week_start[1])),
                                                                 format = '%d-%m-%Y'),
                                     epid_periods_input = epid_periods_input,
                                     reporting_long = reporting_long_1)
      
      epi_data_2 <- produce_epi_data(demography_gp = demography_gp, 
                                     vaccinated_input = vaccinated_input,
                                     subtype_input = epidemic_2,
                                     cm_input = cm_input, 
                                     init_infected = init_infected_2,
                                     transmissibility = transmissibility_2,
                                     susceptibility_vec = susceptibility_vec_2,
                                     start_of_epidemic = as.Date(paste0('01-09-', year(epidemic_dt$week_start[1])),
                                                                 format = '%d-%m-%Y'),
                                     epid_periods_input = epid_periods_input,
                                     reporting_long = reporting_long_2)
      
      weekly_proportions <- epi_data_1[[1]] %>%
        dplyr::rename(c("infections_1" = "infections")) %>%
        dplyr::left_join(epi_data_2[[1]], by = c("week_start")) %>%
        dplyr::rename(c("infections_2" = "infections")) %>%
        mutate(modelled_infections = infections_1 + infections_2) %>%
        mutate(modelled_proportion = infections_1 / modelled_infections) %>%
        dplyr::select(-c("infections_1", "infections_2")) %>%
        dplyr::rename(c("date_formatted" = "week_start")) %>%
        drop_na() %>% 
        ## IF 0, CHANGE TO VERY LOW NUMBER TO AVOID -INF LOG-LIKELIHOOD
        ## IF 1, CHANGE TO VERY HIGH NUMBER TO AVOID -INF LOG-LIKELIHOOD
        mutate(min_modelled_proportion = min(modelled_proportion[modelled_proportion != 0]),
               max_modelled_proportion = max(modelled_proportion[modelled_proportion != 1]),
               modelled_proportion = case_when(modelled_proportion == 0 ~ min_modelled_proportion/1e3,
                                               modelled_proportion == 1 ~ 1 - (1 - max_modelled_proportion)/1e3,
                                               T ~ modelled_proportion)) %>%
        select(!c(min_modelled_proportion, max_modelled_proportion))
      
      weekly_true_proportions <- subtype_season %>%
        filter(subtype == epidemic_1) %>% 
        rename(epidemic_1_positive = value) %>% 
        select(date_formatted, proportion, total_flu, epidemic_1_positive)
      
      weekly_props <- weekly_true_proportions %>% 
        left_join(weekly_proportions, by = 'date_formatted') %>%
        drop_na()
      
      # Merge subtype-specific data and aggregate
      MODELLED_EPI <- rbind(epi_data_1[[2]], epi_data_2[[2]])[, .(expected_cases = sum(expected_cases_gp)),
                                                              by = list(week_start, age_grp, imd_quintile, risk_level, setting)]
      
    } else {
      
      #### RUNNING SINGLE EPIDEMIC ####
      epi_data_1 <- produce_epi_data(demography_gp = demography_gp, 
                                     vaccinated_input = vaccinated_input,
                                     subtype_input = epidemic_1,
                                     cm_input = cm_input, 
                                     init_infected = init_infected_1,
                                     transmissibility = transmissibility_1,
                                     susceptibility_vec = susceptibility_vec_1,
                                     start_of_epidemic = as.Date(paste0('01-09-', year(epidemic_dt$week_start[1])),
                                                                 format = '%d-%m-%Y'),
                                     epid_periods_input = epid_periods_input,
                                     reporting_long = reporting_long_1)
      
      MODELLED_EPI <- epi_data_1[[2]] %>% dplyr::rename(c("expected_cases" = "expected_cases_gp"))
      
    }
    
    # Pivot observed data longer
    observed_data <- melt(epidemic_dt, 
                          measure.vars = c('primary_care', 'secondary_care'),
                          variable.name = 'setting', 
                          value.name = 'observations')
    
    ## Shift observations by the care delays
    observed_data <- rbind(observed_data[setting == 'primary_care'][, observations := shift(observations, n = -care_delays['primary'], fill = 0),
                                                                    by = .(age_grp, imd_quintile, risk_level, index)],
                           observed_data[setting == 'secondary_care'][, observations := shift(observations, n = -care_delays['secondary'], fill = 0),
                                                                      by = .(age_grp, imd_quintile, risk_level, index)])
    
    # Merge with observed data
    time_series_joint <- merge(MODELLED_EPI, 
                               observed_data,
                               by = c('age_grp', 'imd_quintile', 'week_start', 'risk_level', 'setting'), 
                               all.x = F)
    
    # Add 1e-10 where time_series_joint$expected_cases is 0
    time_series_joint[expected_cases == 0, expected_cases := 1e-10]
    
    if (F) {
      time_series_joint %>% ## plot {observations} against {infections x reporting rates}
        ggplot() +
        theme_bw() +
        geom_line(aes(week_start,expected_cases, lty = risk_level, col = setting, group = interaction(setting, risk_level))) +
        geom_point(aes(week_start,observations, shape = risk_level, col = setting, group = interaction(setting, risk_level)),
                   alpha = 0.4) +
        facet_grid(age_grp ~ imd_quintile, scales = 'free')
    }
    
    time_series_primary <- time_series_joint[setting == 'primary_care',]
    time_series_secondary <- time_series_joint[setting == 'secondary_care',]
    
    # Vectorised log likelihood, primary care
    ll_1 <- sum(dpois(x = time_series_primary$observations,
                      lambda = time_series_primary$expected_cases,
                      log  = TRUE), na.rm = TRUE)
    
    # Vectorised log likelihood, secondary care
    ll_2 <- sum(dpois(x = time_series_secondary$observations,
                      lambda = time_series_secondary$expected_cases,
                      log = TRUE), na.rm = TRUE)
    
    if (n_subtypes == 2) {
      
      # Vectorised log-likelihood, subtype proportions
      ll_3 <- sum(dbinom(x = weekly_props$epidemic_1_positive,
                         size = weekly_props$total_flu,
                         p = weekly_props$modelled_proportion,
                         log = TRUE), na.rm = TRUE)
      
    } else { ll_3 <- 0 }
    
    assign(make_key(pars), c(ll_1, ll_2, ll_3), envir = LLcache)
    
    total_ll <- ll_1 + ll_2 + ll_3
    
    if (is.nan(total_ll) | is.infinite(total_ll)) return(-Inf)
    
    return(total_ll)
    
  }
  
  llikelihood <- function(pars) {
    
    if (ll_call_count == 0) { ll_start_time <<- Sys.time() }
    ll_call_count <<- ll_call_count + 1
    if (ll_call_count %% mod_val == 0) {
      now_time <- Sys.time()
      elapsed <- as.numeric(difftime(now_time, ll_start_time, units = 'mins'))
      rate <- ll_call_count / elapsed
      remaining <- (ll_total_calls - ll_call_count) / rate
      writeLines(sprintf("INDEX %s: Generation %d / %d (%.1f%%) | Elapsed: %.1f min \n 
        Est. remaining: %.1f min | Time per 1000 iterations: %.1f min \n  
        Last updated: %s",
                         txt_output, ll_call_count, ll_total_calls,
                         100 * ll_call_count / ll_total_calls, elapsed, remaining, 1000 / rate, as.character(now_time)),
                 txt_out)
    }
    
    if (is.matrix(pars)) {
      
      res <- parallel::mclapply(seq_len(nrow(pars)),
                                function(j) single_ll(pars[j, ]),
                                mc.cores = n_cores)
      
      return(unlist(res))
      
    } else {
      
      return(single_ll(pars))
      
    }
    
  }
  
  
  
  # TODO Should these be lognormal instead?
  
  # set up priors and lower/upper bounds
  {
    # IGPR and IHR centred roughly around the crude estimates (based on attack rate etc.)
    prim_beta  <- beta_pars(initial_parameters[6], 200) 
    sec_beta   <- beta_pars(initial_parameters[12], 200)
    # Susceptibility beta distributed, mean at 0.5
    shape1_susc <- 2
    shape2_susc <- 2
    # R0 centred at 2 
    r0_mean <- 2
    r0_sd   <- 0.4
    # IMD splines normal around 0
    imd_mean <- 0
    imd_sd <- 0.5
    
    ## set bounds
    min_trans <- 0; max_trans <- 1
    min_susc <- 0; max_susc <- 1
    min_r0 <- 1; max_r0 <- 4
    min_log_init_inf <- 0; max_log_init_inf <- 4
    min_reporting <- 0; max_reporting <- 1
    min_spline <- -log(5); max_spline <- log(5) # equivalent to IMD ratios at most
    # 5x higher/lower than each other in IMD 1 vs 3 or IMD 5 vs 3,
    # which would lead to huge ratios between IMD 1 and IMD 5 at the extremes
    
    ## set up sampler
    lower_vals <- c(rep(c(min_trans, rep(min_susc, 3), min_log_init_inf, rep(min_reporting, 12)), n_subtypes),
                    rep(min_spline, 4))
    
    upper_vals <- c(rep(c(max_trans, rep(max_susc, 3), max_log_init_inf, rep(max_reporting, 12)), n_subtypes),
                    rep(max_spline, 4))
    
  }
  
  
  llprior <- function(pars) {
    
    # pars[1] is transmissibility, pars[2:4] are susceptibility (children, adults, older adults), 
    # pars[5] is log of initial infected on 1st september
    
    susceptibility_vec_1 <- fcn_assign_ages(pars[2],
                                            pars[3],
                                            pars[4],
                                            Age_Groups_Names)
    
    R0_1 <- R0_func(susceptibility = susceptibility_vec_1,
                    inf_period = epid_periods_input[2],
                    beta_in = pars[1],
                    cm_in = cm_input,
                    per_capita = T,
                    population_vector = demography_input$population)
    
    reporting_vec_1 <- c(6:n_pars_per_subtype)
    reporting_prim_vec_1 <- c(6:11)
    reporting_sec_vec_1 <- reporting_vec_1[!(reporting_vec_1 %in% reporting_prim_vec_1)]
    
    if (
      
      ## TRANSMISSIBILITY
      pars[1] < min_trans || pars[1] > max_trans ||
      
      ## SUSCEPTIBILITY
      sum(susceptibility_vec_1 < min_susc) > 0 || sum(susceptibility_vec_1 > max_susc) > 0 ||
      
      ## INITIAL INFECTED
      pars[5] < min_log_init_inf || pars[5] > max_log_init_inf ||
      
      ## REPORTING
      any(pars[reporting_vec_1] < min_reporting) || any(pars[reporting_vec_1] > max_reporting) ||
      
      R0_1 < min_r0 || R0_1 > max_r0 
      
    ) {return(-Inf)}
    
    ## SECOND SUBTYPE
    if (n_subtypes == 2){
      
      susceptibility_vec_2 <- fcn_assign_ages(pars[2 + n_pars_per_subtype],
                                              pars[3 + n_pars_per_subtype],
                                              pars[4 + n_pars_per_subtype],
                                              Age_Groups_Names)
      
      R0_2 <- R0_func(susceptibility = susceptibility_vec_2,
                      inf_period = epid_periods_input[2],
                      beta_in = pars[1 + n_pars_per_subtype],
                      cm_in = cm_input,
                      per_capita = T,
                      population_vector = demography_input$population)
      
      reporting_vec_2 <- c(n_pars_per_subtype + 6:n_pars_per_subtype)
      reporting_prim_vec_2 <- c(n_pars_per_subtype + 6:11)
      reporting_sec_vec_2 <- reporting_vec_2[!(reporting_vec_2 %in% reporting_prim_vec_2)]
      
      if (
        
        ## TRANSMISSIBILITY
        pars[1 + n_pars_per_subtype] < min_trans || pars[1 + n_pars_per_subtype] > max_trans ||
        
        ## SUSCEPTIBILITY
        sum(susceptibility_vec_2 < min_susc) > 0 || sum(susceptibility_vec_2 > max_susc) > 0 ||
        
        ## INITIAL INFECTED
        pars[5 + n_pars_per_subtype] < min_log_init_inf || pars[5 + n_pars_per_subtype] > max_log_init_inf ||
        
        ## REPORTING
        any(pars[reporting_vec_2] < min_reporting) || any(pars[reporting_vec_2] > max_reporting) ||
        
        R0_2 < min_r0 || R0_2 > max_r0 
        
      ) {return(-Inf)}
      
    }
    
    ## IMD SPLINES
    if (any(pars[1:4 + n_subtypes*n_pars_per_subtype] < min_spline) ||
        any(pars[1:4 + n_subtypes*n_pars_per_subtype] > max_spline)) {return(-Inf)}
    
    lprob <- 0
    
    # Scaled normal prior on R0 (pnorm(min_r0) is the log normalising constant)
    lprob <- lprob + dnorm(R0_1, mean = r0_mean, sd = r0_sd, log = TRUE) -
      pnorm(min_r0, mean = r0_mean, sd = r0_sd, lower.tail = FALSE, log.p = TRUE)
    # Uniform prior on transmissibility 
    lprob <- lprob + dunif(unname(pars[1]), min = min_trans, max = max_trans, log = TRUE)
    # Beta prior on susceptibility
    lprob <- lprob + sum(dbeta(unname(pars[2:4]), shape1 = shape1_susc, shape2 = shape2_susc, log = TRUE))
    # Uniform prior on initial infected
    lprob <- lprob + dunif(unname(pars[5]), min = min_log_init_inf, max = max_log_init_inf, log = TRUE)
    # Beta prior on primary care reporting rates (pars 5:10), centred at 2%
    lprob <- lprob + sum(dbeta(pars[6:11], 
                               shape1 = prim_beta['a'], 
                               shape2 = prim_beta['b'], log = TRUE))
    # Beta prior on secondary care reporting rates (pars 11:16), centred at 0.5%
    lprob <- lprob + sum(dbeta(pars[12:17], 
                               shape1 = sec_beta['a'], 
                               shape2 = sec_beta['b'], log = TRUE))
    
    if (n_subtypes == 2) {
      
      # priors on second set of epi parameters
      
      lprob <- lprob + dnorm(R0_2, mean = r0_mean, sd = r0_sd, log = TRUE) -
        pnorm(min_r0, mean = r0_mean, sd = r0_sd, lower.tail = FALSE, log.p = TRUE)
      
      lprob <- lprob + dunif(unname(pars[1 + n_pars_per_subtype]), min = min_trans, max = max_trans, log = TRUE)
      
      lprob <- lprob + sum(dbeta(unname(pars[2:4 + n_pars_per_subtype]), shape1 = shape1_susc, shape2 = shape2_susc, log = TRUE))
      
      lprob <- lprob + dunif(unname(pars[5 + n_pars_per_subtype]), min = min_log_init_inf, max = max_log_init_inf, log = TRUE)
      
      lprob <- lprob + sum(dbeta(pars[6:11 + n_pars_per_subtype], 
                                 shape1 = prim_beta['a'], 
                                 shape2 = prim_beta['b'], log = TRUE))
      
      lprob <- lprob + sum(dbeta(pars[12:17 + n_pars_per_subtype], 
                                 shape1 = sec_beta['a'], 
                                 shape2 = sec_beta['b'], log = TRUE))
      
    }
    
    # Normal prior on spline parameters
    lprob <- lprob + sum(dnorm(unname(pars[1:4 + n_subtypes*n_pars_per_subtype]), mean = imd_mean, sd = imd_sd, log = TRUE))
    
    return(lprob)
  }
  
  
  sampler <- function(n = 1) {
    
    out <- matrix(NA, nrow = n, ncol = length(initial_parameters))
    
    for(j in 1:n){
      
      repeat {
        valid <- FALSE
        while (!valid){
          repeat {
            R0_1 <- rnorm(1, mean = r0_mean, sd = r0_sd)
            if (R0_1 >= min_r0 && R0_1 <= max_r0) break
          }
          if (n_subtypes == 2) {
            repeat {
              R0_2 <- rnorm(1, mean = r0_mean, sd = r0_sd)
              if (R0_2 >= min_r0 && R0_2 <= max_r0) break
            }
          }
          
          susc_1 <- rbeta(3, shape1_susc, shape2_susc)
          
          trans_1 <- R0_func(susceptibility = fcn_assign_ages(susc_1[1],
                                                              susc_1[2],
                                                              susc_1[3], 
                                                              Age_Groups_Names),
                             inf_period= epid_periods_input[2],
                             beta_in = 1,
                             cm_in = cm_input,
                             per_capita = TRUE,
                             population_vector = demography_input$population,
                             R0assumed = R0_1,
                             return_beta = TRUE)
          
          if (n_subtypes == 2) {
            
            susc_2 <- rbeta(3, shape1_susc, shape2_susc)
            
            trans_2 <- R0_func(susceptibility = fcn_assign_ages(susc_2[1],
                                                                susc_2[2],
                                                                susc_2[3],  
                                                                Age_Groups_Names),
                               inf_period = epid_periods_input[2],
                               beta_in = 1,
                               cm_in = cm_input,
                               per_capita = TRUE,
                               population_vector = demography_input$population,
                               R0assumed = R0_2,
                               return_beta = TRUE)
            
          }
          
          if (n_subtypes == 2) {
            
            if(trans_1 >= min_trans & trans_1 <= max_trans & 
               trans_2 >= min_trans & trans_2 <= max_trans) valid <- TRUE
            
          } else {
            
            if (trans_1 >= min_trans & trans_1 <= max_trans) valid <- TRUE
            
          }
          
        }
        
        # Check trans is within bounds
        if (trans_1 < min_trans | trans_1 > max_trans) next
        
        if (n_subtypes == 2) {
          if (trans_2 < min_trans | trans_2 > max_trans) next
        }
        
        ## IMD splines normal around 0, SD = 0.5
        imd_spline_samples <- rnorm(4, mean = imd_mean, sd = imd_sd) 
        
        # Check IMD splines are within bounds
        if (sum(imd_spline_samples < min_spline) + sum(imd_spline_samples > max_spline) > 0) next
        
        # Draw reporting rates from Beta priors
        prim_rates_1 <- rbeta(6, shape1 = prim_beta['a'], shape2 = prim_beta['b'])
        sec_rates_1  <- rbeta(6, shape1 = sec_beta['a'],  shape2 = sec_beta['b'])
        
        if (n_subtypes == 2) {
          prim_rates_2 <- rbeta(6, shape1 = prim_beta['a'], shape2 = prim_beta['b'])
          sec_rates_2  <- rbeta(6, shape1 = sec_beta['a'],  shape2 = sec_beta['b'])
        }
        
        out[j, ] <- if(n_subtypes == 2){
          
          c(trans_1, susc_1,
            runif(1, min_log_init_inf, max_log_init_inf),  # log init infected, uniform
            prim_rates_1,                                     # primary care rates
            sec_rates_1,                                      # secondary care rates
            trans_2, susc_2,
            runif(1, min_log_init_inf, max_log_init_inf),  # log init infected, uniform
            prim_rates_2,                                     # primary care rates
            sec_rates_2,                                      # secondary care rates
            imd_spline_samples)       # IMD spline, normal around 0
          
        } else {
          
          c(trans_1, susc_1,
            runif(1, min_log_init_inf, max_log_init_inf),  # log init infected, uniform
            prim_rates_1,                                     # primary care rates
            sec_rates_1,                                  # secondary care rates
            imd_spline_samples)       # IMD spline, normal around 0
          
        }
        
        break
        
      }
      
    }
    
    return(out)
    
  }
  
  ## prior
  prior <- createPrior(density = llprior, 
                       sampler = sampler,
                       lower = lower_vals,
                       upper = upper_vals)
  
  bayesianSetup <- createBayesianSetup(likelihood = llikelihood, 
                                       prior = prior,
                                       parallel = "external")
  
  settings <- list(iterations = (nburn + n_samples*thinning) * n_pop_internal,
                   burnin = nburn,
                   thin = thinning,
                   message = F,
                   nrChains = n_chains,
                   startValue = prior$sampler(n_pop_internal))
  
  out <- runMCMC(bayesianSetup = bayesianSetup, sampler = 'DEzs', settings = settings)
  
  samp <- getSample(out, start = nburn + 1, thin = thinning, coda = TRUE)
  
  nPar <- out$setup$numPars
  
  lookupComponents <- function(chainMat, cache, nPar) {
    parMat <- chainMat[, 1:nPar, drop = FALSE]
    comp <- t(apply(parMat, 1, function(par) {
      val <- mget(make_key(par), envir = cache, ifnotfound = list(c(NA_real_, NA_real_, NA_real_)))[[1]]
      val
    }))
    colnames(comp) <- c("LL1", "LL2","LL3")
    coda::mcmc(cbind(chainMat, comp))
  }
  
  LLcomponents <- if (coda::is.mcmc.list(out$chain)) {
    coda::as.mcmc.list(lapply(out$chain, lookupComponents, cache = LLcache, nPar = nPar))
  } else {
    lookupComponents(out$chain, LLcache, nPar)
  }
  
  for (sample in 1:length(samp)) {
    
    samp[[sample]] <- cbind(samp[[sample]], 
                            LLcomponents[[sample]][nburn + thinning*(1:n_samples), c('LP','LL','LPr','LL1','LL2','LL3')])
    
  }
  
  if (save_all) {
    
    return(out)
    
  } else {
    
    return(samp)
    
  }
  
}





