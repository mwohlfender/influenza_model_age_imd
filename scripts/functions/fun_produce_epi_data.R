##-----------------------##
#### EPIDEMIC FUNCTION ####
##-----------------------##


# subtype_input <- "AH1N1"
# init_infected <- 100
# transmissibility <- 0.2
# susceptibility_vec <- fcn_assign_ages(0.8,
#                                       0.5,
#                                       0.1,
#                                       Age_Groups_Names)
# start_of_epidemic <- lubridate::ymd("2023-09-01")


produce_epi_data <- function(vaccinated_input,
                             subtype_input,
                             cm_input, 
                             init_infected,
                             transmissibility,
                             susceptibility_vec,
                             start_of_epidemic,
                             epid_periods_input,
                             coverage_rates,
                             reporting_long) {
  
  vaccinated_subtype <- vaccinated_input %>% filter(vaccinated_input$subtype == subtype_input)
  
  run_A_shorter <- F
  t_end_in <- if (run_A_shorter) {
    ifelse(subtype_input == 'B', 300, 250)
  } else {300}
  
  # pop <- vaccinated_subtype$pop
  # I0 <- init_infected
  # vacc_cov <- vaccinated_subtype$vaccinated_population
  # ve_inf <- vaccinated_subtype$VE_INF
  # cm <- cm_input
  # trans <- transmissibility
  # susc <- susceptibility_vec
  # lat_per <- epid_periods_input[1]
  # inf_per <- epid_periods_input[2]
  # t_end <- t_end_in
  
  time_series <- run_model(pop = vaccinated_subtype$pop,
                           I0 = init_infected,
                           vacc_cov = vaccinated_subtype$vaccinated_population,
                           ve_inf = vaccinated_subtype$VE_INF,
                           cm = cm_input,
                           trans = transmissibility,
                           susc = susceptibility_vec,
                           lat_per = epid_periods_input[1],
                           inf_per = epid_periods_input[2],
                           t_end = t_end_in)
  
  if (run_A_shorter & subtype_input != 'B') {
    
    add_time_series <- time_series[t %in% 1:50, ]
    add_time_series[, t := t + 250][, infections := 0]
    time_series <- rbind(time_series, add_time_series)
    setorder(time_series, 'age_grp','imd_quintile','risk_level','vaccinated')
    
  }
  
  ## add start date (using first of September throughout)
  time_series <- time_series[, date := start_of_epidemic + t] ## add date
  time_series[, t := NULL]
  
  if (any(is.na(time_series))) {
    
    print("Warning_fun_produce_epi_data_1")
    
    return(-Inf)
    
  }
  
  time_series_fit <- time_series[coverage_rates, on = c("age_grp", "imd_quintile", "risk_level", "vaccinated")]
  time_series_fit[, infections_gp := prop * infections]


  if (any(is.na(time_series_fit))) {
    
    print("Warning_fun_produce_epi_data_2")
    
    return(-Inf)
    
  }
  
  # Aggregate to weekly and join with surveillance data
  time_series_fit[, week_start := last_monday(date)]
  time_series_weekly <- time_series_fit[, .(infections = sum(infections), infections_gp = sum(infections_gp)), 
                                        by = .(week_start, age_grp, imd_quintile, risk_level, vaccinated)]
  
  weekly_infections <- time_series_weekly %>%
    group_by(week_start) %>%
    summarize(infections = sum(infections))
  
  # Join reporting rates (make primary and secondary care versions)
  time_series_long <- merge(rbind(cbind(time_series_weekly, setting = unique(reporting_long$setting)[1]),
                                  cbind(time_series_weekly, setting = unique(reporting_long$setting)[2])), 
                            reporting_long,
                            by = colnames(reporting_long)[colnames(reporting_long) != 'rate'],
                            all.x = TRUE)
  
  # Scale down hospitalisation rate for those ineffectively vaccinated
  # (VE_HOSP must be geq VI_INF)
  time_series_long <- merge(time_series_long, 
                            vaccinated_subtype[start_of_season == year(start_of_epidemic), 
                                               c('age_grp','imd_quintile','risk_level','VE_INF','VE_HOSP')],
                            by = c('age_grp','imd_quintile','risk_level'),
                            all.x = TRUE)
  
  time_series_long[VE_HOSP >= VE_INF & vaccinated & setting == 'secondary_care',
                   rate := rate*((1 - VE_HOSP)/(1 - VE_INF))]
  
  time_series_long[, c('VE_HOSP','VE_INF') := NULL]
  
  time_series_long <- time_series_long %>%
    dplyr::relocate("week_start", "age_grp", "imd_quintile", "risk_level", "vaccinated", "setting", "rate", "infections", "infections_gp")
  
  # Validate reporting rates
  if (any(is.na(time_series_long$rate))) {
    
    print("Warning_fun_produce_epi_data_3")
    
    return(-Inf)
    
  }
  
  if (any(time_series_long$rate <= 0 | time_series_long$rate >= 1)) {
    
    print("Warning_fun_produce_epi_data_4")
    
    return(-Inf)
    
  }
  
  # Aggregate over vaccination status into a mean number of cases per {group x setting}
  # (as the rate varies between vaccinated/unvaccinated, but observations are not vaccination-specific)
  time_series_agg <- time_series_long[, .(expected_visits = sum(infections_gp*rate)),
                                      by = list(week_start, age_grp, imd_quintile, risk_level, setting)]

  return(list(weekly_infections, time_series_agg))
  
}


