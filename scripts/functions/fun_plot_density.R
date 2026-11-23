## FUNCTIONS TO PLOT MCMC SAMPLES ##

plot_density <- function(var, filtered = T){
  
  data <- if(filtered){mcmc_samples_filtered_burned}else{mcmc_samples_filtered}
  
  if(var != 'likelihood'){
    data <- data %>% select(!likelihood)
  }
  
  var_label <- gsub('_rate_','_rate\n', var)
  var_label <- gsub('_spline_','_spline\n', var_label)
  
  var_label <- if(var=='Reff'){'\nReff (calculated after)'}else{
    if(var=='init_infected'){'Initial infected (log10)'}else{var_label}
  }
  
  epid_pars_joining <- epid_pars %>% 
    mutate(true_value = value) %>% 
    filter(name == var) %>% 
    select(epidemic, epidemic_of_season, season, subtype, name, true_value)
  
  if(!grepl('imd_spline', var)){
    
    var_column <- colnames(data)[substr(colnames(data), 1, nchar(colnames(data)) - 7) == var]
    
    data %>%
      select(iteration, epidemic, chain, !!!syms(var_column)) %>%
      pivot_longer(!c(iteration,epidemic,chain)) %>%
      mutate(epidemic_of_season = as.numeric(substr(name, nchar(name), nchar(name))),
             name = substr(name, 1, nchar(name) - 7)) %>% 
      left_join(epid_pars_joining, by = c('epidemic','epidemic_of_season','name')) %>% 
      filter(!is.na(season)) %>% # remove non-real epidemics
      mutate(true_value_NA = case_when(is.na(true_value) ~ T, 
                                       T ~ F)) %>% 
      mutate(true_value = case_when(is.na(true_value) ~ 0, T ~ true_value)) %>% 
      mutate(epidemic = paste0("Season ", epidemic),
             subtype_season = paste0(season, ": ", subtype)) %>% 
      ggplot() +
      geom_density(aes(x = value, fill = as.factor(chain), group = chain), alpha = 0.4) +
      geom_vline(aes(xintercept = true_value, alpha = true_value_NA), lty=2) +
      scale_alpha_manual(values = c(1,0)) +
      theme_lg() + labs(y = ifelse(var=='R0', 'R0 (calculated after)', var)) +
      facet_grid(.~subtype_season, scales = 'free') +
      theme(legend.position = 'none') +
      labs(y = var_label)
    
  }else{
    
    epid_pars_joining <- epid_pars_joining %>% 
      select(!c(subtype, epidemic_of_season)) %>% unique()
    
    data %>% select(iteration, epidemic, chain, !!!syms(var)) %>%
      pivot_longer(!c(iteration,epidemic,chain)) %>%
      left_join(epid_pars_joining, by = c('epidemic','name')) %>% 
      filter(!is.na(season)) %>% # remove non-real epidemics
      mutate(true_value_NA = case_when(is.na(true_value) ~ T, 
                                       T ~ F)) %>% 
      mutate(true_value = case_when(is.na(true_value) ~ 0, T ~ true_value)) %>% 
      mutate(epidemic = paste0("Season ", epidemic)) %>% 
      ggplot() +
      geom_density(aes(x = value, fill = as.factor(chain), group = chain), alpha = 0.4) +
      geom_vline(aes(xintercept = true_value, alpha = true_value_NA), lty=2) +
      scale_alpha_manual(values = c(1,0)) +
      theme_lg() + labs(y = ifelse(var=='R0', 'R0 (calculated after)', var)) +
      facet_grid(.~season, scales = 'free') +
      theme(legend.position = 'none') +
      labs(y = var_label)
    
  }
  
}