

plot_trace <- function(var, filtered = F){
  
  data <- if(filtered){mcmc_samples_filtered_burned}else{mcmc_samples_filtered} 
  
  if(var != 'likelihood'){
    data <- data %>% select(!likelihood)
  }
  
  var_label <- gsub('_rate_','_rate\n', var)
  var_label <- gsub('_spline_','_spline\n', var_label)
  
  var_label <- if(var=='Reff'){'\nReff (calculated after)'}else{
    if(var=='init_infected'){'Initial infected (log10)'}else{var_label}
  }
  if(var=='likelihood'){var_label <- 'Log-likelihood'}
  
  use_colors <- n_distinct(data$job) <= 10
  
  epid_pars_joining <- epid_pars %>% 
    mutate(true_value = value) %>% 
    filter(name == var) %>% 
    select(epidemic, epidemic_of_season, season, subtype, name, true_value)
  
  if(!grepl('imd_spline|likelihood', var)){
    
    var_column <- colnames(data)[substr(colnames(data), 1, nchar(colnames(data)) - 7) == var]
    
    p <- data %>%
      select(iteration, epidemic, chain, job, !!!syms(var_column)) %>% 
      pivot_longer(!c(iteration,epidemic,chain, job)) %>% drop_na() %>% 
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
      geom_hline(aes(yintercept = true_value, alpha = true_value_NA), lty=2) +
      scale_alpha_manual(values = c(1,0)) +
      facet_grid(.~subtype_season) +
      theme_lg() + labs(y = var_label) +
      theme(legend.position = 'none') + labs(x = 'Iteration (1000s)')
    
  }else{
    
    if(var == 'likelihood'){
      
      epid_pars_joining <- epid_pars %>% 
        select(epidemic, season) %>% unique()
      
      p <- data %>% select(iteration, epidemic, chain, job,
                           likelihood) %>%
        pivot_longer(!c(iteration,epidemic,chain,job)) %>%
        left_join(epid_pars_joining, by = c('epidemic')) %>% 
        filter(!is.na(season)) %>% # remove non-real epidemics
        mutate(epidemic = paste0("Season ", epidemic)) %>% 
        ggplot() +
        geom_hline(yintercept = 0, lty=2) +
        facet_grid(.~season, scales = 'free') +
        theme_lg() + labs(y = var_label) + 
        theme(legend.position = 'none') + labs(x = 'Iteration (1000s)')
      
      p_filt <- data %>% select(iteration, epidemic, chain, job,
                                likelihood) %>%
        pivot_longer(!c(iteration,epidemic,chain,job)) %>%
        left_join(epid_pars_joining, by = c('epidemic')) %>% 
        filter(!is.na(season)) %>% # remove non-real epidemics
        mutate(epidemic = paste0("Season ", epidemic)) %>% 
        filter(iteration > max(iteration)/4) %>% 
        ggplot() +
        geom_hline(yintercept = 0, lty=2) +
        facet_grid(.~season, scales = 'free') +
        theme_lg() + labs(y = var_label) + 
        theme(legend.position = 'none') + labs(x = 'Iteration (1000s)')
      
      p2 <- data %>% select(iteration, epidemic, chain,
                            gp_likelihood, hospital_likelihood,
                            subtype_likelihood) %>%
        pivot_longer(!c(iteration,epidemic,chain)) %>%
        left_join(epid_pars_joining, by = c('epidemic')) %>% 
        filter(!is.na(season)) %>% # remove non-real epidemics
        mutate(epidemic = paste0("Season ", epidemic)) %>% 
        group_by(iteration, epidemic, season, name) %>% 
        summarise(value = mean(value, na.rm = T)) %>% 
        filter(value != -Inf) %>% 
        mutate(name = gsub('_',' ', name)) %>% 
        ggplot() +
        geom_bar(aes(iteration/1000, y = value, fill = name),
                 position = 'fill', stat = 'identity', width = 0.5) +
        # geom_line(aes(iteration/1000, y = -value/min(value), 
        # group = name, lty = name)) +
        scale_fill_manual(values = c('#157145','#57A773','#9BD1E5')) +
        facet_grid(.~season, scales = 'free') +
        theme_minimal() + labs(y = 'Proportion of log-likelihood', fill = '') + 
        scale_x_continuous(expand = expansion(c(0,0))) + 
        scale_y_continuous(expand = expansion(c(0,0))) + 
        # scale_linetype_manual(values = c(2,3,1)) + 
        theme(legend.position = 'bottom') +
        labs(x = 'Iteration (1000s)'); suppressWarnings(print(p2))
      
    }else{
      
      epid_pars_joining <- epid_pars_joining %>% 
        select(!c(subtype, epidemic_of_season)) %>% unique()
      
      p <- data %>% select(iteration, epidemic, chain, job, !!!syms(var)) %>%
        pivot_longer(!c(iteration,epidemic,chain, job)) %>%
        left_join(epid_pars_joining, by = c('epidemic','name')) %>% 
        filter(!is.na(season)) %>% # remove non-real epidemics
        mutate(true_value_NA = case_when(is.na(true_value) ~ T, 
                                         T ~ F)) %>% 
        mutate(true_value = case_when(is.na(true_value) ~ 0, T ~ true_value)) %>% 
        mutate(epidemic = paste0("Season ", epidemic)) %>% 
        ggplot() +
        geom_hline(aes(yintercept = true_value, alpha = true_value_NA), lty=2) +
        scale_alpha_manual(values = c(1,0)) +
        facet_grid(.~season) +
        theme_lg() + labs(y = var_label) +
        theme(legend.position = 'none') + labs(x = 'Iteration (1000s)')
      
    }
    
  }
  
  if(var=='init_infected'){
    p <- p + scale_y_log10()
  }
  
  if(use_colors){
    
    p <- p +
      geom_line(aes(x = iteration/1000, y = value, col = as.factor(job), group = chain),
                alpha = 0.4)
    
  }else{
    
    p <- p +
      geom_line(aes(x = iteration/1000, y = value, group = chain), alpha = 0.2)
    
  }
  
  if(var == 'likelihood'){
    
    p_filt <- p_filt +
      geom_line(aes(x = iteration/1000, y = value, col = as.factor(job), group = chain),
                alpha = 0.4)
    
    p + p_filt + suppressWarnings(print(p2)) + plot_layout(nrow = 3)
    
  }else{
    
    p
    
  }
  
}


