## loading in IMD- and age-stratified contact matrices ##

#### SETUP ####
.args <- if (interactive()) c(
  file.path(Dir_Local, Subdir_Data, "contact_matrix", "fitted_matrs_balanced.csv"),
  file.path(Dir_Local, Subdir_Data, "inputs", "contact_matrix.rds")
) else commandArgs(trailingOnly = TRUE)

## read
cm <- read_csv(.args[1])

## summarise
group_vars <- c('p_age_group','c_age_group','p_imd_q','c_imd_q')

summ_matr <- cm %>% 
  group_by(!!!syms(group_vars)) %>% 
  summarise(n = mean(n))


## msw temporary work around
ages_old <- c(0,5,12,18,26,35,50,70,80)
ages_old_names <- paste0(ages_old, '-', lead(ages_old) - 1)
ages_old_names[length(ages_old_names)] <- '80+'

ages_names <- data.table(old = ages_old_names,
                         new = Age_Groups_Names)

summ_matr <- summ_matr %>%
  left_join(ages_names, by = c("p_age_group" = "old")) %>%
  dplyr::rename(c("p_age_group_new" = "new")) %>%
  left_join(ages_names, by = c("c_age_group" = "old")) %>%
  dplyr::rename(c("c_age_group_new" = "new")) %>%
  ungroup()

summ_matr <- summ_matr %>%
  dplyr::select(-c("p_age_group", "c_age_group")) 

summ_matr <- summ_matr %>%
  dplyr::rename(c("p_age_group" = "p_age_group_new", "c_age_group" = "c_age_group_new")) %>%
  relocate(c("p_age_group", "c_age_group"))


## save
write_rds(summ_matr, .args[2])


