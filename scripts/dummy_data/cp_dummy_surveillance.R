


dummy_surveillance <- readRDS(file.path(Dir_Local, Subdir_Data, "dummy_data", "dummy_surveillance.rds"))


dummy_surveillance_2023 <- dummy_surveillance %>%
  filter(index == 1) 

dummy_surveillance_2023_age_imd <- dummy_surveillance_2023 %>%
  group_by(week_start, age_grp, imd_quintile) %>%
  summarize(primary_care = sum(primary_care),
            secondary_care = sum(secondary_care))

dummy_surveillance_2024 <- dummy_surveillance %>%
  filter(index == 2) 

dummy_surveillance_2024_age_imd <- dummy_surveillance_2024 %>%
  group_by(week_start, age_grp, imd_quintile) %>%
  summarize(primary_care = sum(primary_care),
            secondary_care = sum(secondary_care))

dummy_surveillance_2025 <- dummy_surveillance %>%
  filter(index == 3) 

dummy_surveillance_2025_age_imd <- dummy_surveillance_2025 %>%
  group_by(week_start, age_grp, imd_quintile) %>%
  summarize(primary_care = sum(primary_care),
            secondary_care = sum(secondary_care))



plot_dummy_surveillance_2023_age_imd_prim <- ggplot(data = dummy_surveillance_2023_age_imd) +
  geom_line(mapping = aes(x = week_start, y = primary_care)) +
  # scale_y_continuous(breaks = c(0, 1, 2, 3)) +
  facet_grid(age_grp~imd_quintile) +
  theme_bw()

ggsave(filename = paste0(Dir_Local, "/output/figures/plot_dummy_surveillance_2023_age_imd_prim.png"),
       plot = plot_dummy_surveillance_2023_age_imd_prim)


plot_dummy_surveillance_2023_age_imd_sec <- ggplot(data = dummy_surveillance_2023_age_imd) +
  geom_line(mapping = aes(x = week_start, y = secondary_care)) +
  # scale_y_continuous(breaks = c(0, 1, 2, 3)) +
  facet_grid(age_grp~imd_quintile) +
  theme_bw()

ggsave(filename = paste0(Dir_Local, "/output/figures/plot_dummy_surveillance_2023_age_imd_sec.png"),
       plot = plot_dummy_surveillance_2023_age_imd_sec)



plot_dummy_surveillance_2024_age_imd_prim <- ggplot(data = dummy_surveillance_2024_age_imd) +
  geom_line(mapping = aes(x = week_start, y = primary_care)) +
  # scale_y_continuous(breaks = c(0, 1, 2, 3)) +
  facet_grid(age_grp~imd_quintile) +
  theme_bw()

ggsave(filename = paste0(Dir_Local, "/output/figures/plot_dummy_surveillance_2024_age_imd_prim.png"),
       plot = plot_dummy_surveillance_2024_age_imd_prim)


plot_dummy_surveillance_2024_age_imd_sec <- ggplot(data = dummy_surveillance_2024_age_imd) +
  geom_line(mapping = aes(x = week_start, y = secondary_care)) +
  # scale_y_continuous(breaks = c(0, 1, 2, 3)) +
  facet_grid(age_grp~imd_quintile) +
  theme_bw()

ggsave(filename = paste0(Dir_Local, "/output/figures/plot_dummy_surveillance_2024_age_imd_sec.png"),
       plot = plot_dummy_surveillance_2024_age_imd_sec)



plot_dummy_surveillance_2025_age_imd_prim <- ggplot(data = dummy_surveillance_2025_age_imd) +
  geom_line(mapping = aes(x = week_start, y = primary_care)) +
  # scale_y_continuous(breaks = c(0, 1, 2, 3)) +
  facet_grid(age_grp~imd_quintile) +
  theme_bw()

ggsave(filename = paste0(Dir_Local, "/output/figures/plot_dummy_surveillance_2025_age_imd_prim.png"),
       plot = plot_dummy_surveillance_2025_age_imd_prim)


plot_dummy_surveillance_2025_age_imd_sec <- ggplot(data = dummy_surveillance_2025_age_imd) +
  geom_line(mapping = aes(x = week_start, y = secondary_care)) +
  # scale_y_continuous(breaks = c(0, 1, 2, 3)) +
  facet_grid(age_grp~imd_quintile) +
  theme_bw()

ggsave(filename = paste0(Dir_Local, "/output/figures/plot_dummy_surveillance_2025_age_imd_sec.png"),
       plot = plot_dummy_surveillance_2025_age_imd_sec)


