

#### SETUP ####
.args <- #if (interactive()) c(
  c(file.path(Dir_Local, Subdir_Data, "inputs", "subtype_years.rds")
  ) #else commandArgs(trailingOnly = TRUE)


## SUBTYPE-SEASONS
subtype_seasons <- readRDS(.args[1])



subtype_seasons_2023 <- subtype_seasons %>%
  filter(season == "2023/24")

subtype_seasons_2023 %>%
  group_by(subtype) %>%
  summarize(value = sum(value))

subtype_seasons_2023_AH1N1 <- subtype_seasons_2023 %>%
  filter(subtype == "AH1N1")

subtype_seasons_2023_AH3N2 <- subtype_seasons_2023 %>%
  filter(subtype == "AH3N2")

ggplot(data = subtype_seasons_2023_AH1N1) +
  geom_line(mapping = aes(x = date_formatted, y = proportion)) +
  geom_point(data = subtype_seasons_2023_AH1N1 %>% filter(total_flu > 50),
            mapping = aes(x = date_formatted, y = proportion))

ggplot(data = subtype_seasons_2023_AH3N2) +
  geom_line(mapping = aes(x = date_formatted, y = proportion))



subtype_seasons_2024 <- subtype_seasons %>%
  filter(season == "2024/25")

subtype_seasons_2024 %>%
  group_by(subtype) %>%
  summarize(value = sum(value))

subtype_seasons_2024_AH1N1 <- subtype_seasons_2024 %>%
  filter(subtype == "AH1N1")

subtype_seasons_2024_B <- subtype_seasons_2024 %>%
  filter(subtype == "B")

ggplot(data = subtype_seasons_2024_AH1N1) +
  geom_line(mapping = aes(x = date_formatted, y = proportion))

ggplot(data = subtype_seasons_2024_B) +
  geom_line(mapping = aes(x = date_formatted, y = proportion))



subtype_seasons_2025 <- subtype_seasons %>%
  filter(season == "2025/26")

subtype_seasons_2025 %>%
  group_by(subtype) %>%
  summarize(value = sum(value))

subtype_seasons_2025_AH3N2 <- subtype_seasons_2025 %>%
  filter(subtype == "AH3N2")

ggplot(data = subtype_seasons_2025_AH3N2) +
  geom_line(mapping = aes(x = date_formatted, y = proportion))


