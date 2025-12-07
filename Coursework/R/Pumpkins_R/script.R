# load packages
library(dplyr)
library(readr)
library(ggplot2)


# 1) read dataset
pumpkins <- read_csv("pumpkins.csv")

# look at the data
head(pumpkins)


# 2) make a year column from id (first 4 characters)
pumpkins$year <- substr(pumpkins$id, 1, 4)

# find the heaviest pumpkin 
heaviest <- pumpkins %>%
  filter(weight_lbs == max(weight_lbs, na.rm = TRUE))

print(heaviest)


# 3) change lbs to kg
pumpkins <- pumpkins %>%
  mutate(weight_kg = weight_lbs * 0.453592)

head(pumpkins$weight_kg)


# 4) make weight groups
pumpkins <- pumpkins %>%
  mutate(
    weight_class = ifelse(
        weight_kg < 200, "light",
        ifelse(weight_kg < 500, "medium", "heavy")
    )
  )

# check it
table(pumpkins$weight_class)


# 5) plot: estimated vs actual weight
plot1 <- ggplot(pumpkins, 
                aes(x = est_weight, y = weight_kg, color = weight_class)) +
  geom_point(alpha = 0.6) + 
  scale_color_manual(values = c("light" = "yellow", 
                                "medium" = "orange",
                                "heavy" = "red")) +
  labs(
    x = "Estimated weight (lbs)",
    y = "Actual weight (kg)",
    color = "Weight class",
    title = "Estimated vs actual pumpkin weight"
  ) +
  theme_minimal()

plot1

# save plot
ggsave("plot_estimated_vs_actual.png", plot1)


# 6) keep three countries
pumpkins_three <- pumpkins %>%
  filter(country %in% c("Canada", "Japan", "Italy"))

# check it
table(pumpkins_three$country)

# save filtered data 
write_csv(pumpkins_three, "pumpkins_filtered.csv")


# 7) mean weight by country 
mean_country <- pumpkins_three %>%
  group_by(country) %>%
  summarise(mean_weight = mean(weight_kg, na.rm = TRUE))

print(mean_country)

# mean weight by country and variety
mean_variety <- pumpkins_three %>%
  group_by(country, variety) %>%
  summarise(mean_weight = mean(weight_kg, na.rm = TRUE)) %>%
  arrange(mean_weight)

print(mean_variety)


# 8) boxplot for the three countries
box1 <- ggplot(pumpkins_three, aes(x = country, y = weight_kg)) +
  geom_boxplot(fill = "skyblue") +
  labs(
    x = "Country",
    y = "Pumpkin weight (kg)",
    title = "Weight distribution of pumpkins in three countries"
  ) +
  theme_minimal()

box1

# save boxplot
ggsave("plot_boxplot_countries.png", box1,)


# 9) facet plot by pumpkin variety
facet1 <- ggplot(pumpkins_three, aes(x = weight_kg)) +
  geom_histogram(binwidth = 100, fill = "skyblue") +
  facet_wrap(~ variety) +
  labs(
    x = "Pumpkin weight (kg)",
    y = "Count",
    title = "Pumpkin weight distribution by variety"
  ) +
  theme_minimal()

facet1

# save facet plot
ggsave("plot_facet_variety.png", facet1,)