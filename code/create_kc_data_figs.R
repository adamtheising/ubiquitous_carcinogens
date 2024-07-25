# Generate figures on carcinogens data

if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  readr,
  tidyverse,
  data.table, # data wrangling
  dplyr, # data wrangling
  lubridate, # Date object handling
  modelsummary, # regression table generation
  stringr, # string manipulation
  magrittr,
  tidycensus, #Census data
  MASS , #For regressions and modeling
  dplyr,
  readxl,
  janitor
)

getwd()

setwd("C:/Users/tbardot/OneDrive - Environmental Protection Agency (EPA)/Documents/Ubiq. Carcinogens")

kc_data_by_prod <- readRDS("Data/kc_data_by_prod.rds")

kc_prod_count <- count(kc_data_by_prod, major_product_categories)


# Subset to product categories that appear at least twice

prod_hist <- ggplot(kc_prod_count, aes(x = major_product_categories, y = n)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(name = "Major Product Categories", guide = guide_axis(angle=45)) + 
  geom_text(aes(label = n), vjust = -1, family = "serif",
            color="black", check_overlap = TRUE) +
  theme_minimal(base_size = 20, base_family = "serif") +
  ylab("Count")

png(filename = "Plots/majpr_prod_cat.png",
    width = 1080,height = 745, units = 'px')

print(prod_hist)

dev.off()

###########
# Make plot for exposure pathways

kc_data_by_path <- readRDS("Data/kc_data_by_path.rds") 

kc_data_by_path[c('pathway', 'population')] <- 
  str_split_fixed(kc_data_by_path$main_exposure_pathways, "\\(|\\)", 2)

kc_data_by_path <- kc_data_by_path %>% mutate(population = gsub(")", "", population))

kc_path_count <- kc_data_by_path %>%
  mutate(across(where(is.character), str_squish)) %>%
  group_by(pathway, population) %>%
  summarize(n = n())

path_hist <- ggplot(kc_path_count, aes(x = pathway, y = n, fill = population, label = n)) +
  geom_bar(stat = "identity") +
  scale_x_discrete(name = "Main Exposure Pathways", guide = guide_axis(angle=45)) + 
  geom_text(size = 8, position = position_stack(vjust = 0.5)) +
  theme_minimal(base_size =20, base_family = "serif") +
  ylab("Count")
  
  
png(filename = "Plots/main_exp_path.png",
    width = 1080,height = 745, units = 'px')

print(path_hist)

dev.off()


