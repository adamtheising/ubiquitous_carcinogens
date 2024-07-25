## clean known carcinogens data

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

## Read data

kc_data <- read_excel("Copy of exposure_pathway_summary.xlsx", sheet = 2)

## Clean up names and cols

kc_data_cl <- kc_data %>%
  dplyr::select(1:8) %>%
  clean_names() %>%
  mutate(across(everything(), tolower))

saveRDS(kc_data_cl, "Data/kc_data_cl.rds")

## Perform analysis: expand dataset to analyze different exposure pathways and types

kc_data_cl2 <- kc_data_cl %>% 
 separate_longer_delim(cols = major_product_categories, delim = ";") %>%
  mutate(across(where(is.character), str_squish)) %>%
  filter(major_product_categories != "") 

saveRDS(kc_data_cl2, "Data/kc_data_by_prod.rds")

kc_data_cl3 <- kc_data_cl %>% 
  separate_longer_delim(cols = main_exposure_pathways, delim = ";") %>%
  mutate(across(where(is.character), str_squish)) %>%
  filter(main_exposure_pathways != "" & main_exposure_pathways != " ")

saveRDS(kc_data_cl3, "Data/kc_data_by_path.rds")
  
