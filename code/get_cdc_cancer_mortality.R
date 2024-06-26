#### Read in cancer mortality rates from CDC WONDER

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


###################
# Get data from datasets not broken out

#import data for all non-white race categories

nw_rates_us <- read.table(
  "Data/cancer_mortality_non_white2020.txt",
  sep="\t", header=TRUE)

write.csv(nw_rates_us, "Data/cancer_mortality_2020_raw.csv")

cols = c('deaths', 'population', 'crude_rate')

nw_rates_us_clean <- nw_rates_us  %>%
  clean_names %>%
  dplyr::select(-notes) %>%
  mutate(across(all_of(cols), ~gsub("Suppressed", 0, .))) %>%
  mutate(across(all_of(cols), ~gsub("Not Applicable", NA, .))) %>%
  mutate(across(all_of(cols), as.numeric)) %>%
  mutate(race_code = case_when(race_code == "1002-5" ~ "AIAN",
                               race_code == "2054-5" ~ "B-AA",
                               race_code == "2106-3" ~ "WH",
                               race_code == "2131-1" ~ "OTHER",
                               TRUE ~ race_code))

#import dataset with white hispanic/non-hispanic 

wh_rates_us <- read.table(
  "Data/cancer_mortality_wh_eth_2020.txt",
  sep="\t", header=TRUE)

write.csv(wh_rates_us, "Data/cancer_mortality_white_hisp2020_raw.csv")

wh_rates_us_clean <- wh_rates_us  %>%
  clean_names %>%
  dplyr::select(-notes) %>%
  mutate(across(all_of(cols), ~gsub("Suppressed", 0, .))) %>%
  mutate(across(all_of(cols), ~gsub("Not Applicable", NA, .))) %>%
  mutate(across(all_of(cols), as.numeric)) %>%
  mutate(race_code = ifelse(race_code == "2106-3", "WH", race_code)) %>%
  mutate(eth_group = ifelse(ethnicity_code == '2135-2', "H", "NH")) %>% # make hispanic / non-hispanic group
  group_by(leading_cancer_sites, sex, eth_group) %>% 
  mutate(deaths = sum(deaths, na.rm = T)) %>%
  dplyr::select(-ethnicity_code) %>%
  distinct(leading_cancer_sites, sex, eth_group, .keep_all=T)

# Combine datasets

rates_us_comb <- bind_rows(nw_rates_us_clean, wh_rates_us_clean)

write.csv(rates_us_comb, "Data/us_2020_cancer_mortality_comb.csv")


###############################
#data cleaning: getting cancer sites aggregated by certain groups 

assoc_cancer <- read_csv("Data/assoc_cancer.csv") %>%
  dplyr::select(-c(1))

###########################
# combine the rates dataset with the cancer groupings

rates_us_comb2 <- left_join(rates_us_comb, assoc_cancer, relationship = 'many-to-many') %>%
  group_by(associated_cancer_types, sex, race, ethnicity) %>%
  mutate(across(all_of(cols), ~sum(., na.rm=T))) %>%
  mutate(crude_rate = (deaths/population)*100000) %>%
  distinct(associated_cancer_types, sex, race, ethnicity, .keep_all=T) %>%
  dplyr::select(-starts_with('leading_cancer_site')) %>%
  relocate(associated_cancer_types, .before = c(1)) %>%
  arrange(associated_cancer_types, race, sex)

write.csv(rates_us_comb2, "Data/incidence_by_cancer_type.csv")