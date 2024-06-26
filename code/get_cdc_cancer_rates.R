#### Read in cancer incidence rates from CDC WONDER

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



rates_state <- read.table(
  "Data/cancer_incidence_states_rates2020.txt",
  sep="\t", header=TRUE)

rates_us <- read.table(
  "Data/cancer_incidence_us.txt",
  sep="\t", header=TRUE)

us_2020 <- rates_us %>% filter(Year == '2020') %>%
  mutate(Race.Code = case_when(Race.Code == "1002-5" ~ "AIAN",
                               Race.Code == "2054-5" ~ "B-AA",
                               Race.Code == "2106-3" ~ "WH",
                               Race.Code == "2131-1" ~ "OTHER",
                               TRUE ~ Race.Code)) %>%
  clean_names

write.csv(us_2020, "Data/us_2020_cancer_rates_raw.csv")


us_2020_clean <- us_2020 %>%
  dplyr::select(-notes, -year_code) %>%
  mutate(count = ifelse(count == "Suppressed", 0, count),
         count = as.numeric(count),
         population = as.numeric(population),
         crude_rate = ifelse(count == "Suppressed", 0, count)) 

write.csv(us_2020, "Data/us_2020_cancer_rates_disaggregated.csv")

us_2020_agg <- us_2020 %>%
  group_by(race, leading_cancer_sites, sex) %>% 
  mutate(count = ifelse(count == "Suppressed", 0, count),
         count = as.numeric(count),
         count_agg = sum(count, na.rm = T),
         crude_rate = ifelse(count == "Suppressed", 0, count)) %>%
  distinct(race, leading_cancer_sites, sex, .keep_all=T) %>%
  dplyr::select(-starts_with('ethn')) %>%
  ungroup

write.csv(us_2020, "Data/us_2020_cancer_rates_aggregated.csv")

###################
# Get data from datasets not broken out

#import data for all non-white race categories

nw_rates_us <- read.table(
  "Data/lead_cancer_sites_non_white2020.txt",
  sep="\t", header=TRUE)

write.csv(nw_rates_us, "Data/lead_cancer_sites2020_raw.csv")

cols = c('count', 'population', 'crude_rate')

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
  "Data/lead_cancer_sites_white_eth_2020.txt",
  sep="\t", header=TRUE)

write.csv(wh_rates_us, "Data/lead_cancer_sites_white_hisp2020_raw.csv")

wh_rates_us_clean <- wh_rates_us  %>%
  clean_names %>%
  dplyr::select(-notes) %>%
  mutate(across(all_of(cols), ~gsub("Suppressed", 0, .))) %>%
  mutate(across(all_of(cols), ~gsub("Not Applicable", NA, .))) %>%
  mutate(across(all_of(cols), as.numeric)) %>%
  mutate(race_code = ifelse(race_code == "2106-3", "WH", race_code)) %>%
  mutate(eth_group = ifelse(ethnicity_code == '2135-2', "H", "NH")) %>% # make hispanic / non-hispanic group
  group_by(leading_cancer_sites, sex, eth_group) %>% 
  mutate(count = sum(count, na.rm = T)) %>%
  dplyr::select(-ethnicity_code) %>%
  distinct(leading_cancer_sites, sex, eth_group, .keep_all=T)

# Combine datasets

rates_us_comb <- bind_rows(nw_rates_us_clean, wh_rates_us_clean)

write.csv(rates_us_comb, "Data/us_2020_cancer_rates_comb.csv")


###############################
#data cleaning: getting cancer sites aggregated by certain groups 

lead_cancer_sites <- rates_us_comb %>%
  distinct(leading_cancer_sites) %>%
  write.csv("Data/lead_cancer_sites.csv")

assoc_cancer <- read_xlsx("exposure_pathway_summary.xlsx", sheet = 2) %>%
  pivot_longer(cols = c(2:5), names_to = NULL, values_to = "leading_cancer_sites", values_drop_na = T) %>%
  clean_names

write.csv(assoc_cancer, "Data/assoc_cancer.csv")


###########################
# combine the rates dataset with the cancer groupings

rates_us_comb2 <- left_join(rates_us_comb, assoc_cancer, relationship = 'many-to-many') %>%
  group_by(associated_cancer_types, sex, race, ethnicity) %>%
  mutate(across(all_of(cols), ~sum(., na.rm=T))) %>%
  mutate(crude_rate = (count/population)*100000) %>%
  distinct(associated_cancer_types, sex, race, ethnicity, .keep_all=T) %>%
  dplyr::select(-starts_with('leading_cancer_site')) %>%
  relocate(associated_cancer_types, .before = c(1)) %>%
  arrange(associated_cancer_types, race, sex)

write.csv(rates_us_comb2, "Data/incidence_by_cancer_type.csv")
  