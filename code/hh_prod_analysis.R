## Household product purchasing analysis

### Preliminaries
rm(list=ls(all=TRUE))
options(scipen=999)

library(data.table)
library(tidyverse)
library(readxl)
library(openxlsx)

setwd('C:/Users/atheisin/OneDrive - Environmental Protection Agency (EPA)/Documents/GitHub/ubiquitous_carcinogens')

carc.list <- data.table::fread('./output/merged_carcinogen_list.csv')

## Intermediate analysis files from Z. Stanfield (US EPA, ORD) based on his 2021 EHS paper's data work.
all.hhs <- data.table::fread('./raw_data/stanfield_reanalysis/all_carcinogens_freqItemsets_allSingleChems.csv')
afam.hhs <- data.table::fread('./raw_data/stanfield_reanalysis/AfrAmerican_carcinogens_freqItemsets_allSingleChems.csv')
asian.hhs <- data.table::fread('./raw_data/stanfield_reanalysis/Asian_carcinogens_freqItemsets_allSingleChems.csv')
hisp.hhs <- data.table::fread('./raw_data/stanfield_reanalysis/Hispanic_carcinogens_freqItemsets_allSingleChems.csv')
lower.hhs <- data.table::fread('./raw_data/stanfield_reanalysis/lower_carcinogens_freqItemsets_allSingleChems.csv')
midlower.hhs <- data.table::fread('./raw_data/stanfield_reanalysis/mid-lower_carcinogens_freqItemsets_allSingleChems.csv')
cbage.hhs <- data.table::fread('./raw_data/stanfield_reanalysis/CBage_carcinogens_freqItemsets_allSingleChems.csv')

# Merge all together.
together <- merge(all.hhs[, .(itemset, all_support = support)],
                  afam.hhs[, .(itemset, aframer_support = support)],
                  by = 'itemset',
                  all.x = T) %>%
  merge(.,
        asian.hhs[, .(itemset, asian_support = support)],
        by = 'itemset',
        all.x = T) %>%
  merge(.,
        hisp.hhs[, .(itemset, hisp_support = support)],
        by = 'itemset',
        all.x = T) %>%
  merge(.,
        lower.hhs[, .(itemset, lowinc_support = support)],
        by = 'itemset',
        all.x = T) %>%
  merge(.,
        midlower.hhs[, .(itemset, midlowinc_support = support)],
        by = 'itemset',
        all.x = T) %>%
  merge(.,
        cbage.hhs[, .(itemset, fem_cbage_support = support)],
        by = 'itemset',
        all.x = T)

# Rename vars, pull in CAS for future linking.
together.ratios <- merge(together[, .(Chemical = gsub('\\}','',gsub('\\{','',itemset)),
                                `Race/Eth: Non-hispanic black` = aframer_support/all_support,
                                `Race/Eth: Asian` = asian_support/all_support,
                                `Race/Eth: Hispanic` = hisp_support/all_support,
                                `Income: <$15k` = lowinc_support/all_support,
                                `Income: $15k-30k` = midlowinc_support/all_support,
                                `Female, childbearing age` = fem_cbage_support/all_support)
                                ][Chemical == 'Quartz-alpha (SiO2)', Chemical := 'Quartz (SiO2)'
                                ][Chemical == 'Coconut oil acid/Diethanolamine condensate (2:1)',
                                  Chemical := 'Amides, coco, N,N-bis(hydroxyethyl)'
                                ][Chemical == 'Carbaryl', 
                                Chemical := '1-Naphthalenol, 1-(N-methylcarbamate)'],
                         carc.list[, .(CASRN, NAME)],
                         by.x = 'Chemical',
                         by.y = 'NAME')

fwrite(together.ratios, './output/hh_prod_summary.csv')
