## Create a unified list of carcinogenic chemical lists
## February 2024
## Contact: theising.adam@epa.gov
rm(list = ls())
gc()
options(scipen=999)

# Packages
library(tidyverse)
library(data.table)
library(openxlsx)

setwd('C:/Users/atheisin/OneDrive - Environmental Protection Agency (EPA)/Documents/GitHub/ubiquitous_carcinogens/')

# Download comptox ID database if necessary (too large to include on GH)
if('DSSTox_Identifiers_and_CASRN_2021r1.csv' %in% list.files('./raw_data/')){
  comptox.dta <- data.table::fread('./raw_data/DSSTox_Identifiers_and_CASRN_2021r1.csv')
} else {
  # Data download from cluster below and save locally. 
  comptox.dta <- data.table::fread('https://clowder.edap-cluster.com/files/616dd943e4b0a5ca8aeea69d/blob')
  data.table::fwrite(comptox.dta, './raw_data/DSSTox_Identifiers_and_CASRN_2021r1.csv')
}

###############################
# Start by merging IARC data
iarc.clean <- data.table::rbindlist(
  lapply(paste0('./raw_data/', list.files('./raw_data/', pattern = '_IARC')), 
                     function(x) {
                       data.table::fread(x) %>%
                         dplyr::mutate(LIST = gsub('.*Tox_(.+).csv',
                                                   '\\1',
                                                   x))
                       })
)[, .(NAME = `PREFERRED NAME`,
      CASRN,
      LIST,
      DTXSID = gsub('.*details/','',DTXSID))] # 409 (408 unique) chemical entries here.

###########################
# Next, RoC data
roc <- data.table::fread('./raw_data/NTP_ROC15.csv',
                         header = T)[, 1:3
                         ][, LIST := paste0('ROC_', Listing)
                           ][, CASRN := trimws(CASRN)]

# ROC data that is already properly formatted.
roc.match <- merge(roc[!is.na(CASRN)],comptox.dta, 
                   by.x = 'CASRN', by.y = 'casrn'
                   )[, NAME := NULL
                     ][,.(NAME = preferredName,
                          CASRN,
                          LIST,
                          DTXSID = dtxsid)]

# ROC data that needs formatting
roc.nomatch <- roc[!(CASRN %in% roc.match$CASRN)]
roc.crosswalk <- data.table::fread('./raw_data/crosswalks/roc_crosswalk.csv' # This was manually created.
                                   )[is.na(bio_flag) # Drop viruses, etc.
                                     ][is.na(radiation_flag) #
                                       ]#[is.na(alc_tob_flag)] # Drop alcohol & tobacco consumption
                                       
roc.nomatch <- merge(roc.nomatch, roc.crosswalk,
      by.x = c('NAME'),
      by.y = c('NAME')
      )[, .(NAME,
            CASRN = comptox_CASRN,
            LIST,
            DTXSID)]

# Bind cleaned RoC lists together
roc.clean <- rbind(roc.match,
                   merge(roc.nomatch, comptox.dta[, .(dtxsid, preferredName)],
                   by.x = 'DTXSID', by.y = 'dtxsid', all.x = T
                   )[, nametemp := preferredName
                     ][is.na(preferredName), nametemp := NAME
                       ][, .(NAME = nametemp,
                             CASRN,
                             LIST,
                             DTXSID)]
                   )

rm(roc, roc.crosswalk, roc.match, roc.nomatch)

##############################
# Next look at prop 65 data

# data from comptox
p65.comptox <- data.table::fread('./raw_data/CompTox_P65.csv'
)[, .(NAME = `PREFERRED NAME`,
      CASRN,
      DTXSID = gsub('.*details/','',DTXSID))][, LIST := 'Prop65']

# newer data from CA OEHHA website
p65.oehha <- data.table::fread('./raw_data/OEHHA_P65.csv'
)[`Type of Toxicity` %like% 'cancer'
][`CAS No.` %in% c('--','---'), 
  `CAS No.` := ''
]#[, `CAS No.` := gsub('/','-',`CAS No.`)]

# Data that is properly formatted
p65.match <- p65.comptox[CASRN %in% p65.oehha$`CAS No.`]

# Data that needs formatting
p65.nomatch <- p65.oehha[!(`CAS No.` %in% p65.match$CASRN)]
p65.crosswalk <- data.table::fread('./raw_data/crosswalks/p65_crosswalk.csv' # This was manually created.
)[is.na(bio_flag) # Drop viruses, etc.
][is.na(radiation_flag) #
]#[is.na(alc_tob_flag)] # Drop alcohol & tobacco consumption

p65.nomatch <- merge(p65.nomatch, p65.crosswalk,
                     by.x = c('Chemical'),
                     by.y = c('Chemical')
)[, .(NAME = Chemical,
      CASRN = comptox_CAS,
      LIST ='Prop65',
      DTXSID)]

# Bind cleaned RoC lists together
p65.clean <- rbind(p65.match,
                   merge(p65.nomatch, comptox.dta[, .(dtxsid, preferredName)],
                         by.x = 'DTXSID', by.y = 'dtxsid', all.x = T
                   )[, nametemp := preferredName
                   ][is.na(preferredName), nametemp := NAME
                   ][, .(NAME = nametemp,
                         CASRN,
                         LIST,
                         DTXSID)]
)

rm(p65.comptox,p65.crosswalk,p65.match,p65.nomatch,p65.oehha)

# merge and save to csv
ipr <- data.table::dcast(rbind(iarc.clean, roc.clean, p65.clean),
                         NAME + CASRN + DTXSID ~ LIST, length)

data.table::fwrite(ipr, './output/merged_carcinogen_list.csv')
