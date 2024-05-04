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
ipr <- data.table::dcast(rbind(unique(iarc.clean), unique(roc.clean), unique(p65.clean)),
                         NAME + CASRN + DTXSID ~ LIST, length)

# manual correction of duplicates
dupes <- ipr[duplicated(ipr$CASRN)]
all.dupes <- ipr[CASRN %in% dupes$CASRN] #use this list for corrections below

ipr <- ipr[!(NAME %in% 
                c('(+/-)-1,2-Propylene oxide',
                  '4,5-Dihydro-2-mercaptoimidazole',
                  'Coconut oil acid/Diethanolamine condensate (2:1)',
                  'Dibenz(a,h)anthracene',
                  'Benzo(a)pyrene',
                  'Benzo(k)fluoranthene',
                  'Indeno(1,2,3-cd)pyrene',
                  'Quartz-alpha (SiO2)',
                  '4-Vinyl-1-cyclohexene diepoxide',
                  'Polychlorinated biphenyls (containing 60 or more percent chlorine by molecular weight)', #duplicate
                  'Wood dust',
                  'Tobacco Smoking (see Tobacco-Related Exposures)',
                  'Strong inorganic acid mists containing sulfuric acid',
                  'Soots',
                  'Smokeless Tobacco (see Tobacco-Related Exposures)',
                  'Environmental Tobacco Smoke (see Tobacco-Related Exposures)',
                  'Diesel engine exhaust',
                  'Coke oven emissions',
                  'Alcoholic beverages, when associated with alcohol abuse',
                  'Alcoholic beverages'
                  ))
            ][NAME == '1,2-Propylene oxide', IARC2B := 1L
              ][NAME == '1,2-Propylene oxide', Prop65 := 1L
               ][NAME == 'Ethylene thiourea', Prop65 := 1L
                 ][NAME == 'Amides, coco, N,N-bis(hydroxyethyl)', Prop65 := 1L
                   ][NAME == 'Dibenz[a,h]anthracene', ROC_RAHC := 1L
            ][NAME == 'Benzo[a]pyrene', ROC_RAHC := 1L
              ][NAME == 'Benzo[k]fluoranthene', ROC_RAHC := 1L
                ][NAME == 'Indeno[1,2,3-cd]pyrene', ROC_RAHC := 1L
                  ][NAME == 'Quartz (SiO2)', Prop65 := 1L
                    ][NAME == 'Quartz (SiO2)', ROC_Known := 1L
            ][NAME == '4-Vinyl-1-cyclohexene dioxide', IARC2B := 1L
              ][NAME == '4-Vinyl-1-cyclohexene dioxide', Prop65 := 1L
               ][NAME == 'Wood Dust', Prop65 := 1L
                 ][NAME == 'Tobacco smoke', ROC_Known := 1L
                   ][NAME == 'Strong Inorganic Acid Mists Containing Sulfuric Acid', Prop65 := 1L
            ][NAME == 'Soots, tars, and mineral oils (untreated and mildly treated oils and used engine oils)', ROC_Known := 1L
              ][NAME == 'Diesel Exhaust Particulates', Prop65 := 1L
                ][NAME == 'Coke Oven Emissions', Prop65 := 1L
                  ][NAME == 'Alcoholic Beverage Consumption', Prop65 := 1L
                    ][NAME == 'Tobacco, oral use of smokeless products', ROC_Known := 1L]

rm(dupes, all.dupes)

# Fix some encoding issues.
ipr <- ipr[CASRN == '120-80-9', NAME := '1,2-Benzenediol'
    ][CASRN == '13327-32-7', NAME := 'Beryllium hydroxide'
      ][CASRN == '25962-77-0', NAME := 'N,N-Dimethyl-N′-[5-[2-(5-nitro-2-furanyl)ethenyl]-1,3,4-oxadiazol-2-yl]methanimidamide'
        ][CASRN == '1204332-00-2', NAME := 'Trim VX'
          ][NAME == 'Diesel Exhaust Particulates', CASRN := 'DIESEL'
          ]

# Functional uses: pull bulk data from CompTox / Chem expo
# Download comptox ID database if necessary (too large to include on GH)
if('ChemExpo_bulk_functional_uses.csv' %in% list.files('./raw_data/')){
  funct <- data.table::fread('./raw_data/ChemExpo_bulk_functional_uses.csv')
} else {
  # Data download from cluster below and save locally.
  download.file('https://comptox.epa.gov/chemexpo/dl_functional_uses/',
                './raw_data/ChemExpo_bulk_functional_uses.csv')
  funct <- data.table::fread('./raw_data/ChemExpo_bulk_functional_uses.csv')
}

funct.keep <- unique(funct[`Curated CAS` %in% ipr[!(CASRN == '')]$CASRN,
                    .(`Curated CAS`, `Harmonized Functional Use`)
                    ][!(`Harmonized Functional Use` == '')]
                    )[order(`Curated CAS`, `Harmonized Functional Use`)
                      ][, .(lapply(`.SD`, paste, collapse = '; ')),
                        by = 'Curated CAS'
                        ][, V1 := as.character(V1)
                        ][, function_count := stringr::str_count(V1, pattern = ';') + 1]

data.table::setcolorder(funct.keep, c(1,3,2))
names(funct.keep) <- c('CAS', 'Function_count', 'Functions')
rm(funct)

# Product uses: pull bulk data from CompTox / Chem expo
if('ChemExpo_bulk_composition_chemicals.csv' %in% list.files('./raw_data/')){
  comp <- data.table::fread('./raw_data/ChemExpo_bulk_composition_chemicals.csv')
} else {
  # Data download from cluster below and save locally.
  download.file('https://comptox.epa.gov/chemexpo/dl_co_chemicals/',
                './raw_data/ChemExpo_bulk_composition_chems.zip',
                method = 'curl')
  unzip(zipfile = './raw_data/ChemExpo_bulk_composition_chems.zip',
        exdir = './raw_data/')
  comp <- data.table::fread('./raw_data/ChemExpo_bulk_composition_chemicals.csv')
  file.remove('./raw_data/ChemExpo_bulk_composition_chems.zip')
}

comp.keep <- unique(comp[`Curated CAS` %in% ipr[!(CASRN == '')]$CASRN,
                           .(`Curated CAS`, `PUC General Category`)
][!(`PUC General Category` == '')]
)[order(`Curated CAS`, `PUC General Category`)
][, .(lapply(`.SD`, paste, collapse = '; ')),
  by = 'Curated CAS'
][, V1 := as.character(V1)
][, PUC_count := stringr::str_count(V1, pattern = ';') + 1]

data.table::setcolorder(comp.keep, c(1,3,2))
names(comp.keep) <- c('CAS', 'PUC_count', 'PUCs')
rm(comp)

# Merge in the function/PUC
ipr <- merge(ipr,
             funct.keep,
             by.x = 'CASRN',
             by.y = 'CAS',
             all.x = T) %>%
  merge(.,
        comp.keep,
        by.x = 'CASRN',
        by.y = 'CAS',
        all.x = T)

data.table::fwrite(ipr, './output/merged_carcinogen_list.csv')
