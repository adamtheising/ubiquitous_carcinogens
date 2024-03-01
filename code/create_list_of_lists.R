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
                                     ][is.na(alc_tob_flag) # Drop alcohol & tobacco consumption
                                       ][is.na(radiation_flag)] # Drop radiation
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
                       ][, .(NAME,
                             CASRN,
                             LIST,
                             DTXSID)]
                   )

# 
rm(roc, roc.crosswalk, roc.match, roc.nomatch)



# Match on chems with CASRN first
ipr.matched <- data.table::as.data.table(rbind(iarc.p65, roc, fill = T)
)[!is.na(CASRN)
][, DTXSID.mode := Mode(DTXSID),
  by = CASRN
][, NAME.mode := Mode(NAME),
  by = CASRN
][, IUPAC_NAME.mode := Mode(IUPAC_NAME),
  by = CASRN]






##########################
# Next look at prop 65 data
prop65 <- F # Set to true for first pass through.
if (prop65 == T){
  p65.comptox <- data.table::fread('./raw_data/CompTox_P65.csv'
                                   )[, .(NAME = `PREFERRED NAME`,
                                         CASRN,
                                         IUPAC_NAME = `IUPAC NAME`,
                                         DTXSID)][, LIST := 'Prop65']
  p65.oehha <- data.table::fread('./raw_data/OEHHA_P65.csv'
                                 )[`Type of Toxicity` %like% 'cancer'
                                   ][`CAS No.` %in% c('--','---'), 
                                     `CAS No.` := NA_character_
                                     ]#[, `CAS No.` := gsub('/','-',`CAS No.`)]
  
  # First match on CAS Number: 523/632 matches.
  p65.comptox.casmatch <- p65.comptox[CASRN %in% p65.oehha$`CAS No.`]
  p65.oehha.nomatch <- p65.oehha[!(`CAS No.` %in% p65.comptox.casmatch$CASRN)
                                 ][, .(NAME = Chemical,
                                       CASRN = `CAS No.`)]
  
  # Unfortunately need to do some manual cleaning outside R
  # Take this over to excel for a manual match.
  intermediate.match <- rbind(p65.comptox.casmatch, p65.oehha.nomatch, fill = T)
  fwrite(intermediate.match, './raw_data/intermediate_data/p65.csv')
  
  rm(p65.comptox, p65.comptox.casmatch, p65.oehha, p65.oehha.nomatch, intermediate.match)
}

p65.intermediate <- readxl::read_excel('./raw_data/intermediate_data/p65_manualclean.xlsx')

iarc.p65 <- data.table::as.data.table(rbind(iarc.files, p65.intermediate))

Mode <- function(x) {
  ux <- unique(x[!is.na(x)])
  ux[which.max(tabulate(match(x, ux)))]
}


###########################
# Next, RoC data
roc <- data.table::fread('./raw_data/NTP_ROC15.csv',
                         header = T)[, 1:2
                                     ][, LIST := 'ROC']

# Match on chems with CASRN first
ipr.matched <- data.table::as.data.table(rbind(iarc.p65, roc, fill = T)
                                          )[!is.na(CASRN)
                                            ][, DTXSID.mode := Mode(DTXSID),
                                            by = CASRN
                                            ][, NAME.mode := Mode(NAME),
                                              by = CASRN
                                              ][, IUPAC_NAME.mode := Mode(IUPAC_NAME),
                                                by = CASRN]

ipr.wide <- data.table::dcast(ipr.matched[, .(NAME.mode, IUPAC_NAME.mode,
                                              CASRN, DTXSID.mode, LIST)],
                              NAME.mode + IUPAC_NAME.mode + CASRN + DTXSID.mode ~
                                LIST, length)

# Take the file over to excel for some manual edits.
openxlsx::write.xlsx(ipr.wide, './raw_data/intermediate_data/ipr_wide.xlsx')


# Match on non CASRN chems
