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
library(httr)
library(jsonlite)
#library(PubChemR)

setwd('C:/Users/atheisin/OneDrive - Environmental Protection Agency (EPA)/Documents/GitHub/ubiquitous_carcinogens/')

# Download comptox ID database if necessary (too large to include on GH)
if('DSSTox_Identifiers_and_CASRN_2021r1.csv' %in% list.files('./raw_data/')){
  comptox.dta <- data.table::fread('./raw_data/DSSTox_Identifiers_and_CASRN_2021r1.csv')
} else {
  # Data download from cluster below and save locally. 
  comptox.dta <- data.table::fread('https://clowder.edap-cluster.com/files/616dd943e4b0a5ca8aeea69d/blob')
  data.table::fwrite(comptox.dta, './raw_data/DSSTox_Identifiers_and_CASRN_2021r1.csv')
}

# Need a public key to access comptox API
apikey <- data.table::fread('comptox_apikey.csv')$comptox_api

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
                                   )#[is.na(bio_flag) # Drop viruses, etc.
                                     #][is.na(radiation_flag) #
                                      # ]#[is.na(alc_tob_flag)] # Drop alcohol & tobacco consumption
                                       
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
)#[is.na(bio_flag) # Drop viruses, etc.
#][is.na(radiation_flag) #
#]#[is.na(alc_tob_flag)] # Drop alcohol & tobacco consumption

p65.nomatch <- merge(p65.nomatch, p65.crosswalk,
                     by.x = c('Chemical'),
                     by.y = c('Chemical')
)[, .(NAME = Chemical,
      CASRN = comptox_CAS,
      LIST ='Prop65',
      DTXSID)]

# Bind cleaned P65 lists together
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

#######################################
# Pull in EU + UK SVHC lists
eu.svhc <- data.table::as.data.table(
  openxlsx::read.xlsx('./raw_data/EU_SVHC_candidatelist.xlsx',
                      startRow = 4)
)[`Reason.for.inclusion` %like% 'Carcinogenic'
][, .(name = Substance.name,
      cas = `CAS.No.`)
][cas == '10588-01-9, 7789-12-0', cas := '10588-01-9'
][cas == '10124-36-4, 31119-53-6', cas := '10124-36-4'
][cas == '302-01-2, 7803-57-8', cas := '302-01-2'
][cas == '12139-21-8', cas := '1306-19-0']

eu.clean <- rbind(merge(comptox.dta[, .(dtxsid, casrn, preferredName)],
                        eu.svhc[cas != '-'],
                        by.x = 'casrn',
                        by.y = 'cas')[,. (NAME = preferredName,
                                          CASRN = casrn,
                                          LIST ='SVHC_EU',
                                          DTXSID = dtxsid)],
                  eu.svhc[cas == '-'
                  ][, .(NAME = name,
                        CASRN = '',
                        LIST = 'SVHC_EU',
                        DTXSID = '')])

uk.svhc <- data.table::as.data.table(
  openxlsx::read.xlsx('./raw_data/UK_SVHC_candidatelist.xlsx',
                      startRow = 3,
                      sheet = 'Candidate List')
)[`Reason.for.inclusion` %like% 'Carcinogenic'
][, .(name = `Substance.name.Note:.Group.entries.are.split.in.different.rows`,
      cas = `CAS.No.`)
][cas == '10588-01-9, 7789-12-0', cas := '10588-01-9'
][cas == '10124-36-4, 31119-53-6', cas := '10124-36-4'
][cas == '302-01-2, 7803-57-8', cas := '302-01-2'
][cas == '12139-21-8', cas := '1306-19-0']

uk.clean <- rbind(merge(comptox.dta[, .(dtxsid, casrn, preferredName)],
                        uk.svhc[cas != '-'],
                        by.x = 'casrn',
                        by.y = 'cas')[,. (NAME = preferredName,
                                          CASRN = casrn,
                                          LIST ='SVHC_UK',
                                          DTXSID = dtxsid)],
                  uk.svhc[cas == '-'
                  ][, .(NAME = name,
                        CASRN = '',
                        LIST = 'SVHC_UK',
                        DTXSID = '')])

rm(eu.svhc, uk.svhc)

# merge and save to csv
ipr <- data.table::dcast(rbind(unique(iarc.clean), unique(roc.clean), unique(p65.clean),
                         unique(eu.clean), unique(uk.clean)),
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
                  'Alcoholic beverages',
                  'Estrogens, steroidal'
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
                    ][NAME == 'Tobacco, oral use of smokeless products', ROC_Known := 1L
            ][NAME == 'Estrogens, Steroidal', Prop65 := 1L]

rm(dupes, all.dupes)

ipr.noncas <- ipr[CASRN == '']

fwrite(ipr.noncas, './output/noncas_carcinogen_list.csv')


# Fix some encoding issues.
ipr <- ipr[CASRN == '120-80-9', NAME := '1,2-Benzenediol'
    ][CASRN == '13327-32-7', NAME := 'Beryllium hydroxide'
      ][CASRN == '25962-77-0', NAME := 'N,N-Dimethyl-N′-[5-[2-(5-nitro-2-furanyl)ethenyl]-1,3,4-oxadiazol-2-yl]methanimidamide'
        ][CASRN == '1204332-00-2', NAME := 'Trim VX'
          ][NAME == 'Diesel Exhaust Particulates', CASRN := 'DIESEL'
          ][CASRN != ''
          ][!(CASRN %like% 'NOCAS')
          ][!(CASRN == 'DIESEL')]



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


# Merge in supplemental information about ubiquity.

# Available chemical lists on comptox
#all.inv.lists <- httr::GET('https://api-ccte.epa.gov/chemical/list/',
#                           httr::add_headers(`accept` = 'application/json'), 
#                           httr::content_type('application/json'),
#                           httr::add_headers(`x-api-key` = apikey)) %>%
#  httr::content(.) %>%
#  data.table::rbindlist()

inventory.list <- vector(mode = 'list', length = dim(ipr)[1])
pubchemcid.list <- vector(mode = 'list', length = dim(ipr)[1])
for (i in 184:dim(ipr)[1]){
  # API info about list membership
  api.call <- httr::GET(paste0('https://api-ccte.epa.gov/chemical/list/search/by-dtxsid/', ipr$DTXSID[i]),
                        httr::add_headers(`accept` = 'application/json'), 
                        httr::add_headers(`x-api-key` = apikey)) %>%
    httr::content(.) %>%
    data.table::rbindlist(.)
  
  if(dim(api.call)[1] > 0){
    inventory.list[[i]] <- data.table::transpose(api.call[, .(V1, val = 1L)], 
                                                 make.names = 'V1')
  }
  
  #API info about pubchem ID
  api.call2 <- httr::GET(paste0('https://api-ccte.epa.gov/chemical/detail/search/by-dtxsid/', ipr$DTXSID[i]),
                         httr::add_headers(`accept` = 'application/json'), 
                         httr::add_headers(`x-api-key` = apikey)) %>%
    httr::content(.)
  
  if(!(is.null(api.call2$pubchemCid))){
    pubchemcid.list[[i]] <- api.call2$pubchemCid
  }

  if (i %% 10 == 0){print(paste0('Iteration ',i,' complete.'))}
}
rm(api.call, api.call2, i)

pubchemcid.together <- data.table::data.table(PUBCHEMCID = as.integer(purrr::map(pubchemcid.list, 
                                                                      ~ifelse(is.null(.x),NA_integer_,.x))),
                                              .id = 1:dim(ipr)[1]
                                              )[!is.na(PUBCHEMCID)
                                                ][, DTXSID := ipr$DTXSID[.id]
                                                ][, .id := NULL]

inventory.together <- data.table::rbindlist(inventory.list, fill = T, idcol = T
)[, DTXSID := ipr[.id]$DTXSID
][, .id := NULL
][, .(DTXSID, CDR2016, CDR2020,
      TSCA_ACTIVE_NCTI_0823,
      MMDBV1)]
inventory.together[is.na(inventory.together)] <- 0L
# CDR2016,2020 are flags for whether chem was produced/imported into US during those CDR rounds
# TSCA_Active are flags for whether chem is active on the non-confidential TSCA registry
# MMDBV1 are flags for whether chem was detected in the Multimedia Monitoring Database...
#... it is a very rough proxy for ubiquity: https://www.nature.com/articles/s41597-022-01365-8

inventory.together <- merge(inventory.together,
                            pubchemcid.together,
                            all.x = T,
                            by = 'DTXSID')

data.table::setcolorder(inventory.together, 'DTXSID')

rm(inventory.list, pubchemcid.list, pubchemcid.together)

# Merge together
carc.list <- merge(ipr,
                   inventory.together,
                   by = 'DTXSID',
                   all.x = T)

# Pubchem patent/literature info
pubchemids <- unique(carc.list[!is.na(PUBCHEMCID)]$PUBCHEMCID)
pchem.patent.lit <- vector(mode = 'list', length = length(pubchemids))
for (i in 1:length(pubchemids)){
  pchem.patent.lit[[i]] <- data.table::fread(
    paste0('https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/',
    pubchemids[i],
    '/property/PatentCount,LiteratureCount/CSV')
    )
}
patent.lit <- data.table::rbindlist(pchem.patent.lit)

carc.list <- merge(carc.list,
                   patent.lit,
                   by.x = 'PUBCHEMCID',
                   by.y = 'CID',
                   all.x = T)

rm(pubchemids, pchem.patent.lit, patent.lit, i)

# Can pull info for uses in bulk as such, but it's a nested mess.
#temp <- httr::GET(url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/7503/JSON?heading=Uses') %>%
#  httr::content(., as = 'text') %>%
#  jsonlite::fromJSON()

data.table::setcolorder(carc.list,
                        neworder = c("DTXSID","NAME","CASRN","PUBCHEMCID","IARC1",               
                                     "IARC2A","IARC2B","Prop65","ROC_Known",            
                                     "ROC_RAHC","SVHC_EU","SVHC_UK",
                                     "CDR2016","CDR2020","TSCA_ACTIVE_NCTI_0823","MMDBV1",
                                     "Function_count","Functions","PUC_count","PUCs",
                                     "PatentCount","LiteratureCount"
                                      ))

data.table::setorder(carc.list, CASRN)

data.table::fwrite(carc.list, './output/merged_carcinogen_list.csv')
