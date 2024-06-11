# RSEI

### Preliminaries
rm(list=ls(all=TRUE))
options(scipen=999)

library(tidyverse)
library(data.table)
library(readxl)
library(tidycensus)
library(future)
library(furrr)

##############################
### Pull in RSEI microdata ###
setwd('C:/Users/atheisin/OneDrive - Environmental Protection Agency (EPA)/Documents/research/TRI_data')
dta.1 <- data.table::fread('censusmicroblockgroup2022_2022.csv',
               select = c('GeoID','ReleaseNumber','ChemicalNumber',
                          'FacilityNumber','Media','Conc'),
               colClasses = c(GeoID = 'character')
               )[Media < 3] #only air releases

chem.info <- data.table::fread('chemical_data_rsei_v2312.csv'
)[, .(ChemicalNumber, Chemical,
      CASNumber, CASStandard)]

# Merge the chemical info with the concentration info
dta.together <- merge(dta.1,
                      chem.info,
                      by = 'ChemicalNumber'
)[, ChemicalNumber := NULL]

rm(dta.1, chem.info); gc()
               
# Collapse chemical concentration to cbg level
# Implicitly assuming that addition of concentrations from multiple sources
# sums up to a total concentration. This seems a very reasonable assumption,
# note not saying anything explicit about how multiple sources of exposure 
# impact cancer/other risk
dta.together <- dta.together[CASStandard == 'N982', CASStandard := '7440-66-6'
][CASStandard == 'N874', CASStandard := '81-81-2'
][CASStandard == 'N770', CASStandard := '7440-62-2'
][CASStandard == 'N760', CASStandard := '7440-28-0'
][CASStandard == 'N746', CASStandard := '57-24-9'
][CASStandard == 'N740', CASStandard := '7440-22-4'
][CASStandard == 'N725', CASStandard := '7782-49-2'
][CASStandard == 'N590', CASStandard := '130498-29-2'
][CASStandard == 'N583', CASStandard := '108171-26-2'
][CASStandard == 'N575', CASStandard := '59536-65-1'
][CASStandard == 'N535', CASStandard := 'N535'
][CASStandard == 'N530', CASStandard := '25154-52-3'
][CASStandard == 'N511', CASStandard := 'N511'
][CASStandard == 'N503', CASStandard := '54-11-5'
][CASStandard == 'N495', CASStandard := '7440-02-0'
][CASStandard == 'N458', CASStandard := '7439-97-6'
][CASStandard == 'N450', CASStandard := '7439-96-5'
][CASStandard == 'N420', CASStandard := '7439-92-1'
][CASStandard == 'N270', CASStandard := '3194-55-6'
][CASStandard == 'N230', CASStandard := 'N230'
][CASStandard == 'N171', CASStandard := '111-54-6'
][CASStandard == 'N150', CASStandard := '1746-01-6'
][CASStandard == 'N120', CASStandard := '26471-62-5'
][CASStandard == 'N106', CASStandard := '57-12-5'
][CASStandard == 'N100', CASStandard := '7440-50-8'
][CASStandard == 'N096', CASStandard := '7440-48-4'
][CASStandard == 'N090', CASStandard := '7440-47-3'
][CASStandard == 'N084', CASStandard := '88-06-2'
][CASStandard == 'N078', CASStandard := '7440-43-9'
][CASStandard == 'N050', CASStandard := '7440-41-7'
][CASStandard == 'N040', CASStandard := '7440-39-3'
][CASStandard == 'N020', CASStandard := '7440-38-2'
][CASStandard == 'N010', CASStandard := '7440-36-0'
][, .(conc_total = sum(Conc)),
  by = .(GeoID, CASStandard)]

##############################################################
### Get the carcinogen list for trimming down set of chems ###
setwd('C:/Users/atheisin/OneDrive - Environmental Protection Agency (EPA)/Documents/GitHub/ubiquitous_carcinogens/')
carc.list <- data.table::fread('output/merged_carcinogen_list.csv')

dta.trimmed <- dta.together[CASStandard %in% carc.list$CASRN]

rm(dta.together)

##############################################
### Get 2022 ACS data for all block groups ###

parallel.api <- function(st){
  tidycensus::get_acs(
    geography = 'block group',
    variables = c(pop_total = 'B03002_001',
                  pop_white = 'B03002_003',
                  pop_black = 'B03002_004',
                  pop_amerind = 'B03002_005',
                  pop_asian = 'B03002_006',
                  pop_pacisl = 'B03002_007',
                  pop_hisp = 'B03002_012',
                  pov_total = 'C17002_001',
                  pov50 = 'C17002_002',
                  pov99 = 'C17002_003',
                  pov124 = 'C17002_004',
                  pov149 = 'C17002_005',
                  pov184 = 'C17002_006',
                  pov199 = 'C17002_007',
                  med_inc = 'B19013_001'),
    state = st,
    year = 2022,
    geometry = F
  ) %>%
    dplyr::select(GEOID, variable, estimate) %>%
    tidyr::pivot_wider(id_cols = c(GEOID), names_from = variable, values_from = estimate)
}

## State lists
state.list <- c(state.abb, 'DC')

# Loop (in parallel) thru each county/state pair, calling census API
future::plan("multisession", workers = (parallel::detectCores()-3))
cbg.list <- furrr::future_pmap(list(state.list), parallel.api,
                               .options = furrr::furrr_options(seed = NULL))
future::plan('sequential')

# Bind together and clean into variables of interest
cbg.together <- data.table::rbindlist(cbg.list, use.names = T, idcol = 'statelistid')
convert.cols <- names(cbg.together)[-2]
acs.cbg.data <- cbg.together[, (convert.cols) := lapply(.SD, as.integer),
                             .SDcols = convert.cols
][, state := state.list[statelistid]
][, .(GEOID, state, med_inc,
      pop_white,
      pop_black,
      pop_amerind,
      pop_asian,
      pop_pacisl,
      pop_hisp,
      pov100 = pov50+pov99,
      pov200 = pov50+pov99+pov124+pov149+pov184+pov199,
      pop_total)]
rm(convert.cols, cbg.list, cbg.together, parallel.api, state.list)

# Create a vector of pop sums so we can use during weighted mean calculations
convert.cols <- names(acs.cbg.data)[4:12]
demographic.sums <- acs.cbg.data[, lapply(.SD, sum, na.rm = T),
                                 .SDcols = convert.cols]
rm(convert.cols)

# As of 4/29/24, CT is weird with their new CBG codes
ct.weird <- data.table::fread('./raw_data/crosswalks/CT_2022blockcrosswalk.csv',
                               colClasses = 'character'
                               )[, bg2020 := paste0('0', substr(block_fips_2020,1,11))
                                 ][, bg2022 := paste0('0', substr(block_fips_2022,1,11))
                                 ][, .(bg2020, bg2022)
                                   ][, unique(.SD)]

dta.trimmed <- merge(dta.trimmed,
                     ct.weird,
                     by.x = 'GeoID',
                     by.y = 'bg2020',
                     all.x = T)

dta.trimmed[, GeoID := ifelse(is.na(bg2022),
                              GeoID,
                              bg2022)
            ][, bg2022 := NULL]

rm(ct.weird)

# Join demographics with RSEI concentrations:
rsei.dems <- merge(dta.trimmed,
                   acs.cbg.data,
                   by.x = 'GeoID',
                   by.y = 'GEOID')[, pov_above200 := pop_total - pov200]

# Calculate some summary stats
# Calculate average concentration stats by subpopulation. 
rsei.natl <- rsei.dems[, .(mean_concentration_black = weighted.mean(conc_total, pop_black, na.rm = T),
                              mean_concentration_amerind = weighted.mean(conc_total, pop_amerind, na.rm = T),
                              mean_concentration_hispanic = weighted.mean(conc_total, pop_hisp, na.rm = T),
                              mean_concentration_asian = weighted.mean(conc_total, pop_asian, na.rm = T),
                              mean_concentration_pacisl = weighted.mean(conc_total, pop_pacisl, na.rm = T),
                              mean_concentration_nhw = weighted.mean(conc_total, pop_white, na.rm = T),
                              mean_concentration_lowinc = weighted.mean(conc_total, pov200, na.rm = T),
                              mean_concentration_verylowinc = weighted.mean(conc_total, pov100, na.rm = T),
                              mean_concentration_medhighinc = weighted.mean(conc_total, pov_above200, na.rm = T),
                              total_pop_black = sum(pop_black, na.rm = T),
                              total_pop_amerind = sum(pop_amerind, na.rm = T),
                              total_pop_hispanic = sum(pop_hisp, na.rm = T),
                              total_pop_asian = sum(pop_asian, na.rm = T),
                              total_pop_pacisl = sum(pop_pacisl, na.rm = T),
                              total_pop_nhw = sum(pop_white, na.rm = T),
                              total_pop_lowinc = sum(pov200, na.rm = T),
                              total_pop_verylowinc = sum(pov100, na.rm = T),
                              total_pop_medhighinc = sum(pov_above200, na.rm = T),
                              contributing_cbgs = .N),
                          by = .(CASStandard)]

# Next, account for populations living in grid cells with no exposure from TRI facils
rsei.natl <- rsei.natl[
  , .(CASStandard,
      mean_conc_black = mean_concentration_black * (total_pop_black / demographic.sums$pop_black),
      mean_conc_amerind = mean_concentration_amerind * (total_pop_amerind / demographic.sums$pop_amerind),
      mean_conc_hispanic = mean_concentration_hispanic * (total_pop_hispanic / demographic.sums$pop_hisp),
      mean_conc_asian = mean_concentration_asian * (total_pop_asian / demographic.sums$pop_asian),
      mean_conc_pacisl = mean_concentration_pacisl * (total_pop_pacisl / demographic.sums$pop_pacisl),
      mean_conc_nhw = mean_concentration_nhw * (total_pop_nhw / demographic.sums$pop_white),
      mean_conc_lowinc = mean_concentration_lowinc * (total_pop_lowinc / demographic.sums$pov200),
      mean_conc_verylowinc = mean_concentration_verylowinc * (total_pop_verylowinc / demographic.sums$pov100),
      mean_conc_medhighinc = mean_concentration_medhighinc * (total_pop_medhighinc / (demographic.sums$pop_total - demographic.sums$pov200)),
      contributing_cbgs
      )]

rsei.summary <- merge(carc.list[, .(NAME, CASRN, DTXSID)],
                      rsei.natl,
                      by.x = 'CASRN',
                      by.y = 'CASStandard')

data.table::fwrite(rsei.summary, './output/rsei_summary.csv')
