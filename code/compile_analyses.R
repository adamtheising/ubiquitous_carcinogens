## Synthesize multi-dimension disparate exposure analysis.

### Preliminaries
rm(list=ls(all=TRUE))
options(scipen=999)

library(data.table)
library(tidyverse)
library(readxl)
library(openxlsx)

setwd('C:/Users/atheisin/OneDrive - Environmental Protection Agency (EPA)/Documents/GitHub/ubiquitous_carcinogens')

####################
## Carcinogen list
carc.list <- data.table::fread('./output/merged_carcinogen_list.csv')
setorder(carc.list, NAME)

###########################
## NHANES summary tables
nhanes <- T
if (nhanes == T){
  nhanes.analysis <- data.table::fread('./output/nhanes_summary.csv')
  
  # 95th percentile summary table
  nhanes.wide.95ratios <- data.table::dcast(unique(nhanes.analysis[, .(chem, cas, cat, perc_95)]),
                                   chem + cas ~ cat, value.var = 'perc_95'
                                   )[, 4:29 := lapply(.SD, function(x) x/All),
                                     .SDcols = 4:29
                                     ][, .(chem, cas,
                                           `Inc/PL: < 1`,`Inc/PL: between 1 and 2`,
                                           `Race/eth: Hispanic`,`Race/eth: Non-hispanic black`,`Race/eth: Other`,
                                           `Race/gender: Black, female`,`Race/gender: Black, male`,
                                           `Race/gender: Hispanic, female`,`Race/gender: Hispanic, male`,
                                           `Race/gender: Other, female`,`Race/gender: Other, male`,
                                           `Race/income: Hispanic, Inc/PL < 2`,`Race/income: Hispanic, Inc/PL >= 2`,
                                           `Race/income: Black, Inc/PL < 2`, `Race/income: Black, Inc/PL >= 2`,
                                           `Race/income: Other, Inc/PL < 2`, `Race/income: Other, Inc/PL >= 2`
                                           )]
  
  nhanes.wide.95ratios <- nhanes.wide.95ratios[rowSums(is.na(nhanes.wide.95ratios)) <= 16] %>%
    dplyr::select(-chem) %>%
    dplyr::group_by(cas) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), max)) %>% # taking the max of multiple measures
    dplyr::left_join(carc.list[, .(NAME, CASRN)],
                     by = dplyr::join_by(cas == CASRN)) %>%
    data.table::as.data.table()
  
  data.table::setcolorder(nhanes.wide.95ratios,
                          c('NAME', 'cas'))
  # Flag chemicals with 0/1
  nh95.discrete <- data.table::as.data.table(nhanes.wide.95ratios)[, NAME := NULL]
  names(nh95.discrete)[2:18] <- paste0('NHANES (P95) ',names(nh95.discrete)[2:18])
  
  for(j in names(nh95.discrete)[2:18]){
    data.table::set(nh95.discrete,
                    i = which(nh95.discrete[[j]] < 1.5),
                    j = j,
                    value = 0)
    data.table::set(nh95.discrete,
                    i = which(nh95.discrete[[j]] >= 1.5),
                    j = j,
                    value = 1)
  }
  
  
  # Geometric mean summary table
  nhanes.wide.gmratios <- data.table::dcast(unique(nhanes.analysis[, .(chem, cas, cat, mean)]),
                                            chem + cas ~ cat, value.var = 'mean'
  )[, 4:29 := lapply(.SD, function(x) x/All),
    .SDcols = 4:29
  ][, .(chem, cas,
        `Inc/PL: < 1`,`Inc/PL: between 1 and 2`,
        `Race/eth: Hispanic`,`Race/eth: Non-hispanic black`,`Race/eth: Other`,
        `Race/gender: Black, female`,`Race/gender: Black, male`,
        `Race/gender: Hispanic, female`,`Race/gender: Hispanic, male`,
        `Race/gender: Other, female`,`Race/gender: Other, male`,
        `Race/income: Hispanic, Inc/PL < 2`,`Race/income: Hispanic, Inc/PL >= 2`,
        `Race/income: Black, Inc/PL < 2`, `Race/income: Black, Inc/PL >= 2`,
        `Race/income: Other, Inc/PL < 2`, `Race/income: Other, Inc/PL >= 2`
  )]
  
  nhanes.wide.gmratios <- nhanes.wide.gmratios[rowSums(is.na(nhanes.wide.gmratios)) <= 16] %>%
    dplyr::select(-chem) %>%
    dplyr::group_by(cas) %>%
    dplyr::summarise(dplyr::across(dplyr::everything(), max)) %>% # taking the max of multiple measures
    dplyr::left_join(carc.list[, .(NAME, CASRN)],
                     by = dplyr::join_by(cas == CASRN)) %>%
    data.table::as.data.table()
  
  data.table::setcolorder(nhanes.wide.gmratios,
                          c('NAME', 'cas'))
  
  nhgm.discrete <- data.table::as.data.table(nhanes.wide.gmratios)[, NAME := NULL]
  names(nhgm.discrete)[2:18] <- paste0('NHANES (GM) ',names(nhgm.discrete)[2:18])
  
  
  for(j in names(nhgm.discrete)[2:18]){
    data.table::set(nhgm.discrete,
                    i = which(nhgm.discrete[[j]] < 1.5),
                    j = j,
                    value = 0)
    data.table::set(nhgm.discrete,
                    i = which(nhgm.discrete[[j]] >= 1.5),
                    j = j,
                    value = 1)
  }
  
  nhanes.vintage <- unique(nhanes.analysis[, .(cas, year)]
                           )[, start_year := as.numeric(substr(year,1,4))
                           ][, max_year := max(start_year),
                             by = cas
                           ][max_year == start_year
                           ][, .(cas, year)]
}

## HH products summary table
hh_prod <- T
if (hh_prod == T){
  hhprod.sum <- data.table::fread('./output/hh_prod_summary.csv')
  
  hhprod.discrete <- data.table::as.data.table(hhprod.sum)[, Chemical := NULL]
  
  for(j in names(hhprod.discrete)[1:6]){
    data.table::set(hhprod.discrete,
                    i = which(hhprod.discrete[[j]] < 1.5),
                    j = j,
                    value = 0)
    data.table::set(hhprod.discrete,
                    i = which(hhprod.discrete[[j]] >= 1.5),
                    j = j,
                    value = 1)
  }
  
  names(hhprod.discrete)[1:6] <- paste0('HH Prod ', names(hhprod.discrete)[1:6])
}

## RSEI summary table
rsei <- T
if (rsei == T){
  rsei.sum <- data.table::fread('./output/rsei_summary.csv'
                                )[, .(`Race/Eth: Non-hispanic black` = mean_conc_black/mean_conc_nhw,
                                      `Race/Eth: American Indian` = mean_conc_amerind/mean_conc_nhw,
                                      `Race/Eth: Hispanic` = mean_conc_hispanic/mean_conc_nhw,
                                      `Race/Eth: Asian` = mean_conc_asian/mean_conc_nhw,
                                      `Race/Eth: Pacific Islander` = mean_conc_pacisl/mean_conc_nhw,
                                      `Income: Below PL` = mean_conc_lowinc/mean_conc_medhighinc,
                                      `Income: Below 2x PL` = mean_conc_verylowinc/mean_conc_medhighinc),
                                  by = .(NAME, CASRN)]
  
  rsei.discrete <- as.data.table(rsei.sum)[, NAME := NULL]
  
  for(j in names(rsei.discrete)[2:8]){
    data.table::set(rsei.discrete,
                    i = which(rsei.discrete[[j]] < 1.5),
                    j = j,
                    value = 0)
    data.table::set(rsei.discrete,
                    i = which(rsei.discrete[[j]] >= 1.5),
                    j = j,
                    value = 1)
  }
  
  names(rsei.discrete)[2:8] <- paste0('RSEI ', names(rsei.discrete)[2:8])
}

## DW summary table
dw <- T
if (dw == T){
  dw.sum <- data.table::fread('./output/dw_demo_summary.csv'
                                   )[CAS_REGISTRY_NUM != ''
                                     ][, .(`Race/Eth: Non-hispanic black` = mean_concentration_black/mean_concentration_nhw,
                                           `Race/Eth: American Indian` = mean_concentration_amerind/mean_concentration_nhw,
                                           `Race/Eth: Hispanic` = mean_concentration_hispanic/mean_concentration_nhw,
                                           `Race/Eth: Asian` = mean_concentration_asian/mean_concentration_nhw,
                                           `Race/Eth: Pacific Islander` = mean_concentration_pacisl/mean_concentration_nhw,
                                           `Income: Below PL` = mean_concentration_lowinc/mean_concentration_medhighinc,
                                           `Income: Below 2x PL` = mean_concentration_verylowinc/mean_concentration_medhighinc),
                                       by = .(FULL_NAME, CAS_REGISTRY_NUM)]
  
  dw.discrete <- data.table::as.data.table(dw.sum)[, FULL_NAME := NULL]
  
  for(j in names(dw.discrete)[2:8]){
    data.table::set(dw.discrete,
                    i = which(dw.discrete[[j]] < 1.5),
                    j = j,
                    value = 0)
    data.table::set(dw.discrete,
                    i = which(dw.discrete[[j]] >= 1.5),
                    j = j,
                    value = 1)
  }
  
  names(dw.discrete)[2:8] <- paste0('DW ', names(dw.discrete)[2:8])
}

######################################
## Merge summary tables into composite
together <- merge(nhgm.discrete,
                  nh95.discrete,
                  by = 'cas',
                  all = T)

together1 <- merge(together,
                   hhprod.discrete,
                   by.x = 'cas',
                   by.y = 'CASRN',
                   all = T)

together2 <- merge(together1,
                  dw.discrete,
                  by.x = 'cas',
                  by.y = 'CAS_REGISTRY_NUM',
                  all = T)

together3 <- merge(together2,
                   rsei.discrete,
                   by.x = 'cas',
                   by.y = 'CASRN',
                   all = T)

together.all <- merge(carc.list[, .(CASRN, NAME, Functions, PUCs)],
                      together3,
                      by.x = 'CASRN',
                      by.y = 'cas')

together.all <- merge(together.all,
                      nhanes.vintage,
                      by.x = 'CASRN',
                      by.y = 'cas',
                      all.x = T)

rm(together, together1, together2, together3)

together.summary <- 
  together.all[, `NHANES GM Score`:= ifelse(as.numeric(rowSums(!is.na(.SD)) == 0) == 0,
                        rowSums(.SD, na.rm = T),
                        rowSums(.SD)),
               .SDcols = 5:21
               ][, `NHANES P95 Score`:= ifelse(as.numeric(rowSums(!is.na(.SD)) == 0) == 0,
                                              rowSums(.SD, na.rm = T),
                                              rowSums(.SD)),
                 .SDcols = 22:38
               ][, `HH Product Score`:= ifelse(as.numeric(rowSums(!is.na(.SD)) == 0) == 0,
                                                   rowSums(.SD, na.rm = T),
                                                   rowSums(.SD)),
                 .SDcols = 39:44
               ][, `Drinking Water Score`:= ifelse(as.numeric(rowSums(!is.na(.SD)) == 0) == 0,
                                               rowSums(.SD, na.rm = T),
                                               rowSums(.SD)),
                 .SDcols = 45:51
               ][, `RSEI Air Score`:= ifelse(as.numeric(rowSums(!is.na(.SD)) == 0) == 0,
                                               rowSums(.SD, na.rm = T),
                                               rowSums(.SD)),
                 .SDcols = 52:58
               ][,.(CASRN, NAME,
                    `NHANES GM Score`,
                    `NHANES P95 Score`,
                    `HH Product Score`,
                    `Drinking Water Score`,
                    `RSEI Air Score`,
                    Functions,
                    PUCs,
                    `NHANES vintage` = year)
               ][, `Likely disparate flag` := F
               ][`NHANES GM Score` >= 1, `Likely disparate flag` := T
               ][`NHANES P95 Score` >= 1, `Likely disparate flag` := T
               ][(`Drinking Water Score` >= 1) & (`RSEI Air Score` >= 1),
                 `Likely disparate flag` := T
               ][(`Drinking Water Score` >= 1) & (`HH Product Score` >= 1),
                 `Likely disparate flag` := T
               ][(`RSEI Air Score` >= 1) & (`HH Product Score` >= 1),
                 `Likely disparate flag` := T
               ][, `Sum of Scores` := rowSums(.SD, na.rm = T), .SDcols = 3:7
               ][, Functions := gsub(' (EPA)','',Functions, fixed = T)]

data.table::setorder(together.summary, -`Likely disparate flag`, 
                     -`Sum of Scores`, `CASRN`)
data.table::setcolorder(together.summary, 
                        c('CASRN','NAME',
                          'Likely disparate flag',
                          'Sum of Scores'))

###################################
## Write summary tables to xlsx
list_of_datasets <- list('Aggregated Scores' = together.summary,
                         'Combined Discrete Flags' = together.all[, Functions := NULL
                                                                 ][, PUCs := NULL
                                                                 ][, year := NULL
                                                                 ][,1:50],
                         'Carcinogen List' = carc.list,
                         'NHANES GM summary' = nhanes.wide.gmratios,
                         'NHANES P95 summary' = nhanes.wide.95ratios,
                         'HH Products summary' = hhprod.sum,
                         'DW summary' = dw.sum, 
                         'RSEI summary' = rsei.sum)
openxlsx::write.xlsx(list_of_datasets, './output/summary_6_11_2024.xlsx')
