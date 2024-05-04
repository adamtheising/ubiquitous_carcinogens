### Analyze Drinking Water data as exposure proxy

### Preliminaries
rm(list=ls(all=TRUE))
options(scipen=999)

library(tidyverse)
library(data.table)
library(readxl)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

setwd('C:/Users/atheisin/OneDrive - Environmental Protection Agency (EPA)/Documents/GitHub/ubiquitous_carcinogens/')

## Carcinogen list
carc.list <- data.table::fread('output/merged_carcinogen_list.csv')

## Crosswalks
syr.crosswalk <- data.table::as.data.table(
  readxl::read_excel('raw_data/crosswalks/IN_TSAANLYT_TABLE.xlsx')
)

ucmr.crosswalk <- data.table::as.data.table(
  readxl::read_excel('raw_data/crosswalks/UMCR_cas_crosswalk.xlsx',
                     sheet = 'ucmr_chems')
)

haas.crosswalk <- data.table::as.data.table(
  readxl::read_excel('raw_data/crosswalks/UMCR_cas_crosswalk.xlsx',
                     sheet = 'HAAs')
)

## Drinking water data directory: access via W. Austin.
setwd('C:/Users/atheisin/OneDrive - Environmental Protection Agency (EPA)/Data/mdi')

## UCMR summaries
## Code in IF block below creates summary statistics for chem x PWS using UCMR data
clean_summary_ucmr <- F
if (clean_summary_ucmr == T) {
  ucmr.tog <- rbind(
    data.table::fread('ucmr/ucmr_2_main.csv')[, vintage := 'ucmr2'],
    data.table::fread('ucmr/ucmr_3_main.csv')[, vintage := 'ucmr3'],
    data.table::fread('ucmr/ucmr_4_main.csv')[, vintage := 'ucmr4'
          ][, analytical_results_value := analytical_results_value__g_l_
            ][, analytical_results_value__g_l_ := NULL
              ][, ucmr1_sample_type := NULL],
    fill = T
  )
  
  ucmr.tog2 <- merge(ucmr.tog[, contaminant := toupper(contaminant)],
                     ucmr.crosswalk,
                     by.x = 'contaminant',
                     by.y = 'CHEMICAL',
                     all.x = T)
  
  ucmr.carc <- rbind(ucmr.tog2[CAS_NUM %in% carc.list$CASRN],
                     ucmr.tog2[contaminant %in% c('HAA5','HAA6BR','HAA9')]
                     )[!(is.na(mrl) & is.na(analytical_results_value))]
  
  ucmr.dl <- ucmr.carc[, result_dl := ifelse(analytical_results_sign == '=',
                                             analytical_results_value,
                                             mrl/sqrt(2) #if non-detect, set = new DL/sqrt(2)
  )]
  
  # Collapse to chemical x PWS statistics
  ucmr.summary <- ucmr.dl[, .(num_detect = sum(analytical_results_sign == '='),
                              num_total = .N,
                              mean_result = gm_mean(result_dl, na.rm = T),
                              median_result = quantile(result_dl, 0.5, na.rm = T),
                              p95_result = quantile(result_dl, 0.95, na.rm = T)),
                          by = .(pws_id, CAS_REGISTRY_NUM = CAS_NUM, 
                                 FULL_NAME = contaminant, vintage)
  ]
  
  fwrite(ucmr.summary, 'C:/Users/atheisin/OneDrive - Environmental Protection Agency (EPA)/Documents/GitHub/ubiquitous_carcinogens/output/drinkingwater/ucmr_chembypws_summary.csv')
  rm(ucmr.tog, ucmr.tog2, ucmr.carc, ucmr.dl, ucmr.summary)
}

## Code in IF block below trims down six year review data to a managable size
## and saves locally. These files were cleaned and synchronized by W. Austin.
read_syr_data <- F
if (read_syr_data == T) {
  ######################################################
  ## Be careful reading these in, they are >10gb memory
  ## Crashed my shit epa laptop
  
  ## Six Year Review 4: part b of data
  syr4b <- data.table::as.data.table(
    readRDS('syr/syr4_pt2.rds')
    )[, .(analyte_code, pws_id, water_facility_id, sample_id, 
           detection_limit_value, detection_limit_code, detect,
           result_value, unit,
           presence_indicator_code, residual_field_free_chlorine_mg_l,
           residual_field_total_chlorine_mg_l)]
  gc()
  syr4b.carc <- merge(data.table::as.data.table(syr4b)[, analyte_code := as.character(analyte_code)],
                               syr.crosswalk[, .(CODE,CAS_REGISTRY_NUM,FULL_NAME)],
                               by.x = 'analyte_code',
                               by.y = 'CODE')[CAS_REGISTRY_NUM %in% carc.list[CASRN != ""]$CASRN
                               ]
  data.table::fwrite(syr4b.carc,'C:/Users/atheisin/OneDrive - Environmental Protection Agency (EPA)/Documents/GitHub/ubiquitous_carcinogens/output/drinkingwater/syr4b.csv')
  rm(syr4b, syr4b.carc)
  
  ## Six Year Review 4: part a of data
  syr4a <- data.table::as.data.table(
    readRDS('syr/syr4_pt1.rds')
  )[, .(analyte_code, pws_id, water_facility_id, sample_id, 
        detection_limit_value, detection_limit_code, detect,
        result_value, unit,
        presence_indicator_code, residual_field_free_chlorine_mg_l,
        residual_field_total_chlorine_mg_l)]
  gc()
  syr4a.carc <- merge(data.table::as.data.table(syr4a)[, analyte_code := as.character(analyte_code)],
                      syr.crosswalk[, .(CODE,CAS_REGISTRY_NUM,FULL_NAME)],
                      by.x = 'analyte_code',
                      by.y = 'CODE')[CAS_REGISTRY_NUM %in% carc.list[CASRN != ""]$CASRN
                      ]
  data.table::fwrite(syr4a.carc,'C:/Users/atheisin/OneDrive - Environmental Protection Agency (EPA)/Documents/GitHub/ubiquitous_carcinogens/output/drinkingwater/syr4a.csv')
  rm(syr4a, syr4a.carc)
  
  ## Six Year Review 3: part b of data
  syr3b <- data.table::as.data.table(
    readRDS('syr/syr3_pt2.rds')
  )[, .(analyte_id, pws_id, water_facility_id, sampling_point_id, 
        detection_limit_value, detection_limit_code, detect,
        result_value, unit,
        presence_indicator_code, residual_field_free_chlorine_mg_l,
        residual_field_total_chlorine_mg_l)]
  gc()
  syr3b.carc <- merge(data.table::as.data.table(syr3b)[, analyte_id := as.character(analyte_id)],
                      syr.crosswalk[, .(CODE,CAS_REGISTRY_NUM,FULL_NAME)],
                      by.x = 'analyte_id',
                      by.y = 'CODE')[CAS_REGISTRY_NUM %in% carc.list[CASRN != ""]$CASRN
                      ]
  data.table::fwrite(syr3b.carc,'C:/Users/atheisin/OneDrive - Environmental Protection Agency (EPA)/Documents/GitHub/ubiquitous_carcinogens/output/drinkingwater/syr3b.csv')
  rm(syr3b, syr3b.carc)
  
  ## Six Year Review 3: part a of data
  syr3a <- data.table::as.data.table(
    readRDS('syr/syr3_pt1.rds')
  )[, .(analyte_id, pws_id, water_facility_id, sampling_point_id, 
        detection_limit_value, detection_limit_code, detect,
        result_value, unit,
        presence_indicator_code, residual_field_free_chlorine_mg_l,
        residual_field_total_chlorine_mg_l)]
  gc()
  syr3a.carc <- merge(data.table::as.data.table(syr3a)[, analyte_id := as.character(analyte_id)],
                      syr.crosswalk[, .(CODE,CAS_REGISTRY_NUM,FULL_NAME)],
                      by.x = 'analyte_id',
                      by.y = 'CODE')[CAS_REGISTRY_NUM %in% carc.list[CASRN != ""]$CASRN
                      ]
  data.table::fwrite(syr3a.carc,'C:/Users/atheisin/OneDrive - Environmental Protection Agency (EPA)/Documents/GitHub/ubiquitous_carcinogens/output/drinkingwater/syr3a.csv')
  rm(syr3a, syr3a.carc)
  
  ## Six Year Review 2:
  syr2 <- data.table::as.data.table(
    readRDS('syr/syr2.rds')
    )[, .(pws_id, id, sampling_point_id, analyte_id,
                   detect, result_value, unit)]
  
  syr2.carc <- merge(data.table::as.data.table(syr2)[, analyte_id := as.character(analyte_id)],
                     syr.crosswalk[, .(CODE,CAS_REGISTRY_NUM,FULL_NAME)],
                     by.x = 'analyte_id',
                     by.y = 'CODE')[CAS_REGISTRY_NUM %in% carc.list[CASRN != ""]$CASRN
                     ]
  
  data.table::fwrite(syr2.carc,'C:/Users/atheisin/OneDrive - Environmental Protection Agency (EPA)/Documents/GitHub/ubiquitous_carcinogens/output/drinkingwater/syr2.csv')
  
  rm(syr2, syr2.carc)
}

# Directory where intermediate DW data summaries were saved.
setwd('C:/Users/atheisin/OneDrive - Environmental Protection Agency (EPA)/Documents/GitHub/ubiquitous_carcinogens/output/drinkingwater/')

## Code in IF block below cleans and reshapes SYR data
clean_reshape_syr <- F
if (clean_reshape_syr == T){
  #############################################
  ## SYR4 data:
  syr4 <- rbind(data.table::fread('syr4a.csv'),
                data.table::fread('syr4b.csv'))
  
  # Assign missing detection/reporting limits based on lowest nonzero value in class
  min_dl <- syr4[detection_limit_code == 'DL', detection_limit_code := 'MDL'
                 ][detection_limit_code == 'RL', detection_limit_code := 'MRL'
                 ][detection_limit_code == '', detection_limit_code := 'MDL'
                 ][detection_limit_value == 0, detection_limit_value := NA_real_ #replace "DL" of 0 with NA
                 ][, .(min_dl = min(detection_limit_value, na.rm = T)),
                 by = .(CAS_REGISTRY_NUM, detection_limit_code)]
  
  syr4.min_dl <- merge(syr4[detection_limit_code == 'DL', detection_limit_code := 'MDL'
  ][detection_limit_code == 'RL', detection_limit_code := 'MRL'
  ][detection_limit_code == '', detection_limit_code := 'MDL'
  ],
  min_dl,
  by.x = c('CAS_REGISTRY_NUM','detection_limit_code'),
  by.y = c('CAS_REGISTRY_NUM','detection_limit_code')
  )[detection_limit_value == 0, detection_limit_value := NA_real_ #replace "DL" of 0 with NA
  ][is.na(detection_limit_value), detection_limit_value := min_dl #set DL to minimum value if missing
  ][, result_dl := result_value
  ][detect == 0, result_dl := detection_limit_value/sqrt(2) #if non-detect, set = new DL/sqrt(2)
  ]
  
  # Collapse to chemical x PWS statistics
  syr4.summary <- syr4.min_dl[, .(num_detect = sum(detect == 1),
        num_total = .N,
        mean_result = gm_mean(result_dl, na.rm = T),
        median_result = quantile(result_dl, 0.5, na.rm = T),
        p95_result = quantile(result_dl, 0.95, na.rm = T)),
    by = .(pws_id, CAS_REGISTRY_NUM, FULL_NAME)
  ]
  
  data.table::fwrite(syr4.summary, 'syr4_chembypws_summary.csv')
  rm(syr4.min_dl, syr4, min_dl, syr4.summary)
  
  #############################################
  ## SYR3 data:
  syr3 <- rbind(data.table::fread('syr3a.csv'),
                data.table::fread('syr3b.csv'))
  
  # Assign missing detection/reporting limits based on lowest nonzero value in class
  min_dl <- syr3[, .(min_dl = min(detection_limit_value, na.rm = T)),
    by = .(CAS_REGISTRY_NUM, detection_limit_code)]
  
  syr3.min_dl <- merge(syr3,
  min_dl,
  by.x = c('CAS_REGISTRY_NUM','detection_limit_code'),
  by.y = c('CAS_REGISTRY_NUM','detection_limit_code')
  )[is.na(detection_limit_value), detection_limit_value := min_dl #set DL to minimum value if missing
  ][, result_dl := result_value
  ][detect == 0, result_dl := detection_limit_value/sqrt(2) #if non-detect, set = new DL/sqrt(2)
  ]
  
  # Collapse to chemical x PWS statistics
  syr3.summary <- syr3.min_dl[, .(num_detect = sum(detect == 1),
                                  num_total = .N,
                                  mean_result = gm_mean(result_dl, na.rm = T),
                                  median_result = quantile(result_dl, 0.5, na.rm = T),
                                  p95_result = quantile(result_dl, 0.95, na.rm = T)),
                              by = .(pws_id, CAS_REGISTRY_NUM, FULL_NAME)
  ]
  
  data.table::fwrite(syr3.summary, 'syr3_chembypws_summary.csv')
  rm(syr3.min_dl, syr3, min_dl, syr3.summary)
  
  #############################################
  ## SYR2 data:
  syr2 <- data.table::fread('syr2.csv')
  
  syr2.dl <- syr2[, result_dl := result_value
  ][detect == 0, result_dl := result_value/sqrt(2) #if non-detect, set = new DL/sqrt(2)
  ]
  
  # Collapse to chemical x PWS statistics
  syr2.summary <- syr2.dl[, .(num_detect = sum(detect == 1),
                                  num_total = .N,
                                  mean_result = gm_mean(result_dl, na.rm = T),
                                  median_result = quantile(result_dl, 0.5, na.rm = T),
                                  p95_result = quantile(result_dl, 0.95, na.rm = T)),
                              by = .(pws_id, CAS_REGISTRY_NUM, FULL_NAME)
  ]
  
  data.table::fwrite(syr2.summary, 'syr2_chembypws_summary.csv')
  rm(syr2.dl, syr2, syr2.summary)
}

######################################################
## Demographic analysis of SYR & UCMR carcinogen data

## Pull in summary SYR files:
syr4 <- data.table::fread('syr4_chembypws_summary.csv')[, vintage := 'syr4']
syr3 <- data.table::fread(
  'syr3_chembypws_summary.csv'
  )[!(CAS_REGISTRY_NUM %in% unique(syr4$CAS_REGISTRY_NUM))
    ][, vintage := 'syr3']
# note: syr2 is empty-- everything covered in syr3
syr2 <- data.table::fread(
  'syr2_chembypws_summary.csv'
)[!(CAS_REGISTRY_NUM %in% c(unique(syr4$CAS_REGISTRY_NUM), 
                            unique(syr3$CAS_REGISTRY_NUM)))]

## Pull in summary UCMR file:
ucmr <- data.table::fread(
  'ucmr_chembypws_summary.csv'
  )[!(CAS_REGISTRY_NUM %in% c(unique(syr4$CAS_REGISTRY_NUM), 
                              unique(syr3$CAS_REGISTRY_NUM)))]

syr.together <- rbind(syr4,syr3,ucmr)
rm(syr4,syr3,syr2,ucmr)

## Read in DW boundary demographics file.
epic.pops <- data.table::fread(
  'C:/Users/atheisin/OneDrive - Environmental Protection Agency (EPA)/demographics/epic_dems.csv'
  )[, .(pwsid, pws_name, population_served_count,
        service_connections_count, MINORPCT, LOWINCPCT,
        LINGISOPCT, frac_white, frac_black,
        frac_amerind, frac_asian, frac_pacisl,
        frac_hisp, frac_pov99, frac_pov199)
    ][, .(pwsid, pws_name, population_served_count,
                           black_served = frac_black*population_served_count,
                           amerind_served = frac_amerind*population_served_count,
                           hispanic_served = frac_hisp*population_served_count,
                           asian_served = frac_asian*population_served_count,
                           pacisl_served = frac_pacisl*population_served_count,
                           nhw_served = frac_white*population_served_count,
                           lowinc_served = frac_pov199*population_served_count,
                           verylowinc_served = frac_pov99*population_served_count,
                           medhighinc_served = (1-frac_pov199)*population_served_count)]

# Merge six year review w/ pws demos
syr.epic <- merge(syr.together,
                  epic.pops[!is.na(black_served)][!is.na(lowinc_served)],
                  by.x = 'pws_id',
                  by.y = 'pwsid')[CAS_REGISTRY_NUM == '57-74-9',
                                  FULL_NAME := 'CHLORDANE'
                                  ][CAS_REGISTRY_NUM == '118-74-1',
                                    FULL_NAME := 'HEXACHLOROBENZENE'
                                    ][CAS_REGISTRY_NUM == '15972-60-8',
                                      FULL_NAME := 'ALACHOR']

# Calculate average concentration stats by subpopulation. 
syr.natlstats <- syr.epic[, .(mean_concentration_black = weighted.mean(mean_result, black_served, na.rm = T),
                              mean_concentration_amerind = weighted.mean(mean_result, amerind_served, na.rm = T),
                              mean_concentration_hispanic = weighted.mean(mean_result, hispanic_served, na.rm = T),
                              mean_concentration_asian = weighted.mean(mean_result, asian_served, na.rm = T),
                              mean_concentration_pacisl = weighted.mean(mean_result, pacisl_served, na.rm = T),
                              mean_concentration_nhw = weighted.mean(mean_result, nhw_served, na.rm = T),
                              mean_concentration_lowinc = weighted.mean(mean_result, lowinc_served, na.rm = T),
                              mean_concentration_verylowinc = weighted.mean(mean_result, verylowinc_served, na.rm = T),
                              mean_concentration_medhighinc = weighted.mean(mean_result, medhighinc_served, na.rm = T),
                              contributing_pws = .N),
                          by = .(FULL_NAME, CAS_REGISTRY_NUM, vintage)]

data.table::fwrite(syr.natlstats, 'C:/Users/atheisin/OneDrive - Environmental Protection Agency (EPA)/Documents/GitHub/ubiquitous_carcinogens/output/dw_demo_summary.csv')
