## Clean WA State data

rm(list = ls())
options(scipen=999)

library(tidyverse)
library(data.table)
library(openxlsx)

setwd('C:/Users/atheisin/OneDrive - Environmental Protection Agency (EPA)/Documents/GitHub/ubiquitous_carcinogens')

wa.dta <- data.table::as.data.table(
  openxlsx::read.xlsx('./raw_data/state_biomonitoring/WA_2014-2019_data.xlsx', fillMergedCells = T)
)

wa.dta.clean <- data.table::dcast(wa.dta[Study == 'WEBS (2010-2011)'
][Urinary.Sample.Type == 'Creatinine Corrected'
],
Subtopic + Analyte ~ Percentile.Labels,
value.var = 'MeasureVal2'
)[Analyte == '3-phenoxybenzoic acid (3-PBA)',
  CAS := '3739-38-6'
][Analyte == '4-fluoro-3-phenoxybenzoic acid (4F-3PBA)',
  CAS := '77279-89-1'
][Analyte == 'Antimony',
  CAS := '7440-36-0'
][Analyte == 'Arsenic (V) Acid',
  CAS := '7440-38-2'
][Analyte == 'Arsenous (III) Acid',
  CAS := '13464-58-9'
][Analyte == 'Arsenobetaine',
  CAS := '64436-13-1'
][Analyte == 'Arsenocholine',
  CAS := '39895-81-3'
][Analyte == 'Barium',
  CAS := '7440-39-3'
][Analyte == 'Beryllium',
  CAS := '7440-41-7'
][Analyte == 'Cadmium',
  CAS := '7440-43-9'
][Analyte == 'Cesium',
  CAS := '7440-46-2'
][Analyte == 'Cobalt',
  CAS := '7440-48-4'
][Analyte == 'Dimethylarsinic Acid (DMA)',
  CAS := '75-60-5'
][Analyte == 'Lead',
  CAS := '7439-92-1'
][Analyte == 'Molybdenum',
  CAS := '7439-98-7'
][Analyte == 'Monomethylarsonic acid (MMA)',
  CAS := '124-58-3'
][Analyte == 'Platinum',
  CAS := '7440-06-4'
][Analyte == 'TCPy',
  CAS := '6515-38-4'
][Analyte == 'Thallium',
  CAS := '7440-28-0'
][Analyte == 'Total Arsenic',
  CAS := '7440-38-2'
][Analyte == 'Tungsten',
  CAS := '7440-33-7'
][Analyte == 'Uranium',
  CAS := '7440-61-1'
][Analyte == 'cis-3-(2,2-dibromovinyl)-2,2-dimethylcyclopropane carboxylic acid (DBCA)',
  CAS := '63597-73-9'
][Analyte == 'cis-DCCA',
  CAS := '59042-49-8'
][Analyte == 'trans-3-(2,2-dichlorovinyl)-2,2-dimethylcyclopropane carboxylic acid (trans-DCCA)',
  CAS := '55701-05-8'
][, Analyte := paste0(Analyte,' (μg/g CrCl)')
][Subtopic == '$25,000 or less', Subtopic := 'HH Income: $25,000 or less'
][Subtopic == '$25,000 - $49,999', Subtopic := 'HH Income: $25,000 - $49,999'
][Subtopic == '$50,000 - $74,999', Subtopic := 'HH Income: $50,000 - $74,999'
][Subtopic == '$75,000 or more', Subtopic := 'HH Income: $75,000 or more'
][Subtopic == '6-11 years', Subtopic := 'Age: 6-11 years'
][Subtopic == '12-19 years', Subtopic := 'Age: 12-19 years'
][Subtopic == '20 years and older', Subtopic := 'Age: 20 years and older'
][Subtopic == 'White, Non-Hispanic', Subtopic := 'Race/eth: White, Non-Hispanic' 
][Subtopic == 'Hispanic', Subtopic := 'Hispanic'   
][Subtopic == 'Asian, Non-Hispanic', Subtopic := 'Race/eth: Asian, Non-Hispanic' 
][Subtopic == 'Female', Subtopic := 'Gender: Female' 
][Subtopic == 'Male', Subtopic := 'Gender: Male'
][, Year := '2010-2011']

data.table::setcolorder(wa.dta.clean, neworder = c('Analyte','CAS','Subtopic',
                                                   "50th Percentile",
                                                   "75th Percentile",
                                                   "90th Percentile",
                                                   "95th Percentile",
                                                   "Year"))
wa.dta.clean[wa.dta.clean == '<LOD'] <- NA_character_

data.table::fwrite(wa.dta.clean, './output/WA_strat_table.csv')
