## State chemical coverage summary:
rm(list = ls())
options(scipen=999)

library(tidyverse)
library(data.table)
library(readxl)

setwd('C:/Users/atheisin/OneDrive - Environmental Protection Agency (EPA)/Documents/GitHub/ubiquitous_carcinogens')

# Carcinogen working list
carc.list <- data.table::fread('./output/merged_carcinogen_list.csv')

# California chems
ca.dta <- unique(
  data.table::as.data.table(
  readxl::read_excel('./raw_data/state_biomonitoring/BiomonitoringCA_Inventory.xlsx')
)[, .(Chemical = `Chemical measured - abbreviation`,
      CAS = `CAS Number`,
      `Indicates Exposure to`)
  ][Chemical %like% 'PCB', Chemical := 'PCBs'
  ][Chemical %like% 'PCB', CAS := '1336-36-3'
  ][CAS == '14797-73-0', Chemical := 'Perchlorate'
  ][CAS %in% c('1767-28-8','1732-29-2','5522-43-0'), CAS := 'DIESEL' #metabolite
  ][`Indicates Exposure to` == 'Di-2-ethylhexyl phthalate (DEHP)',
    CAS := '117-81-7'
  ][`Indicates Exposure to` == 'Diethyl phthalate (DEP)',
    CAS := '84-66-2'
  ][`Indicates Exposure to` == 'Dicyclohexyl phthalate (DCHP)',
    CAS := '7517-36-4'
  ][`Indicates Exposure to` == 'Benzylbutyl phthalate (BzBP)',
    CAS := '85-68-7'
  ][, `Indicates Exposure to` := NULL
    ][CAS %in% carc.list$CASRN | is.na(CAS)
  ][, Chemical := NULL
  ])[, State := 'California'
    ][, `Monitoring status` := 'Previous & In progress'
      ]

#Washington chems
was.dta <- unique(data.table::fread('./output/WA_strat_table.csv'
)[CAS %in% carc.list$CASRN
][, .(CAS)]
)[, State := 'Washington'
  ][, `Monitoring status` := 'Previous'
    ]

#New Jersey chems
nj.dta <- unique(
  merge(data.table::as.data.table(
  readxl::read_excel('./raw_data/state_biomonitoring/NJHANES_Analytes_2019-2024.xlsx')
  ),
  data.table::as.data.table(
    readxl::read_excel('./raw_data/crosswalks/voc_parent_chems.xlsx') #Adjust for VOC metabolites
  )[, .(CAS, `Metabolite CAS`)][!is.na(`Metabolite CAS`)],
  by.x = 'CAS', by.y = 'Metabolite CAS', all.x = T
  )[!is.na(CAS.y), CAS := CAS.y
][CAS %in% carc.list$CASRN
][,.(CAS)]
)[, State := 'New Jersey'
][, `Monitoring status` := 'Previous & In progress'
]

#New Hampshire chems
nh.2019 <- data.table::as.data.table(
  readxl::read_excel('./raw_data/state_biomonitoring/NH_2019study.xlsx')
)[,2:4][, Year := '2014-2019']

nh.2024 <- data.table::as.data.table(
  readxl::read_excel('./raw_data/state_biomonitoring/NH_2024study.xlsx')
)[,2:4][, Year := '2019-2024']

nh.dta <- unique(rbind(nh.2019, nh.2024
)[CAS %in% carc.list$CASRN
  ][,.(CAS)
  ])[, State := 'New Hampshire'
][, `Monitoring status` := 'Previous & In progress'
]
rm(nh.2019, nh.2024)

#Minnesota chems
mn.dta1 <- unique(data.table::as.data.table(
  readxl::read_excel('./raw_data/state_biomonitoring/MN_healthykids_chems_2019-2024.xlsx')
)[CAS %in% carc.list$CASRN
  ][,.(CAS)
  ])[, State := 'Minnesota'
  ][, `Monitoring status` := 'In progress'
  ]

mn.dta2 <- unique(data.table::as.data.table(
  readxl::read_excel('./raw_data/state_biomonitoring/MN_stratifiedresults_2018-2020.xlsx')
)[CAS %in% carc.list$CASRN
  ][, .(CAS)
  ][, State := 'Minnesota'
  ][, `Monitoring status` := 'Previous'
  ])

mn.dta <- rbind(mn.dta2, mn.dta1
                )[, lapply(.SD, paste0, collapse=" & "), by = .(CAS, State)]

#Iowa chems
ia.dta <- unique(data.table::as.data.table(
  readxl::read_excel('./raw_data/state_biomonitoring/IA_2019-2024.xlsx')
)[CAS %in% carc.list$CASRN
  ][,.(CAS)]
)[, State := 'Iowa'
][, `Monitoring status` := 'In progress'
]

#New York chems
ny.dta <- unique(data.table::as.data.table(
    readxl::read_excel('./raw_data/state_biomonitoring/NY_2019-2024study.xlsx')
  )[CAS %in% carc.list$CASRN
  ][, .(CAS)]
)[, State := 'New York'
  ][, `Monitoring status` := 'In progress'
    ]

#Massachusetts chems
ma.carc <- unique(data.table::as.data.table(
  readxl::read_excel('./raw_data/state_biomonitoring/MA_chemlist_2014-2019.xlsx')
)[CAS %in% carc.list$CASRN
][,.(CAS)]
)[, State := 'Massachusetts'
][, `Monitoring status` := 'Previous'
]

#Michigan chems
mi.carc <- unique(data.table::as.data.table(
  readxl::read_excel('./raw_data/state_biomonitoring/MI_2019-2024_monitorlist.xlsx')
)[CAS %in% carc.list$CASRN
][,.(CAS)]
)[, State := 'Michigan'
][, `Monitoring status` := 'In progress'
]


state.dta <- merge(rbind(ca.dta, ia.dta, mn.dta, nh.dta, nj.dta, ny.dta, was.dta, ma.carc, mi.carc),
                   carc.list[, .(NAME, CASRN)],
                   by.x = 'CAS', by.y = 'CASRN')
names(state.dta)[4] <- 'Carcinogen'

fwrite(state.dta, './output/state_chemlist.csv')
