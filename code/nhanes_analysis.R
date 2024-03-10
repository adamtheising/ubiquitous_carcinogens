# Acquire and clean NHANES data
# Author: A. Theising
# Date: 3/6/2024

rm(list = ls())
options(scipen=999)

library(tidyverse)
library(data.table)
library(survey)

setwd('C:/Users/atheisin/OneDrive - Environmental Protection Agency (EPA)/Documents/GitHub/ubiquitous_carcinogens')

# Carcinogen working list
carcin <- data.table::fread('./output/merged_carcinogen_list.csv')

###############################################################
# Use NHANES 1988 - 2018 compiled data from Nguyen et al:
# https://figshare.com/articles/dataset/NHANES_1988-2018/21743372

if('nhanes_chems.csv' %in% list.files('./raw_data/NHANES/')){
  nhanes.chem <- data.table::fread('./raw_data/NHANES/nhanes_chems.csv')
} else {
  nhanes.chem <- data.table::fread('https://figshare.com/ndownloader/files/41730834')
  data.table::fwrite(nhanes.chem, './raw_data/NHANES/nhanes_chems.csv')
}

if('nhanes_coms.csv' %in% list.files('./raw_data/NHANES/')){
  nhanes.coms <- data.table::fread('./raw_data/NHANES/nhanes_coms.csv')
} else {
  nhanes.coms <- data.table::fread('https://figshare.com/ndownloader/files/41730843')
  data.table::fwrite(nhanes.coms, './raw_data/NHANES/nhanes_coms.csv')
}

if('nhanes_demo.csv' %in% list.files('./raw_data/NHANES/')){
  nhanes.demo <- data.table::fread('./raw_data/NHANES/nhanes_demo.csv')
} else {
  nhanes.demo <- data.table::fread('https://figshare.com/ndownloader/files/41730849')
  data.table::fwrite(nhanes.demo, './raw_data/NHANES/nhanes_demo.csv')
}

if('nhanes_weights.csv' %in% list.files('./raw_data/NHANES/')){
  nhanes.weights <- data.table::fread('./raw_data/NHANES/nhanes_weights.csv')
} else {
  nhanes.weights <- data.table::fread('https://figshare.com/ndownloader/files/41730915')
  data.table::fwrite(nhanes.weights, './raw_data/NHANES/nhanes_weights.csv')
}

if('nhanes_dict.csv' %in% list.files('./raw_data/NHANES/')){
  nhanes.dict <- data.table::fread('./raw_data/NHANES/nhanes_dict.csv')
} else {
  nhanes.dict <- data.table::fread('https://figshare.com/ndownloader/files/41731260')
  data.table::fwrite(nhanes.dict, './raw_data/NHANES/nhanes_dict.csv')
}

################################
# Merge to massive DF
nhanes.merged <- dplyr::full_join(nhanes.demo,
                                  nhanes.chem,
                                  by = c("SEQN",
                                         "SEQN_new", 
                                         "SDDSRVYR")) %>%
  dplyr::full_join(.,
                   nhanes.weights,
                   by = c("SEQN",
                          "SEQN_new", 
                          "SDDSRVYR")) %>%
  dplyr::full_join(.,
                   nhanes.coms,
                   by = c("SEQN",
                          "SEQN_new",
                          "SDDSRVYR"))

rm(nhanes.demo,nhanes.chem,nhanes.weights,nhanes.coms)

###############################################
## List carcinogens that match in NHANES
cas.match <- nhanes.dict[cas_num %in% unique(carcin$CASRN)]
# unique(cas.match$cas_num) #84 unique chems

## Loop through all carcinogens on list.
return.me <- vector(mode = 'list', length = dim(cas.match)[1])
options(survey.lonely.psu="adjust")
for (i in 1:dim(cas.match)[1]){
  ###############################################
  # Subset data for analysis
  
  # Calculating values considering the LOD: Concentrations less than the LOD are 
  # assigned a value equal to the LOD divided by the square root of two for 
  # calculation of geometric means. If the proportion of results below the LOD 
  # was greater than 40%, geometric means were not calculated.
  # Note LOD is attained if comment column == 0, if == 1, then below detection.
  
  subset.x <- nhanes.merged %>% 
    dplyr::filter(!is.na(!!as.symbol(cas.match[i,]$variable_codename_use))) %>%
    dplyr::mutate(maxyear = max(SDDSRVYR)) %>%
    dplyr::filter(SDDSRVYR == maxyear) %>% # Use most recent year
    dplyr::mutate(RIDRETH1 = ifelse(RIDRETH1==2,
                                    1,
                                    RIDRETH1)) %>% # Adjust Mexican-American to Hispanic
    dplyr::mutate(inc.pov = ifelse(INDFMPIR < 2,
                        1, 
                        ifelse(INDFMPIR >= 2,
                               2,
                               3))) %>% #Below or above 2x PL
    dplyr::mutate(inc.pov = tidyr::replace_na(inc.pov,3)) %>%
    dplyr::select(SEQN, SEQN_new, SDDSRVYR,
                  SDMVPSU, SDMVSTRA,
                  RIDAGEYR, RIAGENDR, RIDRETH1, INDFMPIR, inc.pov,
                  dplyr::matches(cas.match[i,]$variable_codename_use),
                  dplyr::all_of(cas.match[i,]$comment_codename_use),
                  maxyear)
  
  ##############################
  ## switch chem var names for ease in using svydesign
  index_chem <- which(colnames(subset.x) == cas.match[i,]$variable_codename_use)
  colnames(subset.x)[index_chem] <- "chem"
  index_weights <- which(colnames(subset.x) == paste0('WT_',cas.match[i,]$variable_codename_use))
  colnames(subset.x)[index_weights] <- "weights"
  index_lodflag <- which(colnames(subset.x) == cas.match[i,]$comment_codename_use)
  colnames(subset.x)[index_lodflag] <- "lodflag"
  rm(index_chem, index_weights,index_lodflag)
  
  ###################################################################
  ## First: quantiles using survey weights
  ## Set chem values to 0, if they land in summary table below, replace with NA
  subset.quan <- subset.x %>%
    dplyr::mutate(chem = ifelse(lodflag == 1,
                                0,
                                chem))
    
  # Set up sampling design:
  dsn.subset <- survey::svydesign(ids = ~SDMVPSU, 
                                  strata = ~SDMVSTRA,
                                  weights = ~weights, 
                                  nest = TRUE, 
                                  data = subset.quan)
  
  ## Merge all subpop quantiles into single df
  together.subset <- rbind(subset.quan %>% #Gender-varying statistics
                             dplyr::reframe(RIAGENDR = survey::svyby(~chem, ~RIAGENDR, dsn.subset, survey::svyquantile, quantiles = 0.5, keep.var=T)[,1],
                                            median = survey::svyby(~chem, ~RIAGENDR, dsn.subset, survey::svyquantile, quantiles = 0.5, keep.var=FALSE)[,2],
                                            perc_75 = survey::svyby(~chem, ~RIAGENDR, dsn.subset, survey::svyquantile, quantiles = 0.75, keep.var=FALSE)[,2],
                                            perc_95 = survey::svyby(~chem, ~RIAGENDR, dsn.subset, survey::svyquantile, quantiles = 0.95, keep.var=FALSE)[,2],
                                            perc_99 = survey::svyby(~chem, ~RIAGENDR, dsn.subset, survey::svyquantile, quantiles = 0.99, keep.var=FALSE)[,2]) %>%
                             dplyr::mutate(cat = ifelse(RIAGENDR == 1,
                                                        "Gender: male", 
                                                        "Gender: female")) %>%
                             dplyr::select(-RIAGENDR) %>%
                             dplyr::relocate(cat),
                           # Race/Ethn varying statistics
                           subset.quan %>%
                             dplyr::reframe(RIDRETH1 = survey::svyby(~chem, ~RIDRETH1, dsn.subset, survey::svyquantile, quantiles = 0.5, keep.var=T)[,1],
                                            median = survey::svyby(~chem, ~RIDRETH1, dsn.subset, survey::svyquantile, quantiles = 0.5, keep.var=T)[,2],
                                            perc_75 = survey::svyby(~chem, ~RIDRETH1, dsn.subset, survey::svyquantile, quantiles = 0.75, keep.var=FALSE)[,2],
                                            perc_95 = survey::svyby(~chem, ~RIDRETH1, dsn.subset, survey::svyquantile, quantiles = 0.95, keep.var=FALSE)[,2],
                                            perc_99 = survey::svyby(~chem, ~RIDRETH1, dsn.subset, survey::svyquantile, quantiles = 0.99, keep.var=FALSE)[,2]) %>%
                             mutate(cat = ifelse(RIDRETH1 %in% c(1,2),
                                                 "Race/eth: Hispanic", 
                                                 ifelse(RIDRETH1 == 3,
                                                        "Race/eth: Non-hispanic white", 
                                                        ifelse(RIDRETH1 == 4,
                                                               "Race/eth: Non-hispanic black", 
                                                               ifelse(RIDRETH1 == 5,
                                                                      "Race/eth: Other",
                                                                      "Race/eth: Missing"))))) %>%
                             dplyr::select(-RIDRETH1) %>%
                             dplyr::relocate(cat),
                           # Income varying statistics
                           subset.quan %>%
                             dplyr::reframe(inc.pov = survey::svyby(~chem, ~inc.pov, dsn.subset, survey::svyquantile, quantiles = 0.5, keep.var=T)[,1],
                                            median = survey::svyby(~chem, ~inc.pov, dsn.subset, survey::svyquantile, quantiles = 0.5, keep.var=FALSE)[,2],
                                            perc_75 = survey::svyby(~chem, ~inc.pov, dsn.subset, survey::svyquantile, quantiles = 0.75, keep.var=FALSE)[,2],
                                            perc_95 = survey::svyby(~chem, ~inc.pov, dsn.subset, survey::svyquantile, quantiles = 0.95, keep.var=FALSE)[,2],
                                            perc_99 = survey::svyby(~chem, ~inc.pov, dsn.subset, survey::svyquantile, quantiles = 0.99, keep.var=FALSE)[,2]) %>%
                             dplyr::mutate(inc.pov = unique(subset.quan$inc.pov)) %>%
                             mutate(cat = ifelse(inc.pov == 1,
                                                 "Income: PIR below 2", 
                                                 ifelse(inc.pov == 2,
                                                        "Income: PIR above 2", 
                                                        "Income: Missing"))) %>%
                             dplyr::select(-inc.pov) %>%
                             dplyr::relocate(cat)
                           ) %>%
    dplyr::mutate(chem = cas.match[i,]$variable_description_use) %>%
    dplyr::mutate(cas = cas.match[i,]$cas_num) %>%
    dplyr::mutate(year = ifelse(unique(subset.quan$maxyear) == 10,
                                '2017-2018',
                                ifelse(unique(subset.quan$maxyear) == 9,
                                       '2015-2016',
                                       ifelse(unique(subset.quan$maxyear) == 8,
                                              '2013-2014',
                                              ifelse(unique(subset.quan$maxyear) == 7,
                                                     '2011-2012',
                                ifelse(unique(subset.quan$maxyear) == 6,
                                       '2009-2010',
                                       ifelse(unique(subset.quan$maxyear) == 5,
                                              '2007-2008',
                                              ifelse(unique(subset.quan$maxyear) == 4,
                                                     '2005-2006',
                                                     ifelse(unique(subset.quan$maxyear) == 3,
                                                            '2003-2004',
                                ifelse(unique(subset.quan$maxyear) == 2,
                                       '2001-2002',
                                       ifelse(unique(subset.quan$maxyear) == 1,
                                              '1999-2000',
                                              '1988-1994')))))))))))  %>% 
    dplyr::mutate(dplyr::across(c(median, perc_99), ~ na_if(., 0)))
  
  rm(subset.quan, dsn.subset)
  
  ###################################################################
  ## Second: geometric mean using survey weights
  ## If >= 60% of samples are detectable, return mean value sans nondetects. Else return NA.
  if (dim(dplyr::filter(subset.x, lodflag == 1))[1] / dim(subset.x)[1] > 0.4){
    together.mean = data.frame(cat = together.subset$cat,
                               mean = rep(NA_real_, dim(together.subset)[1]))
  } else {
    subset.mean <- subset.x %>%
      dplyr::mutate(chem = ifelse(lodflag == 1,
                                  chem / sqrt(2),
                                  chem))
    
    dsn.subset <- survey::svydesign(ids = ~SDMVPSU, 
                                    strata = ~SDMVSTRA,
                                    weights = ~weights, 
                                    nest = TRUE, 
                                    data = subset.mean)
    
    together.mean <- rbind(subset.mean %>% #Gender-varying statistics
      dplyr::reframe(RIAGENDR = survey::svyby(~log(chem), ~RIAGENDR, dsn.subset, 
                                              survey::svymean, na.rm = T)[,1],
                     mean = exp(survey::svyby(~log(chem), ~RIAGENDR, dsn.subset, 
                                              survey::svymean, na.rm = T)[,2])) %>%
      dplyr::mutate(RIAGENDR = unique(subset.mean$RIAGENDR)) %>%
      dplyr::mutate(cat = ifelse(RIAGENDR == 1,
                                 "Gender: male", 
                                 "Gender: female")) %>%
      dplyr::select(-RIAGENDR) %>%
      dplyr::relocate(cat),
      #Race/ethn varying statistics
      subset.mean %>%
        dplyr::reframe(RIDRETH1 = survey::svyby(~log(chem), ~RIDRETH1, dsn.subset, 
                                                survey::svymean, na.rm = T)[,1],
                       mean = exp(survey::svyby(~log(chem), ~RIDRETH1, dsn.subset, 
                                                survey::svymean, na.rm = T)[,2])) %>%
        dplyr::mutate(RIDRETH1 = unique(subset.mean$RIDRETH1)) %>%
        mutate(cat = ifelse(RIDRETH1 %in% c(1,2),
                            "Race/eth: Hispanic", 
                            ifelse(RIDRETH1 == 3,
                                   "Race/eth: Non-hispanic white", 
                                   ifelse(RIDRETH1 == 4,
                                          "Race/eth: Non-hispanic black", 
                                          ifelse(RIDRETH1 == 5,
                                                 "Race/eth: Other",
                                                 "Race/eth: Missing"))))) %>%
        dplyr::select(-RIDRETH1) %>%
        dplyr::relocate(cat),
      # Income varying statistics
      subset.mean %>%
        dplyr::reframe(inc.pov = survey::svyby(~log(chem), ~inc.pov, dsn.subset, 
                                                survey::svymean, na.rm = T)[,1],
                       mean = exp(survey::svyby(~log(chem), ~inc.pov, dsn.subset,
                                                survey::svymean, na.rm = T)[,2])) %>%
        dplyr::mutate(inc.pov = unique(subset.mean$inc.pov)) %>%
        mutate(cat = ifelse(inc.pov == 1,
                            "Income: PIR below 2", 
                            ifelse(inc.pov == 2,
                                   "Income: PIR above 2", 
                                   "Income: Missing"))) %>%
        dplyr::select(-inc.pov) %>%
        dplyr::relocate(cat)
      )
    rm(subset.mean, dsn.subset)
  }
  
  return.me[[i]] <- merge(together.subset, together.mean,
                             by = 'cat') %>%
    dplyr::relocate(chem, cas, cat, mean) %>%
    dplyr::mutate(frac_nd = dim(dplyr::filter(subset.x, lodflag == 1))[1] / dim(subset.x)[1])
  
  rm(together.subset,together.mean) 
  print(paste0('Iteration ',i,' complete.'))
}

