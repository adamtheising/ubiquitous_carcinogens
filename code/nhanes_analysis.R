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
## Metabolites
voc.metab <- data.table::as.data.table(
  readxl::read_excel('./raw_data/crosswalks/voc_parent_chems.xlsx')
)

phth.metab <- data.table::as.data.table(
  readxl::read_excel('./raw_data/crosswalks/phthalate_parent_chems.xlsx')
)

metabs <- rbind(voc.metab, phth.metab)
remove(voc.metab, phth.metab)

###############################################
## List carcinogens that match in NHANES
# cas_num keeps track of metabolites' number in NHANES
# CAS keeps track of parent compound number for presentation, etc.
cas.match <- merge(nhanes.dict[cas_num %in% c(unique(carcin$CASRN),
                                        unique(metabs[!is.na(`Metabolite CAS`)]$`Metabolite CAS`))
                                        ][!(cas_num == '')],
                   metabs,
                   by.x = 'cas_num', by.y = 'Metabolite CAS', all.x = T
                   )[is.na(`Parent Compound`), CAS := cas_num 
                     ][, final_description_use := variable_description_use
                       ][!is.na(`Parent Compound`), 
                       final_description_use := paste0(`Parent Compound`,
                                                          ' metabolite: ',
                                                          variable_description_use)
                       ][CAS %in% c(unique(carcin$CASRN))
                         ][, `Parent Compound` := NULL
                           ][, Metabolite := NULL
                             ][final_description_use == 'Mono-isononyl phthalate (ng/mL)', #mislabeled CAS in source data, I believe.
                               final_description_use := 'Di-isononyl phthalate metabolite: Mono-isononyl phthalate (ng/mL)']

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
                        'Inc/PL: < 2', 
                        ifelse(INDFMPIR >= 2,
                               'Inc/PL: >= 2',
                               'Inc/PL: Missing'))) %>% #Below or above 2x PL
    dplyr::mutate(inc.pov3 = ifelse(INDFMPIR < 1,
                                   'Inc/PL: < 1', 
                                   ifelse(INDFMPIR >= 1 & INDFMPIR < 2,
                                          'Inc/PL: between 1 and 2',
                                          ifelse(INDFMPIR >= 2,
                                                 'Inc/PL: >= 2',
                                                 'Inc/PL: Missing')))) %>% #Income/Pov split 3 ways.
    dplyr::mutate(inc.pov = tidyr::replace_na(inc.pov,'Inc/PL: Missing')) %>%
    dplyr::mutate(inc.pov3 = tidyr::replace_na(inc.pov3,'Inc/PL: Missing')) %>%
    dplyr::mutate(gender = ifelse(RIAGENDR == 1,
                                  "Gender: male", 
                                  "Gender: female")) %>%
    dplyr::mutate(raceeth = ifelse(RIDRETH1 %in% c(1,2),
                                   "Race/eth: Hispanic", 
                                   ifelse(RIDRETH1 == 3,
                                          "Race/eth: Non-hispanic white", 
                                          ifelse(RIDRETH1 == 4,
                                                 "Race/eth: Non-hispanic black", 
                                                 ifelse(RIDRETH1 == 5,
                                                        "Race/eth: Other",
                                                        "Race/eth: Missing"))))) %>%
    dplyr::mutate(race_gender = ifelse(RIDRETH1 %in% c(1,2) & RIAGENDR == 1,
                                       'Race/gender: Hispanic, male',
                                       ifelse(RIDRETH1 %in% c(1,2) & RIAGENDR == 2,
                                              'Race/gender: Hispanic, female',
                                              ifelse(RIDRETH1 == 3 & RIAGENDR == 1,
                                                     'Race/gender: White, male',
                                                     ifelse(RIDRETH1 == 3 & RIAGENDR == 2,
                                                            'Race/gender: White, female',
                                ifelse(RIDRETH1 == 4 & RIAGENDR == 1,
                                                     'Race/gender: Black, male',
                                       ifelse(RIDRETH1 == 4 & RIAGENDR == 2,
                                              'Race/gender: Black, female',
                                              ifelse(RIDRETH1 == 5 & RIAGENDR == 1,
                                                     'Race/gender: Other, male',
                                                     ifelse(RIDRETH1 == 5 & RIAGENDR == 2,
                                                            'Race/gender: Other, female',
                                                            ifelse(is.na(RIDRETH1) & RIAGENDR == 1,
                                                                   'Race/gender: Missing, male',
                                                                   'Race/gender: Missing, female')))))))))) %>%
    dplyr::mutate(race_income = ifelse(RIDRETH1 %in% c(1,2) & INDFMPIR < 2,
                                       'Race/income: Hispanic, Inc/PL < 2',
                                       ifelse(RIDRETH1 %in% c(1,2) & INDFMPIR >= 2,
                                              'Race/income: Hispanic, Inc/PL >= 2',
                                              ifelse(RIDRETH1 == 3 & INDFMPIR < 2,
                                                     'Race/income: White, Inc/PL < 2',
                                                     ifelse(RIDRETH1 == 3 & INDFMPIR >= 2,
                                                            'Race/income: White, Inc/PL >= 2',
                                                            ifelse(RIDRETH1 == 4 & INDFMPIR < 2,
                                                                   'Race/income: Black, Inc/PL < 2',
                                       ifelse(RIDRETH1 == 4 & INDFMPIR >= 2,
                                              'Race/income: Black, Inc/PL >= 2',
                                              ifelse(RIDRETH1 == 5 & INDFMPIR < 2,
                                                     'Race/income: Other, Inc/PL < 2',
                                                     ifelse(RIDRETH1 == 5 & INDFMPIR >= 2,
                                                            'Race/income: Other, Inc/PL >= 2',
                                                            ifelse(is.na(RIDRETH1) & INDFMPIR < 2,
                                                                   'Race/income: Missing, Inc/PL < 2',
                                                                   'Race/income: Missing, Inc/PL >= 2')))))))))) %>%
    dplyr::select(SEQN, SEQN_new, SDDSRVYR,
                  SDMVPSU, SDMVSTRA,
                  RIDAGEYR, raceeth, gender, INDFMPIR, inc.pov, inc.pov3,
                  race_gender, race_income,
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
                                chem)) %>%
    dplyr::mutate(lod = ifelse(lodflag != 1,
                               1L,
                               NA_integer_))
    
  # Set up sampling design:
  dsn.subset <- survey::svydesign(ids = ~SDMVPSU, 
                                  strata = ~SDMVSTRA,
                                  weights = ~weights, 
                                  nest = TRUE, 
                                  data = subset.quan)
  
  ## Merge all subpop quantiles into single df
  together.subset <- rbind(
                            #Everyone
                            subset.quan %>%
                             dplyr::summarise(
                               median = survey::svyquantile(~chem, dsn.subset, quantiles = 0.5)[[1]][1],
                               perc_75 = survey::svyquantile(~chem, dsn.subset, quantiles = 0.75)[[1]][1],
                               perc_95 = survey::svyquantile(~chem, dsn.subset, quantiles = 0.95)[[1]][1],
                               perc_99 = survey::svyquantile(~chem, dsn.subset, quantiles = 0.99)[[1]][1],
                               sample_size = n()
                             ) %>%
                             dplyr::mutate(cat = 'All'),
                           #Gender-varying statistics
                           subset.quan %>% 
                             dplyr::reframe(gender = survey::svyby(~chem, ~gender, dsn.subset, survey::svyquantile, quantiles = 0.5, keep.var=T)[,1],
                                            median = survey::svyby(~chem, ~gender, dsn.subset, survey::svyquantile, quantiles = 0.5, keep.var=FALSE)[,2],
                                            perc_75 = survey::svyby(~chem, ~gender, dsn.subset, survey::svyquantile, quantiles = 0.75, keep.var=FALSE)[,2],
                                            perc_95 = survey::svyby(~chem, ~gender, dsn.subset, survey::svyquantile, quantiles = 0.95, keep.var=FALSE)[,2],
                                            perc_99 = survey::svyby(~chem, ~gender, dsn.subset, survey::svyquantile, quantiles = 0.99, keep.var=FALSE)[,2],
                                            sample_size = survey::svyby(~chem, ~gender, dsn.subset, unwtd.count, keep.var=FALSE)[,2]) %>%
                             dplyr::mutate(cat = gender) %>%
                             dplyr::select(-gender) %>%
                             dplyr::relocate(cat),
                           # Race/Ethn varying statistics
                           subset.quan %>%
                             dplyr::reframe(raceeth = survey::svyby(~chem, ~raceeth, dsn.subset, survey::svyquantile, quantiles = 0.5, keep.var=T)[,1],
                                            median = survey::svyby(~chem, ~raceeth, dsn.subset, survey::svyquantile, quantiles = 0.5, keep.var=FALSE)[,2],
                                            perc_75 = survey::svyby(~chem, ~raceeth, dsn.subset, survey::svyquantile, quantiles = 0.75, keep.var=FALSE)[,2],
                                            perc_95 = survey::svyby(~chem, ~raceeth, dsn.subset, survey::svyquantile, quantiles = 0.95, keep.var=FALSE)[,2],
                                            perc_99 = survey::svyby(~chem, ~raceeth, dsn.subset, survey::svyquantile, quantiles = 0.99, keep.var=FALSE)[,2],
                                            sample_size = survey::svyby(~chem, ~raceeth, dsn.subset, unwtd.count, keep.var=FALSE)[,2]) %>%
                             dplyr::mutate(cat = raceeth) %>%
                             dplyr::select(-raceeth) %>%
                             dplyr::relocate(cat),
                           # Income varying statistics
                           subset.quan %>%
                             dplyr::reframe(inc.pov3 = survey::svyby(~chem, ~inc.pov3, dsn.subset, survey::svyquantile, quantiles = 0.5, keep.var=T)[,1],
                                            median = survey::svyby(~chem, ~inc.pov3, dsn.subset, survey::svyquantile, quantiles = 0.5, keep.var=FALSE)[,2],
                                            perc_75 = survey::svyby(~chem, ~inc.pov3, dsn.subset, survey::svyquantile, quantiles = 0.75, keep.var=FALSE)[,2],
                                            perc_95 = survey::svyby(~chem, ~inc.pov3, dsn.subset, survey::svyquantile, quantiles = 0.95, keep.var=FALSE)[,2],
                                            perc_99 = survey::svyby(~chem, ~inc.pov3, dsn.subset, survey::svyquantile, quantiles = 0.99, keep.var=FALSE)[,2],
                                            sample_size = survey::svyby(~chem, ~inc.pov3, dsn.subset, unwtd.count, keep.var=FALSE)[,2]) %>%
                             dplyr::mutate(cat = inc.pov3) %>%
                             dplyr::select(-inc.pov3) %>%
                             dplyr::relocate(cat),
                           # Race/gender varying statistics
                           subset.quan %>%
                             dplyr::reframe(race_gender = survey::svyby(~chem, ~race_gender, dsn.subset, survey::svyquantile, quantiles = 0.5, keep.var=T)[,1],
                                            median = survey::svyby(~chem, ~race_gender, dsn.subset, survey::svyquantile, quantiles = 0.5, keep.var=FALSE)[,2],
                                            perc_75 = survey::svyby(~chem, ~race_gender, dsn.subset, survey::svyquantile, quantiles = 0.75, keep.var=FALSE)[,2],
                                            perc_95 = survey::svyby(~chem, ~race_gender, dsn.subset, survey::svyquantile, quantiles = 0.95, keep.var=FALSE)[,2],
                                            perc_99 = survey::svyby(~chem, ~race_gender, dsn.subset, survey::svyquantile, quantiles = 0.99, keep.var=FALSE)[,2],
                                            sample_size = survey::svyby(~chem, ~race_gender, dsn.subset, unwtd.count, keep.var=FALSE)[,2]) %>%
                             dplyr::mutate(cat = race_gender) %>%
                             dplyr::select(-race_gender) %>%
                             dplyr::relocate(cat),
                           # Race/income varying statistics
                           subset.quan %>%
                             dplyr::reframe(race_income = survey::svyby(~chem, ~race_income, dsn.subset, survey::svyquantile, quantiles = 0.5, keep.var=T)[,1],
                                            median = survey::svyby(~chem, ~race_income, dsn.subset, survey::svyquantile, quantiles = 0.5, keep.var=FALSE)[,2],
                                            perc_75 = survey::svyby(~chem, ~race_income, dsn.subset, survey::svyquantile, quantiles = 0.75, keep.var=FALSE)[,2],
                                            perc_95 = survey::svyby(~chem, ~race_income, dsn.subset, survey::svyquantile, quantiles = 0.95, keep.var=FALSE)[,2],
                                            perc_99 = survey::svyby(~chem, ~race_income, dsn.subset, survey::svyquantile, quantiles = 0.99, keep.var=FALSE)[,2],
                                            sample_size = survey::svyby(~chem, ~race_income, dsn.subset, unwtd.count, keep.var=FALSE)[,2]) %>%
                             dplyr::mutate(cat = race_income) %>%
                             dplyr::select(-race_income) %>%
                             dplyr::relocate(cat)
                           ) %>%
    dplyr::mutate(chem = cas.match[i,]$final_description_use) %>%
    dplyr::mutate(cas = cas.match[i,]$CAS) %>%
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
    dplyr::mutate(dplyr::across(c(median, perc_75, perc_95, perc_99), ~ na_if(., 0)))
  
  rm(subset.quan, dsn.subset)
  
  ###################################################################
  ## Second: geometric mean using survey weights
  ## If >= 60% of samples are detectable, return mean value sans nondetects. Else return NA.
  subset.mean <- subset.x %>%
    dplyr::mutate(chem = ifelse(lodflag == 1,
                                chem / sqrt(2),
                                chem)) %>%
    dplyr::mutate(lod = ifelse(lodflag != 1,
                               1L,
                               NA_integer_))
  
  dsn.subset <- survey::svydesign(ids = ~SDMVPSU, 
                                  strata = ~SDMVSTRA,
                                  weights = ~weights, 
                                  nest = TRUE, 
                                  data = subset.mean)
  
  together.mean <- rbind(
    #Everyone
    subset.mean %>%
      dplyr::summarise(
        mean = exp(survey::svymean(~log(chem), dsn.subset)[1]),
        detected_samples = sum(ifelse(lod == 1, 1, 0), na.rm = T)
      ) %>%
      dplyr::mutate(cat = 'All'),
    #Gender-varying statistics
    subset.mean %>% 
      dplyr::reframe(gender = survey::svyby(~log(chem), ~gender, dsn.subset, 
                                            survey::svymean, na.rm = T)[,1],
                     mean = exp(survey::svyby(~log(chem), ~gender, dsn.subset, 
                                              survey::svymean, na.rm = T)[,2]),
                     detected_samples = survey::svyby(~lod, ~gender, dsn.subset, 
                                              unwtd.count, keep.var=FALSE)[,2]
      ) %>%
      dplyr::mutate(cat = gender) %>%
      dplyr::select(-gender) %>%
      dplyr::relocate(cat),
    #Race/ethn varying statistics
    subset.mean %>%
      dplyr::reframe(raceeth = survey::svyby(~log(chem), ~raceeth, dsn.subset, 
                                             survey::svymean, na.rm = T)[,1],
                     mean = exp(survey::svyby(~log(chem), ~raceeth, dsn.subset, 
                                              survey::svymean, na.rm = T)[,2]),
                     detected_samples = survey::svyby(~lod, ~raceeth, dsn.subset, 
                                              unwtd.count, keep.var=FALSE)[,2]
      ) %>%
      dplyr::mutate(cat = raceeth) %>%
      dplyr::select(-raceeth) %>%
      dplyr::relocate(cat),
    # Income varying statistics
    subset.mean %>%
      dplyr::reframe(inc.pov3 = survey::svyby(~log(chem), ~inc.pov3, dsn.subset, 
                                              survey::svymean, na.rm = T)[,1],
                     mean = exp(survey::svyby(~log(chem), ~inc.pov3, dsn.subset,
                                              survey::svymean, na.rm = T)[,2]),
                     detected_samples = survey::svyby(~lod, ~inc.pov3, dsn.subset, 
                                              unwtd.count, keep.var=FALSE)[,2]) %>%
      dplyr::mutate(cat = inc.pov3) %>%
      dplyr::select(-inc.pov3) %>%
      dplyr::relocate(cat),
    # Race/gender varying statistics
    subset.mean %>%
      dplyr::reframe(race_gender = survey::svyby(~log(chem), ~race_gender, dsn.subset, 
                                                 survey::svymean, na.rm = T)[,1],
                     mean = exp(survey::svyby(~log(chem), ~race_gender, dsn.subset,
                                              survey::svymean, na.rm = T)[,2]),
                     detected_samples = survey::svyby(~lod, ~race_gender, dsn.subset, 
                                              unwtd.count, keep.var=FALSE)[,2]) %>%
      dplyr::mutate(cat = race_gender) %>%
      dplyr::select(-race_gender) %>%
      dplyr::relocate(cat),
    # Race/income varying statistics
    subset.mean %>%
      dplyr::reframe(race_income = survey::svyby(~log(chem), ~race_income, dsn.subset, 
                                                 survey::svymean, na.rm = T)[,1],
                     mean = exp(survey::svyby(~log(chem), ~race_income, dsn.subset,
                                              survey::svymean, na.rm = T)[,2]),
                     detected_samples = svyby(~lod, ~race_income, dsn.subset, 
                                              unwtd.count, keep.var=FALSE)[,2]) %>%
      dplyr::mutate(cat = race_income) %>%
      dplyr::select(-race_income) %>%
      dplyr::relocate(cat)
  )
  rm(subset.mean, dsn.subset)
  
  
  return.me[[i]] <- merge(together.subset, together.mean,
                             by = 'cat') %>%
    dplyr::relocate(chem, cas, cat, mean) %>%
    dplyr::mutate(frac_nd = detected_samples/sample_size) %>%
    dplyr::mutate(mean = ifelse(frac_nd >= 0.6,
                                mean,
                                NA_real_)) %>%
    dplyr::select(-detected_samples)
  
  rm(together.subset,together.mean,subset.x) 
  print(paste0('Iteration ',i,' complete.'))
}

combined <- data.table::rbindlist(return.me
                                  )[chem != '3-fluoranthene (ng/L)'  #not sure why this made it in.
                                    ][chem == 'N-Acetyl-S-(4-hydroxy-2-butenyl)-L-Cysteine (ng/mL)',
                                      chem := '1,3-Butadiene metabolite: N-Acetyl-S-(4-hydroxy-2-butenyl)-L-Cysteine (ng/mL)'] #the metabolite has no CAS, so name didn't match correctly.

data.table::fwrite(combined, './output/nhanes_summary.csv')
