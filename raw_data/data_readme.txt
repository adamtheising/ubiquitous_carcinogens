Raw chemical list data: sources and metadata

******* IN ROOT RAW_DATA FOLDER *******

"ChemExpo_bulk_composition_chemicals.csv": Version dated 4/30/2024 from https://comptox.epa.gov/chemexpo/get_data/ [Code will download most recent version]

"ChemExpo_bulk_functional_uses.csv": Version dated 4/30/2024 from https://comptox.epa.gov/chemexpo/get_data/ [Code will download most recent version]

"CompTox_IARC1.csv": Version dated 8/13/2021. Downloaded 2/28/2024 from https://comptox.epa.gov/dashboard/chemical-lists/IARC1

"CompTox_IARC2A.csv": Version dated 8/13/2021. Downloaded 2/28/2024 from https://comptox.epa.gov/dashboard/chemical-lists/IARC2A 

"CompTox_IARC2B.csv": Version dated 8/13/2021. Downloaded 2/28/2024 from https://comptox.epa.gov/dashboard/chemical-lists/IARC2A 

"CompTox_P65.csv": Version dated 7/09/2017. Downloaded 2/28/2024 from https://comptox.epa.gov/dashboard/chemical-lists/OEHHA

"DSSTox_Identifiers_and_CASRN_2021r1.csv": Version dated 10/18/2021 from https://gaftp.epa.gov/COMPTOX/Sustainable_Chemistry_Data/Chemistry_Dashboard/2021/DSSTox_Identifiers_and_CASRN_2021r1.csv

"OEHHA_P65.csv": Version dated 12/29/2023. Downloaded 2/28/2024 from https://oehha.ca.gov/media/downloads//p65chemicalslist.csv and file metadata removed from header. Additional metadata from .csv file: In the Listing Mechanism column, "AB" denotes authoritative bodies, "SQE" denotes State's Qualified Experts, "FR" denotes formally required to be labeled or identified, and "LC" denotes Labor Code. For those chemicals for which the basis for listing documentation is available electronically, a hyperlink to the documentation is provided. The identification number indicated in the following list is the Chemical Abstracts Service (CAS) Registry Number. No CAS number is given when several substances are presented as a single listing. The date refers to the initial appearance of the chemical on the list. For those chemicals for which a no significant risk level (NSRL) for carcinogens or maximum allowable dose level (MADL) for reproductive toxicants has been adopted, it is denoted in the column, "NSRL or MADL." For those NSRLs or MADLs for which the risk assessment documentation is available electronically, a hyperlink to the documentation is provided.

"NTP_ROC15.csv": 15th Report on Carcinogens, released 12/21/2021. Downloaded 2/28/2024 from https://ntp.niehs.nih.gov/ntp/roc/content/roc15_casrn_index.xlsx and converted to csv. Additional metadata from .xlsx file: This Excel file contains a list of Chemical Abstracts Service Registry Numbers (CAS #s) of listed substances for which a CAS number is available, and also indicates those substances listed in the Report on Carcinogens for which no CAS # could be identified. The listing status for each substance, i.e., Known to be a Human Carcinogen (Known) or Reasonably Anticipated to be a Human Carcinogen (RAHC), is indicated in the third column.



******* IN CROSSWALKS FOLDER *******

"CT_2022blockcrosswalk.csv": crosswalk between 2020 CBG IDs and 2022 ACS CBG IDs for CT. From https://github.com/CT-Data-Collaborative/2022-block-crosswalk.
"IN_TSAANLYT_TABLE.xlsx": crosswalk between EPA's drinking water chemical IDs and CAS.
"p65_crosswalk.csv": crosswalk between P65 carcinogen names and ChemTox CAS (author-generated).
"phthalate_parent_chems.xlsx": crosswalk of metabolite <-> parent chems (author-generated).
"roc_crosswalk.csv": crosswalk between ROC carcinogen names and ChemTox CAS (author-generated).
"UMCR_cas_crosswalk.xlsx": crosswalk between UMCR chemical names and CAS (author-generated).
"voc_parent_chems.xlsx": crosswalk of metabolite <-> parent chems (author-generated).



******* IN NHANES FOLDER *******

Harmonized NHANES microdata from Nguyen et al. (2023). The code will download 5 files from https://figshare.com/articles/dataset/NHANES_1988-2018/21743372

"nhanes_chems.csv"
"nhanes_coms.csv"
"nhanes_demo.csv"
"nhanes_dict.csv"
nhanes_weights.csv"



******* IN STATE_BIOMONITORING FOLDER *******

"BiomonitoringCA_Tables.xlsx": Inventory of biomonitoring results posted by Biomonitoring CA. Inventory last updated 10/22/2020. Downloaded 3/9/2024 from https://biomonitoring.ca.gov/export/xls/results/chemicals/all?page&_format=xlsx.

"CA_stratifiedresults.xlsx": author-generated file from state reports, presentations and summary tables shared by Nerissa Wu and Kathleen Attfield (CA Biomonitoring) in March 2024.

"IA_2019-2024.xlsx": author-generated file based on chemicals summary shared on IA state hygeience lab website: https://biomonitoring.shl.uiowa.edu/about-our-program. Accessed March 2024.

"MA_chemlist_2014-2019.xlsx": author-generated file based on chemicals summary shared on 2014-2019 MA biomonitoring program from Meg Blanchet (MA DPH).

"MI_2019-2024_monitorlist.xlsx": author-generated file based on chemicals list shared by Rachel Long (MI DHHS).

"MN_healthykids_chems_2019-2024.xlsx": author-generated file based on spreadsheets and reports shared by Jessica Nelson (MN DoH).

"MN_stratifiedresults_2018-2020.xlsx": author-generated file based on spreadsheets and reports shared by Jessica Nelson (MN DoH).

"NH_2019study.xlsx": author-generated file based on information in 2019 NH biomonitoring report available at https://www.dhhs.nh.gov/sites/g/files/ehbemt476/files/inline-documents/2022-02/trace-studyreport-2019.pdf.

"NH_2024study.xlsx": author-generated file based on information received from Melissa Josefiak (Biomonitoring NH).

"NJ_stratifiedresults_2016-2018.xlsx": author-generated file based on papers and summary tables shared by Chang Yu and Tina Fan (NJ DoH).

"NJHANES_Analytes_2019-2024.xlsx": author-generated file based on papers and summary tables shared by Chang Yu and Tina Fan (NJ DoH).

"NY_2019-2024study.xlsx": author-generated file based information available here: https://www.health.ny.gov/environmental/chemicals/chemicals_and_health/biomonitoring.htm

"WA_2014-2019_data.xlsx": Downloaded from here https://doh.wa.gov/data-statistical-reports/washington-tracking-network-wtn/biomonitoring/biomonitoring-dashboards in March 2024.



******* IN STANFIELD_REANALYSIS FOLDER *******

Files in this folder produced by Zach Stanfield's April 2024 reanalysis of data from Stanfield et al. (2021, EHS). Description of files [also see Stanfield (2024, internal EPA memo)]: 

"See the attached zip file.  It contains 24 files.  

The 18 that have a name like “X_carcinogens_freqItemsets_allSingleChems.csv” are the demographic-specific chemical prevalence tables where X in the name is the demographic group.  The file of these that starts with “all” is including all households and therefore all transactions.  This file was used to create the table in slide 2.  Chemicals are ranked by their prevalence/support.  Some of the files don’t have all 44 chemicals.  That means the ones that are missing were not in any products purchased by households of that demographic.  These files were all used together to build the figure on slide 3.

The 6 files that have a name like “MostPrevelantCarcinogenSets_updt_suppDiff_noTies_X_1e-04sup_0.01conf_size1-1.csv” were used to build the tables on slides 4 and 5.  They list the prevalence and rank for each chemical in 1 of 6 underrepresented population groups compared to the chemical’s prevalence and rank across all households.  They are sorted by the rank difference column from largest positive difference (higher potential exposure in that demographic) to largest negative difference (lower potential exposure)."
