##Title: Creation of analysis table 
##Author: Emily Oosterhout
##

setwd("C:/DATA FOOD COMPONENT ANALYSIS/RP2_ChemIBDFood/analysis_table")

#Import libraries
library(readxl)
library(dplyr)
library(tidyverse)
library(data.table)
library(reshape2)
library(ggplot2)
library(ggbreak) 
library(patchwork)
library(rlang)
library(writexl)

##======================================= LOAD DATA FILES ===========================##

mtb_blood <- read_xlsx("blood_mtb_measured.xlsx")
mtb_intake <- read_xlsx('chem_raw_participant_V2.xlsx')
mtb_fecal <- read_xlsx('fecal_full_mtb_measured.xlsx')
disease_activity <- read.table('IBD_diseaseactivity_Arnau.txt', sep = '\t')
metadata <- read_xlsx('metadata_complete.xlsx')

##======================================== CHANGING IDs FROM UMCG# TO UMCGRESEARCH# ======================##

# File with IDs
IDs_LLD <- read.delim("IDs_LLD.txt", header=FALSE)
IDs_IBD <- read.delim("IDs_IBD.txt", header=FALSE)
transform_ID <- read.delim("transform_ID.txt")
ZIC_IBD <- read_xlsx('ZIC_research_ID.xlsx')

#function setting first row as header names
header.true <- function(df) {
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}
#set first row as column names
disease_activity <- header.true(disease_activity)

# Merge to replace the ids Metabolon => UMCG
all_new_ID_raw= merge(ZIC_IBD, disease_activity, by.x = "UMCGnoFromZIC", by.y = "UMCGNoFromZIC") #merge umcg ID with data
all_UMCG_ID <- all_new_ID_raw[, c(2,5:7)]
write_xlsx(all_UMCG_ID, path = 'disease_activity_UMCG.xlsx')

# New df of fecal mtb only IBD
mtb_fec_IBD <- mtb_fecal[1:495,]
names(mtb_fec_IBD)[1] <- 'UMCGIBDDNAID'
all_new_fec <- merge(ZIC_IBD, mtb_fec_IBD, by = "UMCGIBDDNAID", all.y = T)
# Replace IDs with new IDs OR if NA keep the original ID
all_new_fec$UMCGIBDResearchIDorLLDeepID <- ifelse(is.na(all_new_fec$UMCGIBDResearchIDorLLDeepID), all_new_fec$UMCGIBDDNAID, all_new_fec$UMCGIBDResearchIDorLLDeepID)
fec_UMCG_ID <- all_new_fec[,c(3:1695)]

#New df of full fecal metabolites (binding IBD and LLDEEP data)
mtb_fecal <- rbind(fec_UMCG_ID, mtb_fecal[496:750,])


##=================================================== IDENTIFYING SOURCE OF METABOLITES ================================================================##

#clean compound names, intake metabolites
names(mtb_intake) = gsub(pattern = ":", replacement = "_", x = names(mtb_intake))
names(mtb_intake) = gsub(pattern = " ", replacement = "_", x = names(mtb_intake))
names(mtb_intake) = gsub(pattern = "-", replacement = "_", x = names(mtb_intake))
names(mtb_intake) = gsub(pattern = ",", replacement = "_", x = names(mtb_intake))
names(mtb_intake) = gsub(pattern = '"', replacement = "", x = names(mtb_intake))

#select intake, fecal and serum metabolites columns
intake.cols <- as.character(colnames(mtb_intake[,7:1114]))
fecal.cols <- as.character(colnames(mtb_fecal[,2:1685]))
blood.cols <- as.character(colnames(mtb_blood[,2:1184]))

# Modify column names using paste0()
intake <- paste0("int_", intake.cols)
fecal <- paste0("fec_", fecal.cols)
blood <- paste0("ser_", blood.cols)

# Rename selected columns in the data frame
names(mtb_intake)[names(mtb_intake) %in% intake.cols] <- intake
names(mtb_fecal)[names(mtb_fecal) %in% fecal.cols] <- fecal
names(mtb_blood)[names(mtb_blood) %in% blood.cols] <- blood

##=================================================== FILTERING METADATA TABLE ON RELEVANT INFORMATION =======================================##

#New file with relevant metadata for analysis
metadata_relevant <- select(metadata, c('Sex', 'UMCGIBDResearchIDorLLDeepID', 'UMCGIBDDNAID', 'AgeAtFecalSampling', 'FecalCalprotectin', 'AgeDuringVisit', 'BMI', 'Weight', 'DiagnosisCurrent', 'DiagnosisCurrentCDUCOnly', 'Race', 
                                        'AntibioticsWithin3MonthsPriorToSampling', 'MedicationBiologicals', 'MedicationPPI', 'SmokeCurrentSmoker', 'SUMOFKCAL', 'SUMOFEIWITDIER', 'SUMOFEIWITPLANT', 'laxatives', 'antibiotics_merged', 'Irritable_bowel_syndrome', 'IBDinControls', 'BristolStoolScale'))
write_xlsx(metadata_relevant, path = 'metadata_relevant.xlsx')

#Import new file with relevant metadata
metadata_relevant <- read_xlsx('metadata_relevant.xlsx')

#column containing yes/no for calprotectin > 150
metadata_relevant$calprotectin_above150 <- ifelse(metadata_relevant$FecalCalprotectin > 150, 'yes', 'no') 

#column of age
metadata_relevant$age <- rowMeans(metadata_relevant[,c('AgeAtFecalSampling', 'AgeDuringVisit')], na.rm=TRUE)
metadata_relevant$AgeDuringVisit <- NULL
metadata_relevant$AgeAtFecalSampling <- NULL

#Diagnosis column + IBD (yes/no)
metadata_relevant$diagnosis <- ifelse(metadata_relevant$DiagnosisCurrent == 'generalpopulation', 'control', 'IBD')
metadata_relevant$IBD <- ifelse(metadata_relevant$diagnosis == 'IBD' | metadata_relevant$IBDinControls == 'yes', 'yes', 'no')

##=============================================== CREATE FULL ANALYSIS TABLE ==========================##

analysis_table <- merge(metadata_relevant, mtb_intake, by = 'UMCGIBDResearchIDorLLDeepID', all.y = T)
analysis_table <- merge(analysis_table, mtb_fecal, by = 'UMCGIBDResearchIDorLLDeepID', all.x = T)
analysis_table <- merge(analysis_table, mtb_blood, by = 'UMCGIBDResearchIDorLLDeepID', all.x = T)

write_xlsx(analysis_table, path = 'analysis_table.xlsx')


## OPTIONAL ##
#phenotypic data of only fecal metabolites
phen_case_control <- read.delim('phenos_case_control_v1.txt', sep = '\t')
#Change SID column to UMCGIBDResearchIDorLLDeepID
metadata_IBD <- phen_case_control[-c(425:679),]
all_new_phen <- merge(ZIC_IBD, metadata_IBD, by.x = "UMCGIBDDNAID", by.y = 'SID', all.y = T)
# Replace IDs with new IDs OR if NA keep the original ID
all_new_phen$UMCGIBDResearchIDorLLDeepID <- ifelse(is.na(all_new_phen$UMCGIBDResearchIDorLLDeepID), all_new_phen$UMCGIBDDNAID, all_new_phen$UMCGIBDResearchIDorLLDeepID)
phen_UMCG_ID <- all_new_phen[,c(3:206)]

control_recoded <- read.delim('Controls_phenos_recoded.txt', sep = '\t')
IBD_recoded <- read.delim('IBD_phenos_recoded.txt', sep = '\t')
montreal_phenos <- read.delim('montreal_phenos.txt', sep = '\t')
phen_IBD <- read.delim('phenos_IBD_clean_v2.txt', sep = '\t')
phen_LLD <- read.delim('phenos_LLD_clean_v2.txt', sep = '\t')
phenos_montreal <- read.delim('phenos_montreal.txt', sep = '\t')
## == ##
