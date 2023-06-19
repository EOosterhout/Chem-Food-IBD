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

blood_mtb <- read.table("data_1442samples_LLD_baseline_1183plasma_metabolites.txt", sep="\t")
mtb_intake <- read_xlsx('chem_raw_participant_V2.xlsx')
raw_metabolites <- read_xlsx("UVGN-0101-18MLTA+ Client Data Tables CORRECTION METABOLON 220125.xlsx", sheet = "OrigScale", col_names = F)
disease_activity <- read.table('IBD_diseaseactivity_Arnau.txt', sep = '\t')
metadata_relevant <- read_xlsx('metadata_c.xlsx')
phen_case_control <- read.delim('phenos_case_control_v1.txt', sep = '\t')
fecal_phen <- read_xlsx('fecal_phen.xlsx')

##======================================== CHANGING IDs FROM UMCG# TO UMCGRESEARCH# ======================##

# File with IDs
IDs_LLD <- read.delim("IDs_LLD.txt", header=FALSE)
IDs_IBD <- read.delim("IDs_IBD.txt", header=FALSE)
transform_ID <- read.delim("transform_ID.txt")
ZIC_IBD <- read_xlsx('ZIC_research_ID.xlsx')
ID_mtb <- read.delim("key_lld_1183meta_arnau.txt", row.names = 1)
HMDB_ids <- read.delim("HMDB_ids.txt")

#function setting first row as header names
header.true <- function(df) {
  names(df) <- as.character(unlist(df[1,]))
  df[-1,]
}
#set first row as column names
disease_activity <- header.true(disease_activity)

# Merge to replace the IDs to UMCG for disease activity
all_new_ID_raw= merge(ZIC_IBD, disease_activity, by.x = "UMCGnoFromZIC", by.y = "UMCGNoFromZIC") #merge umcg ID with data
all_UMCG_flare <- all_new_ID_raw[, c(2,5:7)]

# add column specifying before a flare (yes/no)
all_UMCG_flare$before_a_flare <- ifelse(all_UMCG_flare$Disease_activity_Categorical == 'before a flare', 'yes', 'no')

write_xlsx(all_UMCG_ID, path = 'disease_activity_UMCG.xlsx')

##blood metabolites
ID_mtb$X=NULL
ID_mtb$X.1=NULL
ID_mtb$X.2=NULL
ID_mtb$X.3=NULL
blood_mtb_id=as.data.frame(t(blood_mtb))
blood_mtb_id=merge(ID_mtb,blood_mtb_id, by="row.names", all.y = T)
row.names(blood_mtb_id)=blood_mtb_id$Row.names
blood_mtb_id$Row.names=NULL
blood_mtb_id_t=as.data.frame(t(blood_mtb_id[,c(17:ncol(blood_mtb_id))]))
blood_mtb_id_t$UMCGIBDResearchIDorLLDeepID <- rownames(blood_mtb_id_t)
mtb_blood <- blood_mtb_id_t %>% relocate(UMCGIBDResearchIDorLLDeepID)


##fecal metabolites
colnames(raw_metabolites)=raw_metabolites[1,]
#Remove unnecessary rows
raw_metabolites=raw_metabolites[-c(1),]
raw_metabolites <- as.data.frame(raw_metabolites)
row.names(raw_metabolites)=make.names(raw_metabolites$BIOCHEMICAL)
colnames(raw_metabolites)=gsub(" ", "_", colnames(raw_metabolites))

#Transform Metabolon HMDB into the primary HMDB id - we can use it later to match with blood metabolites- 
faecal_metabolites=data.frame(Faecal=unique(colnames(raw_metabolites[,14:763])), HMDB_primary=NA)
count=1
for (a in faecal_metabolites$Faecal){
  if (a %in% HMDB_ids$secondary_accession){
    primary_id=subset(HMDB_ids,HMDB_ids$secondary_accession==a)[1,1]
    faecal_metabolites[count,2]=primary_id
  } else{
    for (i in 3:ncol(HMDB_ids)){
      if (a %in% HMDB_ids[,i]){
        primary_id=subset(HMDB_ids,HMDB_ids[,i]==a)[1,i]
        faecal_metabolites[count,2]=primary_id  
      }
    }
  }
  count=count+1
}

#Change Metabolon ids to UMCG ids to connect later to phenotypes 
all_raw=as.data.frame(t(raw_metabolites[,14:ncol(raw_metabolites)]))
colnames(all_raw) <- row.names(raw_metabolites)
IDs_IBD_1=IDs_IBD[,c(1,5)]
IDs_LLD$V2=NULL
IDs=data.frame(rbind(as.matrix(IDs_IBD_1), as.matrix(IDs_LLD)))
rownames(IDs)=IDs$V1
IDs$V1=NULL

# Merge to replace the ids Metabolon => UMCG
row.names(all_raw)=gsub(" ","_", row.names(all_raw))
all_new_ID_raw=merge(IDs,all_raw, by='row.names')

colnames(all_new_ID_raw)[2] <- 'UMCGIBDResearchIDorLLDeepID'
all_new_ID_raw <- all_new_ID_raw[,2:1686]

# Convert "NA" values to NA in the character vector
all_new_ID_raw[is.na(all_new_ID_raw)] <- 0

all_new_ID= as.data.frame(sapply(all_new_ID_raw[,2:1685],as.numeric))
all_new_ID <- cbind(all_new_ID_raw$UMCGIBDResearchIDorLLDeepID, all_new_ID)
colnames(all_new_ID)[1] <- 'UMCGIBDResearchIDorLLDeepID'

#Change UMCG ids to research IDs to connect later to phenotypes 
IDs_IBD_2=IDs_IBD[,c(4,5)]

# New df of fecal mtb only IBD
mtb_fec_IBD <- all_new_ID[1:495,]
all_new_fec <- merge(IDs_IBD_2, mtb_fec_IBD, by.x = 'V5', by.y = 'UMCGIBDResearchIDorLLDeepID')
fec_UMCG_ID <- all_new_fec[,c(2:1686)] # select UMCG ID's + fecal data
fec_UMCG_ID_dt= as.data.frame(sapply(fec_UMCG_ID[,2:1685],as.numeric)) #  change to numeric
fec_UMCG_ID <- cbind(fec_UMCG_ID$V4, fec_UMCG_ID_dt) # bind participant ID's to data
colnames(fec_UMCG_ID)[1] <- 'UMCGIBDResearchIDorLLDeepID'

#New df of full fecal metabolites (binding IBD and LLDEEP data)
mtb_fecal <- rbind(fec_UMCG_ID, all_new_ID[496:750,])
write_xlsx(mtb_fecal, path = 'fecal_mtb_full.xlsx')


##=================================================== IDENTIFYING SOURCE OF METABOLITES ================================================================##

#clean compound names, intake metabolites
names(mtb_intake) = gsub(pattern = ":", replacement = "_", x = names(mtb_intake))
names(mtb_intake) = gsub(pattern = " ", replacement = "_", x = names(mtb_intake))
names(mtb_intake) = gsub(pattern = "-", replacement = "_", x = names(mtb_intake))
names(mtb_intake) = gsub(pattern = ",", replacement = "_", x = names(mtb_intake))
names(mtb_intake) = gsub(pattern = '"', replacement = "", x = names(mtb_intake))
names(mtb_intake) = gsub(pattern = "\\|.*", replacement = "", names(mtb_intake))

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

# Identify duplicated column names
duplicated_columns <- duplicated(names(mtb_intake))
# Remove subsequent duplicated columns
clean_intake <- mtb_intake[, !duplicated_columns]

#full df of all metabolites measured (intake, feces, serum)
metabolites <- merge(clean_intake, mtb_fecal, all.x = T)
metabolites <- merge(metabolites, mtb_blood, all.x = T)

##=================================================== FILTERING METADATA TABLE ON RELEVANT INFORMATION =======================================##

#column containing yes/no for calprotectin > 150
metadata_relevant$calprotectin_above150 <- ifelse(metadata_relevant$FecalCalprotectin > 150, 'yes', 'no') 

#column of age
metadata_relevant$age <- rowMeans(metadata_relevant[,c('AgeAtFecalSampling', 'AgeDuringVisit')], na.rm=TRUE)
metadata_relevant$AgeDuringVisit <- NULL
metadata_relevant$AgeAtFecalSampling <- NULL

#Diagnosis column + IBD (yes/no)
metadata_relevant$diagnosis <- ifelse(metadata_relevant$DiagnosisCurrent == 'generalpopulation', 'control', 'IBD')
metadata_relevant$IBD <- ifelse(metadata_relevant$diagnosis == 'IBD' | metadata_relevant$IBDinControls == 'yes', 'yes', 'no')

#phenotypic data of only fecal metabolites
colnames(phen_case_control)[1] <- 'UMCGIBDResearchIDorLLDeepID'

#Change SID column to UMCGIBDResearchIDorLLDeepID
metadata_IBD <- phen_case_control[-c(425:679),]

all_new_phen <- merge(IDs_IBD, metadata_IBD, by.x = 'V5', by.y = 'UMCGIBDResearchIDorLLDeepID')
phen_UMCG_ID <- all_new_phen[,c(2:205)]
colnames(phen_UMCG_ID)[1] <- 'UMCGIBDResearchIDorLLDeepID'
phen_UMCG_ID <- bind_rows(phen_UMCG_ID, phen_case_control[425:679,])

#Diagnosis column + calprotectin (yes/no)
fecal_phen$diagnosis <- ifelse(fecal_phen$diagnosis == 'Control', 'control', 'IBD')
fecal_phen$calprotectin_above150 <- ifelse(fecal_phen$FecalCalprotectin > 150, 'yes', 'no') 
#merge all phenotypic data
phenos <- bind_rows(metadata_relevant, fecal_phen)
phenos <- merge(phenos, all_UMCG_flare, all = T) #add flare data

##=============================================== CREATE FULL ANALYSIS TABLE ==========================##

analysis_table <- merge(phenos, metabolites, by.y = 'UMCGIBDResearchIDorLLDeepID', all.y = T) #metadata + metabolites
# Remove duplicates in the ID column, keeping the first entry
analysis_table <- analysis_table[!duplicated(analysis_table$UMCGIBDResearchIDorLLDeepID), ]

#clean details analysis table
analysis_table$IBD <- ifelse(analysis_table$diagnosis == 'IBD' | analysis_table$IBDinControls == 'yes', 'yes', 'no')
write_xlsx(analysis_table, path = 'analysis_table.xlsx')
