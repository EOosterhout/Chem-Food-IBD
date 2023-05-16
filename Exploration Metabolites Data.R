##Title: Exploring + Cleaning Metabolite data
##Author: Emily Oosterhout
##

setwd("C:/DATA FOOD COMPONENT ANALYSIS/RP2_ChemIBDFood/Paired_metabolomics")

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


##========================================== LOAD AND CLEAN METABOLITE DATA ===============================##

#blood metabolites
blood_mtb=read.table("data_1442samples_LLD_baseline_1183plasma_metabolites.txt", sep="\t")

#identified metabolites ID's
id_mtb=read.delim("key_lld_1183meta_arnau.txt", row.names = 1)
id_mtb$X=NULL
id_mtb$X.1=NULL
id_mtb$X.2=NULL
id_mtb$X.3=NULL
#Merge metabolite measurements with ID information
blood_mtb_id=as.data.frame(t(blood_mtb))
blood_mtb_id=merge(id_mtb,blood_mtb_id, by="row.names", all.y = T)
row.names(blood_mtb_id)=blood_mtb_id$Row.names
blood_mtb_id$Row.names=NULL
#df only containing measured metabolite in plasma per participant
blood_mtb_id_t=as.data.frame(t(blood_mtb_id[,c(17:ncol(blood_mtb_id))]))
write.table(blood_mtb_id_t, file = 'blood_mtb_measured.txt')

#fecal metabolites
raw_metabolites <- read_excel("UVGN-0101-18MLTA+ Client Data Tables CORRECTION METABOLON 220125.xlsx", sheet = "OrigScale", col_names = T)

#transform table to df
raw_metabolites=as.data.frame(raw_metabolites)
row.names(raw_metabolites)=make.names(raw_metabolites$BIOCHEMICAL) #set measured metabolites as rownames
colnames(raw_metabolites)=gsub(" ", "_", colnames(raw_metabolites))

#File with HMDB metabolites ID
HMDB_ids=read.delim("HMDB_ids.txt")

#Transform Metabolon HMDB into the primary HMDB id --> match with blood metabolites- 
faecal_metabolites=data.frame(Faecal=unique(raw_metabolites$Group_HMDB), HMDB_primary=NA)
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

# File with IDs
IDs_LLD <- read.delim("IDs_LLD.txt", header=FALSE)
IDs_IBD <- read.delim("IDs_IBD.txt", header=FALSE)
transform_ID <- read.delim("transform_ID.txt")
#Change Metabolon ids to UMCG ids to connect later to phenotypes 
all_raw=as.data.frame(t(raw_metabolites[,14:ncol(raw_metabolites)])) #transform raw metabolite df
IDs_IBD=IDs_IBD[,c(1,5)]
IDs_LLD$V2=NULL
IDs=data.frame(rbind(as.matrix(IDs_IBD), as.matrix(IDs_LLD)))
rownames(IDs)=IDs$V1
IDs$V1=NULL

# Merge to replace the ids Metabolon => UMCG
row.names(all_raw)=gsub(" ","_", row.names(all_raw))
all_new_ID_raw=merge(IDs,all_raw, by="row.names") #merge umcg ID with data
rownames(all_new_ID_raw)=all_new_ID_raw$V5
all_new_ID_raw$Row.names=NULL
all_new_ID_raw$V5=NULL
all_new_ID=all_new_ID_raw
all_new_ID=as.data.frame(sapply(all_new_ID,as.numeric))
row.names(all_new_ID)=row.names(all_new_ID_raw)
write.table(all_new_ID, file = 'fecal_mtb_measured.txt')

#SCFA metabolites
SCFA <-  read_excel("TA003-19 UVGN_SCFA_CDT_av.xlsx")
SCFA$ori_Result=SCFA$Result
#Replace values below level quantification (BLOQ) to NA
SCFA$Result[SCFA$Comment...13=="BLOQ"] = "N/Q"
SCFA$Result[SCFA$Result=="N/Q"] = NA
#SCFAsub=SCFA[,c("Unique.Sample.ID...as.labeled.on.tube.","Sample.Amount..gram.","Comment","Analyte","Result")]
SCFAsub=SCFA[,c(2,6,8,9)]
colnames(SCFAsub)=make.names(colnames(SCFAsub))
#Reshape the SCFA table
SCFA_table <- dcast(SCFAsub, ...~Analyte)
#Chance sample names to match metadata
rownames(SCFA_table)=SCFA_table$Unique.Sample.ID....as.labeled.on.tube.
SCFA_table$Unique.Sample.ID....as.labeled.on.tube.=NULL
SCFA_new_ID=merge(IDs,SCFA_table, by="row.names")
SCFA_new_ID$Row.names=NULL
row.names(SCFA_new_ID)=SCFA_new_ID$V5
SCFA_new_ID$V5=NULL
colnames(SCFA_new_ID)=c("Amount_sample_gram", "2_Methylbutyric_acid", "acetic_acid", "butyric_acid", "hexanoic_acid", "isobutyric_acid","isovaleric_acid" ,"propionic_acid", "valeric_acid")
#Create df with UMCG ID's and measured SCFA metabolites
SCFA_ID=SCFA_new_ID 
SCFA_ID=as.data.frame(sapply(SCFA_ID,as.numeric))
row.names(SCFA_ID)=row.names(SCFA_new_ID)
SCFA_ID$Amount_sample_gram=NULL
write.table(SCFA_ID, file = 'SCFA_mtb_measured.txt')

#one file with all metabolites measured in feces
fecal_mtb <- cbind(all_new_ID, SCFA_ID)
write.table(fecal_mtb, file = 'fecal_full_mtb_measured.txt')

##=================================== LOAD INTAKE METABOLITES AND FILTER PAIRED METABOLOMICS ===========================##

#intake data
chem <- read_xlsx('chem_raw_participant_V2.xlsx')

LLD_chem <- chem[chem$UMCGIBDResearchIDorLLDeepID %like% 'LLDeep',] #only selecting LLDEEP samples
row.names(LLD_chem)=LLD_chem$UMCGIBDResearchIDorLLDeepID

#samples with paired metabolomics data
Samples_LLD=row.names(blood_mtb)[(row.names(blood_mtb) %in% row.names(fecal_mtb))]
#samples with intake, fecal and blood metabolomics data
Complete_LLD <- chem$UMCGIBDResearchIDorLLDeepID[chem$UMCGIBDResearchIDorLLDeepID %in% Samples_LLD]
#subset intake, fecal and blood data
LLD_mtbfecal <- subset(all_new_ID, row.names(all_new_ID) %in% Complete_LLD)
write.csv(LLD_mtbfecal, file = 'LLD_pairedmtb_fecal.csv')
LLD_mtbblood <- subset(blood_mtb, row.names(blood_mtb) %in% Complete_LLD)
write.csv(LLD_mtbblood, file = 'LLD_pairedmtb_blood.csv')
LLD_mtbintake <- subset(chem, chem$UMCGIBDResearchIDorLLDeepID %in% Complete_LLD)
write.csv(LLD_mtbintake, file = 'LLD_pairedmtb_intake.csv')



