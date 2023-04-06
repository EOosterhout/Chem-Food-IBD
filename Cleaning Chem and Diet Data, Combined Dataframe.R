##Title: Cleaning Chem and Diet Data, Combined Dataframe
##Author: Emily Oosterhout
## setwd('C:/DATA FOOD COMPONENT ANALYSIS/RP2_ChemIBDFood')
##

##========================= LOAD DATA AND IMPORT LIBRARIES =============================##

#Import libraries
library(readxl)
library(dplyr)
library(tidyverse)
library(data.table)
library(reshape2)
library(ggplot2)
library(ggbreak) 
library(patchwork)
library(writexl)

# Load data
chemdata <- read_xlsx('transposed_ffqchem.xlsx')
ffqdata_raw <- read_xlsx('diet_raw_V2.xlsx')
ffqgroup_explanation <- read_xlsx('ffqgroups_explanation.xlsx')

##===================== FILTERING AND MERGING OF DIET DATA ====================##

#Set all ffqgroup names as character
ffqgroups <- as.character(ffqgroup_explanation$ffq_group)
nodata <- as.character(ffqgroup_explanation$no_chem_data) #ffqgroups that don't have chem data
mergeitem <- as.character(ffqgroup_explanation$merge_item) # merged groups

# new df containing only the portion info from diet data
portion <- ffqdata_raw[, c(ffqgroups)]
portion <- portion[, !(colnames(portion) %in% nodata)]

# new df containing the elements of merged groups
merge_groups <- ffqgroup_explanation[, 4:11]

# Get portion data per merge_group
yogurt_merge <- portion[, (colnames(portion) %in% merge_groups$yogurt_merge)]
coffeecreamer_merge <- portion[, (colnames(portion) %in% merge_groups$coffeecreamer_merge)]
meat_merge <- portion[, (colnames(portion) %in% merge_groups$meat_merge)]
pork_merge <- portion[, (colnames(portion) %in% merge_groups$pork_merge)]
vegetables_fat <- portion[, (colnames(portion) %in% merge_groups$vegetables_fat)]
vegetables_nofat <- portion[, (colnames(portion) %in% merge_groups$vegetables_nofat)]
chocolate <- portion[, (colnames(portion) %in% merge_groups$chocolate)]

# Add median column per merged group to portion data
portion$yogurt_merge <- apply(yogurt_merge, 1, median) # 1 means computing per row, 2 means computing per column
portion$coffeecreamer_merge <- apply(coffeecreamer_merge, 1, median)
portion$meat_merge <- apply(meat_merge, 1, median)
portion$pork_merge <- apply(pork_merge, 1, median)
portion$vegetables_fat <- apply(vegetables_fat, 1, median)
portion$vegetables_nofat <- apply(vegetables_nofat, 1, median)
portion$chocolate <- apply(chocolate, 1, median)

# Remove the columns which are merged
portion <- portion[, !(colnames(portion) %in% ffqgroup_explanation$merge_item)]
portion <- portion[!is.na(portion$chicken),] # remove all NA values from dataframe
ffq_final <- as.character(colnames(portion))

##======================================= TRANSPOSING DIET DATA AND MERGE WITH CHEM DATASET =======================##

#Add LL/UMCG ID's to portion data
portion_ID <- cbind(ffqdata_raw[1:1985, 1], portion)
ID <- as.character(portion_ID$UMCGIBDResearchIDorLLDeepID)

# split df based on ID
split_portion <- split(portion_ID, portion_ID$UMCGIBDResearchIDorLLDeepID)

# write function that drops ID column
remove_column <- function(s) {
  s[!(names(s) %in% 'UMCGIBDResearchIDorLLDeepID')]
} 
# Apply remove_column function
drop_ID <- list()
for (y in ID) {
  m <- remove_column(split_portion[[y]])
  drop_ID[[y]] <- m
}

# Apply transpose function to df
ID_contents <- list()
for (y in ID) {
  n <- transpose(drop_ID[[y]])
  n$ffq_group <- colnames(portion)
  ID_contents[[y]] <- n
}

# Set list element names as colnames for V1 column
col_ID <- list()
for (y in ID) {
  t <- setnames(ID_contents[[y]], old = 'V1', 
           new = y)
  col_ID[[y]] <- t
}

# Melt the list to the final transposed dataframe with first column containing ffq_groups
transposed_portion <- Reduce(full_join, col_ID)

# New df with chem data in mg/g, only ffq groups that are known in diet data
chemdata_mg_g <- chemdata[,2:1191]/100 #chemdata divided by 100 --> mg/g
chemdata_mg_g <- cbind(chemdata$ffq_group, chemdata_mg_g)
colnames(chemdata_mg_g)[1] <- 'ffq_group'

# Filter out chem data for known ffq groups in diet data
ffq_diet <- as.character(transposed_portion$ffq_group)
chemdata_full <- chemdata_mg_g[(chemdata$ffq_group %in% ffqgroup_explanation$final_groups),]


# Bind chem and diet data into one df

chem_diet_data <- full_join(chemdata_full, transposed_portion)
write_xlsx(chem_diet_data, 'chem_diet_combined.xlsx')


