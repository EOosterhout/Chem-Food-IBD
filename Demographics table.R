
###     SCRIPT: SUMMARY STATISTICS 
###     AUTHOR: EMILY OOSTERHOUT
###     DESCRIPTION: SUMMARY STATISTICS OF IBD/LLDEEP PARTICIPANTS

##Contents of this file

## 0. LOAD CLEANED DATA
## 1. CREATE TABLE FUNCTIONS
## 2. CONVERT INTO A TABLE GROUPED BY DIAGNOSIS
## 3. COMPUTE P-VALUES FOR GROUP DIFFERENCES
## 4. ADD P-VALUES TO THE TABLE

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
library(qwraps2) #package for making tables

##=============== LOAD RAW DATA ================##

MD <-read_xlsx("Metadata_raw.xlsx")
MD <- as.data.frame(MD)
chem_part <- read_xlsx('chem_raw_participant.xlsx')


##=============== CLEAN DATA ================##

#List of participants in chem dataset
participants <- as.character(chem_part$UMCGIBDResearchIDorLLDeepID)
MD <- filter(MD, MD$UMCGIBDResearchIDorLLDeepID %in% participants)

# Subset full participant intake to IBD and LLDEEP
IBD <- MD[MD$UMCGIBDResearchIDorLLDeepID %like% "UMCGIBD",]
LLDEEP <- MD[MD$UMCGIBDResearchIDorLLDeepID %like% "LLDeep",]

# Identify IBD participants
IBD_ID <- as.character(IBD$UMCGIBDResearchIDorLLDeepID)
# Create a new column specifying the group for each sample
MD$group <- ifelse(MD$UMCGIBDResearchIDorLLDeepID %in% IBD_ID, "IBD", "LLDEEP")
write_xlsx(MD, path = 'metadata_clean.xlsx')


##=============== CREATE TABLE FUNCTIONS ================##

summary_statistics <-
  list( 
    "Age" = 
      list ("Mean (sd)" = ~ qwraps2::mean_sd(AgeAtFecalSampling, na_rm = TRUE, 
                                               denote_sd = "paren"),
            "Minimum" = ~min(AgeAtFecalSampling, na.rm = TRUE),
            "Maximum" = ~max(AgeAtFecalSampling, na.rm = TRUE)
      ),
    "Gender, n (%)" = 
      list(
        "male" = ~qwraps2::n_perc(Sex == "male", digits = 1, na_rm = TRUE),
        "female" = ~qwraps2::n_perc(Sex == "female", digits = 1, na_rm = TRUE)
      ),
    "BMI" =
      list(
        "Median (Q1, Q3)" = ~qwraps2::median_iqr(BMI, na_rm = TRUE),
        "Minimum" = ~min(BMI, na.rm = TRUE),
        "Maximum" = ~max(BMI, na.rm = TRUE)
      ),
    "Smoking status" =
      list(
        "Smoker" = ~qwraps2::n_perc(SmokeCurrentSmoker == "2", digits = 1, na_rm = TRUE),
        "Non-smoker" = ~qwraps2::n_perc(SmokeCurrentSmoker == "1", digits = 1, na_rm = TRUE)
    ),
    "Medication types used, n (%)" =
      list(
        "Mesalazines" = ~qwraps2::n_perc(MedicationMesalazines == "2", digits = 1, na_rm = TRUE),
        "Steroids" = ~qwraps2::n_perc(MedicationSteroids == "2", digits = 1, na_rm = TRUE),
        "Immunosuppressants" = ~qwraps2::n_perc(MedicationImmu0suppressants == "2", digits = 1, na_rm = TRUE),
        "AntiTNF" = ~qwraps2::n_perc(MedicationAntiTNF == "2", digits = 1, na_rm = TRUE),
        "Antibiotics" = ~qwraps2::n_perc(MedicationAntibiotics == "2", digits = 1, na_rm = TRUE),
        "Vitamins" = ~qwraps2::n_perc(MedicationVitamins == "2", digits = 1, na_rm = TRUE),
        "Minerals" = ~qwraps2::n_perc(MedicationMinerals == "2", digits = 1, na_rm = TRUE),
        "PPI" = ~qwraps2::n_perc(MedicationPPI == "2", digits = 1, na_rm = TRUE)
      )
    )
    

##=============== CONVERT INTO A TABLE IBD vs LLDEEP  ================##

table_by_diagnosis <- summary_table(MD %>% dplyr::group_by(group), summary_statistics)


##=============== COMPUTE P-VALUES FOR GROUP DIFFERENCES  ================##

#Compute the p-values for categorical data (Fisher-Z test)
Gender_pvalue <- frmtp(fisher.test(MD$group, MD$Sex)$p.value)
Smoking_pvalue <- frmtp(fisher.test(MD$group, MD$SmokeCurrentSmoker)$p.value)

Mesalazines_pvalue <- frmtp(fisher.test(MD$group, MD$MedicationMesalazines)$p.value)
Steroids_pvalue <- frmtp(fisher.test(MD$group, MD$MedicationSteroids)$p.value)
Immunosuppressants_pvalue <- frmtp(fisher.test(MD$group, MD$MedicationImmu0suppressants)$p.value)
AntiTNF_pvalue <- frmtp(fisher.test(MD$group, MD$MedicationAntiTNF)$p.value)
Antibiotics_pvalue <- frmtp(fisher.test(MD$group, MD$MedicationAntibiotics)$p.value)
Vitamins_pvalue <- frmtp(fisher.test(MD$group, MD$MedicationVitamins)$p.value)
Minerals_pvalue <- frmtp(fisher.test(MD$group, MD$MedicationMinerals)$p.value)
PPI_pvalue <- frmtp(fisher.test(MD$group, MD$MedicationPPI)$p.value)

#Compute the p-values for normally distributed continuous data (t-test)
Age_pvalue <- frmtp(t.test(AgeAtFecalSampling ~ group, data = MD, paired=FALSE)$p.value)

#Compute the p-values for non-normally distributed continuous data (Wilcoxon Rank Sum test)
BMI_pvalue <- frmtp(wilcox.test(BMI ~ group, data = MD, paired=FALSE)$p.value)


##=============== ADD P-VALUES TO THE TABLE  ================##

#Add an additional column in the table for the p-values
table_by_diagnosis <- cbind(table_by_diagnosis, "P-value" = "")

#add the p-values that should be placed next to the summary statistics
table_by_diagnosis[grepl("Mean..sd.", rownames(table_by_diagnosis)), "P-value"] <- Age_pvalue
table_by_diagnosis[grepl("Mesalazines", rownames(table_by_diagnosis)), "P-value"] <- Mesalazines_pvalue
table_by_diagnosis[grepl("Steroids", rownames(table_by_diagnosis)), "P-value"] <- Steroids_pvalue
table_by_diagnosis[grepl("Immunosuppressants", rownames(table_by_diagnosis)), "P-value"] <- Immunosuppressants_pvalue
table_by_diagnosis[grepl("AntiTNF", rownames(table_by_diagnosis)), "P-value"] <- AntiTNF_pvalue
table_by_diagnosis[grepl("Antibiotics", rownames(table_by_diagnosis)), "P-value"] <- Antibiotics_pvalue
table_by_diagnosis[grepl("Vitamins", rownames(table_by_diagnosis)), "P-value"] <- Vitamins_pvalue
table_by_diagnosis[grepl("Minerals", rownames(table_by_diagnosis)), "P-value"] <- Minerals_pvalue
table_by_diagnosis[grepl("PPI", rownames(table_by_diagnosis)), "P-value"] <- PPI_pvalue


#add the p-values that should be placed next to the group row names
printed_table <- capture.output(print(table_by_diagnosis))
printed_table[grepl("Gender", printed_table)] <-
  sub("&nbsp;&nbsp;\\ \\|$", paste(Gender_pvalue, "|"), printed_table[grepl("Gender", printed_table)])
printed_table[grepl("BMI", printed_table)] <-
  sub("&nbsp;&nbsp;\\ \\|$", paste(BMI_pvalue, "|"), printed_table[grepl("BMI", printed_table)])
printed_table[grepl("Smoking status", printed_table)] <-
  sub("&nbsp;&nbsp;\\ \\|$", paste(Smoking_pvalue, "|"), printed_table[grepl("Smoking status", printed_table)])

#display the table (only works properly when knitting in RMarkdown)
cat(printed_table, sep = "\n")

