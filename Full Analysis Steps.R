##Title: Full Analysis 
##Author: Emily Oosterhout
##

setwd("C:/DATA FOOD COMPONENT ANALYSIS/RP2_ChemIBDFood")

#Import libraries
library(readxl)
library(dplyr)
library(tidyverse)
library(data.table)
library(reshape2)
library(patchwork)
library(rlang)
library(writexl)

#packages for testing/visualizing data distribution
library(ggplot2)
library(ggbreak)
library(ggpubr)
library(FSA)
library(ggsignif)
library(ggrepel)
library(RColorBrewer)
library(gridExtra)
library(jtools)
library(VennDiagram)
library(GGally)
library(MetBrewer)

#packages for linear regression/statistical testing
library(stats)
library(outliers)

#packages for correlation analysis
library(Hmisc)
library(reshape2)
library(RcmdrMisc)
library(psych)
library(corrplot)
library(gdata)
library (plyr)
library(foreach) 
library(ppcor)

## Functions ##

# Function to set outliers as NA within each column using Tukey's fences
remove_outliers <- function(data, multiplier = 1.5) {
  cleaned_data <- data
  
  # Loop through each column in the dataframe
  for (col in names(data)) {
    # Calculate the lower and upper fences
    q1 <- quantile(data[[col]], 0.25, na.rm = TRUE)
    q3 <- quantile(data[[col]], 0.75, na.rm = TRUE)
    iqr <- q3 - q1
    lower_fence <- q1 - multiplier * iqr
    upper_fence <- q3 + multiplier * iqr
    
    # Identify values outside the fences
    outliers <- data[[col]] < lower_fence | data[[col]] > upper_fence
    
    # Set outliers as NA in the dataframe
    cleaned_data[[col]][outliers] <- NA
  }
  
  return(cleaned_data)
}


##=========================================== LOAD DATA, CLEANING NAMES AND SUBSETTING ===============================

data_full <- as.data.frame(read_xlsx('analysis_table.xlsx'))
row.names(data_full) <- data_full$UMCGIBDResearchIDorLLDeepID

# Subset full dataset on plausible intake (sex dependent Willet)

#Males
intake_male <- subset(data_full, data_full$sex == 'male')
intake_male_high <- subset(intake_male, intake_male["SUMOFKCAL"] > 4000)

intake_male <- subset(data_full, data_full$sex == 'male')
intake_male_low <- subset(intake_male, intake_male["SUMOFKCAL"] < 800)

#Females
intake_female <- subset(data_full, data_full$sex == 'female')
intake_female_high <- subset(intake_female, intake_female["SUMOFKCAL"] > 3500)

intake_female <- subset(data_full, data_full$sex == 'female')
intake_female_low <- subset(intake_female, intake_female["SUMOFKCAL"] < 500)

#Merge participants with implausible intake (25 participants)
implausible_intake_male <- rbind(intake_male_high, intake_male_low)
implausible_intake_female <- rbind(intake_female_high, intake_female_low)
implausible_intake <- rbind(implausible_intake_male, implausible_intake_female)

# Remove participants from original dataframe based on presence in imlausible intake
plausible_intake <- data_full[!(rownames(data_full) %in% rownames(implausible_intake)), ]

# Get metabolite data from analysis table using participant list from raw metabolite files
intake <- as.data.frame(read_xlsx('chem_raw_participant_V2.xlsx'))
fecal <- as.data.frame(read_xlsx('fecal_mtb_full.xlsx'))
serum <- as.data.frame(read_xlsx("blood_mtb.xlsx"))

#participant ID's as rownames, for matching correlations
rownames(intake) <- intake$UMCGIBDResearchIDorLLDeepID
rownames(fecal) <- fecal$UMCGIBDResearchIDorLLDeepID
rownames(serum) <- serum$UMCGIBDResearchIDorLLDeepID

#select participant IDs from raw data files
participants_intake <- as.character(intake$UMCGIBDResearchIDorLLDeepID)
participants_fecal <- as.character(fecal$UMCGIBDResearchIDorLLDeepID)
participants_serum <- as.character(serum$UMCGIBDResearchIDorLLDeepID)

#filter full analysis table on metabolite type
intake_mtb <- plausible_intake[plausible_intake$UMCGIBDResearchIDorLLDeepID %in% participants_intake,]
fecal_mtb <- plausible_intake[plausible_intake$UMCGIBDResearchIDorLLDeepID %in% participants_fecal,]
serum_mtb <- plausible_intake[plausible_intake$UMCGIBDResearchIDorLLDeepID %in% participants_serum,]

##========================================== CLEAN INTAKE METABOLITE NAMES TO WORK IN REGRESSION MODEL ===========================================##

intake <- intake[!(rownames(intake) %in% rownames(implausible_intake)), ]

#clean compound names, intake metabolites
names(intake) = gsub(pattern = ":", replacement = "_", x = names(intake))
names(intake) = gsub(pattern = " ", replacement = "_", x = names(intake))
names(intake) = gsub(pattern = "-", replacement = "_", x = names(intake))
names(intake) = gsub(pattern = ",", replacement = "", x = names(intake))
names(intake) = gsub(pattern = '"', replacement = "", x = names(intake))
names(intake) = gsub(pattern = "\\|.*", replacement = "", x = names(intake))
names(intake) = gsub(pattern = "\\(", replacement = "", x = names(intake))
names(intake) = gsub(pattern = "\\)", replacement = "", x = names(intake))
names(intake) = gsub(pattern = "\\+", replacement = "pos", x = names(intake))
names(intake) = gsub(pattern = "\\'", replacement = "", x = names(intake))
names(intake) = gsub(pattern = "\\â±", replacement = "", x = names(intake))

# Get current column names
old_colnames <- colnames(intake)

# Replace numeric column names with letters
new_colnames <- make.names(old_colnames, unique = TRUE)

# Assign new column names to the data frame
colnames(intake) <- new_colnames

row.names(intake) <- intake$UMCGIBDResearchIDorLLDeepID

##===================================== DIFFERENTIAL ABUNDANCE ANALYSIS: INTAKE_IBDvsNon-IBD (LOG TRANSFORMED AND FILTERED DATA) ====================================================

# Columns containing intake metabolites
intake_cols <- grep("^int_", names(intake_mtb), value = TRUE)
metabolites_intake <- intake_mtb[,colnames(intake_mtb) %in% intake_cols]

# Calculate the percentage of non-zero values for each variable
non_zero_pct <- apply(metabolites_intake != 0, 2, mean)
# Filter variables with at least a non-zero value in 20% of the data
filter_20pct <- metabolites_intake[,non_zero_pct >= 0.2]

#add pseudocount to all variables
pseudo <- filter_20pct + 1

# Create a new column specifying IBD (yes/no) for each sample
pseudo_diagnosis <- cbind(intake_mtb$diagnosis, pseudo)
names(pseudo_diagnosis)[1] <- 'diagnosis'

##======================= WILCOXON TEST ======================##
wilcoxon_p <- c() # Initialize empty vector for p-values
# Do "for loop" over selected column names
for (i in 2:1010) {
  
  result <- wilcox.test(pseudo_diagnosis[, i] ~ diagnosis,
                        data = pseudo_diagnosis)
  
  # Stores p-value to the vector with this column name
  wilcoxon_p[[i]]  <- result$p.value
  
}

#store metabolites with raw p-value in new dataframe
wilcoxon_p <- data.frame(metabolites =  names(pseudo_diagnosis[,2:1010]),
                         p_raw = unlist(wilcoxon_p))
wilcoxon_p$p_adjusted <- p.adjust(wilcoxon_p$p_raw, method = "fdr") #add column with p_adj for multiple testing

# prepare a dataframe to plot p values
df <- data.frame(x = c(wilcoxon_p$p_raw, wilcoxon_p$p_adjusted), 
                 type=rep(c("raw", "fdr"),
                          c(length(wilcoxon_p$p_raw),
                            length(wilcoxon_p$p_adjusted))))

# make a histrogram of p values and adjusted p values
wilcoxon_plot <- ggplot(df) +
  geom_histogram(aes(x=x, fill=type)) +
  labs(x = "p-value", y = "Frequency") 

wilcoxon_plot

##============================== VOLCANO PLOT OF DIFFERENTIAL ABUNDANCE ANALYSIS ===============================##

#log2 transformation
intake_log <- log2(pseudo_diagnosis[,2:1010])
intake_log <- cbind(diagnosis = pseudo_diagnosis$diagnosis, intake_log)


#calculate the mean of each metabolite in IBD group
IBD <- filter(intake_log, intake_log$diagnosis == 'IBD')
IBD_m = apply(IBD[,2:1010], 2, mean)

#calcuate the mean of each metabolite in LLDEEP group
control <- filter(intake_log, intake_log$diagnosis == 'control')
control_m = apply(control[,2:1010], 2, mean)

#because the data is already log2 transformed, take the difference between the means.
foldchange <- control_m - IBD_m
hist(foldchange, xlab = "log2 Fold Change (IBD vs LLDEEP)")

#add foldchange to df containing p-values
wilcoxon_p$foldchange <- foldchange

# add a column of NAs
wilcoxon_p$intakedifference <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
wilcoxon_p$intakedifference[wilcoxon_p$foldchange > 0.6 & wilcoxon_p$p_raw < 0.05] <- "YES"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
wilcoxon_p$intakedifference[wilcoxon_p$foldchange < -0.6 & wilcoxon_p$p_raw < 0.05] <- "DOWN"


# Create a new column "label" to df, that will contain the name of metabolites where intake is different between groups (NA in case they are not)
wilcoxon_p$label <- NA
wilcoxon_p$label[wilcoxon_p$intakedifference != "NO"] <- wilcoxon_p$metabolites[wilcoxon_p$intakedifference != "NO"]
wilcoxon_p$label[wilcoxon_p$intakedifference == "NO"] <- wilcoxon_p$metabolites[wilcoxon_p$intakedifference == "NO"]

ggplot(wilcoxon_p, aes(x=foldchange, y=-1*log10(p_raw), col=intakedifference, label=label)) + 
  geom_point() + 
  theme_minimal() +
  theme(legend.position = 'bottom') +
  scale_color_manual(values=c("#999999", "#009E73")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") + 
  geom_text_repel(size = 2) +
  scale_x_continuous(name = 'foldchange (HC - IBD)')

##================================================= SIGNIFICANT INTAKE DIFFERENCE (IBD vs NON-IBD) BASED ON FOLDCHANGE ==============================================##

# Sorts foldchange in decreasing order. Takes 6 first ones. Takes those rows that match
# with foldchange. Takes metabolites. 
highest6 <- wilcoxon_p[wilcoxon_p$foldchange %in% sort(wilcoxon_p$foldchange, decreasing = TRUE)[1:6], ]$metabolites
# From intake table, takes only those metabolites that had highest foldchange
highest6_chem <- intake_log[,colnames(intake_log) %in% highest6]
# Adds colData that includes patient status information
highest6_full <- cbind(pseudo_diagnosis$diagnosis, highest6_chem)
names(highest6_full)[1] <- 'diagnosis'
highest6_full <- na.omit(highest6_full)

# Puts plots in the same picture
gridExtra::grid.arrange(
  
  # Plot 1
  ggplot(highest6_full, aes(x = diagnosis, y = highest6_full[,2], fill = diagnosis)) + 
    geom_boxplot() + 
    geom_signif(comparisons = list(c("IBD", "control")), 
                map_signif_level=TRUE) +
    ylab("Predicted metabolite intake") + # y axis title
    ggtitle(names(highest6_full)[2]) + # main title
    theme_minimal() + 
    theme(title = element_text(size = 7),
          legend.position = 'none',
          axis.text = element_text(size = 7),
          axis.title.x=element_blank()), # makes titles smaller, removes x axis title
  
  # Plot 2
  ggplot(highest6_full, aes(x = diagnosis, y = highest6_full[,3], fill = diagnosis)) + 
    geom_boxplot() + 
    geom_signif(comparisons = list(c("IBD", "control")), 
                map_signif_level=TRUE) +
    ylab("Predicted metabolite intake") + # y axis title
    ggtitle(names(highest6_full)[3]) + # main title
    theme_minimal() + 
    theme(title = element_text(size = 7),
          legend.position = 'none',
          axis.text = element_text(size = 7),
          axis.title.x=element_blank()), # makes titles smaller, removes x axis title
  
  # Plot 3
  ggplot(highest6_full, aes(x = diagnosis, y = highest6_full[,4], fill = diagnosis)) + 
    geom_boxplot() + 
    geom_signif(comparisons = list(c("IBD", "control")), 
                map_signif_level=TRUE) +
    ylab("Predicted metabolite intake") + # y axis title
    ggtitle(names(highest6_full)[4]) + # main title
    theme_minimal() + 
    theme(title = element_text(size = 7),
          legend.position = 'none',
          axis.text = element_text(size = 7),
          axis.title.x=element_blank()), # makes titles smaller, removes x axis title
  
  # Plot 4
  ggplot(highest6_full, aes(x = diagnosis, y = highest6_full[,5], fill = diagnosis)) + 
    geom_boxplot() + 
    geom_signif(comparisons = list(c("IBD", "control")), 
                map_signif_level=TRUE) +
    ylab("Predicted metabolite intake") + # y axis title
    ggtitle(names(highest6_full)[5]) + # main title
    theme_minimal() + 
    theme(title = element_text(size = 7),
          legend.position = 'none',
          axis.text = element_text(size = 7),
          axis.title.x=element_blank()), # makes titles smaller, removes x axis title
  
  # Plot 5
  ggplot(highest6_full, aes(x = diagnosis, y = highest6_full[,6], fill = diagnosis)) + 
    geom_boxplot() + 
    geom_signif(comparisons = list(c("IBD", "control")), 
                map_signif_level=TRUE) +
    ylab("Predicted metabolite intake") + # y axis title
    ggtitle(names(highest6_full)[6]) + # main title
    theme_minimal() + 
    theme(title = element_text(size = 7),
          legend.position = 'none',
          axis.text = element_text(size = 7),
          axis.title.x=element_blank()), # makes titles smaller, removes x axis title
  
  # Plot 6
  ggplot(highest6_full, aes(x = diagnosis, y = highest6_full[,7], fill = diagnosis)) + 
    geom_boxplot() + 
    geom_signif(comparisons = list(c("IBD", "control")), 
                map_signif_level=TRUE) +
    ylab("Predicted metabolite intake") + # y axis title
    ggtitle(names(highest6_full)[7]) + # main title
    theme_minimal() + 
    theme(title = element_text(size = 7),
          legend.position = 'none',
          axis.text = element_text(size = 7),
          axis.title.x=element_blank()), # makes titles smaller, removes x axis title
  
  # 3 columns and 2 rows
  ncol = 3, 
  nrow = 2
)

##===================================== DIFFERENTIAL ABUNDANCE ANALYSIS: INTAKE_calprotectin >150 vs <150 (LOG TRANSFORMED AND FILTERED DATA) ====================================================

# Columns containing intake metabolites
intake_cols <- grep("^int_", names(intake_mtb), value = TRUE)
metabolites_intake <- intake_mtb[,colnames(intake_mtb) %in% intake_cols]

# Calculate the percentage of non-zero values for each variable
non_zero_pct <- apply(metabolites_intake != 0, 2, mean)
# Filter variables with at least a non-zero value in 20% of the data
filter_20pct <- metabolites_intake[,non_zero_pct >= 0.2]

#add pseudocount to all variables
pseudo <- filter_20pct + 1

# Create a new column specifying calprotectin >150 (yes/no) for each sample
pseudo_calprotectin <- cbind(intake_mtb$calprotectin_above150, pseudo)
names(pseudo_calprotectin)[1] <- 'calprotectin_above150'

#pseudo_clean <- remove_outliers(pseudo_diagnosis, 'diagnosis')

##======================= WILCOXON TEST ======================##
wilcoxon_p <- c() # Initialize empty vector for p-values
# Do "for loop" over selected column names
for (i in 2:1010) {
  
  result <- wilcox.test(pseudo_calprotectin[, i] ~ calprotectin_above150,
                        data = pseudo_calprotectin)
  
  # Stores p-value to the vector with this column name
  wilcoxon_p[[i]]  <- result$p.value
  
}

#store metabolites with raw p-value in new dataframe
wilcoxon_p <- data.frame(metabolites =  names(pseudo_calprotectin[,2:1010]),
                         p_raw = unlist(wilcoxon_p))
wilcoxon_p$p_adjusted <- p.adjust(wilcoxon_p$p_raw, method = "fdr") #add column with p_adj for multiple testing

# prepare a dataframe to plot p values
df <- data.frame(x = c(wilcoxon_p$p_raw, wilcoxon_p$p_adjusted), 
                 type=rep(c("raw", "fdr"),
                          c(length(wilcoxon_p$p_raw),
                            length(wilcoxon_p$p_adjusted))))

# make a histrogram of p values and adjusted p values
wilcoxon_plot <- ggplot(df) +
  geom_histogram(aes(x=x, fill=type)) +
  labs(x = "p-value", y = "Frequency") 

wilcoxon_plot


##============================== VOLCANO PLOT OF DIFFERENTIAL ABUNDANCE ANALYSIS ===============================##

#log2 transformation
intake_log <- log2(pseudo_calprotectin[,2:1010])
intake_log <- cbind(calprotectin_above150 = pseudo_calprotectin$calprotectin_above150, intake_log)


#calculate the mean of each metabolite in >150 group
high_calprotectin <- filter(intake_log, intake_log$calprotectin_above150 == 'yes')
high_calprotectin_m = apply(high_calprotectin[,2:1010], 2, mean)

#calcuate the mean of each metabolite in <150 group
low_calprotectin <- filter(intake_log, intake_log$calprotectin_above150 == 'no')
low_calprotectin_m = apply(low_calprotectin[,2:1010], 2, mean)

#because the data is already log2 transformed, take the difference between the means.
foldchange <- low_calprotectin_m - high_calprotectin_m
hist(foldchange, xlab = "log2 Fold Change (<150 vs >150)")

#add foldchange to df containing p-values
wilcoxon_p$foldchange <- foldchange

# add a column of NAs
wilcoxon_p$intakedifference <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
wilcoxon_p$intakedifference[wilcoxon_p$foldchange > 0.6 & wilcoxon_p$p_raw < 0.05] <- "YES"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
wilcoxon_p$intakedifference[wilcoxon_p$foldchange < -0.6 & wilcoxon_p$p_raw < 0.05] <- "DOWN"


# Create a new column "label" to df, that will contain the name of metabolites where intake is different between groups (NA in case they are not)
wilcoxon_p$label <- NA
wilcoxon_p$label[wilcoxon_p$intakedifference != "NO"] <- wilcoxon_p$metabolites[wilcoxon_p$intakedifference != "NO"]
wilcoxon_p$label[wilcoxon_p$intakedifference == "NO"] <- wilcoxon_p$metabolites[wilcoxon_p$intakedifference == "NO"]

ggplot(wilcoxon_p, aes(x=foldchange, y=-1*log10(p_raw), col=intakedifference, label=label)) + 
  geom_point() + 
  theme_minimal() +
  theme(legend.position = 'bottom') +
  scale_color_manual(values=c("#999999", "#009E73")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") + 
  geom_text_repel(size = 2) +
  scale_x_continuous(name = 'foldchange (<150 - >150)')


##================================================= SIGNIFICANT INTAKE DIFFERENCE (IBD vs NON-IBD) BASED ON FOLDCHANGE ==============================================##

# Sorts foldchange in decreasing order. Takes 6 first ones. Takes those rows that match
# with foldchange. Takes metabolites. 
highest6 <- wilcoxon_p[wilcoxon_p$foldchange %in% sort(wilcoxon_p$foldchange, decreasing = TRUE)[1:6], ]$metabolites
# From intake table, takes only those metabolites that had highest foldchange
highest6_chem <- intake_log[,colnames(intake_log) %in% highest6]
# Adds colData that includes patient status information
highest6_full <- cbind(pseudo_calprotectin$calprotectin_above150, highest6_chem)
names(highest6_full)[1] <- 'calprotectin_above150'
highest6_full <- na.omit(highest6_full)

# Puts plots in the same picture
gridExtra::grid.arrange(
  
  # Plot 1
  ggplot(highest6_full, aes(x = calprotectin_above150, y = highest6_full[,2], fill = calprotectin_above150)) + 
    geom_boxplot() + 
    geom_signif(comparisons = list(c("yes", "no")), 
                map_signif_level=TRUE) +
    ylab("Predicted metabolite intake") + # y axis title
    ggtitle(names(highest6_full)[2]) + # main title
    theme_minimal() + 
    theme(title = element_text(size = 7),
          legend.position = 'none',
          axis.text = element_text(size = 7),
          axis.title.x=element_blank()), # makes titles smaller, removes x axis title
  
  # Plot 2
  ggplot(highest6_full, aes(x = calprotectin_above150, y = highest6_full[,3], fill = calprotectin_above150)) + 
    geom_boxplot() + 
    geom_signif(comparisons = list(c("yes", "no")), 
                map_signif_level=TRUE) +
    ylab("Predicted metabolite intake") + # y axis title
    ggtitle(names(highest6_full)[3]) + # main title
    theme_minimal() + 
    theme(title = element_text(size = 7),
          legend.position = 'none',
          axis.text = element_text(size = 7),
          axis.title.x=element_blank()), # makes titles smaller, removes x axis title
  
  # Plot 3
  ggplot(highest6_full, aes(x = calprotectin_above150, y = highest6_full[,4], fill = calprotectin_above150)) + 
    geom_boxplot() + 
    geom_signif(comparisons = list(c("yes", "no")), 
                map_signif_level=TRUE) +
    ylab("Predicted metabolite intake") + # y axis title
    ggtitle(names(highest6_full)[4]) + # main title
    theme_minimal() + 
    theme(title = element_text(size = 7),
          legend.position = 'none',
          axis.text = element_text(size = 7),
          axis.title.x=element_blank()), # makes titles smaller, removes x axis title
  
  # Plot 4
  ggplot(highest6_full, aes(x = calprotectin_above150, y = highest6_full[,5], fill = calprotectin_above150)) + 
    geom_boxplot() + 
    geom_signif(comparisons = list(c("yes", "no")), 
                map_signif_level=TRUE) +
    ylab("Predicted metabolite intake") + # y axis title
    ggtitle(names(highest6_full)[5]) + # main title
    theme_minimal() + 
    theme(title = element_text(size = 7),
          legend.position = 'none',
          axis.text = element_text(size = 7),
          axis.title.x=element_blank()), # makes titles smaller, removes x axis title
  
  # Plot 5
  ggplot(highest6_full, aes(x = calprotectin_above150, y = highest6_full[,6], fill = calprotectin_above150)) + 
    geom_boxplot() + 
    geom_signif(comparisons = list(c("yes", "no")), 
                map_signif_level=TRUE) +
    ylab("Predicted metabolite intake") + # y axis title
    ggtitle(names(highest6_full)[6]) + # main title
    theme_minimal() + 
    theme(title = element_text(size = 7),
          legend.position = 'none',
          axis.text = element_text(size = 7),
          axis.title.x=element_blank()), # makes titles smaller, removes x axis title
  
  # Plot 6
  ggplot(highest6_full, aes(x = calprotectin_above150, y = highest6_full[,7], fill = calprotectin_above150)) + 
    geom_boxplot() + 
    geom_signif(comparisons = list(c("yes", "no")), 
                map_signif_level=TRUE) +
    ylab("Predicted metabolite intake") + # y axis title
    ggtitle(names(highest6_full)[7]) + # main title
    theme_minimal() + 
    theme(title = element_text(size = 7),
          legend.position = 'none',
          axis.text = element_text(size = 7),
          axis.title.x=element_blank()), # makes titles smaller, removes x axis title
  
  # 3 columns and 2 rows
  ncol = 3, 
  nrow = 2
)

##================================================ INTAKE METABOLITES BASED ON INTAKEDIFFERENCE FOUND IN DIAGNOSIS AND calprotectin ====================================================##

#Metabolite names based on intakedifference == YES (only works when NOT running calprotectin part of script)
intakedifference_diagnosis <- as.character(wilcoxon_p$metabolites[wilcoxon_p$intakedifference != "NO"])

#Metabolite names based on intakedifference == YES
intakedifference_calprotectin <- as.character(wilcoxon_p$metabolites[wilcoxon_p$intakedifference != "NO"])

#New df containing only intake metabolites that show difference in intake in 
      ## IBD vs Non-IBD
      ## calprotectin <150 vs calprotectin >150
intake_mtb_1 <- intake_mtb[,colnames(intake_mtb) %in% intakedifference_diagnosis]
intake_mtb_2 <- intake_mtb_1[,colnames(intake_mtb_1) %in% intakedifference_calprotectin]

##========================================== LINEAR REGRESSION: IBD vs NON-IBD ========================================

# Only columns containing metabolite data from original intake df
metabolites_intake <- intake[,7:1114]
colnames(metabolites_intake) <- make.unique(colnames(metabolites_intake), sep = "_") #make sure that column names are all unique

# Calculate the percentage of non-zero values for each variable
non_zero_pct <- apply(metabolites_intake != 0, 2, mean)
# Filter variables with at least a non-zero value in 20% of the data
filter_20pct <- metabolites_intake[,non_zero_pct >= 0.2]

#add pseudocount to all variables
pseudo <- filter_20pct + 1
sum(is.na(pseudo)) #0 NA
# Apply the remove_outliers function to set outliers as NA within each column
pseudo_clean <- remove_outliers(pseudo)
sum(is.na(pseudo_clean)) #106632 NA

#RANK transformation
pseudo_rank <- pseudo_clean %>% mutate_all(~ (length(.) + 1) - rank(.))

## Filtering of analysis table on metabolite type is performed in == LOAD DATA, CLEANING NAMES AND SUBSETTING == #

# Create a new column specifying IBD (yes/no) for each sample
pseudo_diagnosis <- cbind(diagnosis = intake_mtb$diagnosis, pseudo_rank)

# Add covariates to df
pseudo_diagnosis <- cbind(age = intake_mtb$sex, pseudo_diagnosis)
pseudo_diagnosis <- cbind(sex = intake_mtb$age, pseudo_diagnosis)
pseudo_diagnosis <- cbind(BMI = intake_mtb$BMI, pseudo_diagnosis)

#============LINEAR REGRESSION ==============#
#dependent variable: individual diet metabolites
#predictor variables: covariates(age, sex, BMI), IBD(yes/no)

#columns with intake data + predictor variables
metabolite_names <- names(pseudo_diagnosis[,5:996])
predictor_vars <- c('age', 'sex', 'BMI', 'diagnosis')

# Create an empty dataframe to store results
results_df <- data.frame(Intake_Metabolite = character(),
                         Coefficient = numeric(),
                         Estimate = numeric(),
                         PValue = numeric(),
                         RSquared = numeric(),
                         stringsAsFactors = FALSE)

# Perform linear regression over intake_metabolites using a for loop
# Perform linear regression for each dependent variable
for (dep_var in metabolite_names) {
  # Create formula
  formula <- paste(dep_var, paste(predictor_vars, collapse = " + "), sep = " ~ ")
  
  # Perform linear regression
  regression_model <- lm(formula, data = pseudo_diagnosis)
  
  # Extract coefficient estimates, p-values, and R-squared value
  coefficients <- coef(regression_model)
  p_values <- summary(regression_model)$coefficients[, "Pr(>|t|)"]
  r_squared <- summary(regression_model)$r.squared
  
  # Filter significant coefficients (p-value < 0.05)
  significant_coeffs <- coefficients[p_values < 0.05]
  significant_pvalues <- p_values[p_values < 0.05]
  
  # Create a dataframe for the results
  if (length(significant_coeffs) > 0) {
    results <- data.frame(Intake_Metabolite = dep_var,
                          Coefficient = names(significant_coeffs),
                          Estimate = significant_coeffs,
                          PValue = significant_pvalues,
                          RSquared = r_squared,
                          stringsAsFactors = FALSE)
    
    # Append results to the main dataframe
    results_df <- rbind(results_df, results)
  }
}

# Filter results df on significant coefficients 'diagnosis'
linreg_diagnosis <- results_df[results_df$Coefficient == 'diagnosisIBD',]
linreg_diagnosis$p_adjusted <- p.adjust(linreg_diagnosis$PValue, method = "fdr") #add column with p_adj for multiple testing

##============================== VOLCANO PLOT OF LINREG_DIAGNOSIS ===============================##

ggplot(linreg_diagnosis, aes(x=Estimate, y=-1*log10(PValue), label=Intake_Metabolite)) + 
  geom_point() + 
  theme_minimal() +
  theme(legend.position = 'bottom') +
  scale_color_manual(values=c("#999999", "#009E73")) +
  ylim(0, NA) +
  geom_text_repel(size = 2) +
  scale_x_continuous(name = 'Estimate')

##================================================= SIGNIFICANT INTAKE DIFFERENCE (IBD vs NON-IBD) BASED ON FOLDCHANGE ==============================================##

# Sorts pvalue in increasing order. Takes 6 first ones. Takes those rows that match
# with foldchange. Takes metabolites. 
highest6 <- linreg_diagnosis[linreg_diagnosis$RSquared %in% sort(linreg_diagnosis$RSquared, decreasing = T)[1:6], ]$Intake_Metabolite
# From intake table, takes only those metabolites that had highest foldchange
highest6_chem <- pseudo_diagnosis[,colnames(pseudo_diagnosis) %in% highest6]
# Adds colData that includes patient status information
highest6_full <- cbind(pseudo_diagnosis$diagnosis, highest6_chem)
names(highest6_full)[1] <- 'diagnosis'
highest6_full <- na.omit(highest6_full)

# Puts plots in the same picture
gridExtra::grid.arrange(
  
  # Plot 1
  ggplot(highest6_full, aes(x = diagnosis, y = highest6_full[,2], fill = diagnosis)) + 
    geom_boxplot() + 
    geom_signif(comparisons = list(c("IBD", "control")), 
                map_signif_level=TRUE) +
    ylab("Predicted metabolite intake") + # y axis title
    ggtitle(names(highest6_full)[2]) + # main title
    theme_minimal() + 
    theme(title = element_text(size = 7),
          legend.position = 'none',
          axis.text = element_text(size = 7),
          axis.title.x=element_blank()), # makes titles smaller, removes x axis title
  
  # Plot 2
  ggplot(highest6_full, aes(x = diagnosis, y = highest6_full[,3], fill = diagnosis)) + 
    geom_boxplot() + 
    geom_signif(comparisons = list(c("IBD", "control")), 
                map_signif_level=TRUE) +
    ylab("Predicted metabolite intake") + # y axis title
    ggtitle(names(highest6_full)[3]) + # main title
    theme_minimal() + 
    theme(title = element_text(size = 7),
          legend.position = 'none',
          axis.text = element_text(size = 7),
          axis.title.x=element_blank()), # makes titles smaller, removes x axis title
  
  # Plot 3
  ggplot(highest6_full, aes(x = diagnosis, y = highest6_full[,4], fill = diagnosis)) + 
    geom_boxplot() + 
    geom_signif(comparisons = list(c("IBD", "control")), 
                map_signif_level=TRUE) +
    ylab("Predicted metabolite intake") + # y axis title
    ggtitle(names(highest6_full)[4]) + # main title
    theme_minimal() + 
    theme(title = element_text(size = 7),
          legend.position = 'none',
          axis.text = element_text(size = 7),
          axis.title.x=element_blank()), # makes titles smaller, removes x axis title
  
  # Plot 4
  ggplot(highest6_full, aes(x = diagnosis, y = highest6_full[,5], fill = diagnosis)) + 
    geom_boxplot() + 
    geom_signif(comparisons = list(c("IBD", "control")), 
                map_signif_level=TRUE) +
    ylab("Predicted metabolite intake") + # y axis title
    ggtitle(names(highest6_full)[5]) + # main title
    theme_minimal() + 
    theme(title = element_text(size = 7),
          legend.position = 'none',
          axis.text = element_text(size = 7),
          axis.title.x=element_blank()), # makes titles smaller, removes x axis title
  
  # Plot 5
  ggplot(highest6_full, aes(x = diagnosis, y = highest6_full[,6], fill = diagnosis)) + 
    geom_boxplot() + 
    geom_signif(comparisons = list(c("IBD", "control")), 
                map_signif_level=TRUE) +
    ylab("Predicted metabolite intake") + # y axis title
    ggtitle(names(highest6_full)[6]) + # main title
    theme_minimal() + 
    theme(title = element_text(size = 7),
          legend.position = 'none',
          axis.text = element_text(size = 7),
          axis.title.x=element_blank()), # makes titles smaller, removes x axis title
  
  # Plot 6
  ggplot(highest6_full, aes(x = diagnosis, y = highest6_full[,7], fill = diagnosis)) + 
    geom_boxplot() + 
    geom_signif(comparisons = list(c("IBD", "control")), 
                map_signif_level=TRUE) +
    ylab("Predicted metabolite intake") + # y axis title
    ggtitle(names(highest6_full)[7]) + # main title
    theme_minimal() + 
    theme(title = element_text(size = 7),
          legend.position = 'none',
          axis.text = element_text(size = 7),
          axis.title.x=element_blank()), # makes titles smaller, removes x axis title
  
  # 3 columns and 2 rows
  ncol = 3, 
  nrow = 2
)


##========================================== LINEAR REGRESSION: >150 CALPROTECTIN vs <150 CALPROTECTIN ========================================

# Only columns containing metabolite data from original intake df
metabolites_intake <- intake[,7:1114]
colnames(metabolites_intake) <- make.unique(colnames(metabolites_intake), sep = "_") #make sure that column names are all unique

# Calculate the percentage of non-zero values for each variable
non_zero_pct <- apply(metabolites_intake != 0, 2, mean)
# Filter variables with at least a non-zero value in 20% of the data
filter_20pct <- metabolites_intake[,non_zero_pct >= 0.2]

#add pseudocount to all variables
pseudo <- filter_20pct + 1
sum(is.na(pseudo)) #0 NA
# Apply the remove_outliers function to set outliers as NA within each column
pseudo_clean <- remove_outliers(pseudo)
sum(is.na(pseudo_clean)) #106632 NA

#RANK transformation
pseudo_rank <- pseudo_clean %>% mutate_all(~ (length(.) + 1) - rank(.))

## Filtering of analysis table on metabolite (intake, fecal, serum) type is performed in == LOAD DATA, CLEANING NAMES AND SUBSETTING == #

# Create a new column specifying high calprotectin (yes/no) for each sample
pseudo_calprotectin <- cbind(calprotectin_above150 = intake_mtb$calprotectin_above150, pseudo_rank)

# Add covariates to df
pseudo_calprotectin <- cbind(age = intake_mtb$age, pseudo_calprotectin)
pseudo_calprotectin <- cbind(sex = intake_mtb$sex, pseudo_calprotectin)
pseudo_calprotectin <- cbind(BMI = intake_mtb$BMI, pseudo_calprotectin)

#============LINEAR REGRESSION ==============#
#dependent variable: individual diet metabolites
#predictor variables: covariates(age, sex, BMI), calprotectin >150 (yes/no)

#columns with intake data + predictor variables
metabolite_names <- names(pseudo_calprotectin[,5:1013])
predictor_vars <- c('age', 'sex', 'BMI', 'calprotectin_above150')

# Create an empty dataframe to store results
results_df <- data.frame(Intake_Metabolite = character(),
                         Coefficient = numeric(),
                         Estimate = numeric(),
                         PValue = numeric(),
                         RSquared = numeric(),
                         stringsAsFactors = FALSE)

# Perform linear regression over intake_metabolites using a for loop
# Perform linear regression for each dependent variable
for (dep_var in metabolite_names) {
  # Create formula
  formula <- paste(dep_var, paste(predictor_vars, collapse = " + "), sep = " ~ ")
  
  # Perform linear regression
  regression_model <- lm(formula, data = pseudo_calprotectin)
  
  # Extract coefficient estimates, p-values, and R-squared value
  coefficients <- coef(regression_model)
  p_values <- summary(regression_model)$coefficients[, "Pr(>|t|)"]
  r_squared <- summary(regression_model)$r.squared
  
  # Filter significant coefficients (p-value < 0.05)
  significant_coeffs <- coefficients[p_values < 0.05]
  significant_pvalues <- p_values[p_values < 0.05]
  
  # Create a dataframe for the results
  if (length(significant_coeffs) > 0) {
    results <- data.frame(Intake_Metabolite = dep_var,
                          Coefficient = names(significant_coeffs),
                          Estimate = significant_coeffs,
                          PValue = significant_pvalues,
                          RSquared = r_squared,
                          stringsAsFactors = FALSE)
    
    # Append results to the main dataframe
    results_df <- rbind(results_df, results)
  }
}

# Filter results df on significant coefficients 'calprotectin_above150'
linreg_calprotectin <- results_df[results_df$Coefficient == 'calprotectin_above150yes',]
linreg_calprotectin$p_adjusted <- p.adjust(linreg_calprotectin$PValue, method = "fdr") #add column with p_adj for multiple testing

##============================== VOLCANO PLOT OF LINREG_CALPROTECTIN ===============================##

ggplot(linreg_calprotectin, aes(x=Estimate, y=-1*log10(PValue), label=Intake_Metabolite)) + 
  geom_point() + 
  theme_minimal() +
  theme(legend.position = 'bottom') +
  scale_color_manual(values=c("#999999", "#009E73")) +
  ylim(0, NA) +
  geom_text_repel(size = 2) +
  scale_x_continuous(name = 'Estimate')

##================================================= SIGNIFICANT INTAKE DIFFERENCE (CALPROTECTIN) ==============================================##

# Sorts pvalue in increasing order. Takes 6 first ones. Takes those rows that match
# with foldchange. Takes metabolites. 
highest6 <- linreg_calprotectin[linreg_calprotectin$RSquared %in% sort(linreg_calprotectin$RSquared, decreasing = T)[1:6], ]$Intake_Metabolite
# From intake table, takes only those metabolites that had highest foldchange
highest6_chem <- pseudo_calprotectin[,colnames(pseudo_calprotectin) %in% highest6]
# Adds colData that includes patient status information
highest6_full <- cbind(pseudo_calprotectin$calprotectin_above150, highest6_chem)
names(highest6_full)[1] <- 'calprotectin_above150'
highest6_full <- na.omit(highest6_full)

# Puts plots in the same picture
gridExtra::grid.arrange(
  
  # Plot 1
  ggplot(highest6_full, aes(x = calprotectin_above150, y = highest6_full[,2], fill = calprotectin_above150)) + 
    geom_boxplot() + 
    geom_signif(comparisons = list(c("no", "yes")), 
                map_signif_level=TRUE) +
    ylab("Predicted metabolite intake") + # y axis title
    ggtitle(names(highest6_full)[2]) + # main title
    theme_minimal() + 
    theme(title = element_text(size = 7),
          legend.position = 'none',
          axis.text = element_text(size = 7),
          axis.title.x=element_blank()), # makes titles smaller, removes x axis title
  
  # Plot 2
  ggplot(highest6_full, aes(x = calprotectin_above150, y = highest6_full[,3], fill = calprotectin_above150)) + 
    geom_boxplot() + 
    geom_signif(comparisons = list(c("no", "yes")), 
                map_signif_level=TRUE) +
    ylab("Predicted metabolite intake") + # y axis title
    ggtitle(names(highest6_full)[3]) + # main title
    theme_minimal() + 
    theme(title = element_text(size = 7),
          legend.position = 'none',
          axis.text = element_text(size = 7),
          axis.title.x=element_blank()), # makes titles smaller, removes x axis title
  
  # Plot 3
  ggplot(highest6_full, aes(x = calprotectin_above150, y = highest6_full[,4], fill = calprotectin_above150)) + 
    geom_boxplot() + 
    geom_signif(comparisons = list(c("no", "yes")), 
                map_signif_level=TRUE) +
    ylab("Predicted metabolite intake") + # y axis title
    ggtitle(names(highest6_full)[4]) + # main title
    theme_minimal() + 
    theme(title = element_text(size = 7),
          legend.position = 'none',
          axis.text = element_text(size = 7),
          axis.title.x=element_blank()), # makes titles smaller, removes x axis title
  
  # Plot 4
  ggplot(highest6_full, aes(x = calprotectin_above150, y = highest6_full[,5], fill = calprotectin_above150)) + 
    geom_boxplot() + 
    geom_signif(comparisons = list(c("no", "yes")), 
                map_signif_level=TRUE) +
    ylab("Predicted metabolite intake") + # y axis title
    ggtitle(names(highest6_full)[5]) + # main title
    theme_minimal() + 
    theme(title = element_text(size = 7),
          legend.position = 'none',
          axis.text = element_text(size = 7),
          axis.title.x=element_blank()), # makes titles smaller, removes x axis title
  
  # Plot 5
  ggplot(highest6_full, aes(x = calprotectin_above150, y = highest6_full[,6], fill = calprotectin_above150)) + 
    geom_boxplot() + 
    geom_signif(comparisons = list(c("no", "yes")), 
                map_signif_level=TRUE) +
    ylab("Predicted metabolite intake") + # y axis title
    ggtitle(names(highest6_full)[6]) + # main title
    theme_minimal() + 
    theme(title = element_text(size = 7),
          legend.position = 'none',
          axis.text = element_text(size = 7),
          axis.title.x=element_blank()), # makes titles smaller, removes x axis title
  
  # Plot 6
  ggplot(highest6_full, aes(x = calprotectin_above150, y = highest6_full[,7], fill = calprotectin_above150)) + 
    geom_boxplot() + 
    geom_signif(comparisons = list(c("no", "yes")), 
                map_signif_level=TRUE) +
    ylab("Predicted metabolite intake") + # y axis title
    ggtitle(names(highest6_full)[7]) + # main title
    theme_minimal() + 
    theme(title = element_text(size = 7),
          legend.position = 'none',
          axis.text = element_text(size = 7),
          axis.title.x=element_blank()), # makes titles smaller, removes x axis title
  
  # 3 columns and 2 rows
  ncol = 3, 
  nrow = 2
)

##========================================== LINEAR REGRESSION: BEFORE A FLARE vs DURING/AFTER A FLARE ========================================

# Only columns containing metabolite data from original intake df
metabolites_intake <- intake[,7:1114]
colnames(metabolites_intake) <- make.unique(colnames(metabolites_intake), sep = "_") #make sure that column names are all unique

# Calculate the percentage of non-zero values for each variable
non_zero_pct <- apply(metabolites_intake != 0, 2, mean)
# Filter variables with at least a non-zero value in 20% of the data
filter_20pct <- metabolites_intake[,non_zero_pct >= 0.2]

#add pseudocount to all variables
pseudo <- filter_20pct + 1
sum(is.na(pseudo)) #0 NA
# Apply the remove_outliers function to set outliers as NA within each column
pseudo_clean <- remove_outliers(pseudo)
sum(is.na(pseudo_clean)) #106632 NA

#LOG transformation or RANK transformation
pseudo_rank <- pseudo_clean %>% mutate_all(~ (length(.) + 1) - rank(.))

## Filtering of analysis table on metabolite type (intake, fecal, serum) is performed in == LOAD DATA, CLEANING NAMES AND SUBSETTING == #

# Create a new column specifying IBD (yes/no) for each sample
pseudo_flare <- cbind(before_a_flare = intake_mtb$before_a_flare, pseudo_rank)

# Add covariates to df
pseudo_flare <- cbind(age = intake_mtb$age, pseudo_flare)
pseudo_flare <- cbind(sex = intake_mtb$sex, pseudo_flare)
pseudo_flare <- cbind(BMI = intake_mtb$BMI, pseudo_flare)

#============LINEAR REGRESSION ==============#
#dependent variable: individual diet metabolites
#predictor variables: covariates(age, sex, BMI), before a flare (yes/no)

#columns with intake data + predictor variables
metabolite_names <- names(pseudo_flare[,5:1013])
predictor_vars <- c('age', 'sex', 'BMI', 'before_a_flare')

# Create an empty dataframe to store results
results_df <- data.frame(Intake_Metabolite = character(),
                         Coefficient = numeric(),
                         Estimate = numeric(),
                         PValue = numeric(),
                         RSquared = numeric(),
                         stringsAsFactors = FALSE)

# Perform linear regression over intake_metabolites using a for loop
# Perform linear regression for each dependent variable
for (dep_var in metabolite_names) {
  # Create formula
  formula <- paste(dep_var, paste(predictor_vars, collapse = " + "), sep = " ~ ")
  
  # Perform linear regression
  regression_model <- lm(formula, data = pseudo_flare)
  
  # Extract coefficient estimates, p-values, and R-squared value
  coefficients <- coef(regression_model)
  p_values <- summary(regression_model)$coefficients[, "Pr(>|t|)"]
  r_squared <- summary(regression_model)$r.squared
  
  # Filter significant coefficients (p-value < 0.05)
  significant_coeffs <- coefficients[p_values < 0.05]
  significant_pvalues <- p_values[p_values < 0.05]
  
  # Create a dataframe for the results
  if (length(significant_coeffs) > 0) {
    results <- data.frame(Intake_Metabolite = dep_var,
                          Coefficient = names(significant_coeffs),
                          Estimate = significant_coeffs,
                          PValue = significant_pvalues,
                          RSquared = r_squared,
                          stringsAsFactors = FALSE)
    
    # Append results to the main dataframe
    results_df <- rbind(results_df, results)
  }
}

# Filter results df on significant coefficients 'flare'
linreg_flare <- results_df[results_df$Coefficient == 'before_a_flareyes',]
linreg_flare$p_adjusted <- p.adjust(linreg_flare$PValue, method = "fdr") #add column with p_adj for multiple testing

##============================== VOLCANO PLOT OF LINREG_FLARE ===============================##

ggplot(linreg_flare, aes(x=Estimate, y=-1*log10(PValue), label=Intake_Metabolite)) + 
  geom_point() + 
  theme_minimal() +
  theme(legend.position = 'bottom') +
  scale_color_manual(values=c("#999999", "#009E73")) +
  ylim(0, NA) +
  geom_text_repel(size = 2) +
  scale_x_continuous(name = 'Estimate')

##================================================= SIGNIFICANT INTAKE DIFFERENCE (calprotectin) ==============================================##

# Sorts pvalue in increasing order. Takes 6 first ones. Takes those rows that match
# with foldchange. Takes metabolites. 
highest6 <- linreg_flare[linreg_flare$RSquared %in% sort(linreg_flare$RSquared, decreasing = T)[1:6], ]$Intake_Metabolite
# From intake table, takes only those metabolites that had highest foldchange
highest6_chem <- pseudo_flare[,colnames(pseudo_flare) %in% highest6]
# Adds colData that includes patient status information
highest6_full <- cbind(pseudo_flare$before_a_flare, highest6_chem)
names(highest6_full)[1] <- 'before_a_flare'
highest6_full <- na.omit(highest6_full)

# Puts plots in the same picture
gridExtra::grid.arrange(
  
  # Plot 1
  ggplot(highest6_full, aes(x = before_a_flare, y = highest6_full[,2], fill = before_a_flare)) + 
    geom_boxplot() + 
    geom_signif(comparisons = list(c("no", "yes")), 
                map_signif_level=TRUE) +
    ylab("Predicted metabolite intake") + # y axis title
    ggtitle(names(highest6_full)[2]) + # main title
    theme_minimal() + 
    theme(title = element_text(size = 7),
          legend.position = 'none',
          axis.text = element_text(size = 7),
          axis.title.x=element_blank()), # makes titles smaller, removes x axis title
  
  # Plot 2
  ggplot(highest6_full, aes(x = before_a_flare, y = highest6_full[,3], fill = before_a_flare)) + 
    geom_boxplot() + 
    geom_signif(comparisons = list(c("no", "yes")), 
                map_signif_level=TRUE) +
    ylab("Predicted metabolite intake") + # y axis title
    ggtitle(names(highest6_full)[3]) + # main title
    theme_minimal() + 
    theme(title = element_text(size = 7),
          legend.position = 'none',
          axis.text = element_text(size = 7),
          axis.title.x=element_blank()), # makes titles smaller, removes x axis title
  
  # Plot 3
  ggplot(highest6_full, aes(x = before_a_flare, y = highest6_full[,4], fill = before_a_flare)) + 
    geom_boxplot() + 
    geom_signif(comparisons = list(c("no", "yes")), 
                map_signif_level=TRUE) +
    ylab("Predicted metabolite intake") + # y axis title
    ggtitle(names(highest6_full)[4]) + # main title
    theme_minimal() + 
    theme(title = element_text(size = 7),
          legend.position = 'none',
          axis.text = element_text(size = 7),
          axis.title.x=element_blank()), # makes titles smaller, removes x axis title
  
  # Plot 4
  ggplot(highest6_full, aes(x = before_a_flare, y = highest6_full[,5], fill = before_a_flare)) + 
    geom_boxplot() + 
    geom_signif(comparisons = list(c("no", "yes")), 
                map_signif_level=TRUE) +
    ylab("Predicted metabolite intake") + # y axis title
    ggtitle(names(highest6_full)[5]) + # main title
    theme_minimal() + 
    theme(title = element_text(size = 7),
          legend.position = 'none',
          axis.text = element_text(size = 7),
          axis.title.x=element_blank()), # makes titles smaller, removes x axis title
  
  # Plot 5
  ggplot(highest6_full, aes(x = before_a_flare, y = highest6_full[,6], fill = before_a_flare)) + 
    geom_boxplot() + 
    geom_signif(comparisons = list(c("no", "yes")), 
                map_signif_level=TRUE) +
    ylab("Predicted metabolite intake") + # y axis title
    ggtitle(names(highest6_full)[6]) + # main title
    theme_minimal() + 
    theme(title = element_text(size = 7),
          legend.position = 'none',
          axis.text = element_text(size = 7),
          axis.title.x=element_blank()), # makes titles smaller, removes x axis title
  
  # Plot 6
  ggplot(highest6_full, aes(x = before_a_flare, y = highest6_full[,7], fill = before_a_flare)) + 
    geom_boxplot() + 
    geom_signif(comparisons = list(c("no", "yes")), 
                map_signif_level=TRUE) +
    ylab("Predicted metabolite intake") + # y axis title
    ggtitle(names(highest6_full)[7]) + # main title
    theme_minimal() + 
    theme(title = element_text(size = 7),
          legend.position = 'none',
          axis.text = element_text(size = 7),
          axis.title.x=element_blank()), # makes titles smaller, removes x axis title
  
  # 3 columns and 2 rows
  ncol = 3, 
  nrow = 2
)

##================================================ ASSOCIATED DIET METABOLITES WITH CLINICAL OUTCOMES BASED ON LINEAR REGRESSION + VENN DIAGRAM ====================================================

#Associated Diet Metabolite names based on DIAGNOSIS
diagnosis <- as.character(linreg_diagnosis$Intake_Metabolite)

#Associated Diet Metabolite names based on CALPROTECTIN
calprotectin <- as.character(linreg_calprotectin$Intake_Metabolite)

#Associated Diet Metabolite names based on FLARE
flare <- as.character(linreg_flare$Intake_Metabolite)

#New df containing only diet metabolites that are associated with clinical outcomes
## IBD vs Non-IBD
## calprotectin <150 vs calprotectin >150
## After/During calprotectin vs Before calprotectin
intake_mtb_diagnosis <- pseudo_diagnosis[,colnames(pseudo_diagnosis) %in% diagnosis]
intake_mtb_calprotectin <- pseudo_calprotectin[,colnames(pseudo_calprotectin) %in% calprotectin]
intake_mtb_flare <- pseudo_flare[,colnames(pseudo_flare) %in% flare]

# Get the column names from each data frame
col_names_diagnosis <- colnames(intake_mtb_diagnosis)
col_names_calprotectin <- colnames(intake_mtb_calprotectin)
col_names_flare <- colnames(intake_mtb_flare)

# Find the common column names
common_columns <- intersect(intersect(col_names_diagnosis, col_names_calprotectin), col_names_flare) # 2 intake mtb: 
diagnosis_calprotectin <- intersect(col_names_diagnosis, col_names_calprotectin)
diagnosis_flare <- intersect(col_names_diagnosis, col_names_flare)
calprotectin_flare <- intersect(col_names_calprotectin, col_names_flare)

# Convert the metabolite lists to sets
set1 <- unique(diagnosis)
set2 <- unique(calprotectin)
set3 <- unique(flare)

# Create the Venn diagram
venn.plot <- venn.diagram(
  x = list(set1, set2, set3),
  category.names = c("Diagnosis", "Calprotectin", "Flare"),
  fill = c("steelblue", "darkorange1", "forestgreen"),
  alpha = 0.5,
  filename = NULL  
)

# Output the Venn diagram
grid.newpage()
grid.draw(venn.plot)

##=========================== METABOLITES PRESENT FOR INTAKE, FECAL, SERUM ========================

#add pseudocount to all intake variables
pseudo_intake <- intake[,7:1114] + 1
colnames(pseudo_intake) <- make.unique(colnames(pseudo_intake), sep = "_") #make sure that column names are all unique

#add pseudocount to all fecal variables
pseudo_fecal <- fecal[,2:1685] + 1
colnames(pseudo_fecal) <- make.unique(colnames(pseudo_fecal), sep = "_") #make sure that column names are all unique

#add pseudocount to all serum variables
pseudo_serum <- serum[,2:1184] + 1
colnames(pseudo_serum) <- make.unique(colnames(pseudo_serum), sep = "_") #make sure that column names are all unique

# Apply the remove_outliers function to set outliers as NA within each column
pseudo_intake <- remove_outliers(pseudo_intake)
pseudo_fecal <- remove_outliers(pseudo_fecal)
pseudo_serum <- remove_outliers(pseudo_serum)

# Get the common IDs
common_ids <- row.names(pseudo_serum)[(row.names(pseudo_serum)) %in% (row.names(pseudo_fecal))]

# Apply the remove_outliers function to set outliers as NA within each column
pseudo_intake <- remove_outliers(pseudo_intake)
pseudo_fecal <- remove_outliers(pseudo_fecal)
pseudo_serum <- remove_outliers(pseudo_serum)

# Filter data frames based on common IDs
filtered_intake <- pseudo_intake[row.names(pseudo_intake) %in% common_ids, ]
filtered_fecal <- pseudo_fecal[row.names(pseudo_fecal) %in% common_ids, ]
filtered_serum <- pseudo_serum[row.names(pseudo_serum) %in% common_ids, ]

#select intake, fecal and serum metabolites columns
intake.cols <- as.character(colnames(filtered_intake[,1:1108]))
fecal.cols <- as.character(colnames(filtered_fecal[,1:1684]))
blood.cols <- as.character(colnames(filtered_serum[,1:1183]))

# Modify column names using paste0()
intake <- paste0("int_", intake.cols)
fecal <- paste0("fec_", fecal.cols)
blood <- paste0("ser_", blood.cols)

# Rename selected columns in the data frame
names(filtered_intake)[names(filtered_intake) %in% intake.cols] <- intake
names(filtered_fecal)[names(filtered_fecal) %in% fecal.cols] <- fecal
names(filtered_serum)[names(filtered_serum) %in% blood.cols] <- blood

diet <- read_xlsx('chem_diet_combined_V2.xlsx')
diet <- diet[,c(1, 1110:3094)]
diet_t <- as.data.frame(t(diet))
colnames(diet_t) <- diet_t[1,]
diet_t <- diet_t[-1,]

filtered_diet <- diet_t[row.names(diet_t) %in% common_ids, ]

##==================================================== CORRELATION MATRIX: INTAKE METABOLITES ===============================================

##===================================================== INTAKE MEASURED VS DIET METABOLITES ========================

intdiet <- merge(filtered_intake, filtered_diet, by = "row.names", all.x = TRUE)
rownames(intdiet) <- intdiet$Row.names
intdiet$Row.names <- NULL

#RANK transformation
intdiet <- intdiet %>% mutate_all(~ (length(.) + 1) - rank(.))

#filter full set on figure food items
# Find column names containing 'wine'
wine_cols <- grep("wine", colnames(intdiet), value = TRUE)
# Find column names containing 'meat'
meat_cols <- grep("meat", colnames(intdiet), value = TRUE)
# Find column names containing 'fish'
fish_cols <- grep("fish", colnames(intdiet), value = TRUE)
# Find column names containing 'vegetables'
vegetables_cols <- grep("vegetables", colnames(intdiet), value = TRUE)

diet_mtb <- c(wine_cols, meat_cols, fish_cols, vegetables_cols)

#filter on predicted metabolites in food items
# Find column names containing 'resveratrol'
resveratrol_cols <- grep("resveratrol", colnames(intdiet), value = TRUE)
# Find column names containing 'quercetin'
quercetin_cols <- grep("quercetin", colnames(intdiet), value = TRUE)
# Find column names containing 'leucine'
leucine_cols <- grep("leucine", colnames(intdiet), value = TRUE)
# Find column names containing 'EPA'
epa_cols <- grep("20_5", colnames(intdiet), value = TRUE)
# Find column names containing 'DHA'
dha_cols <- grep("22_6", colnames(intdiet), value = TRUE)

pred_mtb <- c(resveratrol_cols, quercetin_cols, leucine_cols, epa_cols, dha_cols)

# Combine diet and microbe cols
selected_cols <- c(diet_mtb, pred_mtb)

# new df with associated diet and microbe metabolites
int_corr <- intdiet[,selected_cols]

#create correlation matrix
cor = rcorr(as.matrix(int_corr), type = "pearson")

# create matrices for R, P and P.adj
cormatrix <- cor$r
pmatrix <- cor$P

# Subset matrix to only show correlations between intake and fecal metabolites
subset_rows <- !grepl("^int_", rownames(cormatrix))
subset_columns <- grepl("^int_", colnames(cormatrix))
subsetted_matrix <- cormatrix[subset_rows, subset_columns]

# Subset p-value to only show p-values for correlations between intake and fecal metabolites
subset_rows <- !grepl("^int_", rownames(pmatrix))
subset_columns <- grepl("^int_", colnames(pmatrix))
subset_p_matrix <- pmatrix[subset_rows, subset_columns]

melt.matrix <- melt(subsetted_matrix)
melt.matrix$p.value <- round(melt(subset_p_matrix)$value, 3) # add p-values to plot dataset
melt.matrix$p.value.adj <- p.adjust(melt.matrix$p.value, method = c('fdr')) # adjust P-values for multiple comparison
melt.matrix$stars <- cut(melt.matrix$p.value, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))  # Create column of significance labels
melt.matrix$stars.adj <- cut(melt.matrix$p.value.adj, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))

melt.matrix <- melt.matrix[!melt.matrix$value == -Inf,] #delete infinite correlation coefficients
melt.matrix <- melt.matrix[!melt.matrix$value == 0,] #delete rows where correlation coefficients = 0

# Sort melt.matrix based on the absolute values of the correlations
sorted_correlations <- melt.matrix[order(abs(melt.matrix$value), decreasing = TRUE), ]
# Select the top 100 correlations
top_100 <- sorted_correlations[1:100, ]

# ========== 2. PLOT MATRIX USING GGPLOT ======== ##

ggplot(sorted_correlations, aes(x = Var1, y = Var2)) +
  geom_tile(colour="grey20", aes(fill=value), size=0.2) + 
  scale_fill_gradient2(name = "Pearson's rho", low="navyblue", high="red", midpoint=0, limits=c(-1,1)) +
  labs(x="",y="") +
  ggtitle("\nCorrelations\nDietary intake and predicted metabolites\n") +
  geom_text(aes(label=stars.adj), position=position_nudge(y=0.15), color="black", size=3) + 
  geom_text(aes(label=sprintf("%1.2f", value)), position=position_nudge(y=-0.1), 
            size=2, colour="black") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle = -45, hjust = 0, size=8)) + 
  theme(axis.text.y=element_text(angle = 0, hjust = 1, size=8)) +
  theme(plot.title = element_text(hjust = 0.5))

##===== CARB/FIBER ========

#RANK transformation
int <- filtered_intake %>% mutate_all(~ (length(.) + 1) - rank(.))

#filter full metabolite set on carbohydrate/fiber + metabolites
# Find column names containing 'fiber'
fiber_cols <- grep("fiber", colnames(int), value = TRUE)
# Find column names containing 'carbohydrate'
carb_cols <- grep("carbohydrate", colnames(int), value = TRUE)

diet_mtb <- c(fiber_cols, carb_cols)

# Microbial metabolites
# Find column names containing 'propionate' (SCFA)
propionate_cols <- grep("propionate", colnames(int), value = TRUE)
# Find column names containing 'acetate' (SCFA)
acetate_cols <- grep("acetate", colnames(int), value = TRUE)
# Find column names containing 'butyrate' (SCFA)
butyrate_cols <- grep("butyrate", colnames(int), value = TRUE)
# Find column names containing 'lactate' (carboxylic acid)
lactate_cols <- grep("lactate", colnames(int), value = TRUE)
# Find column names containing 'succinate' (carboxylic acid)
succinate_cols <- grep("succinate", colnames(int), value = TRUE)

microbe_mtb <- c(propionate_cols, acetate_cols, butyrate_cols, lactate_cols, succinate_cols)


# Combine diet and microbe cols
selected_cols <- c(diet_mtb, microbe_mtb)

# new df with associated diet and microbe metabolites
int_carbfiber <- int[,selected_cols]


#create correlation matrix
cor = rcorr(as.matrix(int_carbfiber), type = "pearson")

# create matrices for R, P
cormatrix <- cor$r
pmatrix <- cor$P

melt.matrix <- melt(cormatrix)
melt.matrix$p.value <- round(melt(pmatrix)$value, 3) # add p-values to plot dataset
melt.matrix$p.value.adj <- p.adjust(melt.matrix$p.value, method = c('fdr')) # adjust P-values for multiple comparison
melt.matrix$stars <- cut(melt.matrix$p.value, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))  # Create column of significance labels
melt.matrix$stars.adj <- cut(melt.matrix$p.value.adj, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))

melt.matrix <- melt.matrix[!melt.matrix$value == -Inf,] #delete infinite correlation coefficients
melt.matrix <- melt.matrix[!melt.matrix$value == 0,] #delete rows where correlation coefficients = 0
melt.matrix <- melt.matrix[!is.na(melt.matrix$Var1),]
# Sort melt.matrix based on the absolute values of the correlations
sorted_correlations <- melt.matrix[order(abs(melt.matrix$value), decreasing = T), ]
# Select the top 100 correlations
top_100 <- sorted_correlations[1:500, ]


# ========== 2. PLOT MATRIX USING GGPLOT ======== ##

ggplot(sorted_correlations, aes(x = Var1, y = Var2)) +
  geom_tile(colour="grey20", aes(fill=value), size=0.2) + 
  scale_fill_gradient2(name = "Pearson's rho", low="navyblue", high="red", midpoint=0, limits=c(-1,1)) +
  labs(x="",y="") +
  ggtitle("\nCorrelations\nDietary carbohydrate/fiber intake metabolites\n") +
  geom_text(aes(label=stars.adj), position=position_nudge(y=0.15), color="black", size=3) + 
  geom_text(aes(label=sprintf("%1.2f", value)), position=position_nudge(y=-0.1), 
            size=2, colour="black") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle = -45, hjust = 0, size=8)) + 
  theme(axis.text.y=element_text(angle = 0, hjust = 1, size=8)) +
  theme(plot.title = element_text(hjust = 0.5))

##============== AMINO ACIDS ===============

#filter full metabolite set on amino acids + metabolites
# Find column names containing 'protein'
protein_cols <- grep("protein", colnames(int), value = TRUE)
# Find column names containing 'alanine'
alanine_cols <- grep("alanine", colnames(int), value = TRUE)
# Find column names containing 'arginine'
arginine_cols <- grep("arginine", colnames(int), value = TRUE)
# Find column names containing 'asparagine'
asparagine_cols <- grep("asparagine", colnames(int), value = TRUE)
# Find column names containing 'aspartic'
aspartic_cols <- grep("aspartic", colnames(int), value = TRUE)
# Find column names containing 'cysteine'
cysteine_cols <- grep("cysteine", colnames(int), value = TRUE)
# Find column names containing 'glutamic'
glutamic_cols <- grep("glutamic", colnames(int), value = TRUE)
# Find column names containing 'glutamine'
glutamine_cols <- grep("glutamine", colnames(int), value = TRUE)
# Find column names containing 'glycine'
glycine_cols <- grep("glycine", colnames(int), value = TRUE)
# Find column names containing 'histidine'
histidine_cols <- grep("histidine", colnames(int), value = TRUE)
# Find column names containing 'isoleucine'
isoleucine_cols <- grep("isoleucine", colnames(int), value = TRUE)
# Find column names containing 'leucine'
leucine_cols <- grep("leucine", colnames(int), value = TRUE)
# Find column names containing 'lysine'
lysine_cols <- grep("lysine", colnames(int), value = TRUE)
# Find column names containing 'methionine'
methionine_cols <- grep("methionine", colnames(int), value = TRUE)
# Find column names containing 'phenylalanine'
phenylalanine_cols <- grep("phenylalanine", colnames(int), value = TRUE)
# Find column names containing 'proline'
proline_cols <- grep("proline", colnames(int), value = TRUE)
# Find column names containing 'serine'
serine_cols <- grep("serine", colnames(int), value = TRUE)
# Find column names containing 'threonine'
threonine_cols <- grep("threonine", colnames(int), value = TRUE)
# Find column names containing 'tryptophan'
tryptophan_cols <- grep("tryptophan", colnames(int), value = TRUE)
# Find column names containing 'tyrosine'
tyrosine_cols <- grep("tyrosine", colnames(int), value = TRUE)
# Find column names containing 'alanine'
valine_cols <- grep("valine", colnames(int), value = TRUE)
# Find column names containing 'selenocysteine'
selenocysteine_cols <- grep("selenocysteine", colnames(int), value = TRUE)

diet_mtb <- c(protein_cols, alanine_cols, arginine_cols, asparagine_cols, aspartic_cols, 
              cysteine_cols, glutamine_cols, glutamic_cols, glycine_cols, histidine_cols, 
              isoleucine_cols, leucine_cols, lysine_cols, methionine_cols, phenylalanine_cols, 
              proline_cols, serine_cols, threonine_cols, tryptophan_cols, tyrosine_cols,
              valine_cols, selenocysteine_cols)

# Microbial metabolites
# Find column names containing 'valine' (BCAA)
valine_cols <- grep("valine", colnames(int), value = TRUE)
# Find column names containing 'isoleucine' (BCAA)
isoleucine_cols <- grep("isoleucine", colnames(int), value = TRUE)
# Find column names containing 'leucine' (BCAA)
leucine_cols <- grep("leucine", colnames(int), value = TRUE)
# Find column names containing 'niacin' 
niacin_cols <- grep("niacin", colnames(int), value = TRUE)
# Find column names containing 'nicotinamide'
nicotinamide_cols <- grep("nicotinamide", colnames(int), value = TRUE)
# Find column names containing 'aminovaleric'
aminovaleric_cols <- grep("aminovaleric", colnames(int), value = TRUE)
# Find column names containing 'dimethylglycine'
dimethylglycine_cols <- grep("dimethylglycine", colnames(int), value = TRUE)
# Find column names containing 'acetylglycine'
acetylglycine_cols <- grep("acetylglycine", colnames(int), value = TRUE)


microbe_mtb <- c(valine_cols, isoleucine_cols, leucine_cols, niacin_cols, nicotinamide_cols,
                 aminovaleric_cols, dimethylglycine_cols, acetylglycine_cols)


# Combine diet and microbe cols
selected_cols <- c(diet_mtb, microbe_mtb)

# new df with associated diet and microbe metabolites
int_aminoacid <- int[,selected_cols]

#create correlation matrix
cor = rcorr(as.matrix(int_aminoacid), type = "pearson")

# create matrices for R, P
cormatrix <- cor$r
pmatrix <- cor$P

melt.matrix <- melt(cormatrix)
melt.matrix$p.value <- round(melt(pmatrix)$value, 3) # add p-values to plot dataset
melt.matrix$p.value.adj <- p.adjust(melt.matrix$p.value, method = c('fdr')) # adjust P-values for multiple comparison
melt.matrix$stars <- cut(melt.matrix$p.value, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))  # Create column of significance labels
melt.matrix$stars.adj <- cut(melt.matrix$p.value.adj, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))

melt.matrix <- melt.matrix[!melt.matrix$value == -Inf,] #delete infinite correlation coefficients
melt.matrix <- melt.matrix[!melt.matrix$value == 0,] #delete rows where correlation coefficients = 0

# Sort melt.matrix based on the absolute values of the correlations
sorted_correlations <- melt.matrix[order(abs(melt.matrix$value), decreasing = TRUE), ]
# Select the top 100 correlations
top_100 <- sorted_correlations[1:500, ]

# ========== 2. PLOT MATRIX USING GGPLOT ======== ##

ggplot(sorted_correlations, aes(x = Var1, y = Var2)) +
  geom_tile(colour="grey20", aes(fill=value), size=0.2) + 
  scale_fill_gradient2(name = "Pearson's rho", low="navyblue", high="red", midpoint=0, limits=c(-1,1)) +
  labs(x="",y="") +
  ggtitle("\nCorrelations\nDietary intake amino acids and metabolites\n") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle = -45, hjust = 0, size=6)) + 
  theme(axis.text.y=element_text(angle = 0, hjust = 1, size=6)) +
  theme(plot.title = element_text(hjust = 0.5))

##======== TRYPTOPHAN =================

#filter full metabolite set on tryptophan + metabolites
# Find column names containing 'tryptophan'
tryptophan_cols <- grep("tryptophan", colnames(int), value = TRUE)

diet_mtb <- c(tryptophan_cols)

# Microbial metabolites
# Find column names containing 'indole' 
indole_cols <- grep("indole", colnames(int), value = TRUE)
# Find column names containing 'indoxylsulfate' (IS) 
indoxylsulfate_cols <- grep("indoxylsulfate", colnames(int), value = TRUE)
# Find column names containing 'tryptamine'
tryptamine_cols <- grep("tryptamine", colnames(int), value = TRUE)


microbe_mtb <- c(indole_cols, indoxylsulfate_cols, tryptamine_cols)

# Combine diet and microbe cols
selected_cols <- c(diet_mtb, microbe_mtb)

# new df with associated diet and microbe metabolites
int_tryptophan <- int[,selected_cols]


#create correlation matrix
cor = rcorr(as.matrix(int_tryptophan), type = "pearson")

# create matrices for R, P
cormatrix <- cor$r
pmatrix <- cor$P

melt.matrix <- melt(cormatrix)
melt.matrix$p.value <- round(melt(pmatrix)$value, 3) # add p-values to plot dataset
melt.matrix$p.value.adj <- p.adjust(melt.matrix$p.value, method = c('fdr')) # adjust P-values for multiple comparison
melt.matrix$stars <- cut(melt.matrix$p.value, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))  # Create column of significance labels
melt.matrix$stars.adj <- cut(melt.matrix$p.value.adj, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))

melt.matrix <- melt.matrix[!melt.matrix$value == -Inf,] #delete infinite correlation coefficients
melt.matrix <- melt.matrix[!melt.matrix$value == 0,] #delete rows where correlation coefficients = 0
melt.matrix <- melt.matrix[!is.na(melt.matrix$value),]
# Sort melt.matrix based on the absolute values of the correlations
sorted_correlations <- melt.matrix[order(abs(melt.matrix$value), decreasing = TRUE), ]
# Select the top 100 correlations
top_100 <- sorted_correlations[1:500, ]


# ========== 2. PLOT MATRIX USING GGPLOT ======== ##

ggplot(sorted_correlations, aes(x = Var1, y = Var2)) +
  geom_tile(colour="grey20", aes(fill=value), size=0.2) + 
  scale_fill_gradient2(name = "Pearson's rho", low="navyblue", high="red", midpoint=0, limits=c(-1,1)) +
  labs(x="",y="") +
  ggtitle("\nCorrelations\nDietary tryptophan and metabolites\n") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle = -45, hjust = 0, size=6)) + 
  theme(axis.text.y=element_text(angle = 0, hjust = 1, size=6)) +
  theme(plot.title = element_text(hjust = 0.5))

##================ POLYPHENOLS ==============

#filter full metabolite set on polyphenols + metabolites
# Find column names containing 'polyphenol'
polyphenol_cols <- grep("polyphenol", colnames(int), value = TRUE)

#flavonoids
# Find column names containing 'epicatechin'
epicatechin_cols <- grep("epicatechin", colnames(int), value = TRUE)
# Find column names containing 'catechin'
catechin_cols <- grep("catechin", colnames(int), value = TRUE)
# Find column names containing 'epigallocatechin'
epigallocatechin_cols <- grep("epigallocatechin", colnames(int), value = TRUE)
# Find column names containing 'hesperatin'
hesperatin_cols <- grep("hesperatin", colnames(int), value = TRUE)
# Find column names containing 'naringenin'
naringenin_cols <- grep("naringenin", colnames(int), value = TRUE)
# Find column names containing 'eriodictyol'
eriodictyol_cols <- grep("eriodictyol", colnames(int), value = TRUE)
# Find column names containing 'apigenin'
apigenin_cols <- grep("apigenin", colnames(int), value = TRUE)
# Find column names containing 'luteolin'
luteolin_cols <- grep("luteolin", colnames(int), value = TRUE)
# Find column names containing 'tangeritin'
tangeritin_cols <- grep("tangeritin", colnames(int), value = TRUE)
# Find column names containing 'chrysin'
chrysin_cols <- grep("chrysin", colnames(int), value = TRUE)
# Find column names containing 'genistein'
genistein_cols <- grep("genistein", colnames(int), value = TRUE)
# Find column names containing 'daidzein'
daidzein_cols <- grep("daidzein", colnames(int), value = TRUE)
# Find column names containing 'kaempferol'
kaempferol_cols <- grep("kaempferol", colnames(int), value = TRUE)
# Find column names containing 'myrestin'
myrestin_cols <- grep("myrestin", colnames(int), value = TRUE)
# Find column names containing 'quercetin'
quercetin_cols <- grep("quercetin", colnames(int), value = TRUE)
# Find column names containing 'cyanidin'
cyanidin_cols <- grep("cyanidin", colnames(int), value = TRUE)
# Find column names containing 'delphinidin'
delphinidin_cols <- grep("delphinidin", colnames(int), value = TRUE)
# Find column names containing 'malvedin'
malvedin_cols <- grep("malvedin", colnames(int), value = TRUE)
# Find column names containing 'pelargonidin'
pelargonidin_cols <- grep("pelargonidin", colnames(int), value = TRUE)

flavonoids <- c(polyphenol_cols, epicatechin_cols, catechin_cols, epigallocatechin_cols,
                hesperatin_cols, naringenin_cols, eriodictyol_cols, apigenin_cols, luteolin_cols,
                tangeritin_cols, chrysin_cols, genistein_cols, daidzein_cols, kaempferol_cols,
                myrestin_cols, quercetin_cols, cyanidin_cols, delphinidin_cols, malvedin_cols, 
                pelargonidin_cols)

#phenolic acids
# Find column names containing 'hydroxybenzoic'
hydroxybenzoic_cols <- grep("hydroxybenzoic", colnames(int), value = TRUE)
# Find column names containing 'protocatechuic'
protocatechuic_cols <- grep("protocatechuic", colnames(int), value = TRUE)
# Find column names containing 'gallic'
gallic_cols <- grep("gallic", colnames(int), value = TRUE)
# Find column names containing 'vanillic'
vanillic_cols <- grep("vanillic", colnames(int), value = TRUE)
# Find column names containing 'ellagic'
ellagic_cols <- grep("ellagic", colnames(int), value = TRUE)
# Find column names containing 'salicyclic'
salicyclic_cols <- grep("salicyclic", colnames(int), value = TRUE)
# Find column names containing 'caffeic'
caffeic_cols <- grep("caffeic", colnames(int), value = TRUE)
# Find column names containing 'ferulic'
ferulic_cols <- grep("ferulic", colnames(int), value = TRUE)
# Find column names containing 'sinapic'
sinapic_cols <- grep("sinapic", colnames(int), value = TRUE)
# Find column names containing 'chlorogenic'
chlorogenic_cols <- grep("chlorogenic", colnames(int), value = TRUE)
# Find column names containing 'coumaric'
coumaric_cols <- grep("coumaric", colnames(int), value = TRUE)
# Find column names containing 'quinic'
quinic_cols <- grep("quinic", colnames(int), value = TRUE)
# Find column names containing 'resveratrol'
resveratrol_cols <- grep("resveratrol", colnames(int), value = TRUE)

phenols <- c(hydroxybenzoic_cols, protocatechuic_cols, gallic_cols, vanillic_cols,
             ellagic_cols, salicyclic_cols, caffeic_cols, ferulic_cols, sinapic_cols,
             chlorogenic_cols, coumaric_cols, quinic_cols, resveratrol_cols)

diet_mtb <- c(flavonoids, phenols)

# Microbial metabolites
# Find column names containing 'quercetin' 
quercetin_cols <- grep("quercetin", colnames(int), value = TRUE)
# Find column names containing 'apigenin'
apigenin_cols <- grep("apigenin", colnames(int), value = TRUE)
# Find column names containing 'naringenin'
naringenin_cols <- grep("naringenin", colnames(int), value = TRUE)

microbe_mtb <- c(quercetin_cols, apigenin_cols, naringenin_cols)

# Combine diet and microbe cols
selected_cols <- c(diet_mtb, microbe_mtb)

# new df with associated diet and microbe metabolites
int_polyphenol <- int[,selected_cols]

#create correlation matrix
cor = rcorr(as.matrix(int_polyphenol), type = "pearson")

# create matrices for R, P 
cormatrix <- cor$r
pmatrix <- cor$P

melt.matrix <- melt(cormatrix)
melt.matrix$p.value <- round(melt(pmatrix)$value, 3) # add p-values to plot dataset
melt.matrix$p.value.adj <- p.adjust(melt.matrix$p.value, method = c('fdr')) # adjust P-values for multiple comparison
melt.matrix$stars <- cut(melt.matrix$p.value, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))  # Create column of significance labels
melt.matrix$stars.adj <- cut(melt.matrix$p.value.adj, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))

melt.matrix <- melt.matrix[!melt.matrix$value == -Inf,] #delete infinite correlation coefficients
melt.matrix <- melt.matrix[!melt.matrix$value == 0,] #delete rows where correlation coefficients = 0
melt.matrix <- melt.matrix[!is.na(melt.matrix$value),]
# Sort melt.matrix based on the absolute values of the correlations
sorted_correlations <- melt.matrix[order(abs(melt.matrix$value), decreasing = TRUE), ]
# Select the top 100 correlations
top_100 <- sorted_correlations[1:5000, ]

# ========== 2. PLOT MATRIX USING GGPLOT ======== ##

ggplot(sorted_correlations, aes(x = Var1, y = Var2)) +
  geom_tile(colour="grey20", aes(fill=value), size=0.2) + 
  scale_fill_gradient2(name = "Pearson's rho", low="navyblue", high="red", midpoint=0, limits=c(-1,1)) +
  labs(x="",y="") +
  ggtitle("\nCorrelations\nDietary tryptophan intake + metabolites\n") +
  geom_text(aes(label=stars.adj), position=position_nudge(y=0.15), color="black", size=3) + 
  geom_text(aes(label=sprintf("%1.2f", value)), position=position_nudge(y=-0.1), 
            size=2, colour="black") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle = -45, hjust = 0, size=8)) + 
  theme(axis.text.y=element_text(angle = 0, hjust = 1, size=8)) +
  theme(plot.title = element_text(hjust = 0.5))


##==================================================== CORRELATION MATRIX: INTAKE VS FECAL METABOLITES ======================================================

#merge intake + fecal metabolites
intfec <- merge(filtered_intake, filtered_fecal, by = "row.names", all.x = TRUE)
rownames(intfec) <- intfec$Row.names
intfec$Row.names <- NULL

#create correlation matrix
cor = rcorr(as.matrix(intfec), type = "pearson")
cor$P.adj <- p.adjust(cor$P, method = c('fdr')) # adjust P-values for multiple comparison
cormatrix <- cor$r

# Subset matrix to only show correlations between intake and fecal metabolites
subset_rows <- grepl("^fec_", rownames(cormatrix))
subset_columns <- grepl("^int_", colnames(cormatrix))
subsetted_matrix <- cormatrix[subset_rows, subset_columns]

melt.matrix <- melt(subsetted_matrix)
melt.matrix <- melt.matrix[!melt.matrix$value == -Inf,] #delete infinite correlation coefficients
melt.matrix <- melt.matrix[!melt.matrix$value == 0,] #delete rows where correlation coefficients = 0

# Sort melt.matrix based on the absolute values of the correlations
sorted_correlations <- melt.matrix[order(abs(melt.matrix$value), decreasing = TRUE), ]
# Select the top 100 correlations
top_100 <- sorted_correlations[1:100, ]


ggscatter(intfec, x = "int_linolenic_acid", y = "fec_X...1483", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "linolenic acid", ylab = "fecal metabolite 1483")

##==================================================== CORRELATION MATRIX: INTAKE VS SERUM METABOLITES ======================================================

#merge intake + serum metabolites
intser <- merge(filtered_intake, filtered_serum, by = "row.names", all.x = TRUE)
rownames(intser) <- intser$Row.names
intser$Row.names <- NULL

#create correlation matrix
cor = rcorr(as.matrix(intser), type = "pearson")
cor$P.adj <- p.adjust(cor$P, method = c('fdr')) # adjust P-values for multiple comparison
cormatrix <- cor$r

# Subset matrix to only show correlations between intake and serum metabolites
subset_rows <- grepl("^ser_", rownames(cormatrix))
subset_columns <- grepl("^int_", colnames(cormatrix))
subsetted_matrix <- cormatrix[subset_rows, subset_columns]

melt.matrix <- melt(subsetted_matrix)
melt.matrix <- melt.matrix[!melt.matrix$value == -Inf,] #delete infinite correlation coefficients
melt.matrix <- melt.matrix[!melt.matrix$value == 0,] #delete rows where correlation coefficients = 0

# Sort melt.matrix based on the absolute values of the correlations
sorted_correlations <- melt.matrix[order(abs(melt.matrix$value), decreasing = TRUE), ]
# Select the top 100 correlations
top_100 <- sorted_correlations[1:100, ]


ggscatter(intser, x = "int_X20_0", y = "ser_meta_994", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Arachidic acid", ylab = "Glycerol triundecanoate")

##================================================ FECES ================================================

##========================== CORRELATION CARBOHYDRATE/FIBER METABOLITES =============================

#merge intake + fecal metabolites
intfec <- merge(filtered_intake, filtered_fecal, by = "row.names", all.x = TRUE)
rownames(intfec) <- intfec$Row.names
intfec$Row.names <- NULL

#RANK transformation
intfec <- intfec %>% mutate_all(~ (length(.) + 1) - rank(.))

#filter full metabolite set on carbohydrate/fiber + fecal metabolites
# Find column names containing 'fiber'
fiber_cols <- grep("fiber", colnames(intfec), value = TRUE)
# Find column names containing 'carbohydrate'
carb_cols <- grep("carbohydrate", colnames(intfec), value = TRUE)

diet_mtb <- c(fiber_cols, carb_cols)

# Microbial metabolites
# Find column names containing 'propionate' (SCFA)
propionate_cols <- grep("propionate", colnames(intfec), value = TRUE)
# Find column names containing 'acetate' (SCFA)
acetate_cols <- grep("acetate", colnames(intfec), value = TRUE)
# Find column names containing 'butyrate' (SCFA)
butyrate_cols <- grep("butyrate", colnames(intfec), value = TRUE)
# Find column names containing 'lactate' (carboxylic acid)
lactate_cols <- grep("lactate", colnames(intfec), value = TRUE)
# Find column names containing 'succinate' (carboxylic acid)
succinate_cols <- grep("succinate", colnames(intfec), value = TRUE)

microbe_mtb <- c(propionate_cols, acetate_cols, butyrate_cols, lactate_cols, succinate_cols)


# Combine diet and microbe cols
selected_cols <- c(diet_mtb, microbe_mtb)

# new df with associated diet and microbe metabolites
intfec_carbfiber <- intfec[,selected_cols]


#create correlation matrix
cor = rcorr(as.matrix(intfec_carbfiber), type = "pearson")

# create matrices for R, P and P.adj
cormatrix <- cor$r
pmatrix <- cor$P

# Subset matrix to only show correlations between intake and fecal metabolites
subset_rows <- grepl("^fec_", rownames(cormatrix))
subset_columns <- grepl("^int_", colnames(cormatrix))
subsetted_matrix <- cormatrix[subset_rows, subset_columns]

# Subset p-value to only show p-values for correlations between intake and fecal metabolites
subset_rows <- grepl("^fec_", rownames(pmatrix))
subset_columns <- grepl("^int_", colnames(pmatrix))
subset_p_matrix <- pmatrix[subset_rows, subset_columns]

melt.matrix <- melt(subsetted_matrix)
melt.matrix$p.value <- round(melt(subset_p_matrix)$value, 3) # add p-values to plot dataset
melt.matrix$p.value.adj <- p.adjust(melt.matrix$p.value, method = c('fdr')) # adjust P-values for multiple comparison
melt.matrix$stars <- cut(melt.matrix$p.value, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))  # Create column of significance labels
melt.matrix$stars.adj <- cut(melt.matrix$p.value.adj, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))

melt.matrix <- melt.matrix[!melt.matrix$value == -Inf,] #delete infinite correlation coefficients
melt.matrix <- melt.matrix[!melt.matrix$value == 0,] #delete rows where correlation coefficients = 0

# Sort melt.matrix based on the absolute values of the correlations
sorted_correlations <- melt.matrix[order(abs(melt.matrix$value), decreasing = TRUE), ]
# Select the top 100 correlations
top_100 <- sorted_correlations[1:100, ]

# ========== 2. PLOT MATRIX USING GGPLOT ======== ##

ggplot(top_100, aes(x = Var1, y = Var2)) +
  geom_tile(colour="grey20", aes(fill=value), size=0.2) + 
  scale_fill_gradient2(name = "Pearson's rho", low="navyblue", high="red", midpoint=0, limits=c(-1,1)) +
  labs(x="",y="") +
  ggtitle("\nCorrelations\nDietary carbohydrate/fiber and microbial metabolites\n") +
  geom_text(aes(label=stars.adj), position=position_nudge(y=0.15), color="black", size=3) + 
  geom_text(aes(label=sprintf("%1.2f", value)), position=position_nudge(y=-0.1), 
            size=2, colour="black") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle = -45, hjust = 0, size=8)) + 
  theme(axis.text.y=element_text(angle = 0, hjust = 1, size=8)) +
  theme(plot.title = element_text(hjust = 0.5))

ggscatter(intfec_carbfiber, x = "int_fiber_dietary", y = "fec_methylsuccinate", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "dietary fiber", ylab = "methylsuccinate")

##========================== CORRELATION AMINO ACID METABOLITES =============================

#merge intake + fecal metabolites
intfec <- merge(filtered_intake, filtered_fecal, by = "row.names", all.x = TRUE)
rownames(intfec) <- intfec$Row.names
intfec$Row.names <- NULL

#RANK transformation
intfec <- intfec %>% mutate_all(~ (length(.) + 1) - rank(.))

#filter full metabolite set on amino acids + fecal metabolites
# Find column names containing 'protein'
protein_cols <- grep("protein", colnames(intfec), value = TRUE)
# Find column names containing 'alanine'
alanine_cols <- grep("alanine", colnames(intfec), value = TRUE)
# Find column names containing 'arginine'
arginine_cols <- grep("arginine", colnames(intfec), value = TRUE)
# Find column names containing 'asparagine'
asparagine_cols <- grep("asparagine", colnames(intfec), value = TRUE)
# Find column names containing 'aspartic'
aspartic_cols <- grep("aspartic", colnames(intfec), value = TRUE)
# Find column names containing 'cysteine'
cysteine_cols <- grep("cysteine", colnames(intfec), value = TRUE)
# Find column names containing 'glutamic'
glutamic_cols <- grep("glutamic", colnames(intfec), value = TRUE)
# Find column names containing 'glutamine'
glutamine_cols <- grep("glutamine", colnames(intfec), value = TRUE)
# Find column names containing 'glycine'
glycine_cols <- grep("glycine", colnames(intfec), value = TRUE)
# Find column names containing 'histidine'
histidine_cols <- grep("histidine", colnames(intfec), value = TRUE)
# Find column names containing 'isoleucine'
isoleucine_cols <- grep("isoleucine", colnames(intfec), value = TRUE)
# Find column names containing 'leucine'
leucine_cols <- grep("leucine", colnames(intfec), value = TRUE)
# Find column names containing 'lysine'
lysine_cols <- grep("lysine", colnames(intfec), value = TRUE)
# Find column names containing 'methionine'
methionine_cols <- grep("methionine", colnames(intfec), value = TRUE)
# Find column names containing 'phenylalanine'
phenylalanine_cols <- grep("phenylalanine", colnames(intfec), value = TRUE)
# Find column names containing 'proline'
proline_cols <- grep("proline", colnames(intfec), value = TRUE)
# Find column names containing 'serine'
serine_cols <- grep("serine", colnames(intfec), value = TRUE)
# Find column names containing 'threonine'
threonine_cols <- grep("threonine", colnames(intfec), value = TRUE)
# Find column names containing 'tryptophan'
tryptophan_cols <- grep("tryptophan", colnames(intfec), value = TRUE)
# Find column names containing 'tyrosine'
tyrosine_cols <- grep("tyrosine", colnames(intfec), value = TRUE)
# Find column names containing 'alanine'
valine_cols <- grep("valine", colnames(intfec), value = TRUE)
# Find column names containing 'selenocysteine'
selenocysteine_cols <- grep("selenocysteine", colnames(intfec), value = TRUE)

diet_mtb <- c(protein_cols, alanine_cols, arginine_cols, asparagine_cols, aspartic_cols, 
              cysteine_cols, glutamine_cols, glutamic_cols, glycine_cols, histidine_cols, 
              isoleucine_cols, leucine_cols, lysine_cols, methionine_cols, phenylalanine_cols, 
              proline_cols, serine_cols, threonine_cols, tryptophan_cols, tyrosine_cols,
              valine_cols, selenocysteine_cols)

# Microbial metabolites
# Find column names containing 'valine' (BCAA)
valine_cols <- grep("valine", colnames(intfec), value = TRUE)
# Find column names containing 'isoleucine' (BCAA)
isoleucine_cols <- grep("isoleucine", colnames(intfec), value = TRUE)
# Find column names containing 'leucine' (BCAA)
leucine_cols <- grep("leucine", colnames(intfec), value = TRUE)
# Find column names containing 'niacin' 
niacin_cols <- grep("niacin", colnames(intfec), value = TRUE)
# Find column names containing 'nicotinamide'
nicotinamide_cols <- grep("nicotinamide", colnames(intfec), value = TRUE)
# Find column names containing 'aminovaleric'
aminovaleric_cols <- grep("aminovaleric", colnames(intfec), value = TRUE)
# Find column names containing 'dimethylglycine'
dimethylglycine_cols <- grep("dimethylglycine", colnames(intfec), value = TRUE)
# Find column names containing 'acetylglycine'
acetylglycine_cols <- grep("acetylglycine", colnames(intfec), value = TRUE)


microbe_mtb <- c(valine_cols, isoleucine_cols, leucine_cols, niacin_cols, nicotinamide_cols,
                 aminovaleric_cols, dimethylglycine_cols, acetylglycine_cols)


# Combine diet and microbe cols
selected_cols <- c(diet_mtb, microbe_mtb)

# new df with associated diet and microbe metabolites
intfec_aminoacid <- intfec[,selected_cols]


#create correlation matrix
cor = rcorr(as.matrix(intfec_aminoacid), type = "pearson")
cor$P.adj <- p.adjust(cor$P, method = c('fdr')) # adjust P-values for multiple comparison

# create matrices for R, P and P.adj
cormatrix <- cor$r
pmatrix <- cor$P
padjmatrix <- cor$P.adj

# Subset matrix to only show correlations between intake and fecal metabolites
subset_rows <- grepl("^fec_", rownames(cormatrix))
subset_columns <- grepl("^int_", colnames(cormatrix))
subsetted_matrix <- cormatrix[subset_rows, subset_columns]

# Subset p-value to only show p-values for correlations between intake and fecal metabolites
subset_rows <- grepl("^fec_", rownames(cormatrix))
subset_columns <- grepl("^int_", colnames(cormatrix))
subset_p_matrix <- pmatrix[subset_rows, subset_columns]

melt.matrix <- melt(subsetted_matrix)
melt.matrix$p.value <- round(melt(subset_p_matrix)$value, 3) # add p-values to plot dataset
melt.matrix$p.value.adj <- p.adjust(melt.matrix$p.value, method = c('fdr')) # adjust P-values for multiple comparison
melt.matrix$stars <- cut(melt.matrix$p.value, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))  # Create column of significance labels
melt.matrix$stars.adj <- cut(melt.matrix$p.value.adj, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))

melt.matrix <- melt.matrix[!melt.matrix$value == -Inf,] #delete infinite correlation coefficients
melt.matrix <- melt.matrix[!melt.matrix$value == 0,] #delete rows where correlation coefficients = 0

# Sort melt.matrix based on the absolute values of the correlations
sorted_correlations <- melt.matrix[order(abs(melt.matrix$value), decreasing = TRUE), ]
# Select the top 100 correlations
top_100 <- sorted_correlations[1:100, ]

# ========== 2. PLOT MATRIX USING GGPLOT ======== ##

ggplot(top_100, aes(x = Var1, y = Var2)) +
  geom_tile(colour="grey20", aes(fill=value), size=0.2) + 
  scale_fill_gradient2(name = "Pearson's rho", low="navyblue", high="red", midpoint=0, limits=c(-1,1)) +
  labs(x="",y="") +
  ggtitle("\nCorrelations\nDietary amino acids and microbial metabolites\n") +
  geom_text(aes(label=stars.adj), position=position_nudge(y=0.15), color="black", size=3) + 
  geom_text(aes(label=sprintf("%1.2f", value)), position=position_nudge(y=-0.1), 
            size=2, colour="black") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle = -45, hjust = 0, size=8)) + 
  theme(axis.text.y=element_text(angle = 0, hjust = 1, size=8)) +
  theme(plot.title = element_text(hjust = 0.5))

ggscatter(intfec_aminoacid, x = "int_alpha_tryptophan", y = "fec_N.lactoyl.leucine", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "dietary alpha tryptophan", ylab = "N-lactoyl-leucine")

##========================== CORRELATION TRYPTOPHAN METABOLITES =============================

#merge intake + fecal metabolites
intfec <- merge(filtered_intake, filtered_fecal, by = "row.names", all.x = TRUE)
rownames(intfec) <- intfec$Row.names
intfec$Row.names <- NULL

#RANK transformation
intfec <- intfec %>% mutate_all(~ (length(.) + 1) - rank(.))

#filter full metabolite set on tryptophan + fecal metabolites
# Find column names containing 'tryptophan'
tryptophan_cols <- grep("tryptophan", colnames(intfec), value = TRUE)

diet_mtb <- c(tryptophan_cols)

# Microbial metabolites
# Find column names containing 'indole' 
indole_cols <- grep("indole", colnames(intfec), value = TRUE)
# Find column names containing 'indoxylsulfate' (IS) 
indoxylsulfate_cols <- grep("indoxylsulfate", colnames(intfec), value = TRUE)
# Find column names containing 'tryptamine'
tryptamine_cols <- grep("tryptamine", colnames(intfec), value = TRUE)


microbe_mtb <- c(indole_cols, indoxylsulfate_cols, tryptamine_cols)

# Combine diet and microbe cols
selected_cols <- c(diet_mtb, microbe_mtb)

# new df with associated diet and microbe metabolites
intfec_tryptophan <- intfec[,selected_cols]


#create correlation matrix
cor = rcorr(as.matrix(intfec_tryptophan), type = "pearson")
cor$P.adj <- p.adjust(cor$P, method = c('fdr')) # adjust P-values for multiple comparison

# create matrices for R, P and P.adj
cormatrix <- cor$r
pmatrix <- cor$P
padjmatrix <- cor$P.adj

# Subset matrix to only show correlations between intake and fecal metabolites
subset_rows <- grepl("^fec_", rownames(cormatrix))
subset_columns <- grepl("^int_", colnames(cormatrix))
subsetted_matrix <- cormatrix[subset_rows, subset_columns]

# Subset p-value to only show p-values for correlations between intake and fecal metabolites
subset_rows <- grepl("^fec_", rownames(cormatrix))
subset_columns <- grepl("^int_", colnames(cormatrix))
subset_p_matrix <- pmatrix[subset_rows, subset_columns]

melt.matrix <- melt(subsetted_matrix)
melt.matrix$p.value <- round(melt(subset_p_matrix)$value, 3) # add p-values to plot dataset
melt.matrix$p.value.adj <- p.adjust(melt.matrix$p.value, method = c('fdr')) # adjust P-values for multiple comparison
melt.matrix$stars <- cut(melt.matrix$p.value, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))  # Create column of significance labels
melt.matrix$stars.adj <- cut(melt.matrix$p.value.adj, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))

melt.matrix <- melt.matrix[!melt.matrix$value == -Inf,] #delete infinite correlation coefficients
melt.matrix <- melt.matrix[!melt.matrix$value == 0,] #delete rows where correlation coefficients = 0
melt.matrix <- melt.matrix[!is.na(melt.matrix$value),]
# Sort melt.matrix based on the absolute values of the correlations
sorted_correlations <- melt.matrix[order(abs(melt.matrix$value), decreasing = TRUE), ]
# Select the top 100 correlations
top_100 <- sorted_correlations[1:72, ]


# ========== 2. PLOT MATRIX USING GGPLOT ======== ##

ggplot(top_100, aes(x = Var1, y = Var2)) +
  geom_tile(colour="grey20", aes(fill=value), size=0.2) + 
  scale_fill_gradient2(name = "Pearson's rho", low="navyblue", high="red", midpoint=0, limits=c(-1,1)) +
  labs(x="",y="") +
  ggtitle("\nCorrelations\nDietary tryptophan and microbial metabolites\n") +
  geom_text(aes(label=stars.adj), position=position_nudge(y=0.15), color="black", size=3) + 
  geom_text(aes(label=sprintf("%1.2f", value)), position=position_nudge(y=-0.1), 
            size=2, colour="black") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle = -45, hjust = 0, size=8)) + 
  theme(axis.text.y=element_text(angle = 0, hjust = 1, size=8)) +
  theme(plot.title = element_text(hjust = 0.5))

ggscatter(intfec_tryptophan, x = "int_d_tryptophan", y = "fec_X5.hydroxyindoleacetate", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "dietary D-tryptophan", ylab = "5-hydroxyindoleacetate")

##========================== CORRELATION POLYPHENOL METABOLITES =============================

#merge intake + fecal metabolites
intfec <- merge(filtered_intake, filtered_fecal, by = "row.names", all.x = TRUE)
rownames(intfec) <- intfec$Row.names
intfec$Row.names <- NULL

#RANK transformation
intfec <- intfec %>% mutate_all(~ (length(.) + 1) - rank(.))

#filter full metabolite set on polyphenols + fecal metabolites
# Find column names containing 'polyphenol'
polyphenol_cols <- grep("polyphenol", colnames(intfec), value = TRUE)

#flavonoids
# Find column names containing 'epicatechin'
epicatechin_cols <- grep("epicatechin", colnames(intfec), value = TRUE)
# Find column names containing 'catechin'
catechin_cols <- grep("catechin", colnames(intfec), value = TRUE)
# Find column names containing 'epigallocatechin'
epigallocatechin_cols <- grep("epigallocatechin", colnames(intfec), value = TRUE)
# Find column names containing 'hesperatin'
hesperatin_cols <- grep("hesperatin", colnames(intfec), value = TRUE)
# Find column names containing 'naringenin'
naringenin_cols <- grep("naringenin", colnames(intfec), value = TRUE)
# Find column names containing 'eriodictyol'
eriodictyol_cols <- grep("eriodictyol", colnames(intfec), value = TRUE)
# Find column names containing 'apigenin'
apigenin_cols <- grep("apigenin", colnames(intfec), value = TRUE)
# Find column names containing 'luteolin'
luteolin_cols <- grep("luteolin", colnames(intfec), value = TRUE)
# Find column names containing 'tangeritin'
tangeritin_cols <- grep("tangeritin", colnames(intfec), value = TRUE)
# Find column names containing 'chrysin'
chrysin_cols <- grep("chrysin", colnames(intfec), value = TRUE)
# Find column names containing 'genistein'
genistein_cols <- grep("genistein", colnames(intfec), value = TRUE)
# Find column names containing 'daidzein'
daidzein_cols <- grep("daidzein", colnames(intfec), value = TRUE)
# Find column names containing 'kaempferol'
kaempferol_cols <- grep("kaempferol", colnames(intfec), value = TRUE)
# Find column names containing 'myrestin'
myrestin_cols <- grep("myrestin", colnames(intfec), value = TRUE)
# Find column names containing 'quercetin'
quercetin_cols <- grep("quercetin", colnames(intfec), value = TRUE)
# Find column names containing 'cyanidin'
cyanidin_cols <- grep("cyanidin", colnames(intfec), value = TRUE)
# Find column names containing 'delphinidin'
delphinidin_cols <- grep("delphinidin", colnames(intfec), value = TRUE)
# Find column names containing 'malvedin'
malvedin_cols <- grep("malvedin", colnames(intfec), value = TRUE)
# Find column names containing 'pelargonidin'
pelargonidin_cols <- grep("pelargonidin", colnames(intfec), value = TRUE)

flavonoids <- c(polyphenol_cols, epicatechin_cols, catechin_cols, epigallocatechin_cols,
                hesperatin_cols, naringenin_cols, eriodictyol_cols, apigenin_cols, luteolin_cols,
                tangeritin_cols, chrysin_cols, genistein_cols, daidzein_cols, kaempferol_cols,
                myrestin_cols, quercetin_cols, cyanidin_cols, delphinidin_cols, malvedin_cols, 
                pelargonidin_cols)

#phenolic acids
# Find column names containing 'hydroxybenzoic'
hydroxybenzoic_cols <- grep("hydroxybenzoic", colnames(intfec), value = TRUE)
# Find column names containing 'protocatechuic'
protocatechuic_cols <- grep("protocatechuic", colnames(intfec), value = TRUE)
# Find column names containing 'gallic'
gallic_cols <- grep("gallic", colnames(intfec), value = TRUE)
# Find column names containing 'vanillic'
vanillic_cols <- grep("vanillic", colnames(intfec), value = TRUE)
# Find column names containing 'ellagic'
ellagic_cols <- grep("ellagic", colnames(intfec), value = TRUE)
# Find column names containing 'salicyclic'
salicyclic_cols <- grep("salicyclic", colnames(intfec), value = TRUE)
# Find column names containing 'caffeic'
caffeic_cols <- grep("caffeic", colnames(intfec), value = TRUE)
# Find column names containing 'ferulic'
ferulic_cols <- grep("ferulic", colnames(intfec), value = TRUE)
# Find column names containing 'sinapic'
sinapic_cols <- grep("sinapic", colnames(intfec), value = TRUE)
# Find column names containing 'chlorogenic'
chlorogenic_cols <- grep("chlorogenic", colnames(intfec), value = TRUE)
# Find column names containing 'coumaric'
coumaric_cols <- grep("coumaric", colnames(intfec), value = TRUE)
# Find column names containing 'quinic'
quinic_cols <- grep("quinic", colnames(intfec), value = TRUE)
# Find column names containing 'resveratrol'
resveratrol_cols <- grep("resveratrol", colnames(intfec), value = TRUE)

phenols <- c(hydroxybenzoic_cols, protocatechuic_cols, gallic_cols, vanillic_cols,
             ellagic_cols, salicyclic_cols, caffeic_cols, ferulic_cols, sinapic_cols,
             chlorogenic_cols, coumaric_cols, quinic_cols, resveratrol_cols)

diet_mtb <- c(flavonoids, phenols)

# Microbial metabolites
# Find column names containing 'quercetin' 
quercetin_cols <- grep("quercetin", colnames(intfec), value = TRUE)
# Find column names containing 'apigenin'
apigenin_cols <- grep("apigenin", colnames(intfec), value = TRUE)
# Find column names containing 'naringenin'
naringenin_cols <- grep("naringenin", colnames(intfec), value = TRUE)

microbe_mtb <- c(quercetin_cols, apigenin_cols, naringenin_cols)

# Combine diet and microbe cols
selected_cols <- c(diet_mtb, microbe_mtb)

# new df with associated diet and microbe metabolites
intfec_polyphenol <- intfec[,selected_cols]

#create correlation matrix
cor = rcorr(as.matrix(intfec_polyphenol), type = "pearson")
cor$P.adj <- p.adjust(cor$P, method = c('fdr')) # adjust P-values for multiple comparison

# create matrices for R, P and P.adj
cormatrix <- cor$r
pmatrix <- cor$P
padjmatrix <- cor$P.adj

# Subset matrix to only show correlations between intake and fecal metabolites
subset_rows <- grepl("^fec_", rownames(cormatrix))
subset_columns <- grepl("^int_", colnames(cormatrix))
subsetted_matrix <- cormatrix[subset_rows, subset_columns]

# Subset p-value to only show p-values for correlations between intake and fecal metabolites
subset_rows <- grepl("^fec_", rownames(cormatrix))
subset_columns <- grepl("^int_", colnames(cormatrix))
subset_p_matrix <- pmatrix[subset_rows, subset_columns]

melt.matrix <- melt(subsetted_matrix)
melt.matrix$p.value <- round(melt(subset_p_matrix)$value, 3) # add p-values to plot dataset
melt.matrix$p.value.adj <- p.adjust(melt.matrix$p.value, method = c('fdr')) # adjust P-values for multiple comparison
melt.matrix$stars <- cut(melt.matrix$p.value, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))  # Create column of significance labels
melt.matrix$stars.adj <- cut(melt.matrix$p.value.adj, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))

melt.matrix <- melt.matrix[!melt.matrix$value == -Inf,] #delete infinite correlation coefficients
melt.matrix <- melt.matrix[!melt.matrix$value == 0,] #delete rows where correlation coefficients = 0
melt.matrix <- melt.matrix[!is.na(melt.matrix$value),]
# Sort melt.matrix based on the absolute values of the correlations
sorted_correlations <- melt.matrix[order(abs(melt.matrix$value), decreasing = TRUE), ]
# Select the top 100 correlations
top_100 <- sorted_correlations[1:100, ]

# ========== 2. PLOT MATRIX USING GGPLOT ======== ##

ggplot(top_100, aes(x = Var1, y = Var2)) +
  geom_tile(colour="grey20", aes(fill=value), size=0.2) + 
  scale_fill_gradient2(name = "Pearson's rho", low="navyblue", high="red", midpoint=0, limits=c(-1,1)) +
  labs(x="",y="") +
  ggtitle("\nCorrelations\nDietary polyphenols and microbial metabolites\n") +
  geom_text(aes(label=stars.adj), position=position_nudge(y=0.2), color="black", size=3) + 
  geom_text(aes(label=sprintf("%1.2f", value)), position=position_nudge(y=-0.1), 
            size=2, colour="black") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle = -45, hjust = 0, size=6)) + 
  theme(axis.text.y=element_text(angle = 0, hjust = 1, size=6)) +
  theme(plot.title = element_text(hjust = 0.5))

ggscatter(intfec_polyphenol, x = "int_alpha_catechin", y = "fec_naringenin", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "dietary alpha-catechin", ylab = "naringenin")

##=================================================== SERUM ========================================

##========================== CORRELATION CARBOHYDRATE/FIBER METABOLITES =============================

#merge intake + seral metabolites
intser <- merge(filtered_intake, filtered_serum, by = "row.names", all.x = TRUE)
rownames(intser) <- intser$Row.names
intser$Row.names <- NULL

#RANK transformation
intser <- intser %>% mutate_all(~ (length(.) + 1) - rank(.))

#filter full metabolite set on carbohydrate/fiber + serum metabolites
# Find column names containing 'fiber'
fiber_cols <- grep("fiber", colnames(intser), value = TRUE)
# Find column names containing 'carbohydrate'
carb_cols <- grep("carbohydrate", colnames(intser), value = TRUE)

diet_mtb <- c(fiber_cols, carb_cols)

# Serum metabolites
# Find column names containing 'propionate' (SCFA)
propionate_cols <- grep(c("meta_394", 'meta_433') , colnames(intser), value = TRUE)
# Find column names containing 'acetate' (SCFA)
acetate_cols <- grep(c("meta_316", 'meta_366', 'meta_374', 'meta_399', 
                       'meta_443', 'meta_465', 'meta_576', 'meta_584',
                       'meta_587', 'meta_697', 'meta_721', 'meta_798',
                       'meta_835', 'meta_927', 'meta_169', 'meta_287'), colnames(intser), value = TRUE)
# Find column names containing 'butyrate' (SCFA)
butyrate_cols <- grep(c('meta_866', 'meta_299'), colnames(intser), value = TRUE)
# Find column names containing 'succinate' (carboxylic acid)
succinate_cols <- grep("meta_421", colnames(intser), value = TRUE)

microbe_mtb <- c(propionate_cols, acetate_cols, butyrate_cols, succinate_cols)

# Combine diet and microbe cols
selected_cols <- c(diet_mtb, microbe_mtb)

# new df with associated diet and serum metabolites
intser_carbfiber <- intser[,selected_cols]

#create correlation matrix
cor = rcorr(as.matrix(intser_carbfiber), type = "pearson")

# create matrices for R, P and P.adj
cormatrix <- cor$r
pmatrix <- cor$P

# Subset matrix to only show correlations between intake and serum metabolites
subset_rows <- grepl("^ser_", rownames(cormatrix))
subset_columns <- grepl("^int_", colnames(cormatrix))
subsetted_matrix <- cormatrix[subset_rows, subset_columns]

# Subset p-value to only show p-values for correlations between intake and serym metabolites
subset_rows <- grepl("^ser_", rownames(pmatrix))
subset_columns <- grepl("^int_", colnames(pmatrix))
subset_p_matrix <- pmatrix[subset_rows, subset_columns]

melt.matrix <- melt(subsetted_matrix)
melt.matrix$p.value <- round(melt(subset_p_matrix)$value, 3) # add p-values to plot dataset
melt.matrix$p.value.adj <- p.adjust(melt.matrix$p.value, method = c('fdr')) # adjust P-values for multiple comparison
melt.matrix$stars <- cut(melt.matrix$p.value, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))  # Create column of significance labels
melt.matrix$stars.adj <- cut(melt.matrix$p.value.adj, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))

melt.matrix <- melt.matrix[!melt.matrix$value == -Inf,] #delete infinite correlation coefficients
melt.matrix <- melt.matrix[!melt.matrix$value == 0,] #delete rows where correlation coefficients = 0

# Sort melt.matrix based on the absolute values of the correlations
sorted_correlations <- melt.matrix[order(abs(melt.matrix$value), decreasing = TRUE), ]
# Select the top 100 correlations
top_100 <- sorted_correlations[1:32, ]

# ========== 2. PLOT MATRIX USING GGPLOT ======== ##

ggplot(top_100, aes(x = Var1, y = Var2)) +
  geom_tile(colour="grey20", aes(fill=value), size=0.2) + 
  scale_fill_gradient2(name = "Pearson's rho", low="navyblue", high="red", midpoint=0, limits=c(-1,1)) +
  labs(x="",y="") +
  ggtitle("\nCorrelations\nDietary carbohydrate/fiber and serum metabolites\n") +
  geom_text(aes(label=stars.adj), position=position_nudge(y=0.15), color="black", size=3) + 
  geom_text(aes(label=sprintf("%1.2f", value)), position=position_nudge(y=-0.1), 
            size=2, colour="black") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle = -45, hjust = 0, size=8)) + 
  theme(axis.text.y=element_text(angle = 0, hjust = 1, size=8)) +
  theme(plot.title = element_text(hjust = 0.5))

ggscatter(intser_carbfiber, x = "int_carbohydrates", y = "ser_meta_421", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "dietary carbohydrates", ylab = "(R)-2-benzylsuccinate")

##========================== CORRELATION AMINO ACID METABOLITES =============================

#merge intake + seral metabolites
intser <- merge(filtered_intake, filtered_serum, by = "row.names", all.x = TRUE)
rownames(intser) <- intser$Row.names
intser$Row.names <- NULL

#RANK transformation
intser <- intser %>% mutate_all(~ (length(.) + 1) - rank(.))

#filter full metabolite set on amino acids + seral metabolites
# Find column names containing 'protein'
protein_cols <- grep("protein", colnames(intser), value = TRUE)
# Find column names containing 'alanine'
alanine_cols <- grep("alanine", colnames(intser), value = TRUE)
# Find column names containing 'arginine'
arginine_cols <- grep("arginine", colnames(intser), value = TRUE)
# Find column names containing 'asparagine'
asparagine_cols <- grep("asparagine", colnames(intser), value = TRUE)
# Find column names containing 'aspartic'
aspartic_cols <- grep("aspartic", colnames(intser), value = TRUE)
# Find column names containing 'cysteine'
cysteine_cols <- grep("cysteine", colnames(intser), value = TRUE)
# Find column names containing 'glutamic'
glutamic_cols <- grep("glutamic", colnames(intser), value = TRUE)
# Find column names containing 'glutamine'
glutamine_cols <- grep("glutamine", colnames(intser), value = TRUE)
# Find column names containing 'glycine'
glycine_cols <- grep("glycine", colnames(intser), value = TRUE)
# Find column names containing 'histidine'
histidine_cols <- grep("histidine", colnames(intser), value = TRUE)
# Find column names containing 'isoleucine'
isoleucine_cols <- grep("isoleucine", colnames(intser), value = TRUE)
# Find column names containing 'leucine'
leucine_cols <- grep("leucine", colnames(intser), value = TRUE)
# Find column names containing 'lysine'
lysine_cols <- grep("lysine", colnames(intser), value = TRUE)
# Find column names containing 'methionine'
methionine_cols <- grep("methionine", colnames(intser), value = TRUE)
# Find column names containing 'phenylalanine'
phenylalanine_cols <- grep("phenylalanine", colnames(intser), value = TRUE)
# Find column names containing 'proline'
proline_cols <- grep("proline", colnames(intser), value = TRUE)
# Find column names containing 'serine'
serine_cols <- grep("serine", colnames(intser), value = TRUE)
# Find column names containing 'threonine'
threonine_cols <- grep("threonine", colnames(intser), value = TRUE)
# Find column names containing 'tryptophan'
tryptophan_cols <- grep("tryptophan", colnames(intser), value = TRUE)
# Find column names containing 'tyrosine'
tyrosine_cols <- grep("tyrosine", colnames(intser), value = TRUE)
# Find column names containing 'alanine'
valine_cols <- grep("valine", colnames(intser), value = TRUE)
# Find column names containing 'selenocysteine'
selenocysteine_cols <- grep("selenocysteine", colnames(intser), value = TRUE)

diet_mtb <- c(protein_cols, alanine_cols, arginine_cols, asparagine_cols, aspartic_cols, 
              cysteine_cols, glutamine_cols, glutamic_cols, glycine_cols, histidine_cols, 
              isoleucine_cols, leucine_cols, lysine_cols, methionine_cols, phenylalanine_cols, 
              proline_cols, serine_cols, threonine_cols, tryptophan_cols, tyrosine_cols,
              valine_cols, selenocysteine_cols)

# Serum metabolites
# Find column names containing 'valine' (BCAA)
valine_cols <- grep(c('meta_613', 'meta_406'), colnames(intser), value = TRUE)
# Find column names containing 'isoleucine' (BCAA)
isoleucine_cols <- grep(c('meta_540', 'meta_544', 'meta_628', 'meta_134', 'meta_397'), colnames(intser), value = TRUE)
# Find column names containing 'leucine' (BCAA)
leucine_cols <- grep(c('meta_549'), colnames(intser), value = TRUE)
# Find column names containing 'niacin' 
niacin_cols <- c('ser_meta_105')
# Find column names containing 'nicotinamide'
nicotinamide_cols <- grep(c('meta_155', 'meta_158'), colnames(intser), value = TRUE)
# Find column names containing 'dimethylglycine'
dimethylglycine_cols <- c('ser_meta_55')
# Find column names containing 'acetylglycine'
acetylglycine_cols <- c('ser_meta_90', 'ser_meta_179')

microbe_mtb <- c(valine_cols, isoleucine_cols, leucine_cols, niacin_cols, nicotinamide_cols, dimethylglycine_cols, acetylglycine_cols)


# Combine diet and microbe cols
selected_cols <- c(diet_mtb, microbe_mtb)

# new df with associated diet and serum metabolites
intser_aminoacid <- intser[,selected_cols]


#create correlation matrix
cor = rcorr(as.matrix(intser_aminoacid), type = "pearson")

# create matrices for R, P and P.adj
cormatrix <- cor$r
pmatrix <- cor$P

# Subset matrix to only show correlations between intake and seral metabolites
subset_rows <- grepl("^ser_", rownames(cormatrix))
subset_columns <- grepl("^int_", colnames(cormatrix))
subsetted_matrix <- cormatrix[subset_rows, subset_columns]

# Subset p-value to only show p-values for correlations between intake and seral metabolites
subset_rows <- grepl("^ser_", rownames(cormatrix))
subset_columns <- grepl("^int_", colnames(cormatrix))
subset_p_matrix <- pmatrix[subset_rows, subset_columns]

melt.matrix <- melt(subsetted_matrix)
melt.matrix$p.value <- round(melt(subset_p_matrix)$value, 3) # add p-values to plot dataset
melt.matrix$p.value.adj <- p.adjust(melt.matrix$p.value, method = c('fdr')) # adjust P-values for multiple comparison
melt.matrix$stars <- cut(melt.matrix$p.value, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))  # Create column of significance labels
melt.matrix$stars.adj <- cut(melt.matrix$p.value.adj, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))

melt.matrix <- melt.matrix[!melt.matrix$value == -Inf,] #delete infinite correlation coefficients
melt.matrix <- melt.matrix[!melt.matrix$value == 0,] #delete rows where correlation coefficients = 0

# Sort melt.matrix based on the absolute values of the correlations
sorted_correlations <- melt.matrix[order(abs(melt.matrix$value), decreasing = TRUE), ]
# Select the top 100 correlations
top_100 <- sorted_correlations[1:100, ]

# ========== 2. PLOT MATRIX USING GGPLOT ======== ##

ggplot(top_100, aes(x = Var1, y = Var2)) +
  geom_tile(colour="grey20", aes(fill=value), size=0.2) + 
  scale_fill_gradient2(name = "Pearson's rho", low="navyblue", high="red", midpoint=0, limits=c(-1,1)) +
  labs(x="",y="") +
  ggtitle("\nCorrelations\nDietary amino acids and serum metabolites\n") +
  geom_text(aes(label=stars.adj), position=position_nudge(y=0.15), color="black", size=3) + 
  geom_text(aes(label=sprintf("%1.2f", value)), position=position_nudge(y=-0.1), 
            size=2, colour="black") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle = -45, hjust = 0, size=8)) + 
  theme(axis.text.y=element_text(angle = 0, hjust = 1, size=8)) +
  theme(plot.title = element_text(hjust = 0.5))

ggscatter(intser_aminoacid, x = "int_protein", y = "ser_meta_105", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "dietary protein", ylab = "niacinamide")

##========================== CORRELATION TRYPTOPHAN METABOLITES =============================

#merge intake + serum metabolites
intser <- merge(filtered_intake, filtered_serum, by = "row.names", all.x = TRUE)
rownames(intser) <- intser$Row.names
intser$Row.names <- NULL

#RANK transformation
intser <- intser %>% mutate_all(~ (length(.) + 1) - rank(.))

#filter full metabolite set on tryptophan + serum metabolites
# Find column names containing 'tryptophan'
tryptophan_cols <- grep("tryptophan", colnames(intser), value = TRUE)

diet_mtb <- c(tryptophan_cols)

# Microbial metabolites
# Find column names containing 'indole-3-propionic acid' (I3PA) 
indoleprop_cols <- c('ser_meta_357')
# Find column names containing 'Indole-3-aldehyde' (IAId) 
indoleald_cols <- c('ser_meta_191')
# Find column names containing 'indoxylsulfate' (IS) 
indoxylsulfate_cols <- c('ser_meta_444')
# Find column names containing 'tryptamine'
tryptamine_cols <- c('ser_meta_454', 'ser_meta_251')


microbe_mtb <- c(indoleprop_cols, indoleald_cols, indoxylsulfate_cols, tryptamine_cols)

# Combine diet and microbe cols
selected_cols <- c(diet_mtb, microbe_mtb)

# new df with associated diet and microbe metabolites
intser_tryptophan <- intser[,selected_cols]


#create correlation matrix
cor = rcorr(as.matrix(intser_tryptophan), type = "pearson")

# create matrices for R, P and P.adj
cormatrix <- cor$r
pmatrix <- cor$P

# Subset matrix to only show correlations between intake and serum  metabolites
subset_rows <- grepl("^ser_", rownames(cormatrix))
subset_columns <- grepl("^int_", colnames(cormatrix))
subsetted_matrix <- cormatrix[subset_rows, subset_columns]

# Subset p-value to only show p-values for correlations between intake and serum metabolites
subset_rows <- grepl("^ser_", rownames(pmatrix))
subset_columns <- grepl("^int_", colnames(pmatrix))
subset_p_matrix <- pmatrix[subset_rows, subset_columns]

melt.matrix <- melt(subsetted_matrix)
melt.matrix$p.value <- round(melt(subset_p_matrix)$value, 3) # add p-values to plot dataset
melt.matrix$p.value.adj <- p.adjust(melt.matrix$p.value, method = c('fdr')) # adjust P-values for multiple comparison
melt.matrix$stars <- cut(melt.matrix$p.value, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))  # Create column of significance labels
melt.matrix$stars.adj <- cut(melt.matrix$p.value.adj, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))

melt.matrix <- melt.matrix[!melt.matrix$value == -Inf,] #delete infinite correlation coefficients
melt.matrix <- melt.matrix[!melt.matrix$value == 0,] #delete rows where correlation coefficients = 0
melt.matrix <- melt.matrix[!is.na(melt.matrix$value),]
# Sort melt.matrix based on the absolute values of the correlations
sorted_correlations <- melt.matrix[order(abs(melt.matrix$value), decreasing = TRUE), ]
# Select the top 100 correlations
top_100 <- sorted_correlations[1:15, ]


# ========== 2. PLOT MATRIX USING GGPLOT ======== ##

ggplot(top_100, aes(x = Var1, y = Var2)) +
  geom_tile(colour="grey20", aes(fill=value), size=0.2) + 
  scale_fill_gradient2(name = "Pearson's rho", low="navyblue", high="red", midpoint=0, limits=c(-1,1)) +
  labs(x="",y="") +
  ggtitle("\nCorrelations\nDietary tryptophan and serum metabolites\n") +
  geom_text(aes(label=stars.adj), position=position_nudge(y=0.15), color="black", size=3) + 
  geom_text(aes(label=sprintf("%1.2f", value)), position=position_nudge(y=-0.1), 
            size=2, colour="black") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle = -45, hjust = 0, size=8)) + 
  theme(axis.text.y=element_text(angle = 0, hjust = 1, size=8)) +
  theme(plot.title = element_text(hjust = 0.5))

ggscatter(intser_tryptophan, x = "int_tryptophan", y = "ser_meta_357", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "dietary tryptophan", ylab = "indole-3-propionic acid")

##========================== CORRELATION POLYPHENOL METABOLITES =============================

#merge intake + serum metabolites
intser <- merge(filtered_intake, filtered_serum, by = "row.names", all.x = TRUE)
rownames(intser) <- intser$Row.names
intser$Row.names <- NULL

#RANK transformation
intser <- intser %>% mutate_all(~ (length(.) + 1) - rank(.))

#filter full metabolite set on polyphenols + serum metabolites
# Find column names containing 'polyphenol'
polyphenol_cols <- grep("polyphenol", colnames(intser), value = TRUE)

#flavonoids
# Find column names containing 'epicatechin'
epicatechin_cols <- grep("epicatechin", colnames(intser), value = TRUE)
# Find column names containing 'catechin'
catechin_cols <- grep("catechin", colnames(intser), value = TRUE)
# Find column names containing 'epigallocatechin'
epigallocatechin_cols <- grep("epigallocatechin", colnames(intser), value = TRUE)
# Find column names containing 'hesperatin'
hesperatin_cols <- grep("hesperatin", colnames(intser), value = TRUE)
# Find column names containing 'naringenin'
naringenin_cols <- grep("naringenin", colnames(intser), value = TRUE)
# Find column names containing 'eriodictyol'
eriodictyol_cols <- grep("eriodictyol", colnames(intser), value = TRUE)
# Find column names containing 'apigenin'
apigenin_cols <- grep("apigenin", colnames(intser), value = TRUE)
# Find column names containing 'luteolin'
luteolin_cols <- grep("luteolin", colnames(intser), value = TRUE)
# Find column names containing 'tangeritin'
tangeritin_cols <- grep("tangeritin", colnames(intser), value = TRUE)
# Find column names containing 'chrysin'
chrysin_cols <- grep("chrysin", colnames(intser), value = TRUE)
# Find column names containing 'genistein'
genistein_cols <- grep("genistein", colnames(intser), value = TRUE)
# Find column names containing 'daidzein'
daidzein_cols <- grep("daidzein", colnames(intser), value = TRUE)
# Find column names containing 'kaempferol'
kaempferol_cols <- grep("kaempferol", colnames(intser), value = TRUE)
# Find column names containing 'myrestin'
myrestin_cols <- grep("myrestin", colnames(intser), value = TRUE)
# Find column names containing 'quercetin'
quercetin_cols <- grep("quercetin", colnames(intser), value = TRUE)
# Find column names containing 'cyanidin'
cyanidin_cols <- grep("cyanidin", colnames(intser), value = TRUE)
# Find column names containing 'delphinidin'
delphinidin_cols <- grep("delphinidin", colnames(intser), value = TRUE)
# Find column names containing 'malvedin'
malvedin_cols <- grep("malvedin", colnames(intser), value = TRUE)
# Find column names containing 'pelargonidin'
pelargonidin_cols <- grep("pelargonidin", colnames(intser), value = TRUE)

flavonoids <- c(polyphenol_cols, epicatechin_cols, catechin_cols, epigallocatechin_cols,
                hesperatin_cols, naringenin_cols, eriodictyol_cols, apigenin_cols, luteolin_cols,
                tangeritin_cols, chrysin_cols, genistein_cols, daidzein_cols, kaempferol_cols,
                myrestin_cols, quercetin_cols, cyanidin_cols, delphinidin_cols, malvedin_cols, 
                pelargonidin_cols)

#phenolic acids
# Find column names containing 'hydroxybenzoic'
hydroxybenzoic_cols <- grep("hydroxybenzoic", colnames(intser), value = TRUE)
# Find column names containing 'protocatechuic'
protocatechuic_cols <- grep("protocatechuic", colnames(intser), value = TRUE)
# Find column names containing 'gallic'
gallic_cols <- grep("gallic", colnames(intser), value = TRUE)
# Find column names containing 'vanillic'
vanillic_cols <- grep("vanillic", colnames(intser), value = TRUE)
# Find column names containing 'ellagic'
ellagic_cols <- grep("ellagic", colnames(intser), value = TRUE)
# Find column names containing 'salicyclic'
salicyclic_cols <- grep("salicyclic", colnames(intser), value = TRUE)
# Find column names containing 'caffeic'
caffeic_cols <- grep("caffeic", colnames(intser), value = TRUE)
# Find column names containing 'ferulic'
ferulic_cols <- grep("ferulic", colnames(intser), value = TRUE)
# Find column names containing 'sinapic'
sinapic_cols <- grep("sinapic", colnames(intser), value = TRUE)
# Find column names containing 'chlorogenic'
chlorogenic_cols <- grep("chlorogenic", colnames(intser), value = TRUE)
# Find column names containing 'coumaric'
coumaric_cols <- grep("coumaric", colnames(intser), value = TRUE)
# Find column names containing 'quinic'
quinic_cols <- grep("quinic", colnames(intser), value = TRUE)
# Find column names containing 'resveratrol'
resveratrol_cols <- grep("resveratrol", colnames(intser), value = TRUE)

phenols <- c(hydroxybenzoic_cols, protocatechuic_cols, gallic_cols, vanillic_cols,
             ellagic_cols, salicyclic_cols, caffeic_cols, ferulic_cols, sinapic_cols,
             chlorogenic_cols, coumaric_cols, quinic_cols, resveratrol_cols)

diet_mtb <- c(flavonoids, phenols)

# serum metabolites
# Find column names containing 'quercetin' 
quercetin_cols <- c('ser_meta_1161', 'ser_meta_1162')


microbe_mtb <- c(quercetin_cols)

# Combine diet and microbe cols
selected_cols <- c(diet_mtb, microbe_mtb)

# new df with associated diet and serum metabolites
intser_polyphenol <- intser[,selected_cols]

#create correlation matrix
cor = rcorr(as.matrix(intser_polyphenol), type = "pearson")

# create matrices for R, P and P.adj
cormatrix <- cor$r
pmatrix <- cor$P

# Subset matrix to only show correlations between intake and serum metabolites
subset_rows <- grepl("^ser_", rownames(cormatrix))
subset_columns <- grepl("^int_", colnames(cormatrix))
subsetted_matrix <- cormatrix[subset_rows, subset_columns]

# Subset p-value to only show p-values for correlations between intake and serum metabolites
subset_rows <- grepl("^ser_", rownames(pmatrix))
subset_columns <- grepl("^int_", colnames(pmatrix))
subset_p_matrix <- pmatrix[subset_rows, subset_columns]

melt.matrix <- melt(subsetted_matrix)
melt.matrix$p.value <- round(melt(subset_p_matrix)$value, 3) # add p-values to plot dataset
melt.matrix$p.value.adj <- p.adjust(melt.matrix$p.value, method = c('fdr')) # adjust P-values for multiple comparison
melt.matrix$stars <- cut(melt.matrix$p.value, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))  # Create column of significance labels
melt.matrix$stars.adj <- cut(melt.matrix$p.value.adj, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))

melt.matrix <- melt.matrix[!melt.matrix$value == -Inf,] #delete infinite correlation coefficients
melt.matrix <- melt.matrix[!melt.matrix$value == 0,] #delete rows where correlation coefficients = 0
melt.matrix <- melt.matrix[!is.na(melt.matrix$value),]
# Sort melt.matrix based on the absolute values of the correlations
sorted_correlations <- melt.matrix[order(abs(melt.matrix$value), decreasing = TRUE), ]
# Select the top 100 correlations
top_100 <- sorted_correlations[1:100, ]


# ========== 2. PLOT MATRIX USING GGPLOT ======== ##

ggplot(top_100, aes(x = Var1, y = Var2)) +
  geom_tile(colour="grey20", aes(fill=value), size=0.2) + 
  scale_fill_gradient2(name = "Pearson's rho", low="navyblue", high="red", midpoint=0, limits=c(-1,1)) +
  labs(x="",y="") +
  ggtitle("\nCorrelations\nDietary polyphenols and serum metabolites\n") +
  geom_text(aes(label=stars), position=position_nudge(y=0.2), color="black", size=3) + 
  geom_text(aes(label=sprintf("%1.2f", value)), position=position_nudge(y=-0.1), 
            size=2, colour="black") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle = -45, hjust = 0, size=6)) + 
  theme(axis.text.y=element_text(angle = 0, hjust = 1, size=5)) +
  theme(plot.title = element_text(hjust = 0.5))

ggscatter(intser_polyphenol, x = "int_protocatechuic_acid.1", y = "ser_meta_1161", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "dietary protocatechuic acid", ylab = "quercetin 3-(6-[4-glucosyl-p-coumaryl]glucosyl)(1->2)-rhamnoside")

