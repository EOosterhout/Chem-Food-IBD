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

#packages for linear regression/statistical testing
library(stats)
library(outliers)

#packages for correlation analysis
library(Hmisc)
library(reshape2)
library(RcmdrMisc)
library(GGally)
library(MetBrewer)
library(psych)
library(corrplot)

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

##===================================== DIFFERENTIAL ABUNDANCE ANALYSIS: INTAKE_flare >150 vs <150 (LOG TRANSFORMED AND FILTERED DATA) ====================================================

# Columns containing intake metabolites
intake_cols <- grep("^int_", names(intake_mtb), value = TRUE)
metabolites_intake <- intake_mtb[,colnames(intake_mtb) %in% intake_cols]

# Calculate the percentage of non-zero values for each variable
non_zero_pct <- apply(metabolites_intake != 0, 2, mean)
# Filter variables with at least a non-zero value in 20% of the data
filter_20pct <- metabolites_intake[,non_zero_pct >= 0.2]

#add pseudocount to all variables
pseudo <- filter_20pct + 1

# Create a new column specifying flare >150 (yes/no) for each sample
pseudo_flare <- cbind(intake_mtb$flare_above150, pseudo)
names(pseudo_flare)[1] <- 'flare_above150'

#pseudo_clean <- remove_outliers(pseudo_diagnosis, 'diagnosis')

##======================= WILCOXON TEST ======================##
wilcoxon_p <- c() # Initialize empty vector for p-values
# Do "for loop" over selected column names
for (i in 2:1010) {
  
  result <- wilcox.test(pseudo_flare[, i] ~ flare_above150,
                        data = pseudo_flare)
  
  # Stores p-value to the vector with this column name
  wilcoxon_p[[i]]  <- result$p.value
  
}

#store metabolites with raw p-value in new dataframe
wilcoxon_p <- data.frame(metabolites =  names(pseudo_flare[,2:1010]),
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
intake_log <- log2(pseudo_flare[,2:1010])
intake_log <- cbind(flare_above150 = pseudo_flare$flare_above150, intake_log)


#calculate the mean of each metabolite in >150 group
high_flare <- filter(intake_log, intake_log$flare_above150 == 'yes')
high_flare_m = apply(high_flare[,2:1010], 2, mean)

#calcuate the mean of each metabolite in <150 group
low_flare <- filter(intake_log, intake_log$flare_above150 == 'no')
low_flare_m = apply(low_flare[,2:1010], 2, mean)

#because the data is already log2 transformed, take the difference between the means.
foldchange <- low_flare_m - high_flare_m
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
highest6_full <- cbind(pseudo_flare$flare_above150, highest6_chem)
names(highest6_full)[1] <- 'flare_above150'
highest6_full <- na.omit(highest6_full)

# Puts plots in the same picture
gridExtra::grid.arrange(
  
  # Plot 1
  ggplot(highest6_full, aes(x = flare_above150, y = highest6_full[,2], fill = flare_above150)) + 
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
  ggplot(highest6_full, aes(x = flare_above150, y = highest6_full[,3], fill = flare_above150)) + 
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
  ggplot(highest6_full, aes(x = flare_above150, y = highest6_full[,4], fill = flare_above150)) + 
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
  ggplot(highest6_full, aes(x = flare_above150, y = highest6_full[,5], fill = flare_above150)) + 
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
  ggplot(highest6_full, aes(x = flare_above150, y = highest6_full[,6], fill = flare_above150)) + 
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
  ggplot(highest6_full, aes(x = flare_above150, y = highest6_full[,7], fill = flare_above150)) + 
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

##================================================ INTAKE METABOLITES BASED ON INTAKEDIFFERENCE FOUND IN DIAGNOSIS AND flare ====================================================##

#Metabolite names based on intakedifference == YES (only works when NOT running flare part of script)
intakedifference_diagnosis <- as.character(wilcoxon_p$metabolites[wilcoxon_p$intakedifference != "NO"])

#Metabolite names based on intakedifference == YES
intakedifference_flare <- as.character(wilcoxon_p$metabolites[wilcoxon_p$intakedifference != "NO"])

#New df containing only intake metabolites that show difference in intake in 
      ## IBD vs Non-IBD
      ## flare <150 vs flare >150
intake_mtb_1 <- intake_mtb[,colnames(intake_mtb) %in% intakedifference_diagnosis]
intake_mtb_2 <- intake_mtb_1[,colnames(intake_mtb_1) %in% intakedifference_flare]

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
sum(is.na(pseudo_clean)) #109332 NA

#LOG transformation or RANK transformation
pseudo_log <- log2(pseudo_clean)
pseudo_rank <- pseudo_clean %>% mutate_all(~ (length(.) + 1) - rank(.))

## Filtering of analysis table on metabolite type is performed in == LOAD DATA, CLEANING NAMES AND SUBSETTING == #

# Create a new column specifying IBD (yes/no) for each sample
pseudo_diagnosis <- cbind(diagnosis = intake_mtb$diagnosis, pseudo_clean)
pseudo_diagnosis <- cbind(diagnosis = intake_mtb$diagnosis, pseudo_log)
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

##============================= SIGNIFICANT COEFFICIENTS ==========================##

#Lactic acid
model_1 <- lm(lactic_acid ~ age + sex + BMI + diagnosis, data = pseudo_diagnosis)
summary(model_1)
plot_summs(model_1)

#Pentanoic acid
model_2 <- lm(pentanoic_acid ~ age + sex + BMI + diagnosis, data = pseudo_diagnosis)
summary(model_2)
plot_summs(model_2)

#Niacin
model_3 <- lm(niacin_total ~ age + sex + BMI + diagnosis, data = pseudo_diagnosis)
summary(model_3)
plot_summs(model_3)

#Oxalic acid
model_4 <- lm(oxalic_acid ~ age + sex + BMI + diagnosis, data = pseudo_diagnosis)
summary(model_4)
plot_summs(model_4)

#Caffeic acid
model_5 <- lm(caffeic_acid.1 ~ age + sex + BMI + diagnosis, data = pseudo_diagnosis)
summary(model_5)
plot_summs(model_5)

#Selenium
model_6 <- lm(selenium_se ~ age + sex + BMI + diagnosis, data = pseudo_diagnosis)
summary(model_6)
plot_summs(model_6)
##============================== VOLCANO PLOT OF LINREG_DIAGNOSIS ===============================##

ggplot(linreg_diagnosis, aes(x=RSquared, y=-1*log10(PValue), label=Intake_Metabolite)) + 
  geom_point() + 
  theme_minimal() +
  theme(legend.position = 'bottom') +
  scale_color_manual(values=c("#999999", "#009E73")) +
  geom_text_repel(size = 2) +
  scale_x_continuous(name = 'R Squared (diagnosis)')

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
sum(is.na(pseudo_clean)) #109332 NA

#LOG transformation or RANK transformation
pseudo_log <- log2(pseudo_clean)
pseudo_rank <- pseudo_clean %>% mutate_all(~ (length(.) + 1) - rank(.))

## Filtering of analysis table on metabolite type is performed in == LOAD DATA, CLEANING NAMES AND SUBSETTING == #

# Create a new column specifying high calprotectin (yes/no) for each sample
pseudo_calprotectin <- cbind(calprotectin_above150 = intake_mtb$calprotectin_above150, pseudo)
pseudo_calprotectin <- cbind(calprotectin_above150 = intake_mtb$calprotectin_above150, pseudo_log)
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

##============================= SIGNIFICANT COEFFICIENTS ==========================##

#Cryptoxanthin beta
model_1 <- lm(cryptoxanthin_beta ~ age + sex + BMI + calprotectin_above150, data = pseudo_calprotectin)
summary(model_1)
plot_summs(model_1)

##============================== VOLCANO PLOT OF LINREG_calprotectin ===============================##

ggplot(linreg_calprotectin, aes(x=RSquared, y=-1*log10(PValue), label=Intake_Metabolite)) + 
  geom_point() + 
  theme_minimal() +
  theme(legend.position = 'bottom') +
  scale_color_manual(values=c("#999999", "#009E73")) +
  geom_text_repel(size = 2) +
  scale_x_continuous(name = 'R Squared (calprotectin)')

##================================================= SIGNIFICANT INTAKE DIFFERENCE (flare) ==============================================##

# Sorts pvalue in increasing order. Takes 6 first ones. Takes those rows that match
# with foldchange. Takes metabolites. 
highest6 <- linreg_flare[linreg_flare$PValue %in% sort(linreg_flare$PValue, decreasing = F)[1:6], ]$Intake_Metabolite
# From intake table, takes only those metabolites that had highest foldchange
highest6_chem <- pseudo_flare[,colnames(pseudo_flare) %in% highest6]
# Adds colData that includes patient status information
highest6_full <- cbind(pseudo_flare$flare_above150, highest6_chem)
names(highest6_full)[1] <- 'flare'
highest6_full <- na.omit(highest6_full)

# Puts plots in the same picture
gridExtra::grid.arrange(
  
  # Plot 1
  ggplot(highest6_full, aes(x = flare, y = highest6_full[,2], fill = flare)) + 
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
  ggplot(highest6_full, aes(x = flare, y = highest6_full[,3], fill = flare)) + 
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
  ggplot(highest6_full, aes(x = flare, y = highest6_full[,4], fill = flare)) + 
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
  ggplot(highest6_full, aes(x = flare, y = highest6_full[,5], fill = flare)) + 
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
  ggplot(highest6_full, aes(x = flare, y = highest6_full[,6], fill = flare)) + 
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
  ggplot(highest6_full, aes(x = flare, y = highest6_full[,7], fill = flare)) + 
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
sum(is.na(pseudo_clean)) #109332 NA

#LOG transformation or RANK transformation
pseudo_log <- log2(pseudo_clean)
pseudo_rank <- pseudo_clean %>% mutate_all(~ (length(.) + 1) - rank(.))

## Filtering of analysis table on metabolite type is performed in == LOAD DATA, CLEANING NAMES AND SUBSETTING == #

# Create a new column specifying IBD (yes/no) for each sample
pseudo_flare <- cbind(before_a_flare = intake_mtb$before_a_flare, pseudo)
pseudo_flare <- cbind(before_a_flare = intake_mtb$before_a_flare, pseudo_log)
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
                          PValue = significant_pvalues,
                          RSquared = r_squared,
                          stringsAsFactors = FALSE)
    
    # Append results to the main dataframe
    results_df <- rbind(results_df, results)
  }
}

# Filter results df on significant coefficients 'diagnosis'
linreg_flare <- results_df[results_df$Coefficient == 'before_a_flareyes',]
linreg_flare$p_adjusted <- p.adjust(linreg_flare$PValue, method = "fdr") #add column with p_adj for multiple testing

##============================= SIGNIFICANT COEFFICIENTS ==========================##

#Cryptoxanthin beta
model_1 <- lm(daidzin ~ age + sex + BMI + calprotectin_above150, data = pseudo_calprotectin)
summary(model_1)
plot_summs(model_1)

##============================== VOLCANO PLOT OF LINREG_DIAGNOSIS ===============================##

ggplot(linreg_flare, aes(x=RSquared, y=-1*log10(PValue), label=Intake_Metabolite)) + 
  geom_point() + 
  theme_minimal() +
  theme(legend.position = 'bottom') +
  scale_color_manual(values=c("#999999", "#009E73")) +
  geom_text_repel(size = 2) +
  scale_x_continuous(name = 'R Squared (flare)')

##================================================= SIGNIFICANT INTAKE DIFFERENCE (FLARE) ==============================================##

# Sorts pvalue in increasing order. Takes 6 first ones. Takes those rows that match
# with foldchange. Takes metabolites. 
highest6 <- linreg_flare[linreg_flare$p_adjusted %in% sort(linreg_flare$p_adjusted, decreasing = F)[1:6], ]$Intake_Metabolite
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

#Associated Diet Metabolite names based on flare
flare <- as.character(linreg_flare$Intake_Metabolite)

#Associated Diet Metabolite names based on flare
flare <- as.character(linreg_flare$Intake_Metabolite)

#New df containing only diet metabolites that are associated with clinical outcomes
## IBD vs Non-IBD
## flare <150 vs flare >150
## After/During Flare vs Before Flare
intake_mtb_diagnosis <- pseudo_diagnosis[,colnames(pseudo_diagnosis) %in% diagnosis]
intake_mtb_flare <- pseudo_flare[,colnames(pseudo_flare) %in% flare]
intake_mtb_flare <- pseudo_flare[,colnames(pseudo_flare) %in% flare]

# Get the column names from each data frame
col_names_diagnosis <- colnames(intake_mtb_diagnosis)
col_names_flare <- colnames(intake_mtb_flare)
col_names_flare <- colnames(intake_mtb_flare)

# Find the common column names
common_columns <- intersect(intersect(col_names_diagnosis, col_names_flare), col_names_flare) # 2 intake mtb: 
diagnosis_flare <- intersect(col_names_diagnosis, col_names_flare)
diagnosis_flare <- intersect(col_names_diagnosis, col_names_flare)
flare_flare <- intersect(col_names_flare, col_names_flare)

# Convert the metabolite lists to sets
set1 <- unique(diagnosis)
set2 <- unique(flare)
set3 <- unique(flare)

# Create the Venn diagram
venn.plot <- venn.diagram(
  x = list(set1, set2, set3),
  category.names = c("Diagnosis", "flare", "Flare"),
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
subsetted_matrix <- filtered_cor[subset_rows, subset_columns]

melt.matrix <- melt(subsetted_matrix)
melt.matrix <- melt.matrix[!melt.matrix$value == -Inf,] #delete infinite correlation coefficients
melt.matrix <- melt.matrix[!melt.matrix$value == 0,] #delete rows where correlation coefficients = 0

# Sort melt.matrix based on the absolute values of the correlations
sorted_correlations <- melt.matrix[order(abs(melt.matrix$value), decreasing = TRUE), ]
# Select the top 100 correlations
top_100 <- sorted_correlations[1:100, ]


ggscatter(intfec, x = "int_X35_dicaffeoylquinic_acid", y = "fec_cysteine.s.sulfate", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "3,5-dicaffeoylquinic acid", ylab = "cysteine-s-sulfate")

##==================================================== CORRELATION MATRIX: INTAKE VS FECAL METABOLITES ======================================================

#merge intake + fecal metabolites
intser <- merge(filtered_intake, filtered_serum, by = "row.names", all.x = TRUE)
rownames(intser) <- intser$Row.names
intser$Row.names <- NULL

#create correlation matrix
cor = rcorr(as.matrix(intser), type = "pearson")
cor$P.adj <- p.adjust(cor$P, method = c('fdr')) # adjust P-values for multiple comparison
cormatrix <- cor$r

# Subset matrix to only show correlations between intake and fecal metabolites
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

