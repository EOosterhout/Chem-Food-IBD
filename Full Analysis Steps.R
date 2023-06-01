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

#packages for linear regression/statistical testing
library(stats)

#packages for correlation analysis
library(Hmisc)
library(reshape2)
library(RcmdrMisc)
library(GGally)
library(MetBrewer)
library(psych)
library(corrplot)

## Functions ##


##=========================================== LOAD DATA, CLEANING NAMES AND SUBSETTING ===============================

data_full <- read_xlsx('analysis_table.xlsx')
plausible_intake <- data_full %>% filter(between(SUMOFKCAL, 800, 5000)) # filter on plausible intake (800 < intake < 5000)

# Get metabolite data from analysis table using participant list from raw metabolite files
intake <- read_xlsx('chem_raw_participant_V2.xlsx')
fecal <- read_xlsx('fecal_mtb.xlsx')
serum <- read_xlsx("blood_mtb_measured.xlsx")

participants_intake <- as.character(intake$UMCGIBDResearchIDorLLDeepID)
participants_fecal <- as.character(fecal$UMCGIBDResearchIDorLLDeepID)
participants_serum <- as.character(serum$UMCGIBDResearchIDorLLDeepID)

intake_mtb <- plausible_intake[plausible_intake$UMCGIBDResearchIDorLLDeepID %in% participants_intake,]
fecal_mtb <- plausible_intake[plausible_intake$UMCGIBDResearchIDorLLDeepID %in% participants_fecal,]
serum_mtb <- plausible_intake[plausible_intake$UMCGIBDResearchIDorLLDeepID %in% participants_serum,]

#clean compound names, intake metabolites
names(intake_mtb) = gsub(pattern = ":", replacement = "_", x = names(intake_mtb))
names(intake_mtb) = gsub(pattern = " ", replacement = "_", x = names(intake_mtb))
names(intake_mtb) = gsub(pattern = "-", replacement = "_", x = names(intake_mtb))
names(intake_mtb) = gsub(pattern = ",", replacement = "", x = names(intake_mtb))
names(intake_mtb) = gsub(pattern = '"', replacement = "", x = names(intake_mtb))
names(intake_mtb) = gsub(pattern = "\\|.*", replacement = "", x = names(intake_mtb))
names(intake_mtb) = gsub(pattern = "\\(", replacement = "", x = names(intake_mtb))
names(intake_mtb) = gsub(pattern = "\\)", replacement = "", x = names(intake_mtb))
names(intake_mtb) = gsub(pattern = "\\+", replacement = "pos", x = names(intake_mtb))
names(intake_mtb) = gsub(pattern = "\\'", replacement = "", x = names(intake_mtb))
names(intake_mtb) = gsub(pattern = "\\â±", replacement = "", x = names(intake_mtb))


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

#pseudo_clean <- remove_outliers(pseudo_diagnosis, 'diagnosis')

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

##===================================== DIFFERENTIAL ABUNDANCE ANALYSIS: INTAKE_CALPROTECTIN >150 vs <150 (LOG TRANSFORMED AND FILTERED DATA) ====================================================

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

##================================================ INTAKE METABOLITES BASED ON INTAKEDIFFERENCE FOUND IN DIAGNOSIS AND CALPROTECTIN ====================================================##

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
pseudo_diagnosis <- cbind(diagnosis = intake_mtb$diagnosis, pseudo)

# Add covariates to df
pseudo_diagnosis <- cbind(age = intake_mtb$age, pseudo_diagnosis)
pseudo_diagnosis <- cbind(sex = intake_mtb$sex, pseudo_diagnosis)
pseudo_diagnosis <- cbind(BMI = intake_mtb$BMI, pseudo_diagnosis)

#============LINEAR REGRESSION ==============#
#dependent variable: individual diet metabolites
#predictor variables: covariates(age, sex, BMI), IBD(yes/no)

#columns with intake data + predictor variables
metabolite_names <- names(pseudo_diagnosis[,5:1013])
metabolite_names
predictor_vars <- c('age', 'sex', 'BMI', 'diagnosis')

# Create an empty dataframe to store results
results_df <- data.frame(Intake_Metabolite = character(),
                         Coefficient = numeric(),
                         PValue = numeric(),
                         RSquared = numeric(),
                         stringsAsFactors = FALSE)

# Perform linear regression over intake_metabolites using a for loop
for (dep_var in metabolite_names) {
  # Create formula
  formula <- paste(dep_var, paste(predictor_vars, collapse = " + "), sep = " ~ ")
  
  tryCatch({
    # Perform linear regression
    regression_model <- lm(formula, data = pseudo_diagnosis)
    
    # Extract significant coefficients and p-values
    coefficients <- coef(regression_model)
    p_values <- summary(regression_model)$coefficients[, 4]
    
    # Extract R-squared value
    r_squared <- summary(regression_model)$r.squared
    
    # Filter significant coefficients (p-value < 0.05)
    significant_coeffs <- coefficients[p_values < 0.05]
    
    # Check if there are any significant coefficients
    if (length(significant_coeffs) > 0) {
      # Create a dataframe for the results
      results <- data.frame(Intake_Metabolite = dep_var,
                            Coefficient = names(significant_coeffs),
                            PValue = significant_coeffs,
                            RSquared = r_squared,
                            stringsAsFactors = FALSE)
      
      # Append results to the main dataframe
      results_df <- rbind(results_df, results)
    }
  }, error = function(e) {
    # Skip the dependent variable and continue to the next iteration
    message(paste("Skipping", dep_var, "due to an error:", conditionMessage(e)))
  })
}

# Filter results df on significant coefficients 'diagnosis'
linreg_diagnosis <- results_df[results_df$Coefficient == 'diagnosisIBD',]
linreg_diagnosis$p_adjusted <- p.adjust(linreg_diagnosis$PValue, method = "fdr") #add column with p_adj for multiple testing

##========================================== LINEAR REGRESSION: >150 CALPROTECTIN vs <150 CALPROTECTIN ========================================

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
pseudo_calprotectin <- cbind(calprotectin_above150 = intake_mtb$calprotectin_above150, pseudo)

# Add covariates to df
pseudo_calprotectin <- cbind(age = intake_mtb$age, pseudo_calprotectin)
pseudo_calprotectin <- cbind(sex = intake_mtb$sex, pseudo_calprotectin)
pseudo_calprotectin <- cbind(BMI = intake_mtb$BMI, pseudo_calprotectin)

#============LINEAR REGRESSION ==============#
#dependent variable: individual diet metabolites
#predictor variables: covariates(age, sex, BMI), calprotectin >150 (yes/no)

#columns with intake data + predictor variables
metabolite_names <- names(pseudo_calprotectin[,5:1013])
metabolite_names
predictor_vars <- c('age', 'sex', 'BMI', 'calprotectin_above150')

# Create an empty dataframe to store results
results_df <- data.frame(Intake_Metabolite = character(),
                         Coefficient = numeric(),
                         PValue = numeric(),
                         RSquared = numeric(),
                         stringsAsFactors = FALSE)

# Perform linear regression over intake_metabolites using a for loop
for (dep_var in metabolite_names) {
  # Create formula
  formula <- paste(dep_var, paste(predictor_vars, collapse = " + "), sep = " ~ ")
  
  tryCatch({
    # Perform linear regression
    regression_model <- lm(formula, data = pseudo_calprotectin)
    
    # Extract significant coefficients and p-values
    coefficients <- coef(regression_model)
    p_values <- summary(regression_model)$coefficients[, 4]
    
    # Extract R-squared value
    r_squared <- summary(regression_model)$r.squared
    
    # Filter significant coefficients (p-value < 0.05)
    significant_coeffs <- coefficients[p_values < 0.05]
    
    # Check if there are any significant coefficients
    if (length(significant_coeffs) > 0) {
      # Create a dataframe for the results
      results <- data.frame(Intake_Metabolite = dep_var,
                            Coefficient = names(significant_coeffs),
                            PValue = significant_coeffs,
                            RSquared = r_squared,
                            stringsAsFactors = FALSE)
      
      # Append results to the main dataframe
      results_df <- rbind(results_df, results)
    }
  }, error = function(e) {
    # Skip the dependent variable and continue to the next iteration
    message(paste("Skipping", dep_var, "due to an error:", conditionMessage(e)))
  })
}

# Filter results df on significant coefficients 'calprotectin_above150'
linreg_calprotectin <- results_df[results_df$Coefficient == 'calprotectin_above150yes',]
linreg_calprotectin$p_adjusted <- p.adjust(linreg_calprotectin$PValue, method = "fdr") #add column with p_adj for multiple testing

##========================================== LINEAR REGRESSION: BEFORE A FLARE vs DURING/AFTER A FLARE ========================================

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
pseudo_flare <- cbind(before_a_flare = intake_mtb$before_a_flare, pseudo)

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
for (dep_var in metabolite_names) {
  # Create formula
  formula <- paste(dep_var, paste(predictor_vars, collapse = " + "), sep = " ~ ")
  
  tryCatch({
    # Perform linear regression
    regression_model <- lm(formula, data = pseudo_flare)
    
    # Extract significant coefficients and p-values
    coefficients <- coef(regression_model)
    p_values <- summary(regression_model)$coefficients[, 4]
    
    # Extract R-squared value
    r_squared <- summary(regression_model)$r.squared
    
    # Filter significant coefficients (p-value < 0.05)
    significant_coeffs <- coefficients[p_values < 0.05]
    
    # Check if there are any significant coefficients
    if (length(significant_coeffs) > 0) {
      # Create a dataframe for the results
      results <- data.frame(Intake_Metabolite = dep_var,
                            Coefficient = names(significant_coeffs),
                            PValue = significant_coeffs,
                            RSquared = r_squared,
                            stringsAsFactors = FALSE)
      
      # Append results to the main dataframe
      results_df <- rbind(results_df, results)
    }
  }, error = function(e) {
    # Skip the dependent variable and continue to the next iteration
    message(paste("Skipping", dep_var, "due to an error:", conditionMessage(e)))
  })
}

# Filter results df on significant coefficients 'diagnosis'
linreg_flare <- results_df[results_df$Coefficient == 'before_a_flareyes',]
linreg_flare$p_adjusted <- p.adjust(linreg_flare$PValue, method = "fdr") #add column with p_adj for multiple testing
