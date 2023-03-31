##Title: Conversion Algorithm and exploratory plots
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

# Load data
chemdata <- read_xlsx('transposed_ffqchem.xlsx')
ffqdata_raw <- read_xlsx('diet_raw_V2.xlsx')
ffqgroup_explanation <- read_xlsx('ffqgroups_full.xlsx')

##===================== FILTERING AND MERGING OF DIET DATA =====================================##

#Set all ffqgroup names as character
ffqgroups <- as.character(colnames(ffqdata_raw[, 24:133]))
nodata <- as.character(ffqgroup_explanation$no_chem_data) #ffqgroups that don't have chem data
mergeitem <- as.character(ffqgroup_explanation$merge_item) # merged groups
# new df containing only the portion info from diet data
portion <- ffqdata_raw[, c(ffqgroups)]
portion <- portion[, !(colnames(portion) %in% nodata)]

# new df containing the elements of merged groups
merge_groups <- ffqgroup_explanation[, 4:9]

# Get portion data per merge_group
yogurt_merge <- portion[, (colnames(portion) %in% merge_groups$yogurt_merge)]
coffeecreamer_merge <- portion[, (colnames(portion) %in% merge_groups$coffeecreamer_merge)]
meat_merge <- portion[, (colnames(portion) %in% merge_groups$meat_merge)]
pork_merge <- portion[, (colnames(portion) %in% merge_groups$pork_merge)]
vegetables_fat <- portion[, (colnames(portion) %in% merge_groups$vegetables_fat)]
vegetables_nofat <- portion[, (colnames(portion) %in% merge_groups$vegetables_nofat)]

# Add median column per merged group to portion data
portion$yogurt_merge <- apply(yogurt_merge, 1, median) # 1 means computing per row, 2 means computing per column
portion$coffeecreamer_merge <- apply(coffeecreamer_merge, 1, median)
portion$meat_merge <- apply(meat_merge, 1, median)
portion$pork_merge <- apply(pork_merge, 1, median)
portion$vegetables_fat <- apply(vegetables_fat, 1, median)
portion$vegetables_nofat <- apply(vegetables_nofat, 1, median)

# Remove the columns which are merged
portion <- portion[, !(colnames(portion) %in% mergeitem)]

##======================================= CONVERSION DIET DATA TO CHEM DATA =======================##

# New df with filtered/merged portion data
ffqdata <- cbind(ffqdata_raw[, 1:6], portion)
ffqdata <- cbind(ffqdata, ffqdata_raw[, 161:180])
ffqdata <- cbind(ffqdata, ffqdata_raw[, 7:23])

##===================== EXPLORATORY PLOTS: PLOT CHEM PER FFQ ELEMENT ==============================##

# Create new df with percentages of components
sum_row_chemdata <- rowSums(chemdata[ ,2:1191], na.rm = TRUE)
sum_row_chemdata <- as.data.frame(sum_row_chemdata)
colnames(sum_row_chemdata)[1] <- 'row_sum'

# Calculate percentages per chem_element and add to percentage df
chem_elements <- as.character(colnames(chemdata[ ,2:1191]))
ffqgroup <- as.character(chemdata$ffq_group)
for (x in chem_elements) {
 y <- chemdata[,2:1191]/sum_row_chemdata$row_sum
}
percentage_chemdata <- cbind(chemdata$ffq_group, y) # new df with percentages of chem_elements per ffq_group
colnames(percentage_chemdata)[1] <- 'ffq_group'

#the function melt reshapes df from wide to long
df <- melt(percentage_chemdata, id.vars = 'ffq_group') #percentage
df <- df[complete.cases(df),]
dt <- melt(chemdata, id.vars = 'ffq_group') #quantity
dt <- dt[complete.cases(dt),]

# Test a function on one piece to develop graph
ggplot(subset(df, ffq_group == "fish_fatty"), aes(x = variable, y = value)) +
    geom_bar(stat = "identity") +
    labs(title = subset(df, ffq_group == 'fish_fatty')$ffq_group) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


# Plot chem_element distribution for each ffqgroup using a for loop
for (group in ffqgroup) { 
  Filename <- paste("plot_", group, ".png", sep="")
  myPlot <-ggplot(subset(dt, ffq_group == group), aes(x = variable, y = value)) +
  geom_bar(stat = "identity") +
  labs(title = subset(dt, ffq_group == group)$ffq_group) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(filename = Filename, myPlot, width=40, height=15)
}


