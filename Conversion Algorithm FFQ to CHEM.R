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
library(ggbreak) 
library(patchwork)

# Load data
chemdata <- read_xlsx('transposed_ffqchem.xlsx')
ffqdata_raw <- read_xlsx('diet_raw_V2.xlsx')
ffqgroup_explanation <- read_xlsx('ffqgroups_explanation.xlsx')

##===================== EXPLORATORY PLOTS: PLOT CHEM PER FFQ ELEMENT ====================##

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
ggplot(subset(df, ffq_group == "fish_fatty"), aes(x = variable, y = value)) + # Per ffq group
  geom_bar(stat = "identity") +
  labs(title = subset(df, ffq_group == 'fish_fatty')$ffq_group) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

ggplot(subset(df, ffq_group == "fish_fatty"), aes(x = variable, y = value)) + # Per element
  geom_bar(stat = "identity", width = 0.5) + 
  labs(title = subset(df, ffq_group == 'fish_fatty')$ffq_group) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5)) + 
  scale_y_break(c(0.05, 0.1), scales = 'free')

# Plot chem_element distribution for each ffqgroup using a for loop
for (group in ffqgroup) { 
  Filename <- paste("plot_", group, ".png", sep="")
  myPlot <-ggplot(subset(dt, ffq_group == group), aes(x = variable, y = value)) +
    geom_bar(stat = "identity") +
    labs(title = subset(dt, ffq_group == group)$ffq_group) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(filename = Filename, myPlot, width=40, height=15)
}

for (chem in chem_elements) { 
  Filename <- paste("plot_", chem, ".png", sep="")
  myPlot <-ggplot(subset(df, variable == chem), aes(x = ffq_group, y = value)) +
    geom_bar(stat = "identity") +
    labs(title = subset(df, variable == chem)$variable) + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  ggsave(filename = Filename, myPlot, width=30, height=15)
}

##===================== FILTERING AND MERGING OF DIET DATA =================================##

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
portion <- portion[, !(colnames(portion) %in% ffqgroup_explanation$merge_item)]
portion <- portion[!is.na(portion$chicken),] # remove all NA values from dataframe
ffq_final <- as.character(colnames(portion))

##============================= TRANSPOSING DIET DATA AND MERGE WITH CHEM DATASET =======================##

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
transposed_portion <- bind_cols(col_ID)
transposed_portion <-transposed_portion %>% select(-contains('ffq_group'))
transposed_portion$ffq_group <- ffq_final
transposed_portion <- transposed_portion %>%
  select(ffq_group, everything())

# New df with chem data in mg/g
chemdata_mg_g <- chemdata[,2:1191]/100 #chemdata divided by 100 --> mg/g
chemdata_mg_g <- cbind(chemdata$ffq_group, chemdata_mg_g)
colnames(chemdata_mg_g)[1] <- 'ffq_group'

# Bind chem and diet data into one df
chem_diet_data <- cbind(chemdata_mg_g, transposed_portion)
write_xlsx(chem_diet_data, 'chem_diet_combined.xlsx')

##======================================= CONVERSION FFQ_GROUP/DAY TO CHEM/DAY =======================================##
