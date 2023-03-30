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

##============= PLOT CHEM PER FFQ ELEMENT ============##

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

 
