##Title: Exploratory plots chem_data
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

##===================== COMPUTE PERCENTAGES (RELATIVE), MELT DF'S ============##

# Create new df with percentages of components
sum_row_chemdata <- rowSums(chemdata[ ,2:1191], na.rm = TRUE)
sum_row_chemdata <- as.data.frame(sum_row_chemdata)
colnames(sum_row_chemdata)[1] <- 'row_sum'

# Calculate percentages per chem_element and add to percentage df
chem_elements <- as.character(colnames(chemdata[ ,2:1191]))
ffqgroup <- as.character(chemdata$ffq_group)
for (x in chem_elements) {
  y <- (chemdata[,2:1191]/sum_row_chemdata$row_sum)*100
}
percentage_chemdata <- cbind(chemdata$ffq_group, y) # new df with percentages of chem_elements per ffq_group
colnames(percentage_chemdata)[1] <- 'ffq_group'

#the function melt reshapes df from wide to long
d_relative <- melt(percentage_chemdata, id.vars = 'ffq_group') #relative
d_relative <- d_relative[complete.cases(d_relative),]
d_absolute <- melt(chemdata, id.vars = 'ffq_group') #absolute
d_absolute <- d_absolute[complete.cases(d_absolute),]

##========================================= EXPLORATORY PLOTS ========================================##

# Select data to plot and perform log transformation
dat <- subset(d_absolute, ffq_group == "fish_fatty" & value != 0) # select data to use for graph
dat$log_abs_value <- log10(dat$value + 1)

# Plot per ffq_group
dat %>% 
ggplot(aes(x = reorder(variable,-value), y = log_abs_value), reorder(variable, -log_abs_value, median)) + # Per ffq group
  geom_bar(stat = "identity", fill = '#FF6666') +
  labs(title = dat$ffq_group, x = 'chemical element', y = 'log(absolute content)' ) +
  theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

# Plot all groups at once, stacked content(lcfa)
# Stacked + percent
variable_lcfa <- d_relative[d_relative$variable %like% ":", ] # selecting only lcfa content
dat <- subset(variable_lcfa, value != 0)

ggplot(dat, aes(fill= variable, y=value, x=ffq_group)) + 
geom_bar(position= position_fill(reverse = TRUE), stat="identity", color = 'white') +
labs(x = '', y = 'relative content' ) +
coord_flip() +
theme_minimal() +
theme(legend.key.size = unit(0.3, 'cm'), #change legend key size
      legend.key.height = unit(0.3, 'cm'), #change legend key height
      legend.key.width = unit(0.3, 'cm'), #change legend key width
      legend.title = element_text(size=10), #change legend title font size
      legend.text = element_text(size=8),
      legend.position = 'none',
      axis.text.x = element_text(vjust=0.5, size = 10)) 

##============================================== LOOP FOR MAKING PLOTS ==========================##

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