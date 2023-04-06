##Title: Conversion Algorithm FFQ/day to CHEM/day
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
chem_diet <- read_xlsx('chem_diet_combined.xlsx')

##======================================= CONVERSION FFQ_GROUP/DAY TO CHEM/DAY =======================================##
chemdata <- as.character(colnames(chem_diet[,2:1191]))

# Calculation of chemical content per day per ffq_group
diet_conv <- list()
for (x in chemdata) {
  output <- chem_diet[,x]*chem_diet$LLDeep_0688
  diet_conv[[x]] <- output
  diet_conv[[x]] <- cbind(chem_diet$ffq_group, diet_conv[[x]])
}

# Chem intake per ffq element for one participant
LLDeep_0688 <- Reduce(full_join, diet_conv) # Transforming list to dataframe, keeping only one column containing the ffqgroups
colnames(LLDeep_0688) [1] <- 'ffq_group'


# Sum intake per food compound --> total intake of food compound/day
participant <- as.data.frame(colSums(LLDeep_0688[,2:1191], na.rm = TRUE))
colnames(participant)[1] <- 'value'

##===================== COMPUTE PERCENTAGES (RELATIVE), MELT DF'S ============##

# Create new df with percentages of components
sum_row_chemdata <- rowSums(LLDeep_0688[ ,2:1191], na.rm = TRUE)
sum_row_chemdata <- as.data.frame(sum_row_chemdata)
colnames(sum_row_chemdata)[1] <- 'row_sum'

# Calculate percentages per chem_element and add to percentage df
chemdata <- as.character(colnames(chem_diet[,2:1191]))
ffqgroup <- as.character(chem_diet$ffq_group)

for (x in chemdata) {
  y <- (LLDeep_0688[,2:1191]/sum_row_chemdata$row_sum)*100
}

percentage_chemdata <- cbind(chem_diet$ffq_group, y) # new df with percentages of chem_elements per ffq_group
colnames(percentage_chemdata)[1] <- 'ffq_group'

#the function melt reshapes df from wide to long, needed for plotting data
d_relative <- melt(percentage_chemdata, id.vars = 'ffq_group') #relative
d_relative <- d_relative[complete.cases(d_relative),]
d_absolute <- melt(LLDeep_0688, id.vars = 'ffq_group') #absolute
d_absolute <- d_absolute[complete.cases(d_absolute),]

collapse <- function(df){
  df %>% 
    group_by(variable) %>%
    summarise(value = sum(value, na.rm = TRUE),
              ffq_group = NULL)
}

participant_intake <- collapse(d_absolute)
##=========================================== EXPLORATORY PLOTS ===================================##

# Select data to plot and perform log transformation
variable_lcfa <- participant_intake[participant_intake$variable %like% ":", ]
dat <- subset(variable_lcfa, value != 0)
dat$log_abs_value <- log10(dat$value + 1)

# Select only intake of 50 most consumed food compounds
dat_high <- dat[with(dat,order(-log_abs_value)),]
dat_high <- dat_high[1:50,]

# Plot total intake of food compound per day for one participant
dat_high %>% 
  ggplot(aes(x = reorder(variable, -value), y = log_abs_value)) + # Per compound
  geom_bar(stat = "identity", fill = 'steelblue') +
  coord_flip() +
  labs(x = '', y = 'log(absolute content)') +
  theme_minimal() + theme(axis.text.x = element_text(hjust = 1, size = 7))


# Plot all intake at once, stacked content(lcfa)
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
        legend.position = 'bottom',
        axis.text.x = element_text(vjust=0.5, size = 10)) 
