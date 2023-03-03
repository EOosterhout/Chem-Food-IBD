##Title: Request compound information in FoodB
##Author: Emily Oosterhout
##

#Import libraries
library(readxl)
library(dplyr)
library(tidyverse)

#Load data
compounds <- read_xlsx("Compound.xlsx")
food <- read.csv('Food.csv')
content <- read.csv('Content.csv')
id <- read_xlsx('idlist.xlsx')

#Filter for foods asked in FFQ and add food_id column to idlist
ffq <- food[match(id$public_id, food$public_id), ]
id <- cbind(id, ffq$id) 
colnames(id)[5] <- 'food_id'
write.csv(id, file = "IDffq.csv")

#Filter food item compounds
FFQ_ID <- id$food_id

ffq_comp <- filter(compounds, compounds$food_id == FFQ_ID)
ffq_comp <- as.data.frame(ffq_comp)
ffq_comp <- subset(ffq_comp, select = -c(orig_food_id))

ffq_cont <- filter(content, content$food_id == FFQ_ID)
ffq_cont <- subset(ffq_cont, select = -c(orig_food_id))

#Adding data from compound and content file together
full_chem <- ffq_comp %>% full_join(ffq_cont)
full_chem <- subset(full_chem, select = -c(id, source_id, orig_citation, creator_id, updater_id, created_at, updated_at))
write.table(full_chem, file = 'full_chem.txt')

#Add ffq_group to food_id
write.table(full_chem$food_id, file = 'ffqfood.txt')

foodffq <- read_xlsx('ffqfood.xlsx')

full_chem <- full_chem[order(full_chem$food_id), ] #descending order
chem_ffq <- cbind(full_chem, foodffq)


# Splitting dataframe into the different ffq elements 
# expect 62 elements minus two elements not found in compound/content datafiles --> 60
splitchem_1 <- split.data.frame(chem_ffq, chem_ffq$ffq_group_1)
splitchem_2 <- split.data.frame(chem_ffq, chem_ffq$ffq_group_2)
splitchem_3 <- split.data.frame(chem_ffq, chem_ffq$ffq_group_3)

splitchem_ffq <- c(splitchem_1, splitchem_2, splitchem_3)
splitchem_ffq[["NA"]] <- NULL

#Saving each element of list output as .txt file
sapply(names(splitchem_ffq), 
       function (x) write.table(splitchem_ffq[[x]], file=paste(x, "txt", sep=".") ) 
       )
