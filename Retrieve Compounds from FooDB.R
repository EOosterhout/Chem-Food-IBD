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

#Filter for foods asked in FFQ
ffq <- food[match(id$public_id, food$public_id), ]
id <- cbind(id, ffq$id) 
colnames(id)[2] <- 'food_id'
write.csv(id, file = "IDffq.csv")

#Filter food item compounds
FFQ_ID <- id$food_id

ffq_comp <- filter(compounds, compounds$food_id == FFQ_ID)
ffq_comp <- as.data.frame(ffq_comp)

ffq_cont <- filter(content, content$food_id == FFQ_ID)

#Remove NA values from dataframe (optional)
x <- complete.cases(ffq_comp$orig_content)
ffq_comp <- filter(ffq_comp, x)

y <- complete.cases(ffq_cont$orig_content)
ffq_cont <- filter(ffq_cont, y)

#Splitting dataframe into the different fooditems
splitcomp <- split.data.frame(ffq_comp, ffq_comp$food_id)
splitcont <- split.data.frame(ffq_cont, ffq_cont$food_id)

#Saving each element of listoutput as .csv file
sapply(names(splitcomp), 
       function (x) write.table(splitcomp[[x]], file=paste(x, "txt", sep=".") ) 
       )
sapply(names(splitcont), 
       function (x) write.table(splitcont[[x]], file=paste(x, "txt", sep=".") ) 
       )  
