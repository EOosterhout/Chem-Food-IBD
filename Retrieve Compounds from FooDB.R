##Title: Request compound information in FoodB
##Author: Emily Oosterhout
## setwd('C:/DATA FOOD COMPONENT ANALYSIS/RP2_ChemIBDFood')
##

##========================= LOAD DATA AND IMPORT LIBRARIES =============================##

#Import libraries
library(readxl)
library(dplyr)
library(tidyverse)
library(data.table)

# Load data
food <- read.csv('Food.csv')
content <- read.csv('Content.csv')
compounds <- read.csv('Compound.csv')
nutrients <- read.csv('Nutrient.csv')
id <- read_xlsx('idlist.xlsx')

##========================= CLEAN DATA AND FILL/REMOVE NA'S ===========================##

# Filter for foods asked in FFQ and add food_id column to idlist
ffq <- food[match(id$public_id, food$public_id), ]
id <- cbind(id, ffq$id) 
colnames(id)[8] <- 'food_id'

#Get ids, subset table 'content'
foodid <- as.character(id$food_id)
ffq_content <- content[content$food_id %in% foodid, ] #filter content table based on food_id in character

# Filter and merge content with ffq elements
chem_ffq <- merge(ffq_content, id, by = 'food_id', ) # add ffq elements to chem data
chem_ffq[chem_ffq == ""] <- NA # set all blank space to NA
chem_ffq <- chem_ffq %>% drop_na(orig_content) # drop all rows that don't contain chem orig_content data

# set column names in compounds dataframe to colnames in content dataframe
colnames(compounds)[3] <- 'orig_source_name'
colnames(compounds)[1] <- 'source_id'

# Merge content source name with compound source name
chem_ffq_comp <- left_join(chem_ffq, compounds, by ="source_id", ) # joins based on values present in dataframe x
chem_ffq$orig_source_name.x <- as.character(chem_ffq_comp$orig_source_name.x) # set both columns in different dataframes to the same format
chem_ffq_comp$orig_source_name <- coalesce(chem_ffq_comp$orig_source_name.x, chem_ffq_comp$orig_source_name.y) # join source_name in one column
chem_ffq_comp <- subset(chem_ffq_comp, select = -c(orig_source_name.x, orig_source_name.y)) # remove the merged columns

# set column names in nutrients dataframe to colnames in content dataframe
colnames(nutrients)[5] <- 'orig_source_name'
colnames(nutrients)[1] <- 'source_id'

# Merge content source name with nutrients source name
chem_ffq_comp_nut <- left_join(chem_ffq_comp, nutrients, by ="source_id", ) # joins based on values present in dataframe x
chem_ffq_comp_nut$orig_source_name <- coalesce(chem_ffq_comp_nut$orig_source_name.x, chem_ffq_comp_nut$orig_source_name.y) # join source_name in one column
chem_ffq_comp_nut <- subset(chem_ffq_comp_nut, select = -c(orig_source_name.x, orig_source_name.y)) # remove the merged columns

# Fill rest of NA values using fill()
chem_ffq_comp_nut_fill <- chem_ffq_comp_nut %>%
  group_by(source_id) %>%
  fill(orig_source_name, .direction = 'updown') %>%
  ungroup

# check if number of NA in orig_source_name is changed and what the NA's in the merged dataframe are.
sum(is.na(chem_ffq$orig_source_name)) #47271 NA Values
sum(is.na(chem_ffq_comp$orig_source_name)) #17160 NA values
sum(is.na(chem_ffq_comp_nut$orig_source_name)) #12062 NA values
sum(is.na(chem_ffq_comp_nut_fill$orig_source_name)) #7376 NA values

##======================== EXPLORING REMAINING NA's ===============================##

## Checking NA's + source_id's 
df <- select(chem_ffq_comp_nut
             , c('source_id','orig_source_name', 'orig_content', 'orig_unit'))
sum(df$source_id == 0) #5894 items with source_id equal to 0
df2 <- df[is.na(df$orig_source_name),]
sum(df2$orig_content == 0) 
# 5241 items have no source name, where orig_content is equal to 0
sum(df2$orig_content == 0 & df2$source_id == 0) 
# 2160 items have source_id equal to 0 & and orig_content equal to 0

split_source_id <- split(chem_ffq_comp_nut, chem_ffq_comp_nut$source_id)

first_sourceid <- lapply(split_source_id,'[',1,) #show first row of each source_id
first_sourceid <- rbindlist(first_sourceid, fill = TRUE)
first_NA <- first_sourceid[is.na(first_sourceid$orig_source_name),] #40 NA as first row
#Check if NA in first row contains a name in another row
sourceid <- as.character(first_NA$source_id)
NA_sourceid <- chem_ffq_comp_nut[chem_ffq_comp_nut$source_id %in% sourceid, ]
NA_sourceid_name <- subset(NA_sourceid, select = c(source_id, orig_source_name))

##================== SPLITTING DATA ON FOOD_ID AND COLLAPSING BASED ON SOURCE_ID ==========##

#Split whole chem df on food_id
split_food_id <- split(chem_ffq_comp_nut_fill, chem_ffq_comp_nut$food_id)

# Write function that groups dataframes in the list by source id, takes mean of orig_content
collapse <- function(df) { df %>% 
    group_by(source_id) %>%
    summarise( id=first(id),
               source_type=first(source_type),
               food_id=first(food_id),
               orig_food_id=first(orig_food_id),
               orig_food_common_name=first(orig_food_common_name),
               orig_food_scientific_name=first(orig_food_scientific_name),
               orig_food_part= first(orig_food_part),
               orig_source_id=first(orig_source_id),
               orig_source_name= first(orig_source_name),
               orig_content = mean(orig_content, na.rm = TRUE),
               orig_min = mean(orig_min, na.rm = TRUE),
               orig_max = mean(orig_max, na.rm = TRUE),
               orig_unit = first(orig_unit),
               citation_type=first(citation_type),
               standard_content = mean(standard_content, na.rm = TRUE),
               preparation_type=first(preparation_type))
}

# Loop the collapse function over the original split dataframe and store in list
split_foodid_collapse <- list() #empty list to store output of for loop in
for (x in foodid) {
  output <- collapse(split_food_id[[x]])
  split_foodid_collapse[[x]] <- output
}

# Creating a new dataframe with mean orig_content per compound 
chem_collapse <- rbindlist(split_foodid_collapse, fill = TRUE) #--> 9694 obs in total
chem_ffq_collapse <- merge(chem_collapse, id, by = 'food_id', ) # adding ffq elements to dataframe

#Check NA's in source_name
sum(is.na(chem_ffq_collapse$orig_source_name)) # 280 NA items in orig_source_name
NA_sourcename <- chem_ffq_collapse[is.na(chem_ffq_collapse$orig_source_name),]

##================= SPLITTING COLLAPSED DATA (SOURCE_ID) ON FFQ_GROUP AND COLLAPSING AGAIN BASED ON SOURCE_ID ============##

# Split full dataframe based on ffq_element
split_1 <- split(chem_ffq_collapse, chem_ffq_collapse$ffq_group_1)
split_2 <- split(chem_ffq_collapse, chem_ffq_collapse$ffq_group_2)
split_3 <- split(chem_ffq_collapse, chem_ffq_collapse$ffq_group_3)
split_4 <- split(chem_ffq_collapse, chem_ffq_collapse$ffq_group_4)
split_5 <- split(chem_ffq_collapse, chem_ffq_collapse$ffq_group_5)
split_6 <- split(chem_ffq_collapse, chem_ffq_collapse$ffq_group_6)
ffq_chem_list <- c(split_1, split_2, split_3, split_4, split_5, split_6)

# Loop the collapse function over the ffq split dataframe and store in list
ffqelement <- read_xlsx('ffq_group.xlsx') # data containing ffq_groups
ffqgroup <- as.character(ffqelement$ffq_group) # setting ffq_groups as character

split_ffq_collapse <- list() #empty list to store output of for loop in
for (y in ffqgroup) {
  out <- collapse(ffq_chem_list[[y]])
  split_ffq_collapse[[y]] <- out
}

# Creating a new dataframe with mean orig_content per compound --> 9837 obs in total
ffq_chem_collapse <- rbindlist(split_ffq_collapse, fill = TRUE)

# Check NA's in Source_Name
sum(is.na(ffq_chem_collapse$orig_source_name)) # 285 NA items in orig_source_name
NA_name <- ffq_chem_collapse[is.na(ffq_chem_collapse$orig_source_name),]

# Clean DF by removing remaining NA's in orig_source_name
ffq_chem_collapse_clean <- ffq_chem_collapse[!is.na(ffq_chem_collapse$orig_source_name),]



















