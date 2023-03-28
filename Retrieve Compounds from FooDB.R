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
library(reshape2)

# Load data
food <- read.csv('Food.csv')
content <- read.csv('Content.csv')
compounds <- read.csv('Compound.csv')
nutrients <- read.csv('Nutrient.csv')
id <- read_xlsx('idlist.xlsx')
ffq <- read_xlsx('foodid_ffqgroup.xlsx')

##========================= CLEAN DATA AND FILL/REMOVE NA'S ===========================##

# Filter for foods asked in FFQ and add food_id column to idlist
ffq_match <- food[match(id$public_id, food$public_id), ]
id <- cbind(id, ffq_match$id) 
colnames(id)[8] <- 'food_id'
write.csv2(id, file = 'ID_foodid.csv')

#Get ids, subset table 'content'
foodid <- as.character(id$food_id)
ffq_content <- content[content$food_id %in% foodid, ] #filter content table based on food_id in character

# Filter
chem_ffq <- ffq_content
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

##================== SPLITTING DATA ON FOOD_ID AND COLLAPSING BASED ON SOURCE_ID ==========##

#Split whole chem df on food_id
split_food_id <- split(chem_ffq_comp_nut_fill, chem_ffq_comp_nut_fill$food_id)

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
chem_ffq_collapse <- merge(chem_collapse, ffq, by = 'food_id', allow.cartesian = TRUE ) # adding ffq elements to dataframe

#Check NA's in source_name
sum(is.na(chem_ffq_collapse$orig_source_name)) # 430 NA items in orig_source_name
NA_sourcename <- chem_ffq_collapse[is.na(chem_ffq_collapse$orig_source_name),]

##================= SPLITTING COLLAPSED DATA (SOURCE_ID) ON FFQ_GROUP AND COLLAPSING AGAIN BASED ON SOURCE_ID ============##

# Split full dataframe based on ffq_element
ffq_chem_list <- split(chem_ffq_collapse, chem_ffq_collapse$ffq_group)

# Loop the collapse function over the ffq split dataframe and store in list
ffqgroup <- as.character(ffq$ffq_group) # setting ffq_groups as character

split_ffq_collapse <- list() #empty list to store output of for loop in
for (y in ffqgroup) {
  out <- collapse(ffq_chem_list[[y]])
  split_ffq_collapse[[y]] <- out
}

##======================== NEW COLLAPSED DATAFRAME  ============================##

# Melt list into one dataframe and merge the ffq_groups to one column
ffq_chem_melt <- melt(split_ffq_collapse, measure.vars = 'food_id', level = 'ffq_group') 
colnames(ffq_chem_melt) [19] <- 'ffq_group'
colnames(ffq_chem_melt) [18] <- 'food_id'

# Check NA's in Source_Name
sum(is.na(ffq_chem_melt$orig_source_name)) # 285 NA items in orig_source_name
NA_name <- ffq_chem_melt[is.na(ffq_chem_melt$orig_source_name),]

# Make seperate dataframe for energy per food item
energy <- ffq_chem_melt[ffq_chem_melt$orig_source_name == 'Energy', ]
  write.csv2(energy, file = 'energy.csv')

##======================= FINAL DATAFRAME, ENERGY IN SEPERATE FILE, CONVERTING UNITS TO MG/100G =============================##

# Create final df containing ffq_group, source_name, orig_content, orig_unit
final_ffqchem <- subset(ffq_chem_melt, select = c(ffq_group, food_id, orig_source_name, orig_content, orig_unit, source_id))

# Clean DF by removing remaining NA's in orig_source_name and setting all where unit is mg/100g to the same format
final_ffqchem$orig_unit <- gsub('mg/100 g', 'mg/100g', final_ffqchem$orig_unit)

# Set all rows where unit is mg/kg '' to mg/kg
final_ffqchem$orig_unit <- gsub('mg/kg puree', 'mg/kg', final_ffqchem$orig_unit)
final_ffqchem$orig_unit <- gsub('mg/kg fresh sample', 'mg/kg', final_ffqchem$orig_unit)
final_ffqchem$orig_unit <- gsub('mg/kg fresh weight', 'mg/kg', final_ffqchem$orig_unit)
# Set all rows where unit is µg/kg '' to µg/kg
final_ffqchem$orig_unit <- gsub('ug/kg fresh weight', 'µg/kg', final_ffqchem$orig_unit)
final_ffqchem$orig_unit <- gsub('ug/kg dry weight', 'µg/kg', final_ffqchem$orig_unit)
# Set all rows where unit is Î+--TE to mg-TE
final_ffqchem$orig_unit <- gsub('Î±-TE', 'mg-TE', final_ffqchem$orig_unit)

# Conversions from measured unit (mg/kg or µg/kg) to mg/100g
final_ffqchem$new_content <- ifelse(final_ffqchem$orig_unit == 'mg/kg', final_ffqchem$orig_content/10, final_ffqchem$orig_content)
final_ffqchem$new_unit <- gsub('mg/kg', 'mg/100g', final_ffqchem$orig_unit)
final_ffqchem$new_content <- ifelse(final_ffqchem$orig_unit == 'µg/kg', final_ffqchem$orig_content/1000, final_ffqchem$new_content)
final_ffqchem$new_unit <- gsub('µg/kg', 'mg/100g', final_ffqchem$new_unit)

# Vitamin D and A measured in IU to mg/100g
final_ffqchem$new_content <- ifelse(final_ffqchem$orig_source_name == 'Vitamin D', final_ffqchem$orig_content/0.000025, final_ffqchem$new_content)
final_ffqchem$new_content <- ifelse(final_ffqchem$food_id == '721' & final_ffqchem$orig_source_name == 'Vitamin A, IU', final_ffqchem$orig_content/0.0006, final_ffqchem$new_content)
final_ffqchem$new_content <- ifelse(final_ffqchem$food_id != '721' & final_ffqchem$orig_source_name == 'Vitamin A, IU', final_ffqchem$orig_content/0.0003, final_ffqchem$new_content)
final_ffqchem$new_unit <- gsub('IU', 'mg/100g', final_ffqchem$new_unit)

# Vitamin A measured in RE to mg/100g
plant_A <- c(649, 65, 939, 277, 98, 16, 125, 175, 22, 274, 24, 630, 47, 788, 105, 605, 59, 836, 268)
plant_A <- as.character(plant_A)
for (a in plant_A) {
  z <- ifelse(final_ffqchem$orig_source_name == 'Vitamin A, total', final_ffqchem$orig_content/0.006, final_ffqchem$new_content)
  final_ffqchem$new_content <- z
}
animal_A <- c(667, 714, 709, 710, 711, 605, 633, 873, 506, 761, 549, 669, 634)
animal_A <- as.character(animal_A)
for (b in plant_A) {
  t <- ifelse(final_ffqchem$orig_source_name == 'Vitamin A, total', final_ffqchem$orig_content/0.001, final_ffqchem$new_content)
  final_ffqchem$new_content <- t
}
final_ffqchem$new_unit <- gsub('RE', 'mg/100g', final_ffqchem$new_unit)

# For Vitamin E and B3 the content stays the same as 1 NE or mg-TE = 1 mg
final_ffqchem$new_unit <- gsub('NE', 'mg/100g', final_ffqchem$new_unit)
final_ffqchem$new_unit <- gsub('mg-TE', 'mg/100g', final_ffqchem$new_unit)
