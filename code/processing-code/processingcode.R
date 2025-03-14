###############################
# processing script
#
#this script loads the raw data, processes and cleans it 
#and saves it as Rds file in the processed-data folder
#
# Note the ## ---- name ---- notation
# This is done so one can pull in the chunks of code into the Quarto document
# see here: https://bookdown.org/yihui/rmarkdown-cookbook/read-chunk.html
###############################


## ---- packages --------
#Load required packages
library(readxl) #for loading Excel files
library(dplyr) #for data processing/cleaning
library(tidyr) #for data processing/cleaning
library(skimr) #for nice visualization of data 
library(here) #to set paths


## ---- load-data --------
#Set data paths
full_data <- here::here("data","raw-data","creek_1-3day_avg_weather_2021-2023.xlsx")

#Load data from individual sheets in full data
#Note: Weather data includes daily weather measurements, three day averages from the date of sampling, and monthly averages
#I will look at just weather on the day of sampling

#CRISPR Sero-Seq data (Salmonella serovar distributions) from each sampling date/location 
css_raw <- readxl::read_excel(full_data, sheet = "CREEK1")

#Weather data for system A
weatherA_raw <- readxl::read_excel(full_data, sheet = "GAINES 1 DAY")

#Weather data for system B
weatherB_raw <- readxl::read_excel(full_data, sheet = "WATHORT 1 DAY ")

#Weather data for systems C and D
#Note: The weather tower that was nearest systems C and D was the same
weatherCD_raw <- readxl::read_excel(full_data, sheet = "WATUGA 1 DAY ")

## ---- exploredata --------
#Take a look at the data
dplyr::glimpse(css_raw)

#Summary
summary(css_raw)

#Skim
skimr::skim(css_raw)

#Weather A
glimpse(weatherA_raw)
summary(weatherA_raw)
skim(weatherA_raw)

#Weather B
glimpse(weatherB_raw)
summary(weatherB_raw)
skim(weatherB_raw)

#Weather CD
glimpse(weatherCD_raw)
summary(weatherCD_raw)
skim(weatherCD_raw)

## ---- clean-data-1 --------
#With the css data, I would only like to focus my analysis on Salmonella serovars found in all systems

#Replace all NAs with 0
css <- css_raw %>%
  mutate_all(~replace(., is.na(.), 0))

#Identify columns containing serovar information
serovar_cols <- colnames(css)[-(1:5)]

#Change all character vectors to numeric vectors
css <- css %>%
  mutate(across(all_of(serovar_cols), ~ if(is.character(.)) as.numeric(.) else .))

#Find serovars present in each system
# serovars_in_systems <- css %>%
#   group_by(System) %>%
#   summarise(across(all_of(serovar_cols), ~ any(. > 0), .names = "present_{.col}")) %>%
#   summarise(across(starts_with("present_"), all))  #Keep only serovars found in all systems
# 
# #Extract serovar names present in all 4 systems
# common_serovars <- names(serovars_in_systems)[which(serovars_in_systems[1, ]==TRUE)]
# common_serovars <- gsub("present_", "", common_serovars) #Remove prefix to get original name
# common_serovars #Display common serovar names
# 
# #Filter based on common serovars
# css <- css %>%
#   select(1:5, all_of(common_serovars)

#Count occurances of each serovar in each system
serovars_in_systems <- css %>%
  group_by(System) %>%
  summarize(across(all_of(serovar_cols), ~ sum(. > 0), .names = ".count_{.col}"))

#Find serovars present at least 5 times in each system
serovars_meeting_criteria <- serovars_in_systems %>%
  summarize(across(starts_with("count_"), ~ all(. >= 5))) %>%
  select(where(~ . == TRUE))





## ---- cleandata2 --------
# Now we see that there is one person with a height of 6. 
# that could be a typo, or someone mistakenly entered their height in feet.
# If we don't know, we might need to remove this person.
# But let's assume that we somehow know that this is meant to be 6 feet, so we can convert to centimeters.
d2 <- d1 %>% dplyr::mutate( Height = replace(Height, Height=="6",round(6*30.48,0)) )
#height values seem ok now
skimr::skim(d2)


## ---- cleandata3 --------
# now let's look at weight
# there is a person with weight of 7000, which is impossible,
# and one person with missing weight.
# Note that the original data had an empty cell. 
# The codebook says that's not allowed, it should have been NA.
# R automatically converts empty values to NA.
# If you don't want that, you can adjust it when you load the data.
# to be able to analyze the data, we'll remove those individuals as well.
# Note: Some analysis methods can deal with missing values, so it's not always necessary to remove them. 
# This should be adjusted based on your planned analysis approach. 
d3 <- d2 %>%  dplyr::filter(Weight != 7000) %>% tidyr::drop_na()
skimr::skim(d3)


## ---- cleandata4 --------
# We also want to have Gender coded as a categorical/factor variable
# we can do that with simple base R code to mix things up
d3$Gender <- as.factor(d3$Gender)  
skimr::skim(d3)


## ---- cleandata5 --------
#now we see that there is another NA, but it's not "NA" from R 
#instead it was loaded as character and is now considered as a category.
#There is also an individual coded as "N" which is not allowed.
#This could be mistyped M or a mistyped NA. If we have a good guess, we could adjust.
#If we don't we might need to remove that individual.
#well proceed here by removing both the NA and N individuals
#since this keeps an empty category, I'm also using droplevels() to get rid of it
d4 <- d3 %>% dplyr::filter( !(Gender %in% c("NA","N")) ) %>% droplevels()
skimr::skim(d4)


## ---- savedata --------
# all done, data is clean now. 
# Let's assign at the end to some final variable
# makes it easier to add steps above
processeddata <- d4
# location to save file
save_data_location <- here::here("data","processed-data","processeddata.rds")
saveRDS(processeddata, file = save_data_location)



## ---- notes --------
# anything you don't want loaded into the Quarto file but 
# keep in the R file, just give it its own label and then don't include that label
# in the Quarto file

# Dealing with NA or "bad" data:
# removing anyone who had "faulty" or missing data is one approach.
# it's often not the best. based on your question and your analysis approach,
# you might want to do cleaning differently (e.g. keep individuals with some missing information)

# Saving data as RDS:
# I suggest you save your processed and cleaned data as RDS or RDA/Rdata files. 
# This preserves coding like factors, characters, numeric, etc. 
# If you save as CSV, that information would get lost.
# However, CSV is better for sharing with others since it's plain text. 
# If you do CSV, you might want to write down somewhere what each variable is.
# See here for some suggestions on how to store your processed data:
# http://www.sthda.com/english/wiki/saving-data-into-r-data-format-rds-and-rdata
