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
library(lubridate) #for date formatting
library(psych) #for correlation analysis


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

## ---- clean-css-data --------
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

#Compute correlation matrix for serovars to find co-occurring serovars
sero_only <- as.matrix(select(css, -c(ID, Month, Site, System, Season))) #Only serovar proportions
sero_corr <- psych::corr.test(sero_only, method = "spearman")

#Converting results of correlation matrix into more accessible format
ut <- upper.tri(sero_corr$r) #Denote upper half of correlation matrix, will remove redundancy
sero_corr_format <- data.frame(
  row = rownames(sero_corr$r)[row(sero_corr$r)[ut]],
  column = rownames(sero_corr$r)[col(sero_corr$r)[ut]],
  cor = sero_corr$r[ut],
  p = sero_corr$p[ut]
)
#Filter correlations to only include p-values below 0.05
sero_corr_format <- sero_corr_format %>%
  filter(p <= 0.05)

#Count occurrences of each serovar in each system
serovars_in_systems <- css %>%
  group_by(System) %>%
  summarize(across(all_of(serovar_cols), ~ sum(. > 0), .names = "count_{.col}"))

#Find serovars present at least 5 times in each system
serovars_meeting_criteria <- serovars_in_systems %>%
  summarize(across(starts_with("count_"), ~ all(. >= 5)))

#Extract the serovar names that meet the criteria
serovars_to_keep <- names(serovars_meeting_criteria)[which(serovars_meeting_criteria[1, ] == TRUE)]
serovars_to_keep <- gsub("count_", "", serovars_to_keep) # Remove prefix to get original names
serovars_to_keep

#Cross-reference serovars that meet criteria with serovars that were correlated together
sero_corr_format <- sero_corr_format %>%
  filter(row %in% serovars_to_keep | column %in% serovars_to_keep)
#Of the four serovars that occur 5+ times in each system, Give I and Muenchen I are correlated
#with Infantis. Muenchen I is also correlated with Aqua/Inverness I will add these three 
#serovars to the processed data. Of note, Give I is also correlated with Rubislaw, which also
#occurs 5+ times in each system

#Get unique row and column names that are not already in the list of serovars to keep
rows_not_in_list <- setdiff(unique(sero_corr_format$row), serovars_to_keep)
cols_not_in_list <- setdiff(unique(sero_corr_format$column), serovars_to_keep)
all_not_in_list <- union(rows_not_in_list, cols_not_in_list) #Combine

#Filter the original data based on the identified serovars
filtered_css <- css %>%
  select(1:5, all_of(c(serovars_to_keep, all_not_in_list)))

#Add binary variables for each serovar that indicate presence/absence
filtered_css <- filtered_css %>%
  mutate(`Give I Prev` = ifelse(`Give I` > 0, 1, 0)) %>%
  mutate(`Muenchen I Prev` = ifelse(`Muenchen I` > 0, 1, 0)) %>%
  mutate(`Rubislaw Prev` = ifelse(Rubislaw > 0, 1, 0)) %>%
  mutate(`Typhimurium Prev` = ifelse(Typhimurium > 0, 1, 0)) %>%
  mutate(`Aqua/Inverness Prev` = ifelse(`Aqua/Inverness` > 0, 1, 0)) %>%
  mutate(`Infantis Prev` = ifelse(Infantis > 0, 1, 0))

## ---- clean-weather-data --------
#Goal: incorporate weather data into css data

#Define start date of the study
start_date <- as.Date("2021-11-01")

#Change date format in the css data (month of study -> yyyy/mm)
filtered_css <- filtered_css %>%
  mutate(Date = start_date %m+% months(filtered_css$Month)) %>%
  mutate(YearMonth = format(Date, "%Y-%m"))

#Define new data frames of the weather data for manipulation
weatherA <- weatherA_raw
weatherB <- weatherB_raw
weatherCD <- weatherCD_raw

#Change the date format of each of the weather datasets (year and day of year -> yyyy/mm)
weatherA <- weatherA %>%
  mutate(Date = as.Date(`Julian Day` - 1, origin = paste0(Year, "-01-01"))) %>%
  mutate(YearMonth = format(Date, "%Y-%m"))
weatherB <- weatherB %>%
  mutate(Date = as.Date(`Julian Day` - 1, origin = paste0(Year, "-01-01"))) %>%
  mutate(YearMonth = format(Date, "%Y-%m"))
weatherCD <- weatherCD %>%
  mutate(Date = as.Date(`Julian Day` - 1, origin = paste0(Year, "-01-01"))) %>%
  mutate(YearMonth = format(Date, "%Y-%m"))

#Separate css data by system
cssA <- filter(filtered_css, System == "A")
cssB <- filter(filtered_css, System == "B")
cssCD <- filter(filtered_css, System == "C" | System == "D") #sites c and d kept together bc of shared weather station

#Join the weather data with the css data
siteA <- full_join(weatherA, cssA, by = "YearMonth") %>% na.omit()
siteB <- full_join(weatherB, cssB, by = "YearMonth") %>% na.omit()
siteCD <- full_join(weatherCD, cssCD, by = "YearMonth") %>% na.omit()

#Check summaries of each combined site dataset
summary(siteA)
summary(siteB)
summary(siteCD)

#Remove redundant variables from site data
siteA <- siteA %>% select(-c(...2, Date.x, Date.y))
siteB <- siteB %>% select(-c(Date.x, Date.y, `Site ID`))
siteCD <- siteCD %>% select(-c(Date.x, Date.y))

## ---- savedata --------
#Combine individual site datasets
processed_data <- rbind(siteA, siteB, siteCD)

#Location to save file
save_data_location <- here::here("data","processed-data","processeddata.rds")
saveRDS(processed_data, file = save_data_location)

