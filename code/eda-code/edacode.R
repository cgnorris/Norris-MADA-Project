## ---- packages --------
#Load needed packages
library(here) #for data loading/saving
library(tidyverse)
library(skimr)
library(ggplot2)
library(corrplot)

## ---- load-data --------
#Set data path of processed data
data_path <- here("data", "processed-data", "processeddata.rds")

#Load data
data <- readRDS(data_path)

## ---- overall-data-exploration ----
#Get an overhead view of the processed data
str(data) #Get structure
summary(data) #Get summary statistics for all variables in the data
skimr::skim(data) #Get more detailed summary statistics 
dim(data) #Get dimenstions of the data frame (240 x 27)
names(data) #Get all column names
sapply(data, class) #Check data types of each column
colSums(is.na(data)) #Check each column for NA values

## ---- numeric-data-exploration ----
#Make a data-frame only containing numeric variables
numeric_vars <- data %>% select(where(is.numeric))

#Histograms for numeric variables
numeric_vars %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = value)) +
  geom_histogram(bins = 30, fill = "grey70", color = "black") +
  facet_wrap(~ variable, scales = "free") +
  labs(title = "Distributions of Numeric Variables") +
  theme_minimal()

#Boxplots for numeric variables
numeric_vars %>%
  select(-c(`Julian Day`, Year)) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  ggplot(aes(x = variable, y = value)) +
  geom_boxplot(fill = "orange", color = "black") +
  coord_flip() +
  labs(title = "Boxplots of Numeric Variables") +
  theme_minimal()

## ---- categorical-data-exploration ----
#Create data frame for categorical variables
cat_vars <- data %>% select(where(is.character))

#Bar graphs for categorical variables
cat_vars %>%
  select(-ID) %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  group_by(variable, value) %>%
  summarise(count = n(), .groups = 'drop') %>%
  ggplot(aes(x = value, y = count, fill = variable)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ variable, scales = "free") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Categorical Variable Distributions", x = "Value", y = "Count")
#Note: The highest number of positive samples here were in spring, system C, and in June of 2022,
#August of 2023, and October of 2023

## ---- correlations ----