## ---- packages --------
#Load needed packages
library(here) #for data loading/saving
library(tidyverse)
library(skimr)
library(ggplot2)
library(corrplot)
library(car)

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

#Save bar-graphs for categorical variables as Figure 1
ggsave("results/figures/Figure1.png", last_plot())

## ---- correlations ----
#Create correlation matrix using continuous data
corr <- cor(numeric_vars)
corr
#Create correlation plot from matrix
corrplot <- {
  corrplot(corr, 
           title = "Correlation Matrix of Weather Variables and Serovars",
           mar = c(0, 0, 1, 0),
           number.cex = 0.5,
           number.digits = 2);
  recordPlot()
}

#Save correlation plot as Figure 2
ggsave(filename = "results/figures/Figure2.png", plot = replayPlot(corrplot))

#Strong positive correlation between month and year - use year-month engineered variable
#Strong positive correlation between min and max air temp - include interaction term
#Strong positive correlation between total rain and avg humidity - include interaction term
#Strong positive correlation between avg wind speed and max wind speed - use avg wind speed
#Strong negative correlation between total rain and total solar radiation - include interaction term

## ---- scatterplots ----
#Create subset of numeric variables to exclude prevalence data
numeric_vars_no_prev <- numeric_vars %>%
  select(-c(`Give I Prev`, `Infantis Prev`, `Muenchen I Prev`, `Rubislaw Prev`, 
            `Typhimurium Prev`, `Aqua/Inverness Prev`))

#Create scatterplot matrix from other variables
scatterplotMatrix(numeric_vars_no_prev)
#Complex relationships between month/day and weather variables (expected due to relationship of seasons) - use season as proxy


#Create individual plots for variables that seem to have interactions from the matrix
#Max vs. Min Temp
ggplot(processed_data, aes(x = `Min Air Temperature(F)`, y = `Max Air Temperature(F)`)) +
  geom_point() +
  geom_smooth(method = lm) +
  theme_minimal()

#Humidity vs. Rain
ggplot(processed_data, aes(x = `Avg Relative Humidity(%)`, y = `Total Rain(in)`)) +
  geom_point() +
  geom_smooth(method = lm) +
  theme_minimal()
#Many zeros for rain - interaction may not be significant

#Rain and solar radiation
ggplot(processed_data, aes(x = `Total Solar Radiation(MJ/m^2)`, y = `Total Rain(in)`)) +
  geom_point() +
  geom_smooth(method = lm) +
  theme_minimal()