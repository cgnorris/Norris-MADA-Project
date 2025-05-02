#Analysis script

#Load needed packages
library(ggplot2)
library(broom) 
library(here) 
library(tidymodels)
library(tidyverse)
library(betareg)
library(ranger)


#Load data
df <- readRDS("data/processed-data/processeddata.rds")

# ---- extra-data-processing ----
#Remove variables from original processed data that wont be used in modelling
df <- df %>%
  select(-c(Year, `Julian Day`, `Max Wind Speed (m/s)`, ID, Month))

# ---- model-definition ----

#Beta regresssion (used when response is bounded 0-1, as with proportions)
#Beta regression models are not natively supported in the tinymodels framework, but it can be added as a custom model
#Define a new model
set_new_model("beta_reg")

#Set the model mode
set_model_mode("beta_reg", mode = "regression")

#Set the model engine
set_model_engine("beta_reg", mode = "regression", eng = "betareg")

#Set the model fitting function
set_fit(
  model = "beta_reg",
  eng = "betareg",
  mode = "regression",
  value = list(
    interface = "formula",
    protect = c("formula", "data"),
    func = c(pkg = "betareg", fun = "betareg"),
    defaults = list()
  )
)

#Set how to extract predictions
set_pred(
  model = "beta_reg",
  eng = "betareg",
  mode = "regression",
  type = "numeric",
  value = list(
    pre = NULL,
    post = NULL,
    func = c(pkg = "stats", fun = "predict"),
    args = list(object = expr(object), newdata = expr(new_data), type = "response")
  )
)

#Create model specification
beta_reg <- function(mode = "regression") {
  parsnip::new_model_spec("beta_reg", args = list(), mode = mode, eng_args = NULL, method = NULL, engine = NULL)
}

beta_reg <- beta_reg() %>%
  set_engine("betareg")

#Beta models won't work with exact 0s and 1s - transform proportion data accordingly
epsilon <- 1e-4

#Values < 0.0001 become 0.0001
#Values > 0.9999 become 0.9999
df$`Aqua/Inverness` <- pmin(pmax(df$`Aqua/Inverness`, epsilon), 1 - epsilon)
df$`Give I` <- pmin(pmax(df$`Give I`, epsilon), 1 - epsilon)
df$Infantis <- pmin(pmax(df$Infantis, epsilon), 1 - epsilon)
df$`Muenchen I` <- pmin(pmax(df$`Muenchen I`, epsilon), 1 - epsilon)
df$Typhimurium <- pmin(pmax(df$Typhimurium, epsilon), 1 - epsilon)
df$Rubislaw <- pmin(pmax(df$Rubislaw, epsilon), 1 - epsilon)

#Define logistic regression model
log_reg <- logistic_reg(mode = "classification") %>%
  set_engine("glm")

#Define random forest models
#Regression - proportions
rf_reg_model <- rand_forest(mode = "regression", trees = 1000) %>%
  set_engine("ranger")

#Classification - prevalence
rf_class_model <- rand_forest(mode = "classification", trees = 1000) %>%
  set_engine("ranger")


