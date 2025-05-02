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

# ---- data-splitting ----
#Split data into 80% training data and 20% testing data
#Set random seed
rngseed = 1234
set.seed(rngseed)

#Data splitting (80% training, 20% testing)
data_split <- initial_split(df, prop = 0.8)
train <- training(data_split)
test <- testing(data_split)

#Cross-validation: 5-fold CV
set.seed(rngseed)
cv_folds <- vfold_cv(train, v = 5)


# ---- recipe-creation ----
#Proportion model for Give I
give_prop <- recipe(data = train, `Give I` ~ `Muenchen I` + Rubislaw + Typhimurium + `Aqua/Inverness` +
                      Infantis + `Max Air Temperature(F)` + `Min Air Temperature(F)` + 
                      `Avg Relative Humidity(%)` + `Avg Wind Speed(mph)` + `Total Solar Radiation(MJ/m^2)` +
                      `Total Rain(in)` + Site + System + Season) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_interact(terms = ~ `Max Air Temperature(F)`:`Min Air Temperature(F)` + 
                  `Avg Relative Humidity(%)`:`Total Rain(in)` + 
                  `Total Rain(in)`:`Total Solar Radiation(MJ/m^2)`)


#Prevalence model for Give I
give_prev <- recipe(data = train, `Give I Prev` ~ `Muenchen I Prev` + `Rubislaw Prev` + 
                      `Typhimurium Prev` + `Aqua/Inverness Prev` +
                      `Infantis Prev` + `Max Air Temperature(F)` + `Min Air Temperature(F)` + 
                      `Avg Relative Humidity(%)` + `Avg Wind Speed(mph)` + `Total Solar Radiation(MJ/m^2)` +
                      `Total Rain(in)` + Site + System + Season) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_interact(terms = ~ `Max Air Temperature(F)`:`Min Air Temperature(F)` + 
                  `Avg Relative Humidity(%)`:`Total Rain(in)` + 
                  `Total Rain(in)`:`Total Solar Radiation(MJ/m^2)`)

#Proportion model for Muenchen I
muenchen_prop <- recipe(data = train, `Muenchen I` ~ `Give I` + Rubislaw + Typhimurium + `Aqua/Inverness` +
                      Infantis + `Max Air Temperature(F)` + `Min Air Temperature(F)` + 
                      `Avg Relative Humidity(%)` + `Avg Wind Speed(mph)` + `Total Solar Radiation(MJ/m^2)` +
                      `Total Rain(in)` + Site + System + Season) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_interact(terms = ~ `Max Air Temperature(F)`:`Min Air Temperature(F)` + 
                  `Avg Relative Humidity(%)`:`Total Rain(in)` + 
                  `Total Rain(in)`:`Total Solar Radiation(MJ/m^2)`)

#Prevalence model for Muenchen I
muenchen_prev <- recipe(data = train, `Muenchen I Prev` ~ `Give I Prev` + `Rubislaw Prev` + 
                      `Typhimurium Prev` + `Aqua/Inverness Prev` +
                      `Infantis Prev` + `Max Air Temperature(F)` + `Min Air Temperature(F)` + 
                      `Avg Relative Humidity(%)` + `Avg Wind Speed(mph)` + `Total Solar Radiation(MJ/m^2)` +
                      `Total Rain(in)` + Site + System + Season) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_interact(terms = ~ `Max Air Temperature(F)`:`Min Air Temperature(F)` + 
                  `Avg Relative Humidity(%)`:`Total Rain(in)` + 
                  `Total Rain(in)`:`Total Solar Radiation(MJ/m^2)`)

#Proportion model for Rubislaw
rubislaw_prop <- recipe(data = train, `Rubislaw` ~ `Give I` + `Muenchen I` + Typhimurium + `Aqua/Inverness` +
                      Infantis + `Max Air Temperature(F)` + `Min Air Temperature(F)` + 
                      `Avg Relative Humidity(%)` + `Avg Wind Speed(mph)` + `Total Solar Radiation(MJ/m^2)` +
                      `Total Rain(in)` + Site + System + Season) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_interact(terms = ~ `Max Air Temperature(F)`:`Min Air Temperature(F)` + 
                  `Avg Relative Humidity(%)`:`Total Rain(in)` + 
                  `Total Rain(in)`:`Total Solar Radiation(MJ/m^2)`)

#Prevalence model for Rubislaw
rubislaw_prev <- recipe(data = train, `Rubislaw Prev` ~ `Give I Prev` + `Muenchen I Prev` + 
                      `Typhimurium Prev` + `Aqua/Inverness Prev` +
                      `Infantis Prev` + `Max Air Temperature(F)` + `Min Air Temperature(F)` + 
                      `Avg Relative Humidity(%)` + `Avg Wind Speed(mph)` + `Total Solar Radiation(MJ/m^2)` +
                      `Total Rain(in)` + Site + System + Season) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_interact(terms = ~ `Max Air Temperature(F)`:`Min Air Temperature(F)` + 
                  `Avg Relative Humidity(%)`:`Total Rain(in)` + 
                  `Total Rain(in)`:`Total Solar Radiation(MJ/m^2)`)

#Proportion model for Typhimurium
typhimurium_prop <- recipe(data = train, Typhimurium ~ `Give I` + `Muenchen I` + Rubislaw + `Aqua/Inverness` +
                      Infantis + `Max Air Temperature(F)` + `Min Air Temperature(F)` + 
                      `Avg Relative Humidity(%)` + `Avg Wind Speed(mph)` + `Total Solar Radiation(MJ/m^2)` +
                      `Total Rain(in)` + Site + System + Season) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_interact(terms = ~ `Max Air Temperature(F)`:`Min Air Temperature(F)` + 
                  `Avg Relative Humidity(%)`:`Total Rain(in)` + 
                  `Total Rain(in)`:`Total Solar Radiation(MJ/m^2)`)

#Prevalence model for Typhimurium
typhimurium_prev <- recipe(data = train, `Typhimurium Prev` ~ `Give I Prev` + `Muenchen I Prev` + 
                      `Rubislaw Prev` + `Aqua/Inverness Prev` +
                      `Infantis Prev` + `Max Air Temperature(F)` + `Min Air Temperature(F)` + 
                      `Avg Relative Humidity(%)` + `Avg Wind Speed(mph)` + `Total Solar Radiation(MJ/m^2)` +
                      `Total Rain(in)` + Site + System + Season) %>%
  step_normalize(all_numeric_predictors()) %>%
  step_dummy(all_nominal_predictors()) %>%
  step_interact(terms = ~ `Max Air Temperature(F)`:`Min Air Temperature(F)` + 
                  `Avg Relative Humidity(%)`:`Total Rain(in)` + 
                  `Total Rain(in)`:`Total Solar Radiation(MJ/m^2)`)

# ---- workflow creation ----
#Define metrics for evaluation
reg_metrics <- metric_set(rmse, rsq)

#Workflows for all models
#Give I beta
give_beta_wf <- workflow() %>%
  add_recipe(give_prop) %>%
  add_model(beta_reg)

#Give I logistic
give_log_wf <- workflow() %>%
  add_recipe(give_prev) %>%
  add_model(log_reg)

#Give I regression rf
give_regrf_wf <- workflow() %>%
  add_recipe(give_prop) %>%
  add_model(rf_reg_model)

#Give I classification rf
give_classrf_wf <- workflow() %>%
  add_recipe(give_prev) %>%
  add_model(rf_class_model)

#Muenchen I beta
muenchen_beta_wf <- workflow() %>%
  add_recipe(muenchen_prop) %>%
  add_model(beta_reg)

#Muenchen I logistic
muenchen_log_wf <- workflow() %>%
  add_recipe(muenchen_prev) %>%
  add_model(log_reg)

#Muenchen I regression rf
muenchen_regrf_wf <- workflow() %>%
  add_recipe(muenchen_prop) %>%
  add_model(rf_reg_model)

#Muenchen I classification rf
muenchen_classrf_wf <- workflow() %>%
  add_recipe(muenchen_prev) %>%
  add_model(rf_class_model)
            
#Rubislaw beta
rubislaw_beta_wf <- workflow() %>%
  add_recipe(rubislaw_prop) %>%
  add_model(beta_reg)

#Rubislaw logistic
rubislaw_log_wf <- workflow() %>%
  add_recipe(rubislaw_prev) %>%
  add_model(log_reg)

#Rubislaw I regression rf
rubislaw_regrf_wf <- workflow() %>%
  add_recipe(rubislaw_prop) %>%
  add_model(rf_reg_model)

#Rubislaw I classification rf
rubislaw_classrf_wf <- workflow() %>%
  add_recipe(rubislaw_prev) %>%
  add_model(rf_class_model)

#typhimurium beta
typhimurium_beta_wf <- workflow() %>%
  add_recipe(typhimurium_prop) %>%
  add_model(beta_reg)

#typhimurium logistic
typhimurium_log_wf <- workflow() %>%
  add_recipe(typhimurium_prev) %>%
  add_model(log_reg)

#typhimurium I regression rf
typhimurium_regrf_wf <- workflow() %>%
  add_recipe(typhimurium_prop) %>%
  add_model(rf_reg_model)

#typhimurium I classification rf
typhimurium_classrf_wf <- workflow() %>%
  add_recipe(typhimurium_prev) %>%
  add_model(rf_class_model)

## ---- beta-regression-fitting ----
#Cross-validation for Give beta regression model
give_beta_cv_results <- fit_resamples(
  give_beta_wf,
  resamples = cv_folds,
  metrics = reg_metrics
)

#Cross-validation for Muenchen beta regression model
muenchen_beta_cv_results <- fit_resamples(
  muenchen_beta_wf,
  resamples = cv_folds,
  metrics = reg_metrics
)

#Cross-validation for Rubislaw beta regression model
rubislaw_beta_cv_results <- fit_resamples(
  rubislaw_beta_wf,
  resamples = cv_folds,
  metrics = reg_metrics
)

#Cross-validation for Typhimurium beta regression model
typhimurium_beta_cv_results <- fit_resamples(
  typhimurium_beta_wf,
  resamples = cv_folds,
  metrics = reg_metrics
)