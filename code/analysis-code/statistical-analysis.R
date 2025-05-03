#Analysis script

#Load needed packages
library(ggplot2)
library(broom) 
library(here) 
library(tidymodels)
library(tidyverse)
library(betareg)
library(ranger)
library(knitr)
library(car)
library(vip)
library(DALEX)
library(DALEXtra)

#Load data
df <- readRDS("data/processed-data/processeddata.rds")

# ---- extra-data-processing ----
#Remove variables from original processed data that wont be used in modelling
df <- df %>%
  select(-c(Year, `Julian Day`, `Max Wind Speed (m/s)`, ID, Month))

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

#Convert prevalence data to factors
df$`Give I Prev` <- as.factor(df$`Give I Prev`)
df$`Aqua/Inverness Prev` <- as.factor(df$`Aqua/Inverness Prev`)
df$`Rubislaw Prev` <- as.factor(df$`Rubislaw Prev`)
df$`Typhimurium Prev` <- as.factor(df$`Typhimurium Prev`)
df$`Infantis Prev` <- as.factor(df$`Infantis Prev`)
df$`Muenchen I Prev` <- as.factor(df$`Muenchen I Prev`)

# ---- model-definition ----
#Define logistic regression model
log_reg <- logistic_reg(mode = "classification") %>%
  set_engine("glm")

#Define random forest model
#Regression - proportions
rf_reg_model <- rand_forest(mode = "regression", trees = 1000) %>%
  set_engine("ranger", importance = "permutation")

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
  step_dummy(all_nominal_predictors()) %>%
  step_interact(terms = ~ `Max Air Temperature(F)`:`Min Air Temperature(F)` + 
                  `Avg Relative Humidity(%)`:`Total Rain(in)` + 
                  `Total Rain(in)`:`Total Solar Radiation(MJ/m^2)`)

#Proportion model for Muenchen I
muenchen_prop <- recipe(data = train, `Muenchen I` ~ `Give I` + Rubislaw + Typhimurium + `Aqua/Inverness` +
                      Infantis + `Max Air Temperature(F)` + `Min Air Temperature(F)` + 
                      `Avg Relative Humidity(%)` + `Avg Wind Speed(mph)` + `Total Solar Radiation(MJ/m^2)` +
                      `Total Rain(in)` + Site + System + Season) %>%
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
  step_dummy(all_nominal_predictors()) %>%
  step_interact(terms = ~ `Max Air Temperature(F)`:`Min Air Temperature(F)` + 
                  `Avg Relative Humidity(%)`:`Total Rain(in)` + 
                  `Total Rain(in)`:`Total Solar Radiation(MJ/m^2)`)

#Proportion model for Rubislaw
rubislaw_prop <- recipe(data = train, `Rubislaw` ~ `Give I` + `Muenchen I` + Typhimurium + `Aqua/Inverness` +
                      Infantis + `Max Air Temperature(F)` + `Min Air Temperature(F)` + 
                      `Avg Relative Humidity(%)` + `Avg Wind Speed(mph)` + `Total Solar Radiation(MJ/m^2)` +
                      `Total Rain(in)` + Site + System + Season) %>%
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
  step_dummy(all_nominal_predictors()) %>%
  step_interact(terms = ~ `Max Air Temperature(F)`:`Min Air Temperature(F)` + 
                  `Avg Relative Humidity(%)`:`Total Rain(in)` + 
                  `Total Rain(in)`:`Total Solar Radiation(MJ/m^2)`)

#Proportion model for Typhimurium
typhimurium_prop <- recipe(data = train, Typhimurium ~ `Give I` + `Muenchen I` + Rubislaw + `Aqua/Inverness` +
                      Infantis + `Max Air Temperature(F)` + `Min Air Temperature(F)` + 
                      `Avg Relative Humidity(%)` + `Avg Wind Speed(mph)` + `Total Solar Radiation(MJ/m^2)` +
                      `Total Rain(in)` + Site + System + Season) %>%
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
  step_dummy(all_nominal_predictors()) %>%
  step_interact(terms = ~ `Max Air Temperature(F)`:`Min Air Temperature(F)` + 
                  `Avg Relative Humidity(%)`:`Total Rain(in)` + 
                  `Total Rain(in)`:`Total Solar Radiation(MJ/m^2)`)

# ---- workflow creation ----
#Define metrics for evaluation
reg_metrics <- metric_set(yardstick::rmse, yardstick::rsq)

#Define metrics for classification evaluation
class_metrics <- metric_set(accuracy, roc_auc)

#Workflows for all models
#Give I logistic
give_log_wf <- workflow() %>%
  add_recipe(give_prev) %>%
  add_model(log_reg)

#Give I rf
give_regrf_wf <- workflow() %>%
  add_recipe(give_prop) %>%
  add_model(rf_reg_model)

#Muenchen I logistic
muenchen_log_wf <- workflow() %>%
  add_recipe(muenchen_prev) %>%
  add_model(log_reg)

#Muenchen I rf
muenchen_regrf_wf <- workflow() %>%
  add_recipe(muenchen_prop) %>%
  add_model(rf_reg_model)

#Rubislaw logistic
rubislaw_log_wf <- workflow() %>%
  add_recipe(rubislaw_prev) %>%
  add_model(log_reg)

#Rubislaw I rf
rubislaw_regrf_wf <- workflow() %>%
  add_recipe(rubislaw_prop) %>%
  add_model(rf_reg_model)

#typhimurium logistic
typhimurium_log_wf <- workflow() %>%
  add_recipe(typhimurium_prev) %>%
  add_model(log_reg)

#typhimurium I rf
typhimurium_regrf_wf <- workflow() %>%
  add_recipe(typhimurium_prop) %>%
  add_model(rf_reg_model)

## ---- beta-regression-CV ----
#Standardize the predictors involved in interaction
train$`Max Air Temperature(F)_std` <- as.numeric(scale(train$`Max Air Temperature(F)`))
train$`Min Air Temperature(F)_std` <- as.numeric(scale(train$`Min Air Temperature(F)`))
train$`Avg Relative Humidity(%)_std` <- as.numeric(scale(train$`Avg Relative Humidity(%)`))
train$`Total Solar Radiation(MJ/m^2)_std` <- as.numeric(scale(train$`Total Solar Radiation(MJ/m^2)`))
train$`Total Rain(in)_std` <- as.numeric(scale(train$`Total Rain(in)`))

#Recalculate interaction terms using standardized predictors
train$Interaction_Temp <- train$`Max Air Temperature(F)_std` * train$`Min Air Temperature(F)_std`
train$Interaction_Humidity_Rain <- train$`Avg Relative Humidity(%)_std` * train$`Total Rain(in)_std`
train$Interaction_Solar_Rain <- train$`Total Solar Radiation(MJ/m^2)_std` * train$`Total Rain(in)_std`

#Re-run cross validation
set.seed(rngseed)
cv_folds <- vfold_cv(train, v = 5)

#Define a list of response variables and their corresponding formulas
response_formulas <- list(
  `Give I` = `Give I` ~ `Muenchen I` + Rubislaw + Typhimurium + `Aqua/Inverness` +
    Infantis + `Max Air Temperature(F)_std` + `Min Air Temperature(F)_std` +
    `Avg Relative Humidity(%)_std` + `Avg Wind Speed(mph)` +
    `Total Solar Radiation(MJ/m^2)_std` + `Total Rain(in)_std` + Site + System + Season +
    `Interaction_Temp` + `Interaction_Humidity_Rain` + `Interaction_Solar_Rain`,
  
  `Muenchen I` = `Muenchen I` ~ Rubislaw + Typhimurium + `Aqua/Inverness` + Infantis + `Give I` +
    `Max Air Temperature(F)_std` + `Min Air Temperature(F)_std` +
    `Avg Relative Humidity(%)_std` + `Avg Wind Speed(mph)` +
    `Total Solar Radiation(MJ/m^2)_std` + `Total Rain(in)_std` + Site + System + Season +
    `Interaction_Temp` + `Interaction_Humidity_Rain` + `Interaction_Solar_Rain`,
    
  `Rubislaw` = Rubislaw ~ `Muenchen I` + Typhimurium + `Aqua/Inverness` + Infantis + `Give I` +
  `Max Air Temperature(F)_std` + `Min Air Temperature(F)_std` +
    `Avg Relative Humidity(%)_std` + `Avg Wind Speed(mph)` +
    `Total Solar Radiation(MJ/m^2)_std` + `Total Rain(in)_std` + Site + System + Season +
    `Interaction_Temp` + `Interaction_Humidity_Rain` + `Interaction_Solar_Rain`,
  
  `Typhimurium` = Typhimurium ~ `Muenchen I` + Rubislaw + `Aqua/Inverness` + Infantis + `Give I` +
    `Max Air Temperature(F)_std` + `Min Air Temperature(F)_std` +
    `Avg Relative Humidity(%)_std` + `Avg Wind Speed(mph)` +
    `Total Solar Radiation(MJ/m^2)_std` + `Total Rain(in)_std` + Site + System + Season +
    `Interaction_Temp` + `Interaction_Humidity_Rain` + `Interaction_Solar_Rain`
)

#Initialize a data frame to store summary results
summary_table <- data.frame(
  Response = character(),
  Average_RMSE = numeric(),
  SE_RMSE = numeric(),
  stringsAsFactors = FALSE
)

#Perform manual cross-validation for each response variable
for (response_name in names(response_formulas)) {
  #Get the formula for the current response variable
  current_formula <- response_formulas[[response_name]]
  
  #Initialize a vector to store RMSE for each fold for the current response
  rmse_results <- numeric(length(cv_folds$splits))
  
  #Loop through each fold
  for (i in seq_along(cv_folds$splits)) {
    #Get the training and validation sets for the current fold
    train_fold <- analysis(cv_folds$splits[[i]])
    val_fold <- assessment(cv_folds$splits[[i]])
    
    #Fit the beta regression model on the training fold
    beta_model <- betareg(current_formula, data = train_fold)
    
    #Make predictions on the validation fold
    val_preds <- predict(beta_model, newdata = val_fold, type = "response")
    
    #Calculate RMSE for the current fold
    rmse_results[i] <- sqrt(mean((val_fold[[response_name]] - val_preds)^2))
  }
  
  #Calculate the average RMSE and SE across all folds for the current response
  avg_rmse <- mean(rmse_results)
  se_rmse <- sd(rmse_results)
  
  #Add the results to the summary table
  summary_table <- rbind(summary_table, data.frame(
    Response = response_name,
    Average_RMSE = avg_rmse,
    SE_RMSE = se_rmse
  ))
}

#Create a publication-ready summary table
kable(summary_table, format = "markdown", digits = 3, caption = "Cross-Validation Results")

# ---- beta-fitting ----

#Fit Give model
give_beta_model <- betareg(response_formulas$`Give I`, data = train)
summary(give_beta_model)
vif(give_beta_model)

#Get summary of Give model coefficients
give_beta_summary <- tidy(give_beta_model) %>%
  select(term, estimate, std.error, p.value) %>%
  rename(Coefficent = estimate,
         Std_Error = std.error,
         p_value = p.value) %>%
  mutate(Response = "Give I")

#Fit Muenchen model
muenchen_beta_model <- betareg(response_formulas$`Muenchen I`, data = train)
summary(muenchen_beta_model)
vif(muenchen_beta_model)

#Get summary of Muenchen coefficients
muenchen_beta_summary <- tidy(muenchen_beta_model) %>%
  select(term, estimate, std.error, p.value) %>%
  rename(Coefficent = estimate,
         Std_Error = std.error,
         p_value = p.value) %>%
  mutate(Response = "Muenchen I")

#Tried to fit a beta model for Rubislaw, model failed to run without error

#Fit Typhimurium model
typhimurium_beta_model <- betareg(response_formulas$Typhimurium, data = train)
summary(typhimurium_beta_model)
vif(typhimurium_beta_model)

#Get summary of Typhimurium coefficients
typhimurium_beta_summary <- tidy(typhimurium_beta_model) %>%
  select(term, estimate, std.error, p.value) %>%
  rename(Coefficent = estimate,
         Std_Error = std.error,
         p_value = p.value) %>%
  mutate(Response = "Typhimurium")

#Combine summaries
beta_model_results <- rbind(give_beta_summary, muenchen_beta_summary, typhimurium_beta_summary) %>%
  select(Response, term, Coefficent, Std_Error, p_value) %>%
  filter(p_value <= 0.05)

beta_model_results

#Fit null models
give_beta_null <- betareg(`Give I` ~ 1, data = train)
muenchen_beta_null <- betareg(`Muenchen I` ~ 1, data = train)
typhimurium_beta_null <- betareg(Typhimurium ~ 1, data = train)


#Extract RMSE and R-squared, compare
#Initalize summary table
summary_table <- data.frame(
  Model = character(),
  RMSE = numeric(),
  McFadden_R2 = numeric(),
  stringsAsFactors = FALSE
)

#Make list of beta models
beta_models <- list(give_beta_model = give_beta_model,
                    give_beta_null  = give_beta_null,
                    muenchen_beta_model = muenchen_beta_model,
                    muenchen_beta_null = muenchen_beta_null,
                    typhimurium_beta_model = typhimurium_beta_model,
                    typhimurium_beta_null = typhimurium_beta_null
                    )

#Loop through the models to calculate metrics
for (model_name in names(beta_models)) {
  #Get the current model
  model <- beta_models[[model_name]]
  print(beta)
  
  # --- Calculate RMSE ---
  residuals <- residuals(model, type = "response")
  rmse <- sqrt(mean(residuals^2))
  
  # --- Extract R-squared values ---
  pseudo_r_squared <- summary(model)$pseudo.r.squared
  mcfadden_r2 <- pseudo_r_squared[1]  # McFadden's R-squared
  
  # --- Add metrics to the summary table ---
  summary_table <- rbind(
    summary_table,
    data.frame(
      Model = model_name,
      RMSE = round(rmse, 3),
      McFadden_R2 = round(mcfadden_r2, 3)
    )
  )
}

kable(summary_table, format = "markdown", caption = "Summary of Beta Regression Models")
#Models only performed slightly better than the null models

#Make same transformations to testing data
#Standardize the predictors involved in interaction
test$`Max Air Temperature(F)_std` <- as.numeric(scale(test$`Max Air Temperature(F)`))
test$`Min Air Temperature(F)_std` <- as.numeric(scale(test$`Min Air Temperature(F)`))
test$`Avg Relative Humidity(%)_std` <- as.numeric(scale(test$`Avg Relative Humidity(%)`))
test$`Total Solar Radiation(MJ/m^2)_std` <- as.numeric(scale(test$`Total Solar Radiation(MJ/m^2)`))
test$`Total Rain(in)_std` <- as.numeric(scale(test$`Total Rain(in)`))

#Recalculate interaction terms using standardized predictors
test$Interaction_Temp <- test$`Max Air Temperature(F)_std` * test$`Min Air Temperature(F)_std`
test$Interaction_Humidity_Rain <- test$`Avg Relative Humidity(%)_std` * test$`Total Rain(in)_std`
test$Interaction_Solar_Rain <- test$`Total Solar Radiation(MJ/m^2)_std` * test$`Total Rain(in)_std`

#Make predictions using testing data
give_beta_preds <- predict(give_beta_model, newdata = test, type = "response")
muenchen_beta_preds <- predict(muenchen_beta_model, newdata = test, type = "response")
typhimurium_beta_preds <- predict(typhimurium_beta_model, newdata = test, type = "response")

#Combine results with testing data
give_beta_results <- test %>%
  mutate(predicted = give_beta_preds, actual = `Give I`)

muenchen_beta_results <- test %>%
  mutate(predicted = muenchen_beta_preds, actual = `Muenchen I`)

typhimurium_beta_results <- test %>%
  mutate(predicted = typhimurium_beta_preds, actual = Typhimurium)

#Get metrics
give_beta_metrics <- give_beta_results %>%
  metrics(truth = actual, estimate = predicted)

muenchen_beta_metrics <- muenchen_beta_results %>%
  metrics(truth = actual, estimate = predicted)

typhimurium_beta_metrics <- typhimurium_beta_results %>%
  metrics(truth = actual, estimate = predicted)

#Print metrics
give_beta_metrics
muenchen_beta_metrics
typhimurium_beta_metrics

# ---- prop-random-forest-fitting ----
#Cross-validation
give_prop_rf_results <- fit_resamples(give_regrf_wf, resamples = cv_folds, metrics = reg_metrics)
muenchen_prop_rf_results <- fit_resamples(muenchen_regrf_wf, resamples = cv_folds, metrics = reg_metrics)
rubislaw_prop_rf_results <- fit_resamples(rubislaw_regrf_wf, resamples = cv_folds, metrics = reg_metrics)
typhimurium_prop_rf_results <- fit_resamples(typhimurium_regrf_wf, resamples = cv_folds, metrics = reg_metrics)

#Collect metrics
give_prop_rf_metrics <- collect_metrics(give_prop_rf_results)
muenchen_prop_rf_metrics <- collect_metrics(muenchen_prop_rf_results)
rubislaw_prop_rf_metrics <- collect_metrics(rubislaw_prop_rf_results)
typhimurium_prop_rf_metrics <- collect_metrics(typhimurium_prop_rf_results)

#Compare metrics
give_prop_rf_metrics 
muenchen_prop_rf_metrics 
rubislaw_prop_rf_metrics
typhimurium_prop_rf_metrics

#Compare to null model 
null_mod <- null_model() %>%
  set_engine("parsnip") %>%
  set_mode("regression")

#Make recipes for the null model
give_null_recipe <- recipe(`Give I` ~ 1, data = train)
muenchen_null_recipe <- recipe(`Muenchen I` ~ 1, data = train)
rubislaw_null_recipe <- recipe(`Rubislaw` ~ 1, data = train)
typhimurium_null_recipe <- recipe(`Typhimurium` ~ 1, data = train)

#Create a workflow for the null model
give_null_workflow <- workflow() %>%
  add_recipe(give_null_recipe) %>%
  add_model(null_mod)

muenchen_null_workflow <- workflow() %>%
  add_recipe(muenchen_null_recipe) %>%
  add_model(null_mod)

rubislaw_null_workflow <- workflow() %>%
  add_recipe(rubislaw_null_recipe) %>%
  add_model(null_mod)

typhimurium_null_workflow <- workflow() %>%
  add_recipe(typhimurium_null_recipe) %>%
  add_model(null_mod)

#Fit null models with cross validation
set.seed(rngseed)
give_null_results <- fit_resamples(
  give_null_workflow,
  resamples = cv_folds,
  metrics = reg_metrics
)

muenchen_null_results <- fit_resamples(
  muenchen_null_workflow,
  resamples = cv_folds,
  metrics = reg_metrics
)

rubislaw_null_results <- fit_resamples(
  rubislaw_null_workflow,
  resamples = cv_folds,
  metrics = reg_metrics
)

typhimurium_null_results <- fit_resamples(
  typhimurium_null_workflow,
  resamples = cv_folds,
  metrics = reg_metrics
)

#Show metrics
collect_metrics(give_null_results)
collect_metrics(muenchen_null_results)
collect_metrics(rubislaw_null_results)
collect_metrics(typhimurium_null_results)

#Fit the random forest workflows to training data
give_regrf_fit <- fit(give_regrf_wf, data = train)
muenchen_regrf_fit <- fit(muenchen_regrf_wf, data = train)
rubislaw_regrf_fit <- fit(rubislaw_regrf_wf, data = train)
typhimurium_regrf_fit <- fit(typhimurium_regrf_wf, data = train)

# ---- logistic-regression-fitting ----
#Cross-validation
give_log_results <- fit_resamples(give_log_wf, resamples = cv_folds, metrics = class_metrics)
muenchen_log_results <- fit_resamples(muenchen_log_wf, resamples = cv_folds, metrics = class_metrics)
rubislaw_log_results <- fit_resamples(rubislaw_log_wf, resamples = cv_folds, metrics = class_metrics)
typhimurium_log_results <- fit_resamples(typhimurium_log_wf, resamples = cv_folds, metrics = class_metrics)

#Collect metrics
give_log_metrics <- collect_metrics(give_log_results)
muenchen_log_metrics <- collect_metrics(muenchen_log_results)
rubislaw_log_metrics <- collect_metrics(rubislaw_log_results)
typhimurium_log_metrics <- collect_metrics(typhimurium_log_results)

#Compare metrics
give_log_metrics 
muenchen_log_metrics 
rubislaw_log_metrics
typhimurium_log_metrics

#Compare to null model 
#Define null model
null_class_mod <- null_model() %>%
  set_engine("parsnip") %>%
  set_mode("classification")

#Make recipes for the null model
give_null_class_recipe <- recipe(`Give I Prev` ~ 1, data = train)
muenchen_null_class_recipe <- recipe(`Muenchen I Prev` ~ 1, data = train)
rubislaw_null_class_recipe <- recipe(`Rubislaw Prev` ~ 1, data = train)
typhimurium_null_class_recipe <- recipe(`Typhimurium Prev` ~ 1, data = train)

#Create a workflow for the null model
give_null_class_workflow <- workflow() %>%
  add_recipe(give_null_class_recipe) %>%
  add_model(null_class_mod)

muenchen_null_class_workflow <- workflow() %>%
  add_recipe(muenchen_null_class_recipe) %>%
  add_model(null_class_mod)

rubislaw_null_class_workflow <- workflow() %>%
  add_recipe(rubislaw_null_class_recipe) %>%
  add_model(null_class_mod)

typhimurium_null_class_workflow <- workflow() %>%
  add_recipe(typhimurium_null_class_recipe) %>%
  add_model(null_class_mod)

#Fit null models with cross validation
set.seed(rngseed)
give_null_class_results <- fit_resamples(
  give_null_class_workflow,
  resamples = cv_folds,
  metrics = class_metrics
)

muenchen_null_class_results <- fit_resamples(
  muenchen_null_class_workflow,
  resamples = cv_folds,
  metrics = class_metrics
)

rubislaw_null_class_results <- fit_resamples(
  rubislaw_null_class_workflow,
  resamples = cv_folds,
  metrics = class_metrics
)

typhimurium_null_class_results <- fit_resamples(
  typhimurium_null_class_workflow,
  resamples = cv_folds,
  metrics = class_metrics
)

#Show metrics
collect_metrics(give_null_class_results)
give_log_metrics 
collect_metrics(muenchen_null_class_results)
muenchen_log_metrics 
collect_metrics(rubislaw_null_class_results)
rubislaw_log_metrics
collect_metrics(typhimurium_null_class_results)
typhimurium_log_metrics

#All models performed better than the null models

#Fit logistic regression models to the training data
give_log_fit <- fit(give_log_wf, data = train)
muenchen_log_fit <- fit(muenchen_log_wf, data = train)
rubislaw_log_fit <- fit(rubislaw_log_wf, data = train)
typhimurium_log_fit <- fit(typhimurium_log_wf, data = train)

#Extract and summarize the fitted models
summary(pull_workflow_fit(give_log_fit)$fit)
summary(pull_workflow_fit(muenchen_log_fit)$fit)
summary(pull_workflow_fit(rubislaw_log_fit)$fit)
summary(pull_workflow_fit(typhimurium_log_fit)$fit)
 
#Make predictions from each model using the training data
give_log_train_preds <- predict(give_log_fit, new_data = train, type = "class")
give_log_train_prob_preds <- predict(give_log_fit, new_data = train, type = "prob")
muenchen_log_train_preds <- predict(muenchen_log_fit, new_data = train, type = "class")
muenchen_log_train_prob_preds <- predict(muenchen_log_fit, new_data = train, type = "prob")
rubislaw_log_train_preds <- predict(rubislaw_log_fit, new_data = train, type = "class")
rubislaw_log_train_prob_preds <- predict(rubislaw_log_fit, new_data = train, type = "prob")
typhimurium_log_train_preds <- predict(typhimurium_log_fit, new_data = train, type = "class")
typhimurium_log_train_prob_preds <- predict(typhimurium_log_fit, new_data = train, type = "prob")

#Combine predictions with actual values
give_results <- bind_cols(
  train,
  predictions = give_log_train_preds,
  prob_1 = give_log_train_prob_preds$.pred_1
)
muenchen_results <- bind_cols(
  train,
  predictions = muenchen_log_train_preds,
  prob_1 = muenchen_log_train_prob_preds$.pred_1
)
rubislaw_results <- bind_cols(
  train,
  predictions = rubislaw_log_train_preds,
  prob_1 = rubislaw_log_train_prob_preds$.pred_1
)
typhimurium_results <- bind_cols(
  train,
  predictions = typhimurium_log_train_preds,
  prob_1 = typhimurium_log_train_prob_preds$.pred_1
)

#Calculate accuracy and roc-auc
give_accuracy <- accuracy(give_results, truth = "Give I Prev", estimate = .pred_class)
give_roc_auc <- roc_auc(give_results, truth = "Give I Prev", prob_1)
muenchen_accuracy <- accuracy(muenchen_results, truth = "Muenchen I Prev", estimate = .pred_class)
muenchen_roc_auc <- roc_auc(muenchen_results, truth = "Muenchen I Prev", prob_1)
rubislaw_accuracy <- accuracy(rubislaw_results, truth = "Rubislaw Prev", estimate = .pred_class)
rubislaw_roc_auc <- roc_auc(rubislaw_results, truth = "Rubislaw Prev", prob_1)
typhimurium_accuracy <- accuracy(typhimurium_results, truth = "Typhimurium Prev", estimate = .pred_class)
typhimurium_roc_auc <- roc_auc(typhimurium_results, truth = "Typhimurium Prev", prob_1)

#Combine into df
log_training_results <- data.frame(
  Response = c("Give I", "Muenchen I", "Rubislaw", "Typhimurium"),
  Accuracy = c(give_accuracy$.estimate, muenchen_accuracy$.estimate, rubislaw_accuracy$.estimate, typhimurium_accuracy$.estimate),
  ROC_AUC = c(give_roc_auc$.estimate, muenchen_roc_auc$.estimate, rubislaw_roc_auc$.estimate, typhimurium_roc_auc$.estimate)
)

#Repeat for testing data
give_log_test_preds <- predict(give_log_fit, new_data = test, type = "class")
give_log_test_prob_preds <- predict(give_log_fit, new_data = test, type = "prob")
muenchen_log_test_preds <- predict(muenchen_log_fit, new_data = test, type = "class")
muenchen_log_test_prob_preds <- predict(muenchen_log_fit, new_data = test, type = "prob")
rubislaw_log_test_preds <- predict(rubislaw_log_fit, new_data = test, type = "class")
rubislaw_log_test_prob_preds <- predict(rubislaw_log_fit, new_data = test, type = "prob")
typhimurium_log_test_preds <- predict(typhimurium_log_fit, new_data = test, type = "class")
typhimurium_log_test_prob_preds <- predict(typhimurium_log_fit, new_data = test, type = "prob")

#Combine predictions with actual values
give_results_test <- bind_cols(
  test,
  predictions = give_log_test_preds,
  prob_1 = give_log_test_prob_preds$.pred_1
)
muenchen_results_test <- bind_cols(
  test,
  predictions = muenchen_log_test_preds,
  prob_1 = muenchen_log_test_prob_preds$.pred_1
)
rubislaw_results_test <- bind_cols(
  test,
  predictions = rubislaw_log_test_preds,
  prob_1 = rubislaw_log_test_prob_preds$.pred_1
)
typhimurium_results_test <- bind_cols(
  test,
  predictions = typhimurium_log_test_preds,
  prob_1 = typhimurium_log_test_prob_preds$.pred_1
)

#Calculate accuracy and roc-auc
give_accuracy_test <- accuracy(give_results_test, truth = "Give I Prev", estimate = .pred_class)
give_roc_auc_test <- roc_auc(give_results_test, truth = "Give I Prev", prob_1)
muenchen_accuracy_test <- accuracy(muenchen_results_test, truth = "Muenchen I Prev", estimate = .pred_class)
muenchen_roc_auc_test <- roc_auc(muenchen_results_test, truth = "Muenchen I Prev", prob_1)
rubislaw_accuracy_test <- accuracy(rubislaw_results_test, truth = "Rubislaw Prev", estimate = .pred_class)
rubislaw_roc_auc_test <- roc_auc(rubislaw_results_test, truth = "Rubislaw Prev", prob_1)
typhimurium_accuracy_test <- accuracy(typhimurium_results_test, truth = "Typhimurium Prev", estimate = .pred_class)
typhimurium_roc_auc_test <- roc_auc(typhimurium_results_test, truth = "Typhimurium Prev", prob_1)

#Combine into df
log_testing_results <- data.frame(
  Response = c("Give I", "Muenchen I", "Rubislaw", "Typhimurium"),
  Accuracy = c(give_accuracy_test$.estimate, muenchen_accuracy_test$.estimate, rubislaw_accuracy_test$.estimate, typhimurium_accuracy_test$.estimate),
  ROC_AUC = c(give_roc_auc_test$.estimate, muenchen_roc_auc_test$.estimate, rubislaw_roc_auc_test$.estimate, typhimurium_roc_auc_test$.estimate)
)

# ---- model-interpretation ----
# -- Logistic Regression Interpretation --
#Create explainers for logistic regression models
give_log_explainer <- explain_tidymodels(
  give_log_fit,
  data = train %>% select(-`Give I Prev`), #Exclude response variable from data
  y = as.numeric(train$`Give I Prev`), #Provide the response variable
  label = "Give I Logistic Regression"
)

muenchen_log_explainer <- explain_tidymodels(
  muenchen_log_fit,
  data = train %>% select(-`Muenchen I Prev`),
  y = as.numeric(train$`Muenchen I Prev`),
  label = "Muenchen I Logistic Regression"
)

rubislaw_log_explainer <- explain_tidymodels(
  rubislaw_log_fit,
  data = train %>% select(-`Rubislaw Prev`),
  y = as.numeric(train$`Rubislaw Prev`),
  label = "Rubislaw Logistic Regression"
)

typhimurium_log_explainer <- explain_tidymodels(
  typhimurium_log_fit,
  data = train %>% select(-`Typhimurium Prev`),
  y = as.numeric(train$`Typhimurium Prev`),
  label = "Typhimurium Logistic Regression"
)

# -- Beta Regression Interpretation --
# Create explainers for beta regression models
give_beta_explainer <- explain(
  model = give_beta_model,
  data = train %>% select(-`Give I`),  #Exclude response variable from data
  y = train$`Give I`,                 #Provide the response variable
  label = "Give I Beta Regression"
)

muenchen_beta_explainer <- explain(
  model = muenchen_beta_model,
  data = train %>% select(-`Muenchen I`),
  y = train$`Muenchen I`,
  label = "Muenchen I Beta Regression"
)

typhimurium_beta_explainer <- explain(
  model = typhimurium_beta_model,
  data = train %>% select(-Typhimurium),
  y = train$Typhimurium,
  label = "Typhimurium Beta Regression"
)

# -- Random Forest Interpretation --
# Create explainers for random forest models
give_rf_explainer <- explain_tidymodels(
  give_regrf_fit,
  data = train %>% select(-`Give I`),
  y = train$`Give I`,
  label = "Give I Random Forest"
)

muenchen_rf_explainer <- explain_tidymodels(
  muenchen_regrf_fit,
  data = train %>% select(-`Muenchen I`),
  y = train$`Muenchen I`,
  label = "Muenchen I Random Forest"
)

rubislaw_rf_explainer <- explain_tidymodels(
  rubislaw_regrf_fit,
  data = train %>% select(-Rubislaw),
  y = train$Rubislaw,
  label = "Rubislaw Random Forest"
)

typhimurium_rf_explainer <- explain_tidymodels(
  typhimurium_regrf_fit,
  data = train %>% select(-Typhimurium),
  y = train$Typhimurium,
  label = "Typhimurium Random Forest"
)

# -- Feature Importance --
# Logistic Regression
vip(give_log_fit, geom = "point") + ggtitle("Variable Importance: Give I Logistic Regression")
vip(muenchen_log_fit, geom = "point") + ggtitle("Variable Importance: Muenchen I Logistic Regression")
vip(rubislaw_log_fit, geom = "point") + ggtitle("Variable Importance: Rubislaw Logistic Regression")
vip(typhimurium_log_fit, geom = "point") + ggtitle("Variable Importance: Typhimurium Logistic Regression")

# Beta Regression
vip(give_beta_model, geom = "point") + ggtitle("Variable Importance: Give I Beta Regression")
vip(muenchen_beta_model, geom = "point") + ggtitle("Variable Importance: Muenchen I Beta Regression")
vip(typhimurium_beta_model, geom = "point") + ggtitle("Variable Importance: Typhimurium Beta Regression")

# Random Forest
vip(give_regrf_fit, geom = "point") + ggtitle("Variable Importance: Give I Random Forest")
vip(muenchen_regrf_fit, geom = "point") + ggtitle("Variable Importance: Muenchen I Random Forest")
vip(rubislaw_regrf_fit, geom = "point") + ggtitle("Variable Importance: Rubislaw Random Forest")
vip(typhimurium_regrf_fit, geom = "point") + ggtitle("Variable Importance: Typhimurium Random Forest")

# -- Residual Analysis --
# Logistic Regression
plot(model_performance(give_log_explainer)) + ggtitle("Residuals: Give I Logistic Regression")
plot(model_performance(muenchen_log_explainer)) + ggtitle("Residuals: Muenchen I Logistic Regression")
plot(model_performance(rubislaw_log_explainer)) + ggtitle("Residuals: Rubislaw Logistic Regression")
plot(model_performance(typhimurium_log_explainer)) + ggtitle("Residuals: Typhimurium Logistic Regression")

# Beta Regression
plot(model_performance(give_beta_explainer)) + ggtitle("Residuals: Give I Beta Regression")
plot(model_performance(muenchen_beta_explainer)) + ggtitle("Residuals: Muenchen I Beta Regression")
plot(model_performance(typhimurium_beta_explainer)) + ggtitle("Residuals: Typhimurium Beta Regression")

# Random Forest
plot(model_performance(give_rf_explainer)) + ggtitle("Residuals: Give I Random Forest")
plot(model_performance(muenchen_rf_explainer)) + ggtitle("Residuals: Muenchen I Random Forest")
plot(model_performance(rubislaw_rf_explainer)) + ggtitle("Residuals: Rubislaw Random Forest")
plot(model_performance(typhimurium_rf_explainer)) + ggtitle("Residuals: Typhimurium Random Forest")