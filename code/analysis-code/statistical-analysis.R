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
reg_metrics <- metric_set(rmse, rsq)

#Define metrics for classification evaluation
class_metrics <- metric_set(accuracy, roc_auc)

#Workflows for all models
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
#Initialize a data frame to store summary results for all models
model_results <- data.frame(
  Response = character(),
  Term = character(),
  Coefficient = numeric(),
  Std_Error = numeric(),
  P_Value = numeric(),
  stringsAsFactors = FALSE
)

#Fit models and extract results
for (response_name in names(response_formulas)) {
  #Get the formula for the current response variable
  current_formula <- response_formulas[[response_name]]
  
  #Fit the beta regression model on the training data
  beta_model <- betareg(current_formula, data = train)
  
  #Extract model results (coefficients, standard errors, and p-values)
  tidy_model <- tidy(beta_model) %>%
    select(term, estimate, std.error, p.value) %>%
    rename(
      Coefficient = estimate,
      Std_Error = std.error,
      P_Value = p.value
    )
  
  #Add the response variable name to the results
  tidy_model <- tidy_model %>%
    mutate(Response = response_name)
  
  #Append the results to the main results data frame
  model_results <- bind_rows(model_results, tidy_model)
}

#Create a publication-ready summary table
kable(model_results, format = "markdown", digits = 3, caption = "Model Results")

give_beta_model <- betareg(response_formulas$`Give I`, data = train)
summary(give_beta_model)
vif(give_beta_model)

give_beta_summary <- tidy(give_beta_model) %>%
  select(term, estimate, std.error, p.value) %>%
  rename(Coefficent = estimate,
         Std_Error = std.error,
         p_value = p.value) %>%
  mutate(Response = "Give I")

muenchen_beta_model <- betareg(response_formulas$`Muenchen I`, data = train)
summary(muenchen_beta_model)
vif(muenchen_beta_model)

muenchen_beta_summary <- tidy(muenchen_beta_model) %>%
  select(term, estimate, std.error, p.value) %>%
  rename(Coefficent = estimate,
         Std_Error = std.error,
         p_value = p.value) %>%
  mutate(Response = "Muenchen I")

typhimurium_beta_model <- betareg(response_formulas$Typhimurium, data = train)
summary(typhimurium_beta_model)
vif(typhimurium_beta_model)

typhimurium_beta_summary <- tidy(typhimurium_beta_model) %>%
  select(term, estimate, std.error, p.value) %>%
  rename(Coefficent = estimate,
         Std_Error = std.error,
         p_value = p.value) %>%
  mutate(Response = "Typhimurium")

beta_model_results <- rbind(give_beta_summary, muenchen_beta_summary, typhimurium_beta_summary) %>%
  select(Response, term, Coefficent, Std_Error, p_value) %>%
  filter(p_value <= 0.05)

beta_model_results

#Refit each model using significant predictors
final_give_beta_model <- betareg(`Give I` ~ `Total Rain(in)_std` + System + `Avg Relative Humidity(%)_std` +
                                  `Total Solar Radiation(MJ/m^2)_std` + Interaction_Humidity_Rain + 
                                   Interaction_Solar_Rain, data = train)

final_muenchen_beta_model <- betareg(`Muenchen I` ~ 1, data = train) #Just a null model :/

final_typhimurium_model <- betareg(Typhimurium ~ `Total Rain(in)_std`, data = train)

#Re-do summaries
give_beta_summary <- tidy(give_beta_model) %>%
  select(term, estimate, std.error, p.value) %>%
  rename(Coefficent = estimate,
         Std_Error = std.error,
         p_value = p.value) %>%
  mutate(Response = "Give I")

muenchen_beta_summary <- tidy(muenchen_beta_model) %>%
  select(term, estimate, std.error, p.value) %>%
  rename(Coefficent = estimate,
         Std_Error = std.error,
         p_value = p.value) %>%
  mutate(Response = "Muenchen I")

typhimurium_beta_summary <- tidy(typhimurium_beta_model) %>%
  select(term, estimate, std.error, p.value) %>%
  rename(Coefficent = estimate,
         Std_Error = std.error,
         p_value = p.value) %>%
  mutate(Response = "Typhimurium")

beta_model_results <- rbind(give_beta_summary, muenchen_beta_summary, typhimurium_beta_summary) %>%
  select(Response, term, Coefficent, Std_Error, p_value) %>%
  filter(p_value <= 0.05)

beta_model_results

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
                    typhimurium_beta_model = typhimurium_beta_model, 
                    muenchen_beta_model = muenchen_beta_model)

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

#Models did not perform better than null models

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

