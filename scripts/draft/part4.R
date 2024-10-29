setwd('/Users/henrylouie/module1-group-11/data')
library(tidyverse)
library(randomForest)
library(tidymodels)
library(ggplot2)
library(yardstick)
library(glmnet)
library(modelr)
library(rsample)

set.seed(102022)

## LASSO:
#set up data into training and testing
biomarker <- biomarker_clean %>% 
  select(-ados) %>% 
  mutate(class = as.numeric(group == 'ASD'))

partitions <- biomarker %>% 
  initial_split(prop = 0.8)

bio_x_train <- training(partitions) %>% 
  select(-c(group, class)) %>% 
  as.matrix()

bio_y_train <- training(partitions) %>% 
  pull(class)

# using glmnet to find optimal lambda for LASSO
cv_out <- cv.glmnet(bio_x_train,bio_y_train,
                    family = 'binomial',
                    nfolds = 5,
                    type.measure = 'deviance')
cvout_df <- tidy(cv_out)
plot(cv_out)

fit <- glmnet(bio_x_train, bio_y_train, family = 'binomial')
fit_df <- tidy(fit)
plot(fit)

#choosing lambda for final model
lambda <- exp(-1.8)
best_model <- glmnet(bio_x_train,bio_y_train, family = 'binomial', lambda = lambda)
best_model_df <- tidy(best_model)

#selecting proteins for panel
protein_select <- best_model_df %>% 
  filter(term != '(Intercept)') %>% 
  pull(term)
protein_select

#rebuild training and testing with new panel
biomarker_model <- biomarker_clean %>% 
  select(group, any_of(protein_select)) %>% 
  mutate(class = as.numeric(group == 'ASD')) %>% 
  select(-group)

set.seed(46340)
bio_partitions <- biomarker_model %>% 
  initial_split(prop = 0.8)

bio_train <- training(bio_partitions)

final_fit <- glm(class ~., data = bio_train,
                 family = 'binomial')

#Using yardstick package to find accuracy of model with the 5 protein panel
class_metric <- metric_set(accuracy)

pred_df <- testing(bio_partitions) %>% 
  add_predictions(final_fit, type = 'response') %>% 
  mutate(pred.class = (pred > 0.5),
         group = factor(class, labels = c('TD', 'ASD')),
         pred.group = factor(pred.class, labels = c('TD', 'ASD'))) 

pred_df %>% class_metric(truth = group, estimate = pred.group, event_level = 'second')
#Result is 0.774 accuracy


## Elastic Net
set.seed(795843)

en_partitions <- biomarker %>% 
  initial_split(prop = 0.8)
en_x <- training(en_partitions) %>% 
  select(-c(group,class)) %>% 
  as.matrix()
en_y <- training(en_partitions) %>% 
  pull(class)

en_cv <- cv.glmnet(en_x, en_y, alpha = 0.75, family = 'binomial')
en_lam = en_cv$lambda.min

plot(en_cv$glmnet.fit, 'lambda', label = F)

en_model <- glmnet(en_x,en_y, family = 'binomial', lambda = en_lam)
en_model_df <- tidy(en_model)

en_prot_select <- en_model_df %>% 
  filter(term != '(Intercept)') %>% 
  pull(term)
en_prot_select

en_biomarker_model <- biomarker_clean %>% 
  select(group, any_of(protein_select)) %>% 
  mutate(class = as.numeric(group == 'ASD')) %>% 
  select(-group)
en_final_part <- en_biomarker_model %>% 
  initial_split(prop = 0.8)

en_fit <- glm(class~., data = training(en_final_part), family = 'binomial')

en_pred_df <- testing(en_final_part) %>% 
  add_predictions(en_fit, type = 'response') %>% 
  mutate(pred.class = (pred > 0.5),
         group = factor(class, labels = c('TD', 'ASD')),
         pred.group = factor(pred.class, labels = c('TD', 'ASD'))) 

en_pred_df %>% class_metric(truth = group, estimate = pred.group, event_level = 'second')

# accuracy of 0.677