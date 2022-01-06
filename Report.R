### This script contains all the R code used to create the analyses presented in the final report ---
## The report was submitted as the Capstone project for the Harvard X data science professional certificate course series
## Author: Tero Jalkanen

## packages --------------


## Library for plotting small molecule chemical structure
if(!require(ChemmineR)){
  if(!require(BiocManager))  install.packages("BiocManager", repos = "http://cran.us.r-project.org")
  BiocManager::install("ChemmineR", force = TRUE)}

## Data wrangling and visualization
if(!require(tidyverse)) install.packages("tidyverse", repos = "http://cran.us.r-project.org")
#### Includes: ggplot2 (data visualisation), dplyr (data manipulation), tidyr (data tidying), readr (data import), purrr (functional programming), tibble (tibbles, a modern re-imagining of data frames), stringr (strings), forcats (factors), readxl (.xls and .xlsx sheets), rvest (web scraping), lubridate (dates and date-times), broom (turning models into tidy data)
if(!require(gridExtra)) install.packages("gridExtra", repos = "http://cran.us.r-project.org") # For multiple plots
#if(!require(cowplot)) install.packages("cowplot", repos = "http://cran.us.r-project.org") # For multiple plots
if(!require(corrplot)) install.packages("corrplot", repos = "http://cran.us.r-project.org") # Visualizing correlation
# DiagrammeR only works for html output
#if(!require(DiagrammeR)) install.packages("DiagrammeR", repos = "http://cran.us.r-project.org") # Building diagrams with R

## Modelling
if(!require(caret)) install.packages("caret", repos = "http://cran.us.r-project.org")
if(!require(tidymodels)) install.packages("tidymodels", repos = "http://cran.us.r-project.org")
if(!require(MASS)) install.packages("MASS", repos = "http://cran.us.r-project.org")
if(!require(discrim)) install.packages("discrim", repos = "http://cran.us.r-project.org") #LDA and QDA
if(!require(ranger)) install.packages("ranger", repos = "http://cran.us.r-project.org") # Random forest
if(!require(glmnet)) install.packages("glmnet", repos = "http://cran.us.r-project.org") #Ridge and Lasso
if(!require(rpart)) install.packages("rpart", repos = "http://cran.us.r-project.org") #CART
if(!require(vip)) install.packages("vip", repos = "http://cran.us.r-project.org") #variable importance for trees
if(!require(xgboost)) install.packages("xgboost", repos = "http://cran.us.r-project.org") # boosted trees
#if(!require()) install.packages("", repos = "http://cran.us.r-project.org")

## Typography
if(!require(bookdown)) install.packages("bookdown", repos = "http://cran.us.r-project.org") #Fig. references etc.
#if(!require(kableExtra)) install.packages("kableExtra", repos = "http://cran.us.r-project.org") # Nicer table formatting
if(!require(equatiomatic)) install.packages("equatiomatic", repos = "http://cran.us.r-project.org") # automatic way to show model equations

## Set graphical theme
theme_set(theme_bw())

## Download data -------
#download modified data, the final dataset
df <- read.csv(file = "data/final_dataset.csv")

# One example molecule
molecule <- ChemmineR::pubchemCidToSDF(as.numeric(df$pubchem[3]))
plot(molecule[1], print = FALSE, print_cid = "")

### EDA --------

#Histogram of the retention times
hist1 <- df %>% ggplot(aes(x = rt)) +
  geom_histogram(fill = "gray", color = "black") +
  ggtitle("a)") +
  labs(x = "Retention time (s)", y = "Number of molecules #") +
  geom_vline(xintercept = 300, color = "red", lty = 2) #not retained molecules below this line

#Histogram with logarithmic y-axis
hist2 <- df %>% 
  ggplot(aes(x = rt)) +
  geom_histogram(aes(fill = rt<300), color = "black") + # Fill color by retainment status
  ggtitle("b)") +
  labs(x = "Retention time (s)", y = "Number of molecules #") +
  geom_vline(xintercept = 300, color = "red", lty = 2) + #not retained molecules below this line
  scale_y_log10()

grid.arrange(hist1, 
             hist2, 
             ncol = 2,
             widths = c(2,3))

rm(hist1, hist2)

# Code to create data splitting diagram
# output only works with html but not pdf
DiagrammeR::grViz("digraph {
  graph [layout = dot, rankdir = TB]

  node [shape = rectangle]
  rec1 [label = 'Original data']
  rec3 [label =  'Train data']
  rec2 [label = 'Test data']
  rec4 [label = 'Fold 1']
  rec5 [label = 'Fold 2 ... Fold k-1']
  rec6 [label = 'Fold k']

  # edge definitions with the node IDs
  rec1 -> rec2
  rec1 -> rec3 -> rec4
  rec3 -> rec5
  rec3 -> rec6
  }",
  height = 300)

## Predictors/features can be divided into subcategories
### below we define groups of different kinds of predictors

# Number of chemical elements in each molecule
elements <- c("C", "H", "N", "S", "Cl", "O", "F", "I", "Br", "Si", "P")

# Number of chemical groups and special atoms in each molecule
chem_groups <- c("RNH2", "R2NH", "R3N", "ROPO3", "ROH", "RCHO", "RCOR", "RCOOH", "RCOOR", "ROR", "RCCH", "RCN", "RINGS", "AROMATIC")

#Values for physico chemical properties
phys_chem <- c("MW", "Ncharges", "ALogP", "ALogp2", "AMR", "apol", "naAromAtom", "nAromBond", "TopoPSA", "fragC", "nHBAcc", "nHBDon", "nRotB", "VABC", "Zagreb", "ECCEN", "WPATH", "WPOL")

# Boxplot of chemical element distributions
# Not sure if this is the best choise for discrete value variables, but can't think of a better choise
box1 <- df %>% dplyr::select(tidyselect::all_of(elements)) %>% 
  pivot_longer(cols = tidyselect::all_of(elements), names_to = "element_name", values_to = "atom_counts") %>% 
  ggplot(aes(x = element_name, y = atom_counts)) +
  geom_boxplot() +
  labs(x = "Chemical element", y = "Number of atoms in a molecule") +
  ggtitle("a)")

# Violin plot to display carbon and hydrogen
violin1 <- df %>% dplyr::select(C, H) %>% 
  pivot_longer(cols = c(C, H), names_to = "element_name", values_to = "atom_counts") %>% 
  ggplot(aes(x = element_name, y = atom_counts)) +
  geom_violin(fill = "gray") +
  labs(x = "Chemical element", y = "Number of atoms in a molecule") +
  ggtitle("b)")

grid.arrange(box1, 
             violin1, 
             ncol = 2,
             widths = c(3,1))


rm("box1", "violin1")


#The outcome and ID variables
df %>% dplyr::select(rt, pubchem, MF) %>% 
  head() %>%
  rename("Retention time (s)" = rt, "PubChem ID" = pubchem, "Molecular Formula" = MF) %>% 
  knitr::kable(caption = "The outcome variable and the two ID variables") 

## Boxplot not very good since values are counts
# df %>% select(all_of(chem_groups)) %>% 
#   pivot_longer(cols = all_of(chem_groups), names_to = "group_name", values_to = "group_counts") %>% 
#   ggplot(aes(x = group_name, y = group_counts)) +
#   geom_boxplot() +
#   labs(x = "Chemical group", y = "Number of groups in a molecule") +
#   ggtitle("Distributions of chemical groups in the dataset") +
#   coord_flip()

#Dot plot, where the size is proportional to number of observations
df %>% dplyr::select(all_of(chem_groups)) %>% # get chemical groups
  pivot_longer(cols = tidyselect::all_of(chem_groups), names_to = "group_name", values_to = "group_counts") %>% 
  group_by(group_name, group_counts) %>% # Group by name of the group and amount in a molecule
  summarise(n = n()) %>%  # Summarize number of observations
  ggplot(aes(x = group_name, y = group_counts)) +
  geom_point(aes(size = n, color = n)) +
  labs(x = "Chemical group", y = "Number of groups in a molecule") +
  coord_flip()

#Boxplot of phys chem properties

df %>% dplyr::select(tidyselect::all_of(phys_chem)) %>% 
  pivot_longer(cols = tidyselect::all_of(phys_chem), names_to = "prop_name", values_to = "prop_value") %>% 
  # Lets take properties with large values and put them to a separate facet
  mutate(value_type = if_else(prop_name %in% c("WPATH", "fragC", "ECCEN"), "Large values", "Small values")) %>% 
  ggplot(aes(x = prop_name, y = prop_value)) +
  geom_boxplot() + coord_flip() +
  # Free scales to allow for different values of property value scale
  facet_grid(. ~value_type, scales = "free") +
  labs(y= "Property value", x = "")

# Count values
pc_count_vars <- c("nRotB", "nHBDon", "nHBAcc", "Ncharges", "nAromBond", "naAromAtom")

df %>% dplyr::select(tidyselect::all_of(pc_count_vars)) %>% 
  pivot_longer(cols = tidyselect::all_of(pc_count_vars), 
               names_to = "prop_name", 
               values_to = "prop_value") %>%
  ggplot(aes(x = prop_value)) +
  geom_histogram(fill = "red", color = "black", binwidth = 1) +
  # Facet each count variable into individual pane
  facet_wrap(facets = "prop_name") +
  scale_y_log10() +
  labs(x = "Property value", y = "Number of observations") +
  # Print the largest observation value as text on each facet
  geom_text(data = df %>% dplyr::select(tidyselect::all_of(pc_count_vars)) %>% 
              pivot_longer(cols = tidyselect::all_of(pc_count_vars), 
                           names_to = "prop_name", 
                           values_to = "prop_value") %>%
              group_by(prop_name) %>% 
              summarise(range_high = range(prop_value)[2]), #high end of range
            aes(x = 20, y = 60000, label = paste("Largest value:", range_high))
  )

#Create a pretty correlogram here
#Build correlogram
M <- df %>% dplyr::select(-MF, -pubchem) %>% cor()

# rename rows and columns for better readability
original_rownames <- row.names(M)
#Let's change column names into numbers to make the graph less cluttered
colnames(M) <- c("rt", 2:dim(M)[2])
# Add numbers to row names as well
row.names(M) <- paste0(row.names(M), paste(c("", rep(",", times = 43)), c("", 2:44)))

corrplot(M, method = 'ellipse', tl.cex = 0.7, type = "upper")

### Visualizing not_retained molecules -----------

#Box-plot logP by class
box_logp <- df %>% 
  mutate(not_retained = (rt<300)) %>%
  ggplot(aes(y = ALogP, x = not_retained)) +
  geom_boxplot() +
  ggtitle("a)") +
  labs(x = "Molecule not retained")

#Box-plot MW by class
box_MW <- df %>% 
  mutate(not_retained = (rt<300)) %>%
  ggplot(aes(y = MW, x = not_retained)) +
  geom_boxplot() +
  ggtitle("b)") +
  labs(x = "Molecule not retained", y = "Molecular weight")

#Plots side-by-side
grid.arrange(box_logp, 
             box_MW, 
             ncol = 2
)

#remove plot objects
rm("box_logp", "box_MW")

### PCA fit -----------------

# Add classification variable
df2 <- df %>%
  mutate(not_retained = rt<300) %>%
  dplyr::select(-rt, -pubchem, - MF) # Remove outcome and ID variables

#Fit PCA for predictors
pca_fit <- prcomp(x = dplyr::select(df2, -not_retained), scale = TRUE)

#plot(pca_fit)

# Biplot with the two classes
biplot_1 <-
  ggplot(data = data.frame(PC1 = pca_fit$x[,1], PC2 = pca_fit$x[,2], not_retained = df2$not_retained),
       aes(x = PC1, y = PC2, color = not_retained)) +
  geom_point(alpha = 0.5) +
  labs(x = "Principal component 1", y = "Principal component 2") +
  theme(legend.position = "bottom") +
  ggtitle("a)")

#Proportion of variance explained
pr.var <- pca_fit$sdev^2
#pr.var
#proportion of variance explained
pve <- pr.var / sum(pr.var)

#Cumulative proportion of variance explained
var_plot <-
  ggplot(data = data.frame(cum_var = cumsum(pve), PC = 1:length(pve)),
       aes(x = PC, y = cum_var)) +
  geom_line() +
  geom_point(size = 2, color = "red") +
  labs(x = "Principal component", y = "Variance explained") +
  ggtitle("b)")

#Plots side-by-side
pca_plots <- grid.arrange(biplot_1,
             var_plot,
             ncol = 2
             )

ggsave("pca_plots.png", plot = pca_plots, dpi = 150, device = "png", width = 7, height = 3, units = "in")

rm("df2", "pca_fit", "biplot_1", "var_plot", "pve", "pr.var", "pca_plots")

## Visualizing not_retained molecules -----------

scatter1 <- df %>% mutate(not_retained = (rt<300)) %>%
  ggplot(aes(y = rt, x = ALogP)) +
  geom_point(aes(color = not_retained)) +
  labs(y = "Retention time (s)")

ggsave("scatterplot.png", plot = scatter1, dpi = 150, device = "png", width = 4, height = 2.5, units = "in")

rm("scatter1")

#### Classification -----------------------

### Test/Train split

set.seed(1234)
#split to test train 25:75
classification_split <- initial_split(data = df %>% mutate(not_retained = (rt<300) %>% as.factor()), 
                                      prop = 3/4, 
                                      strata = not_retained)


#Create test and train
df_test <- testing(classification_split)
df_train <- training(classification_split)

### Create 10-fold cross-validation dataset for train

train_folds <- vfold_cv(df_train, v = 10, strata = not_retained)


# How are the retained vs non-retained molecules divided?
test_train_class <- df_test %>%  
  count(not_retained) %>% 
  mutate(prop = n/sum(n), set = "test") %>% 
  bind_rows(
    df_train %>%  
      count(not_retained) %>% 
      mutate(prop = n/sum(n), set = "train")
  )

# Print amount of molecules as table
test_train_class %>% 
  knitr::kable(caption = "The amount and proportion of observations in different classes for the test and train sets")

rm("test_train_class")

# all predictors


# recipe for classification drop retention time as "ID variable"
class_recipe_all <-
  recipe(not_retained ~ ., data = df_train) %>% #retention is the outcome we want to predict
  update_role(rt, pubchem, MF, new_role = "ID") %>%  # rt, pubchem id, and molecular formula as ID variables (not predictors)
  step_zv(all_predictors()) # remove zero variance variables (if any)

## Logistic regression model specification
lr_model_spec <-
  logistic_reg() %>%
  set_engine("glm") %>%
  set_mode("classification")

# Workflow
class_wflow <-
  workflow() %>%
  add_model(lr_model_spec) %>%
  add_recipe(class_recipe_all)




#fit logistic regression using tidymodels all predictors
class_lr_fit1 <-
  class_wflow %>%
  fit_resamples(resamples = train_folds)

# Let's collect metrics from our 10-fold cross-validation
# collect_metrics(class_lr_fit1) %>%
#   knitr::kable(caption = "Metrics for Logistic regression model with all predictors")


## compare models in tabular format
model_comp <-
  collect_metrics(class_lr_fit1) %>% mutate(model = "Logistic", complexity = "all vars")


#Function for adding new model metrics to comparison table
add_model_comp <- function(models = model_comp, new_fit, model, complexity){
  new_comparison <- models %>% 
    bind_rows(collect_metrics(new_fit) %>% 
                mutate(model = model, complexity = complexity)
    )
  return(new_comparison)
}

# Simple model with just two predictors
# New simple recipe:

class_recipe_simple <-
  recipe(not_retained ~ MW + ALogP, data = df_train)

class_wflow2 <-
  workflow() %>%
  add_model(lr_model_spec) %>%
  add_recipe(class_recipe_simple)


#10-fold cv with simple LogReg-model
class_lr_fit2 <-
  class_wflow2 %>%
  fit_resamples(resamples = train_folds)

# Let's collect metrics from our 10-fold cross-validation
model_comp <- add_model_comp(models = model_comp,
                             new_fit =class_lr_fit2,
                             model = "Logistic",
                             complexity = "simple"
                            )

# # Halfway model (takeaway the atom count variables)
# 
# Number of chemical elements in each molecule
elements <- c("C", "H", "N", "S", "Cl", "O", "F", "I", "Br", "Si", "P")
# Variables which are highly correlated with other variables
high_corr <- c("ALogp2", "apol", "ECCEN", "WPOL")

#New model recipe
class_recipe_intermed <-
  recipe(not_retained ~ ., data = df_train) %>% #retention is the outcome we want to predict
  update_role(rt, pubchem, MF, new_role = "ID") %>%  # rt, pubchem id, and molecular formula as ID variables (not predictors)
  update_role(all_of(elements), new_role = "ID") %>% #Make chemical elements just ID variables
  update_role(all_of(high_corr), new_role = "ID") %>% # Make some highly correlated variables ID variables
  step_zv(all_predictors()) # remove zero variance variables (if any)

class_wflow3 <-
  workflow() %>%
  add_model(lr_model_spec) %>%
  add_recipe(class_recipe_intermed)


#10-fold cv with simple LogReg-model
class_lr_fit3 <-
  class_wflow3 %>%
  fit_resamples(resamples = train_folds)

# Let's collect metrics from our 10-fold cross-validation
model_comp <- add_model_comp(models = model_comp,
                             new_fit =class_lr_fit3,
                             model = "Logistic",
                             complexity = "intermediate"
                            )

#Let's write the metrics to a csv-file, so that we save time when knitting several times
write.csv(model_comp, file = "classificationmetrics1.csv", row.names = FALSE)


# Read cross-validation results from a file
model_comp <- read.csv("classificationmetrics1.csv") %>% as.tibble()

model_comp %>% 
  dplyr::select(-.config) %>% 
  knitr::kable(caption = "Metrics for different logistic regression models were calculated based on 10-fold cross-validation", digits = 3)

set.seed(111)

### Linear discriminant analysis ----------------------------

## model specification

lda_model_spec <- discrim_linear() %>%
  set_mode("classification") %>%
  set_engine("MASS")

# Workflow for simple LDA model
lda_wf_1 <-
  workflow() %>%
  add_model(lda_model_spec) %>%
  add_recipe(class_recipe_simple)

# Fit simple model:
#10-fold cv with simple LDA-model
class_lda_fit1 <-
  lda_wf_1 %>%
  fit_resamples(resamples = train_folds)

# Add metrics to model comparison
model_comp <- add_model_comp(new_fit = class_lda_fit1, model = "LDA", complexity = "simple")

# Workflow for intermediate LDA model
lda_wf_2 <-
  workflow() %>%
  add_model(lda_model_spec) %>%
  add_recipe(class_recipe_intermed)

# Fit intermediate model:
#10-fold cv with LDA-model
class_lda_fit2 <-
  lda_wf_2 %>%
  fit_resamples(resamples = train_folds)

# Add metrics to model comparison
model_comp <- add_model_comp(new_fit = class_lda_fit2, model = "LDA", complexity = "intermediate")


# Workflow for LDA model with all predictors
lda_wf_3 <-
  workflow() %>%
  add_model(lda_model_spec) %>%
  add_recipe(class_recipe_all)

# Fit full model:
#10-fold cv with full LDA-model
class_lda_fit3 <-
  lda_wf_3 %>%
  fit_resamples(resamples = train_folds)

# Add metrics to model comparison
model_comp <- add_model_comp(new_fit = class_lda_fit3, model = "LDA", complexity = "all vars")

## Quadratic discriminant analysis -------------------

## model specification
qda_model_spec <- discrim_quad() %>%
  set_mode("classification") %>%
  set_engine("MASS")

# Workflow for simple QDA model
qda_wf_1 <-
  workflow() %>%
  add_model(qda_model_spec) %>%
  add_recipe(class_recipe_simple)

# Fit simple model:
#10-fold cv with simple QDA-model
class_qda_fit1 <-
  qda_wf_1 %>%
  fit_resamples(resamples = train_folds)

# Add metrics to model comparison
model_comp <- add_model_comp(new_fit = class_qda_fit1, model = "QDA", complexity = "simple")

# Workflow for intermediate QDA model
qda_wf_2 <-
  workflow() %>%
  add_model(qda_model_spec) %>%
  add_recipe(class_recipe_intermed)

# Fit intermediate model:
#10-fold cv with QDA-model
class_qda_fit2 <-
  qda_wf_2 %>%
  fit_resamples(resamples = train_folds) #Fails for two folds due to rank deficiency

# Add metrics to model comparison
model_comp <- add_model_comp(new_fit = class_qda_fit2, model = "QDA", complexity = "intermediate")

## Random forest -----------------

# Model specification
rf_model_spec <-
  rand_forest(trees = 1000) %>%
  set_engine("ranger") %>%
  set_mode("classification")

## Workflow for random forest with all variables

rf_wf <- workflow() %>%
  add_model(rf_model_spec) %>%
  add_recipe(class_recipe_all)


# #Simple model test
# rf_fit1 <- rf_wf %>% fit(df_train)
#
# #Predictions
# rf_preds <- rf_fit1$fit$fit$fit$predictions %>% as.data.frame()
# names(rf_preds) <- c("p_retained", "p_not_retained")
#
# rf_preds <- rf_preds %>% mutate(pred = if_else(p_not_retained >= 0.5, TRUE, FALSE) %>% as.factor())
#
# caret::confusionMatrix(rf_preds$pred, df_train$not_retained)

#10-fold cross-validation
rf_fit <-
  rf_wf %>%
  fit_resamples(resamples = train_folds)

# Add metrics to model comparison
model_comp <- add_model_comp(new_fit = rf_fit, model = "Random Forest", complexity = "all vars")

#Write metrics to a csv-file
write.csv(model_comp, file = "classificationmetrics2.csv", row.names = FALSE)


#Read metrics from file
model_comp <- read.csv("classificationmetrics2.csv") %>% as.tibble()

## Compare classification models -----------------


# Baseline prediction performance for accuracy = all molecules not retained
bl_pred <- sum(df$rt>300)/length(df$rt)
## Compare models

#accuracy
model_comp %>% mutate(model_t = paste(model, complexity, sep = " ")) %>% 
  filter(.metric == "accuracy") %>% 
  ggplot(aes(x = model_t, y = mean, color = (mean > bl_pred))) +
  geom_pointrange(aes(ymin = mean - std_err, ymax = mean + std_err)) +
  geom_hline(yintercept = bl_pred, color = "red", lty = 2) +
  labs(y = "Accuracy", x = "Model & Complexity") +
  scale_color_discrete(name="Accuracy above baseline:") +
  theme(legend.position = "top") +
  coord_flip()

#ROC auc
model_comp %>% mutate(model_t = paste(model, complexity, sep = " ")) %>% 
  filter(.metric == "roc_auc") %>% 
  ggplot(aes(x = model_t, y = mean)) +
  geom_pointrange(aes(ymin = mean - std_err, ymax = mean + std_err)) +
  labs(y = "ROC AUC", x = "Model & Complexity") +
  coord_flip()

#Tabular comparison for two best models
model_comp %>% 
  dplyr::filter(model %in% c("Logistic", "Random Forest")) %>% 
  dplyr::filter(complexity == "all vars") %>% 
  dplyr::select(-".config") %>% 
  knitr::kable(
    caption = "Accuracy and ROC AUC for the two best models",
    digits = 4
  )

### Compare random forest and logistic regression ------------
# Let's compare the two best models by splitting the training data into two parts

#split train data 50:50
classification_split <- initial_split(data = df_train %>% 
                                        mutate(not_retained = (rt<300) %>% 
                                                 as.factor()), 
                                      prop = 1/2, 
                                      strata = not_retained)


#Create two training folds
train1 <- testing(classification_split)
train2 <- training(classification_split)


# Metrics which we will check
multi_metric <- metric_set(accuracy, yardstick::sensitivity, yardstick::specificity, f_meas)

## Train models with first half of training data, i.e. train1
#fit logistic regression using tidymodels all predictors
LR_fit_train1 <- 
  class_wflow %>% 
  fit(data = train1)

#Random forest on first fold of train data
RF_fit_train1 <- 
  rf_wf %>% 
  fit(data = train1)

##compare models
models <- list("Logistic regression, fold 1" = LR_fit_train1,
               "Random forest, fold 1" = RF_fit_train1)

preds1 <- imap_dfr(models, augment, 
                   new_data = train2, .id = "model")

## Train models with second half of training data, i.e. train2

#Logistic regression
LR_fit_train2 <- 
  class_wflow %>% 
  fit(data = train2)

#Random forest fold2
RF_fit_train2 <- 
  rf_wf %>% 
  fit(data = train2)

##compare models
models <- list("Logistic regression, fold 2" = LR_fit_train2,
               "Random forest, fold 2" = RF_fit_train2)

preds2 <- imap_dfr(models, augment, 
                   new_data = train1, .id = "model")

## Comparison of both folds ----

fold1 <- preds1 %>% group_by(model) %>% roc_curve(truth = not_retained, estimate = .pred_FALSE)
fold2 <- preds2 %>% group_by(model) %>% roc_curve(truth = not_retained, estimate = .pred_FALSE)

bind_rows(fold1, fold2) %>% 
  separate(col = model, into = c("model", "fold"), sep = ", ") %>% 
  ggplot(aes(x = 1-specificity, y = sensitivity, color = model, lty = model)) +
  geom_line(size = 1.2, alpha = 0.9) +
  geom_abline(lty = 3) +
  facet_grid(. ~ fold) +
  theme(legend.position = "top")

## Training the final model with logistic regression -----
final_fit <- 
  class_wflow %>% 
  fit(data = df_train)

#Regression coefficients with p < 0.05
final_fit %>% 
  tidy %>% 
  filter(p.value<0.05) %>% 
  filter(abs(estimate) > 1) %>% 
  knitr::kable(caption = "Regression coefficients for terms with a p-value below 0.05 and absolute estimated value above one for the final logistic regression model")


#Confusion matrix
augment(final_fit, new_data = df_test) %>%
  conf_mat(truth = not_retained, estimate = .pred_class) %>% 
  autoplot(type = "heatmap")



#metrics for final classification model
augment(final_fit, new_data = df_test) %>%
  multi_metric(truth = not_retained, estimate = .pred_class) %>% 
  knitr::kable(caption = "Accuracy, sensitivity, specificity, and F1 score as metrics for the final classification model performance")


### Predicting retention time ------------------------------

## Let's use the mean value  of retention time as the "null model"
# This is the baseline for prediction accuracy we should be able to beat

# Mean rt for retained molecules
rt_mean <- df %>% dplyr::filter(rt>=300) %>% pull(rt) %>% mean()

# Function for calculating RMSE
calculate_rmse <- function(prediction, truth){
  RMSE <- sqrt(mean((truth - prediction)^2))
  return(RMSE)
}

#RMSE baseline from the null model
baseline_prediction_error <- calculate_rmse(prediction = rt_mean, 
                                            truth = filter(df, rt >=300) %>% pull(rt))

# Creating the data for regression modelling
### Test/Train split

## filter out non-retained molecules
df <- df %>% dplyr::filter(rt>=300)

set.seed(1234)
#split to test train 25:75
regression_split <- initial_split(data = df, 
                                  prop = 3/4#, strata = rt
)

#Create test and train
df_test <- testing(regression_split)
df_train <- training(regression_split)


### Create 10-fold cross-validation dataset for train

train_folds <- vfold_cv(df_train, v = 10)

# is retention time similarly distributed in test and train?

hist1 <- df_train %>% 
  ggplot(aes(x = rt)) + 
  geom_histogram(fill = "red", color = "black") + 
  labs(x = "Retention time (s)") + 
  ggtitle("a) train set")

hist2 <- df_test %>% 
  ggplot(aes(x = rt)) + 
  geom_histogram(fill = "red", color = "black") + 
  labs(x = "Retention time (s)") + 
  ggtitle("b) test set")

grid.arrange(hist1, 
             hist2, 
             ncol = 2
)

rm(hist1, hist2)

## Linear regression

#Simple model
lm_fit <- lm(rt ~ MW + ALogP, data = df_train)
# summary(lm_fit)
equatiomatic::extract_eq(lm_fit)

## Tidymodels fit with 10-fold cross-validation

## Linear regression model specification
lm_spec <- 
  linear_reg() %>% 
  set_mode("regression") %>% 
  set_engine("lm")


### Simple lm model -----
recipe_simple <- 
  recipe(rt ~ MW + ALogP, data = df_train) #retention time is the outcome we want to predict

# Workflow

lm_wflow <- 
  workflow() %>% 
  add_model(lm_spec) %>% 
  add_recipe(recipe_simple)


#fit linear regression using tidymodels
lm_fit1 <- 
  lm_wflow %>% #fit(df_train)
  fit_resamples(resamples = train_folds)

# Let's collect metrics from our 10-fold cross-validation
collect_metrics(lm_fit1) %>% 
  dplyr::select(-".config") %>% 
  knitr::kable(caption = "RMSE and R-squared metrics for the simple linear regression model",
               digits = 4)

## Which variables have near zero variance?
low_variance <- caret::nearZeroVar(df_train)

## Model with low variance variables removed
recipe_intermed <- 
  recipe(rt ~ ., data = df_train) %>%  #retention time is the outcome we want to predict
  update_role(pubchem, MF, new_role = "ID") %>%   #pubchem id, and molecular formula as ID variables (not predictors)
  update_role(all_of(low_variance), new_role = "ID")

# Workflow
lm_wflow <- 
  workflow() %>% 
  add_model(lm_spec) %>% 
  add_recipe(recipe_intermed)

#fit linear regression using tidymodels
lm_fit2 <- 
  lm_wflow %>% #fit(df_train)
  fit_resamples(resamples = train_folds)

# all predictors
recipe_all <- 
  recipe(rt ~ ., data = df_train) %>% #retention time is the outcome we want to predict
  update_role(pubchem, MF, new_role = "ID")

# Workflow

lm_wflow <- 
  workflow() %>% 
  add_model(lm_spec) %>% 
  add_recipe(recipe_all)

#fit linear regression using tidymodels all predictors
lm_fit3 <-
  lm_wflow %>%
  fit_resamples(resamples = train_folds)


# Let's collect metrics from our 10-fold cross-validation
collect_metrics(lm_fit2) %>% 
  mutate(model = "intermediate") %>% 
  bind_rows(collect_metrics(lm_fit3) %>% mutate(model = "all vars")) %>% 
  dplyr::select(-".config") %>% 
  knitr::kable(caption = "RMSE and R-squared metrics for the linear regression models with more predictors",
               digits = 4)

#Let's collect the most important features for the model with all variables
#fit model with the entire training data
lm_fit <- workflow() %>%
  add_model(lm_spec) %>% 
  add_recipe(recipe_all) %>% 
  fit(df_train)

#10 Most statistically significant predictors
lm_fit %>% tidy() %>% 
  filter(p.value < 0.05) %>% 
  arrange(p.value, desc(abs(estimate))) %>% 
  head(6) %>% 
  knitr::kable(caption = "Model terms with the smallest p-values")

#Ridge regression

#specification
ridge_spec <- linear_reg(mixture = 0, penalty = tune()) %>%
  set_mode("regression") %>%
  set_engine("glmnet")

#ridge recipe
ridge_recipe <- 
  recipe(formula = rt ~ ., data = df_train) %>% 
  update_role(pubchem, MF, new_role = "ID") %>% 
  step_zv(all_predictors()) %>% 
  step_normalize(all_predictors()) #ridge is scale sensitive, hence the normalization of vars

#ridge workflow
ridge_workflow <- workflow() %>% 
  add_recipe(ridge_recipe) %>% 
  add_model(ridge_spec)

#tuning parameter values
penalty_grid <- grid_regular(penalty(range = c(-5, 5)), levels = 50)

## Fit 10-fold CV to find best value for regularization
#Note: the fitting will give a warning if/when the regularization eliminates all predictors
tune_res <- tune_grid(
  ridge_workflow,
  resamples = train_folds, 
  grid = penalty_grid
)

#plot metrics
autoplot(tune_res)

#best RMSE
best_penalty <- select_best(tune_res, metric = "rmse")
#best_penalty

#Performance with the best regularization parameter value
#Not as good as linear regression
tune_res %>% 
  collect_metrics() %>% 
  filter(penalty == best_penalty$penalty) %>% 
  dplyr::select(-".config") %>% 
  knitr::kable(digits = 5, caption = "RMSE and R-squared metrics for the best performing ridge regression model")

# Lasso regression

#Lasso recipe
lasso_recipe <- ridge_recipe

#Model specs
lasso_spec <- 
  linear_reg(penalty = tune(), mixture = 1) %>% 
  set_mode("regression") %>% 
  set_engine("glmnet") 

#Lasso workflow
lasso_workflow <- workflow() %>% 
  add_recipe(lasso_recipe) %>% 
  add_model(lasso_spec)

#Fit 10-fold cv to tune penalty parameter
#Note: the fitting will give a warning if/when the regularization eliminates all predictors
tune_res <- tune_grid(
  lasso_workflow,
  resamples = train_folds, 
  grid = penalty_grid
)

autoplot(tune_res)

#best RMSE
best_penalty <- select_best(tune_res, metric = "rmse")

#Performance with the best regularization parameter value
#Not as good as linear regression
tune_res %>% 
  collect_metrics() %>% 
  filter(penalty == best_penalty$penalty) %>% 
  dplyr::select(-".config") %>% 
  knitr::kable(digits = 4, caption = "RMSE and R-squared metrics for the best performing lasso model")

#recipe for polynomial regression 
rec_poly <- recipe(rt ~ ., data = df_train) %>%
  update_role(MF, pubchem, new_role = "ID") %>% 
  step_poly(all_of(phys_chem), degree = 3) # include terms up to third power for terms related to phys_chem properties

# workflow for polynomial regression fit
poly_wf <- 
  workflow() %>%
  add_model(lm_spec) %>% #specification for linear regression
  add_recipe(rec_poly)

# Fit polynomial regression with 10-fold CV
poly_fit <- 
  poly_wf %>% 
  fit_resamples(resamples = train_folds)

poly_fit %>% 
  collect_metrics() %>% 
  dplyr::select(-".config") %>% 
  knitr::kable(digits = 4, caption = "RMSE and R-squared metrics for the cubic regression model")

## regression trees ----

reg_tree_spec <- decision_tree() %>% 
  set_engine("rpart") %>% 
  set_mode("regression")

reg_tree_wf <- workflow() %>% 
  add_recipe(recipe_all) %>% 
  add_model(reg_tree_spec)

set.seed(2134)

# Fit regression trees with 10-fold CV
reg_tree_fit <- 
  reg_tree_wf %>% 
  fit_resamples(resamples = train_folds)

reg_tree_fit %>% 
  collect_metrics() %>% 
  dplyr::select(-".config") %>% 
  knitr::kable(digits = 4, caption = "RMSE and R-squared metrics for a regression tree model")

#Let's try tuning cost complexity

reg_tree_wf <- workflow() %>%
  add_recipe(recipe_all) %>% 
  add_model(reg_tree_spec %>% set_args(cost_complexity = tune()))

param_grid <- grid_regular(cost_complexity(range = c(-4, -1)), levels = 5)

#Find best value for cost_compexity
tune_res <- tune_grid(
  reg_tree_wf, 
  resamples = train_folds, 
  grid = param_grid
)

#Plot tuning results
autoplot(tune_res)

#best model
best_tree <- select_best(tune_res, metric = "rmse")

tune_res %>% 
  collect_metrics() %>% 
  filter(cost_complexity == best_tree$cost_complexity) %>% 
  dplyr::select(-".config") %>% 
  knitr::kable(digits = 4, caption = "RMSE and R-squared metrics for the best performing regression tree model")

#Best model
reg_tree_final <- finalize_workflow(reg_tree_wf, best_tree)

#Fit tree on entire train data to visualize variable importance
reg_tree_final_fit <- fit(reg_tree_final, data = df_train)

#plot variable importance
reg_tree_final_fit %>%
  extract_fit_parsnip() %>%
  vip()

## Boosted trees with xgboost -----

#Boost specification
boost_spec <- boost_tree(trees = 1000, tree_depth = 2) %>%
  set_engine("xgboost") %>%
  set_mode("regression")

boost_wf <- workflow() %>% 
  add_recipe(recipe_all) %>% 
  add_model(boost_spec)

set.seed(2134)

# Fit boosted trees with 10-fold CV
boost_fit <- 
  boost_wf %>% 
  fit_resamples(resamples = train_folds)

boost_fit %>% 
  collect_metrics() %>% 
  dplyr::select(-".config") %>% 
  knitr::kable(digits = 4, caption = "RMSE and R-squared metrics for a boosted tree model containing 1000 trees")

#Let's try to tune the model parameters
### Let's try tuning with grid search only two parameters -------------------------
#Boost specification
boost_spec <- boost_tree(trees = tune(),
                         tree_depth = tune()
                         ) %>%
  set_engine("xgboost") %>%
  set_mode("regression")

#parameters to be tuned
boost_parameters <- parameters(boost_spec)

#Tuning range
boost_parameters %>% pull_dials_object("trees")
boost_parameters %>% pull_dials_object("tree_depth")

#Let's use a smaller tree depth
boost_parameters <- boost_parameters %>% update(tree_depth = tree_depth(c(1,4)))
#Let's update the number of trees
boost_parameters <- boost_parameters %>% update(trees = trees(c(500,1500)))

boost_parameters %>% pull_dials_object("tree_depth")
boost_parameters %>% pull_dials_object("trees")
# A regular grid with three levels for trees and 4 for tree_depth
boost_grid <- grid_regular(boost_parameters, levels = c(trees = 3, tree_depth = 4))

#workflow for assessing the grid
boost_wf <-
  workflow() %>%
  add_model(boost_spec) %>%
  add_recipe(recipe_all)

set.seed(1234)
## grid search with rmse as the main metric
# trees 500, 1000, 1500 & tree_depth 1, 2, 3, 4
boost_tune <-
  boost_wf %>%
  tune_grid(
    train_folds,
    grid = boost_grid,
    metrics = metric_set(rmse)
  )

#save tuning results
xgb_tuning_results <- collect_metrics(boost_tune)

#Let's see if tree depth of 5 improves results
#Let's update tree depth
boost_parameters <- boost_parameters %>% update(tree_depth = tree_depth(c(5,5)))
#Let's update the number of trees
boost_parameters <- boost_parameters %>% update(trees = trees(c(1000,1500)))

boost_parameters %>% pull_dials_object("tree_depth")
boost_parameters %>% pull_dials_object("trees")
# A regular grid with three levels for trees and 4 for tree_depth
boost_grid2 <- grid_regular(boost_parameters, levels = c(trees = 2, tree_depth = 1))

set.seed(1234)
## grid search with rmse as the main metric
# tree depth 5 and trees 1000 & 1500
boost_tune2 <-
  boost_wf %>%
  tune_grid(
    train_folds,
    grid = boost_grid2,
    metrics = metric_set(rmse)
  )

#Let's combine the new results with the first round of tuning
xgb_tuning_results_all <-
  xgb_tuning_results %>%
  full_join(collect_metrics(boost_tune2))

# Write results to a file
write_csv(xgb_tuning_results_all, file = "regressionmetrics.csv")

#load tuning results from a csv-file
xgb_results <- read_csv(file = "regressionmetrics.csv")

#plotting the tuning results
xgb_results %>% 
  mutate(trees = as.factor(trees)) %>% 
  ggplot(aes(x = tree_depth, y = mean, color = trees)) +
  geom_line(aes(lty = trees)) +
  geom_point(size = 2) +
  theme(legend.position = "top") +
  labs(x = "Tree depth", y = "RMSE", color = "# Trees", lty = "# Trees")

# Best models as table
#show_best(boost_tune) %>% select(-".estimator", -".config")
xgb_results %>% 
  dplyr::select(-".estimator", - ".config") %>% 
  arrange(mean) %>% 
  head(2) %>% 
  knitr::kable(digits = 4, caption = "RMSE for the two best models found during parameter tuning")

### Final regression model training and evaluation -------------------

#selecting the best
#select_best(boost_tune, metric = "rmse")

boost_parameters <- tibble(
  trees = 1000,
  tree_depth = 5
)

#updating the workflow with best parameters
final_boost_wf <- 
  boost_wf %>% 
  finalize_workflow(boost_parameters)

#Final regression model fit
final_boost_fit <- 
  final_boost_wf %>% 
  fit(df_train)

## Fit test data
df_test_preds <- augment(final_boost_fit, df_test)

## Metrics for test set
regression_metrics <- metric_set(rmse, rsq, mae)
regression_metrics(df_test_preds, truth = rt, estimate = .pred) %>% 
  knitr::kable(digits = 4, caption = "RMSE, R-squared, and MAE of the final model evaluated on the test data")

## Visualizing residuals

df_test_preds %>% 
  dplyr::mutate(residuals = rt - .pred) %>% 
  ggplot(aes(x = residuals)) +
  geom_histogram(fill = "red", color = "black", binwidth = 10) +
  geom_vline(xintercept = c(-36, 36), lty = 2) #observed mean variability in retention time measurements

#how_many are good predictions
df_test_preds %>% 
  mutate(residuals = rt - .pred) %>% 
  mutate(good_pred = abs(residuals) < 36) %>% 
  group_by(good_pred) %>% 
  summarise(n = n()) %>% 
  mutate(prop = n/sum(n)) %>% 
  knitr::kable(digits = 3, caption = "Number and proportion of predictions within and outside 36 s of the true value")

## Visualizing predictions to true values
mod_perf.pl <-
  df_test_preds %>%
  ggplot(aes(x = rt, y = .pred)) +
  geom_point(alpha = 0.2) +
  geom_abline(slope = 1) +
  labs(x = "True values (s)", y = "Predicted values (s)")

ggsave("finalmodelplot.png", plot = mod_perf.pl, dpi = 200, device = "png", width = 5, height = 3, units = "in")
