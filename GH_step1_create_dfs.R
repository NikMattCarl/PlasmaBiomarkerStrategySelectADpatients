rm(list = ls())

#libraries that are used somewhere in the analyses
library(dplyr)
library(tidyr)
library(caret)
library(gridExtra)
library(cutpointr)

# PREPARE DATA

#load the data (contains demographics, biomarkers, tau PET data)
db=read.csv(file="your_data.csv")

#define plasma_markers_log10 as a list of biomarkers to evaluate

#split into training and test sets
set.seed(123) #for reproducibility
train_indices=sample(seq_len(nrow(db)), size = 0.8 * nrow(db))
test_set=db[train_indices,] 
train_set=db[-train_indices,]

#test_set and train_set are saved for future use