#Train and test models to predict binary AB-status

temp_train=na.omit(train_set[,c("sid","age","gender","apoe4", "AB_status",plasma_markers_log10,"mubada_group","mubada","mubada_high")])
temp_test=na.omit(test_set[,c("sid","age","gender","apoe4","AB_status",plasma_markers_log10,"mubada_group","mubada","mubada_high")])

biomarkers_roc=biomarkers_roc_test=m=list()

#Data frame to store results for models to predict AB status
results_AB=data.frame(biomarker=character(),
                       AUC_cv10_train=numeric(),
                       train_lw95=numeric(),
                       train_hi95=numeric(),
                       AUC_test=numeric(),
                       test_lw95=numeric(),
                       test_hi95=numeric())

#Set train control
trControl <- trainControl(method = "cv", number = 10, classProbs = TRUE,savePredictions = TRUE,
                          summaryFunction = twoClassSummary)

#Loop through individual biomarkers
for (i in plasma_markers_log10) {
  temp_train$biomarker=temp_train[,i]
  temp_test$biomarker=temp_test[,i]

  #logistic regression models
  m[[i]] <- train(AB_status~., data = temp_train[,c("AB_status","biomarker","age","gender","apoe4")], 
                                        method = "glm", trControl=trControl,tuneLength = 10,
                                        metric="ROC", preProcess=c("center","scale"))

  #Performance in training set
  AUC_cv10= m[[i]]$results$ROC
  ROCSD=m[[i]]$results$ROCSD
  
  #ROC curve in training set
  temp_train$preds_caret = as.numeric(unlist(predict(m[[i]], temp_train, type="prob")["AB_Pos"]))
  biomarkers_roc[[i]]=pROC::roc(temp_train$AB_status~temp_train$preds_caret, ci=TRUE)

  #Performance in test set
  temp_test$preds_caret <- as.numeric(unlist(predict(m[[i]], temp_test, type="prob")["AB_Pos"]))
  biomarkers_roc_test[[plasma_markers_log10[i]]]=pROC::roc(temp_test$AB_status~temp_test$preds_caret, ci=TRUE)
  glm_test_auc=as.numeric(pROC::auc(biomarkers_roc_test[[plasma_markers_log10[i]]]))
  glm_test_ci=as.numeric(pROC::ci(biomarkers_roc_test[[plasma_markers_log10[i]]]))
  
  #store all results in results_AB
}

#Train random forest model
rf1 <- train(AB_status~., data = temp_train[,c("AB_status",plasma_markers_log10,"age","gender","apoe4")], 
             method = "rf",   trControl=trControl, tuneLength = 10, metric="ROC")
final_mtry=rf1$finalModel$mtry #tuning of mytry provided by caret

#tuning of n_tree
model_list = list()
tunegrid = expand.grid(.mtry = final_mtry)
control = trainControl(method="repeatedcv", number=3, repeats=10, search="grid")

for (n_tree in c(100, 500, 1000, 1500)) {
  set.seed(3333)
  fit = train(AB_status~., data=temp_train[,c("AB_status",plasma_markers_log10,"age","gender","apoe4")], 
               method="rf", metric="ROC", 
               tuneGrid=tunegrid, trControl= control, ntree=n_tree,
               preProcess = c("center","scale"))
  key = toString(n_tree)
  model_list[[key]] = fit
}
results = resamples(model_list)

#refit with tuned hyperparameters
rf1 = train(AB_status~., data=temp_train[,c("AB_status",plasma_markers_log10,"age","gender","apoe4")], 
            method="rf", metric="ROC", 
            tuneGrid=tunegrid, trControl= trControl, ntree=500) # #max ROC at ntree=500
           
#performance in training data
auc_rf1=rf1$results[rf1$results$mtry==final_mtry,"ROC"]
sd_rf1=rf1$results[rf1$results$mtry==final_mtry,"ROCSD"]

#check performance in test data (as shown above for individual biomarkers)
#store results in results_AB

#extreme gradient boosting
xgb1 <- train(AB_status~., data =  temp_train[,c("AB_status",plasma_markers_log10,"age","gender","apoe4")],  
              method = "xgbTree",    trControl=trControl,verbosity=0, metric="ROC",
              preProcess = c("center","scale"))
besttune=xgb1$bestTune #by tuning

#performance in training data
xgb_auc_cv10=xgb1$results[xgb1$results$eta==besttune$eta &
                            xgb1$results$max_depth==besttune$max_depth &
                            xgb1$results$nrounds==besttune$nrounds &
                            xgb1$results$gamma==besttune$gamma &
                            xgb1$results$colsample_bytree==besttune$colsample_bytree &
                            xgb1$results$min_child_weight==besttune$min_child_weight &
                            xgb1$results$subsample==besttune$subsample,"ROC"]
xgb_auc_sd=xgb1$results[xgb1$results$eta==besttune$eta &
                          xgb1$results$max_depth==besttune$max_depth &
                          xgb1$results$nrounds==besttune$nrounds &
                          xgb1$results$gamma==besttune$gamma &
                          xgb1$results$colsample_bytree==besttune$colsample_bytree &
                          xgb1$results$min_child_weight==besttune$min_child_weight &
                          xgb1$results$subsample==besttune$subsample,"ROCSD"]

#performance in test data (as shown above)
#store results in results_AB

# SVM
svm_model <- train(AB_status~., data = temp_train[,c("AB_status",plasma_markers_log10,"age","gender","apoe4")],  
                   method = "svmRadial",trControl=trControl,verboseIter=0, metric="ROC",
                   preProcess = c("center","scale"))
besttune=svm_model$bestTune
#performance in training data
svm_auc_cv10=svm_model$results[svm_model$results$sigma==besttune$sigma &
                                 svm_model$results$C==besttune$C ,"ROC"]
svm_auc_sd=svm_model$results[svm_model$results$sigma==besttune$sigma &
                               svm_model$results$C==besttune$C ,"ROCSD"]

#performance in test data (as shown above)
#store results in results_AB

#Compare each model to the p-tau217 model
pvals=NULL
for (i in results_AB$biomarker) {
  roc_compare=roc.test(biomarkers_roc_test[["PL_ptau217"]],
           biomarkers_roc_test[[i]],
           methods="delong",
           paired=TRUE)  
  pvals=c(pvals,roc_compare$p.value) #pvalues for comparison of ROC curves
}
#store the p-values together with results_AB

#Compare models with a sparse P-tau217 model, without covariates
m_basic <- train(AB_status~PL_ptau217,
                 data = temp_train[,c("AB_status","PL_ptau217","age","gender","apoe4")], #ensure that the same dataset is used
                 method = "glm", trControl=trControl,tuneLength = 10, metric="ROC")

#Performance in training set
temp_train$preds_caret <- as.numeric(unlist(predict(m_basic, temp_train, type="prob")["AB_Pos"]))
biomarkers_roc_basic=pROC::roc(temp_train$AB_status~temp_train$preds_caret, ci=TRUE)
#check performance in test data (as shown above)

#Comparison of test AUC to p-tau217 alone (without covariatese)
pvals=NULL
for (i in results_AB$biomarker) {
  roc_compare=roc.test(biomarkers_roc_test_basic,
              biomarkers_roc_test[[i]],
              methods="delong",
              paired=TRUE)  
  pvals=c(pvals,roc_compare$p.value) #pvalues for comparison of ROC curves
}
#store the p-values together with results_AB


#
#Performance of P-tau217 to rule out AB-negatives
#

#refit model on whole train_set
m=glm(AB_status~PL_ptau217, data=train_set, family="binomial")
train_set$preds=predict(m, newdata=train_set, type="response")

#start with a "one cutpoint" approach - Youden

opt_cut_youden <- cutpointr(train_set, PL_ptau217, AB_status,
                        direction = ">=", pos_class = "AB_Pos",
                        neg_class = "AB_Neg", method = maximize_metric,
                        metric = youden,
                        boot_runs=1000)
youden_cp=opt_cut_youden$optimal_cutpoint
boot_ci(opt_cut_youden, sensitivity) #get 95% sens
boot_ci(opt_cut_youden, specificity) #get 95% spec

#performance of youden cutpoint in test set
test_set$Predicted_AB_class=ifelse(test_set$PL_ptau217>youden_cp,"AB_Pos",
                                    ifelse(test_set$PL_ptau217<=youden_cp,"AB_Neg",NA))    
sens_test=table(test_set$Predicted_AB_class, test_set$AB_status)["AB_Pos","AB_Pos"]/sum(table(test_set$Predicted_AB_class, test_set$AB_status)[,"AB_Pos"])
spec_test=table(test_set$Predicted_AB_class, test_set$AB_status)["AB_Neg","AB_Neg"]/sum(table(test_set$Predicted_AB_class, test_set$AB_status)[,"AB_Neg"])

#bootstrap sens and spec in test set
boot_results_test=data.frame(boot=numeric(),
                             sens=numeric(),
                             spec=numeric())
for (i in 1:1000) {
  boot_id=data.frame(sid=sample(test_set$sid, replace=TRUE), boot_sid=seq(1,nrow(test_set),1))
  boot_id=left_join(boot_id, test_set, by="sid")
  boot_id$Predicted_AB_class=ifelse(boot_id$PL_ptau217>youden_cp,"AB_Pos",
                                    ifelse(boot_id$PL_ptau217<=youden_cp,"AB_Neg",NA))
  sens_boot=table(boot_id$Predicted_AB_class,boot_id$AB_status)["AB_Pos","AB_Pos"]/(sum(table(boot_id$Predicted_AB_class,boot_id$AB_status)[,"AB_Pos"]))
  spec_boot=table(boot_id$Predicted_AB_class,boot_id$AB_status)["AB_Neg","AB_Neg"]/(sum(table(boot_id$Predicted_AB_class,boot_id$AB_status)[,"AB_Neg"]))
 
  sub=data.frame(boot=i,
                 sens=sens_boot,
                 spec=spec_boot)
  boot_results_test=rbind(boot_results_test,sub)
}
boot_sens_ci_test=quantile(boot_results_test$sens, probs=c(0.025,0.975))
boot_spec_ci_test=quantile(boot_results_test$spec, probs=c(0.025,0.975))

#
#gray zone approach, 2 cutpoints
#

#evaluate gray zone for plasma ptau217 vs amyloid PET
temp_train=na.omit(train_set[,c("sid","AB_status","amyloid_PET","PL_ptau217")])
temp_test=na.omit(test_set[,c("sid","AB_status","amyloid_PET","PL_ptau217")])

m=glm(AB_status~PL_ptau217, family="binomial", data=temp_train)
temp_train$preds=predict(m, temp_train, type="response")
temp_test$preds=predict(m, temp_test, type="response")

#procedure to select cut-points for 2-cutpoint approach, based on finding cutpoints which minimizes gray zone while keeping sens and spec>90%
results=data.frame(cutoff=numeric(),
                   pos_prob=numeric(),
                   neg_prob=numeric())
temp_train=temp_train[order(temp_train$PL_ptau217, decreasing=FALSE),]
for (i in 1:nrow(temp_train)) {
  cutp=temp_train[i,"PL_ptau217"]
  meanposprob=mean(temp_train[temp_train$PL_ptau217>cutp,"preds"], na.rm=TRUE)
  meannegprob=mean(temp_train[temp_train$PL_ptau217<cutp,"preds"], na.rm=TRUE)
  
  sub=data.frame(cutoff=cutp,
                 pos_prob=meanposprob,
                 neg_prob=meannegprob)
  results=rbind(results,sub)
}

results_cupoint=data.frame(pos_neg_prob=numeric(),
                           high_cutpoint=numeric(),
                           low_cutpoint=numeric(),
                           sens_train=numeric(),
                           sens_train_lw95=numeric(),
                           sens_train_hi95=numeric(),
                           
                           spec_train=numeric(),
                           spec_train_lw95=numeric(),
                           spec_train_hi95=numeric(),
                           
                           range_gz_train=numeric(),
                           range_train_lw95=numeric(),
                           range_train_hi95=numeric(),
                           
                           sens_test=numeric(),
                           sens_test_lw95=numeric(),
                           sens_test_hi95=numeric(),
                           
                           spec_test=numeric(),
                           spec_test_lw95=numeric(),
                           spec_test_hi95=numeric(),
                           
                           range_gz_test=numeric(),
                           range_test_lw95=numeric(),
                           range_test_hi95=numeric())
probs_try=seq(from=0.85,to=0.95, by=0.01) #iterate through possible levels of probability
for (k in 1:length(probs_try)) {
  
  current_prob=probs_try[k]
  
  #lowest cutpoint for at least current_prob pos_prob
  pos_group=results[results$pos_prob>current_prob,]
  minpos_group=pos_group[order(pos_group$cutoff, decreasing=FALSE),]
  hi_cutpoint=minpos_group[1,"cutoff"]
  
  #highest cutpoint for at least 1-current_prob pos_prob
  neg_group=results[results$neg_prob<(1-current_prob),]
  maxneg_group=neg_group[order(neg_group$cutoff, decreasing=TRUE),]
  low_cutpoint=maxneg_group[1,"cutoff"]
  
  #use optimal cutpoints to predict AB class, for carrying forward a filtered data set
  train_set$Predicted_AB_class=ifelse(train_set$PL_ptau217>=hi_cutpoint,"AB_Pos",
                                      ifelse(train_set$PL_ptau217<=low_cutpoint,"AB_Neg",
                                             as.character(train_set$AB_status)))    
  sens_train=table(train_set$Predicted_AB_class, train_set$AB_status)["AB_Pos","AB_Pos"]/sum(table(train_set$Predicted_AB_class, train_set$AB_status)[,"AB_Pos"])
  spec_train=table(train_set$Predicted_AB_class, train_set$AB_status)["AB_Neg","AB_Neg"]/sum(table(train_set$Predicted_AB_class, train_set$AB_status)[,"AB_Neg"])
  range_gz_train=prop.table(table(train_set$PL_ptau217>low_cutpoint & train_set$PL_ptau217<hi_cutpoint))["TRUE"]
  
  #bootstrap sens and spec as shown above 
  
  #check performance in the test set
  test_set$Predicted_AB_class=ifelse(test_set$PL_ptau217>=hi_cutpoint,"AB_Pos",
                                     ifelse(test_set$PL_ptau217<=low_cutpoint,"AB_Neg",
                                            as.character(test_set$AB_status)))    
  test_table=table(test_set$Predicted_AB_class,test_set$AB_status)
  sens_test=test_table["AB_Pos","AB_Pos"]/sum(test_table[,"AB_Pos"])
  spec_test=test_table["AB_Neg","AB_Neg"]/sum(test_table[,"AB_Neg"])
  range_gz_test=prop.table(table(test_set$PL_ptau217>low_cutpoint & test_set$PL_ptau217<hi_cutpoint))["TRUE"]
  
  #bootstrap sens and spec in test set as shown above 

  sub_prob=data.frame(pos_neg_prob=current_prob,
                      high_cutpoint=hi_cutpoint,
                      low_cutpoint=low_cutpoint,
                      
                      sens_train=sens_train,
                      sens_train_lw95=boot_sens_ci[1],
                      sens_train_hi95=boot_sens_ci[2],
                      
                      spec_train=spec_train,
                      spec_train_lw95=boot_spec_ci[1],
                      spec_train_hi95=boot_spec_ci[2],
                      
                      range_gz_train=range_gz_train,
                      range_train_lw95=boot_gz[1],
                      range_train_hi95=boot_gz[2],
                      
                      sens_test=sens_test,
                      sens_test_lw95=boot_sens_ci_test[1],
                      sens_test_hi95=boot_sens_ci_test[2],
                      
                      spec_test=spec_test,
                      spec_test_lw95=boot_spec_ci_test[1],
                      spec_test_hi95=boot_spec_ci_test[2],
                      
                      range_gz_test=range_gz_test,
                      range_test_lw95=boot_gz_test[1],
                      range_test_hi95=boot_gz_test[2])
  results_cupoint=rbind(results_cupoint,sub_prob)
}

#read results_cupoint to identify probability cutpoint with at least 90% sens, 90% spec and minimal gray zone
#in our case probability=0.88
#use this to identify cutpoints on plasma p-tau217 scale
results=data.frame(cutoff=numeric(),
                   pos_prob=numeric(),
                   neg_prob=numeric())
temp_train=temp_train[order(temp_train$PL_ptau217, decreasing=FALSE),]
for (i in 1:nrow(temp_train)) {
  cutp=temp_train[i,"PL_ptau217"]
  meanposprob=mean(temp_train[temp_train$PL_ptau217>cutp,"preds"], na.rm=TRUE)
  meannegprob=mean(temp_train[temp_train$PL_ptau217<cutp,"preds"], na.rm=TRUE)
  sub=data.frame(cutoff=cutp,
                 pos_prob=meanposprob,
                 neg_prob=meannegprob)
  results=rbind(results,sub)
}

#final hi cutpoint 
pos_range=results[results$pos_prob>0.88,]
minpos=pos_range[order(pos_range$cutoff, decreasing=FALSE),]
hi_cutpoint=minpos[1,"cutoff"]

#final low cutpoint
neg_range=results[results$neg_prob<0.12,]
maxneg=neg_range[order(neg_range$cutoff, decreasing=TRUE),]
low_cutpoint=maxneg[1,"cutoff"]

#use optimal final cutpoints for gray zone to predict AB class
train_set$Predicted_AB_class=ifelse(train_set$PL_ptau217>=hi_cutpoint,"AB_Pos",
                                    ifelse(train_set$PL_ptau217<=low_cutpoint,"AB_Neg",
                                           as.character(train_set$AB_status)))    
sens=table(train_set$Predicted_AB_class, train_set$AB_status)["AB_Pos","AB_Pos"]/sum(table(train_set$Predicted_AB_class, train_set$AB_status)[,"AB_Pos"])
spec=table(train_set$Predicted_AB_class, train_set$AB_status)["AB_Neg","AB_Neg"]/sum(table(train_set$Predicted_AB_class, train_set$AB_status)[,"AB_Neg"])
range_gz_train=prop.table(table(train_set$PL_ptau217>low_cutpoint & train_set$PL_ptau217<hi_cutpoint))["TRUE"]

#bootstrap sens and spec in training set as shown above

#check performance in the test set
test_set$Predicted_AB_class=ifelse(test_set$PL_ptau217>=hi_cutpoint,"AB_Pos",
                                   ifelse(test_set$PL_ptau217<=low_cutpoint,"AB_Neg",
                                          as.character(test_set$AB_status)))    
test_table=table(test_set$Predicted_AB_class,test_set$AB_status)
sens_test=test_table["AB_Pos","AB_Pos"]/sum(test_table[,"AB_Pos"])
spec_test=test_table["AB_Neg","AB_Neg"]/sum(test_table[,"AB_Neg"])
range_gz_test=prop.table(table(test_set$PL_ptau217>low_cutpoint & test_set$PL_ptau217<hi_cutpoint))["TRUE"]

#bootstrap sens and spec in test set as shown above



