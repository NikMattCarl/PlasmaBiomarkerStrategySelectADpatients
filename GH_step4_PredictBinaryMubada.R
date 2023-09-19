#Train and test models to predict binary MUBADA status

m=list()
results_mubada=data.frame(biomarker=character(),
                      AUC_cv10_train=numeric(),
                      train_lw95=numeric(),
                      train_hi95=numeric(),
                      AUC_test=numeric(),
                      test_lw95=numeric(),
                      test_hi95=numeric())

biomarkers_roc=biomarkers_roc_test=list()

trControl <- trainControl(method = "cv",
                          number = 10,
                          classProbs = TRUE,
                          savePredictions = TRUE,
                          summaryFunction = twoClassSummary)

for (i in plasma_markers_log10) {
  train_set_df$biomarker=train_set_df[,i]
  test_set_df$biomarker=test_set_df[,i]

  m[[i]] <- train(mubada_high~., data = train_set_df[,c("mubada_high","biomarker","age","gender","apoe4")], 
                                        method = "glm",
                                        trControl=trControl,
                                        tuneLength = 10,
                                        metric="ROC")
  
  m_glm=m[[i]] 
  AUC_cv10=m_glm$results$ROC
  ROCSD=m_glm$results$ROCSD

  #get ROC cuvrves
  train_set_df$preds_caret <- as.numeric(unlist(predict(m_glm, train_set_df, type="prob")["high_mubada"]))
  biomarkers_roc[[i]]=pROC::roc(train_set_df$mubada_high~train_set_df$preds_caret, ci=TRUE)

  #get performance in test set
  temp_test=na.omit(test_set_df[,c("biomarker","mubada_high","age","gender","apoe4")])
  temp_test$preds_caret <- as.numeric(unlist(predict(m_glm, temp_test, type="prob")["high_mubada"]))
  biomarkers_roc_test[[i]]=pROC::roc(temp_test$mubada_high~temp_test$preds_caret, ci=TRUE)
  glm_test_auc=as.numeric(pROC::auc(biomarkers_roc_test[[i]]))
  glm_test_ci=as.numeric(pROC::ci(biomarkers_roc_test[[i]]))
  
  #store results in results_mubada
}

#Random forest
trControl <- trainControl(method = "cv",  number = 10,  classProbs = TRUE,     savePredictions = TRUE,
                          summaryFunction = twoClassSummary)

rf1 = train(mubada_high~., data = train_set_df[,c("mubada_high",plasma_markers_log10,"age","gender","apoe4")], 
             method = "rf",
             trControl=trControl,
             tuneLength = 10,
             metric="ROC")
final_mtry=rf1$finalModel$mtry
auc_rf1=rf1$results[rf1$results$mtry==final_mtry,"ROC"]
sd_rf1=rf1$results[rf1$results$mtry==final_mtry,"ROCSD"]

#performance in test data (as shown above)
#store results in results_mubada

#XGB
xgb1 = train(mubada_high~., data = train_set_df[,c("mubada_high",plasma_markers_log10,"age","gender_baseline_variable","apoe4")], 
              method = "xgbTree",
              trControl=trControl,
              verbosity=0)

besttune=xgb1$bestTune
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
#store results in results_mubada

#  SVR
svm_model = train(mubada_high~., data = train_set_df[,c("mubada_high",plasma_markers_log10,"age","gender","apoe4")], 
              method = "svmRadial",
              trControl=trControl,
              verbosity=0)

besttune=svm_model$bestTune
svm_auc_cv10=svm_model$results[svm_model$results$sigma==besttune$sigma &
                            svm_model$results$C==besttune$C ,"ROC"]
svm_auc_sd=svm_model$results[svm_model$results$sigma==besttune$sigma &
                                 svm_model$results$C==besttune$C ,"ROCSD"]

#performance in test data (as shown above)
#store results in results_mubada

#now compare all test ROC curves to p-tau217 alone
pvals=NULL
for (i in results_mubada$biomarker) {
  rt=roc.test(biomarkers_roc_test[["PL_ptau217"]],
              biomarkers_roc_test[[i]],
              methods="delong",
              paired=TRUE)  
  pvals=c(pvals,rt$p.value)
}
results_mubada$pval_comp_ptau217_test=pvals


#now calculate classification metrics for different levels of the plasma biomarker
temp=na.omit(train_set_predABpos[,c("PL_ptau217","mubada_high")])
temp$biomarker=temp$PL_ptau217
m=glm(mubada_high~biomarker, data=temp, family="binomial")
temp$preds=predict(m, newdata=temp, type="response")

#use cutpointr to calcaulte tn, fn, tp, fp for all possible cutpoints  
opt_cut_tp = cutpointr(temp, biomarker, mubada_high,
                          direction = ">=", pos_class = "high_mubada",
                          neg_class = "low_med_mubada", method = maximize_metric,
                          metric = tp)
curve_cutpoints=as.data.frame(opt_cut_tp$roc_curve[[1]])
plot_df=data.frame(cutpoint=curve_cutpoints$x.sorted)
plot_df["Saved Scans"]=(plot_df$tn+plot_df$fn)/(plot_df$tp+plot_df$fp+plot_df$tn+plot_df$fn)
plot_df$'True Positive'=plot_df$tp/(plot_df$tp+plot_df$fp+plot_df$tn+plot_df$fn)
plot_df$'False Positive'=plot_df$fp/(plot_df$tp+plot_df$fp+plot_df$tn+plot_df$fn)
plot_df$'True Negative'=plot_df$tn/(plot_df$tp+plot_df$fp+plot_df$tn+plot_df$fn)
plot_df$'False Negative'=plot_df$fn/(plot_df$tp+plot_df$fp+plot_df$tn+plot_df$fn)
plot_df$'False Negative Rate'=plot_df$fn/(plot_df$tp+plot_df$fn)
plot_df$sens=plot_df$tp/(plot_df$tp+plot_df$fn)
plot_df$spec=plot_df$tn/(plot_df$tn+plot_df$fp)
  
#evaluate cutpoints for specific FNR in the test set
results_train_test_saved=data.frame(group=character(),
                                    cutpoint_choice=character(),
                                    FN=numeric(),
                                    FP=numeric(),
                                    TN=numeric(),
                                    TP=numeric(),
                                    FNR=numeric(),
                                    Saved.scans=numeric(),
                                    cp=numeric(),
                                    sens=numeric(),
                                    spec=numeric())

fnr_vector=c(0.05,0.1,0.2)
for (i in fnr_vector) {
  t=plot_df[plot_df$'False.Negative.Rate'<=i,]
  tmax=t[t$Saved.Scans==max(t$Saved.Scans),]
  
  results_train_test_saved_sub=data.frame(group="train",
                                          cutpoint_choice=paste("FNR<=",i,sep=""),
                                          FN=tmax$'False.Negative',
                                          FP=tmax$'False.Positive',
                                          TN=tmax$'True.Negative',
                                          TP=tmax$'True.Positive',
                                          FNR=tmax$'False.Negative.Rate',
                                          Saved.scans=tmax$Saved.Scans,
                                          cp=tmax$cutpoint,
                                          sens=tmax$sens,
                                          spec=tmax$spec)
  
  #store results
}

#evaluate cut-points in test set
for (i in results_train_test_saved$cp) {
  temp=test_set[test_set$Predicted_AB_class=="AB_Pos", ]
  test_table=table(temp$mubada_high,temp$biomarker>= as.numeric(i))
  test_fn=test_table["high_mubada","FALSE"]
  test_fp=test_table["low_med_mubada","TRUE"]
  test_tn=test_table["low_med_mubada","FALSE"]
  test_tp=test_table["high_mubada","TRUE"]
  test_FN=test_fn/(sum(test_table))
  test_FP=test_fp/(sum(test_table))
  test_TN=test_tn/(sum(test_table))
  test_TP=test_tp/(sum(test_table))
  test_FNR=test_fn/(test_fn+test_tp)
  test_saved_scans=(test_fn+test_tn)/sum(test_table)
  sens=test_tp/( test_tp+ test_fn)
  spec=test_tn/( test_tn+ test_fp)
  #store results
}

