#Train and test models to predict continuous MUBADA

# fuction to calculate r2 
lm_r2=function(y_predicted, y_observed) {
  residuals = y_observed - y_predicted
  RMSE = sqrt(mean(residuals^2))
  y_observed_mean = mean(y_observed)
  # Calculate total sum of squares
  tss =  sum((y_observed - y_observed_mean)^2 )
  # Calculate residual sum of squares
  rss =  sum(residuals^2)
  # Calculate R-squared
  rsq  =  1 - (rss/tss)
  return(rsq)
}

m=list()
results=data.frame(biomarker=character(),
                   R2_cv10_train=numeric(),
                   R2_test=numeric())

trControl <- trainControl(method = "cv",number = 10,search = "grid")

for (i in plasma_markers_log10) {

  train_set_predABpos$biomarker=train_set_predABpos[,i]
  test_set_predABpos$biomarker=test_set_predABpos[,i]
  
  m[[i]] <- train(mubada~., data = train_set_predABpos[,c("mubada","biomarker","age","gender","apoe4")], 
                                        method = "lm",
                                        trControl=trControl,
                                        metric="Rsquared", preProcess=c("center","scale"))
  m_lm=m[[i]]
  r2_cv10=m_lm$results$Rsquared
  
  train_set_predABpos$predicted = predict(m_lm, train_set_predABpos)  
  test_set_predABpos$predicted = predict(m_lm, test_set_predABpos) #for plot

  lm_train_r2=lm_r2(y_predicted=train_set_predABpos$predicted,y_observed=train_set_predABpos$mubada)
  lm_test_r2=lm_r2(y_predicted=test_set_predABpos$predicted,y_observed=test_set_predABpos$mubada)
  
  #store results
}

#Random forest
rf1 = train(mubada~., data = train_set_predABpos[,c("mubada",plasma_markers_log10,"age","gender","apoe4")], 
             method = "rf", trControl=trControl,   metric="Rsquared", preProcess=c("center","scale"))
final_mtry=rf1$finalModel$mtry

#manual tuning for n_tree
model_list = list()
tunegrid = expand.grid(.mtry = final_mtry)
control = trainControl(method="repeatedcv", number=3, repeats=10, search="grid")

for (n_tree in c(100, 500, 1000, 1500)) {
  set.seed(3333)
  fit = train(mubada~., data=train_set_predABpos[,c("mubada",plasma_markers_log10,"age","gender","apoe4")], 
                 method="rf", metric="Rsquared", 
               tuneGrid=tunegrid, trControl= trControl, ntree=n_tree,
               preProcess = c("center","scale"))
  
  key = toString(n_tree)
  model_list[[key]] = fit
}
resultstemp <- resamples(model_list)

rf1 = train(mubada~., data=train_set_predABpos[,c("mubada",plasma_markers_log10,"age","gender","apoe4")], 
             method="rf", metric="Rsquared", 
             tuneGrid=tunegrid, trControl= trControl, ntree=1500, #by manual tuning
             preProcess = c("center","scale"))

r2_cv10=rf1$results[rf1$results$mtry==rf1$finalModel$mtry,"Rsquared"]#
test_set_predABpos$predicted = predict(rf1, test_set_predABpos)
rf1_test_r2=lm_r2(y_predicted=test_set_predABpos$predicted, y_observed=test_set_predABpos$mubada)
test_set_predABpos$predicted = predict(rf1, test_set_predABpos) #for plot
#store results

#XGB
xgb1 = train(mubada~., data = train_set_predABpos[,c("mubada",plasma_markers_log10,"age","gender","apoe4")], 
              method = "xgbTree",
              trControl=trControl,
              verbosity=0,
              metric="Rsquared",
              preProcess = c("center","scale"))
              
besttune=xgb1$bestTune
xgb_r2_cv10=xgb1$results[xgb1$results$eta==besttune$eta &
                           xgb1$results$max_depth==besttune$max_depth &
                           xgb1$results$nrounds==besttune$nrounds &
                           xgb1$results$gamma==besttune$gamma &
                           xgb1$results$colsample_bytree==besttune$colsample_bytree &
                           xgb1$results$min_child_weight==besttune$min_child_weight &
                           xgb1$results$subsample==besttune$subsample,"Rsquared"]

test_set_predABpos$predicted = predict(xgb1, test_set_predABpos)
xgb_test_r2=lm_r2(y_predicted=test_set_predABpos$predicted, y_observed=test_set_predABpos$mubada)
test_set_predABpos$predicted = predict(rf1, test_set_predABpos) # for plot
#store results

