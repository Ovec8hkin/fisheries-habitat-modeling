library(tidyverse)  
library(mgcv)
library(DHARMa)
library(MASS)

simulateNegbin <- function(modfit, nsims=250, offsetval=1){
  muval = predict(modfit, type = "response")*offsetval  #Get the mean with offset
  nObs = length(muval)  
  thetaval = modfit$family$getTheta(trans=TRUE)  #Get theta not log transforme
  sim = replicate(nsims,rnbinom(nObs,mu=muval, size=thetaval))  #Simulate negative binomial data
  return(sim)
}

get_dharma_residuals <- function(gam, data){
  pred_vals = mgcv::predict.gam(gam, type="response")
  
  # Subset simulations for better residual vs predicted dharma plot
  simvals=simulateNegbin(gam,250,offsetval=1)
  
  DHARMaRes = createDHARMa(simulatedResponse = simvals, 
                           observedResponse = data$catch, 
                           fittedPredictedResponse = pred_vals
  )
  return(DHARMaRes)
}

compute_aic_scores <- function(...){
  aic_scores = AIC(...)
  aic_scores = aic_scores[order(aic_scores$AIC),]
  delta_aic = aic_scores$AIC - min(aic_scores$AIC)
  aic_scores = cbind(aic_scores, delta_aic)
  aic_scores = round(aic_scores, 3)
  return(aic_scores)
}

check_model <- function(model, data){
  residuals = get_dharma_residuals(model, data)
  resids = plot(residuals)
  dispersion = testDispersion(residuals)
  uniformity = testUniformity(residuals)
  zeroinflation = testZeroInflation(residuals)
  outliers = testOutliers(residuals)
  
  model_check_stats = c(dispersion$statistic, uniformity$statistic, zeroinflation$statistic)
  model_check_pval = c(dispersion$p.value, uniformity$p.value, zeroinflation$p.value)
  
  model_check = data.frame(model_check_stats, model_check_pval)
  colnames(model_check) = c("Statistic", "p-value")
  row.names(model_check) = c("Dispersion", "Uniformity", "Zero-Inflation")
  model_check = round(model_check, 3)
  return(model_check)
}

repredict_data <- function(model, data){
  re_predicted = predict(model, type="response")
  re_predicted_data = cbind(data, re_predicted)
  diff = re_predicted_data$catch - re_predicted_data$re_predicted
  re_predicted_data = cbind(re_predicted_data, diff)
  
  return(re_predicted_data)
}

predict_data <- function(model, data){
  predicted = predict(model, data, type="response")
  pred_new = cbind(data, predicted)
  pred_new$floored = floor(predicted)
  return(pred_new)
}

dredge_models <- function(global.model){
  d = MuMIn::dredge(global.model, extra = list(
    "*" = function(x) {
      s <- summary(x)
      c(Rsq = s$r.sq, dev = s$dev.exp)
    })
  )
  return(d)
}

format_model_table <- function(models){
  models.table = models[, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 12)]
  models.table[, c(1, 7, 8, 10, 11)] = round(models.table[, c(1, 7, 8, 10, 11)], 3)
  names(models.table) = c("Int", "b_salt", "b_temp", "depth", "s_salt", "s_temp", "Rsq", "dev", "df", "AICc", "∆")
  return(models.table)
}






