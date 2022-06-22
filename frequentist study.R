library(rugarch)
library(parallel)
library(rmgarch)
library(zoo)
library(vars)

setwd("~/Downloads/Praca magisterska/Kod")
set.seed(12345)

data <- read.csv('study_data.csv', header=TRUE, sep=',')

var_coeffs <- data.frame()
var_coeffs_ext <- data.frame()
varf_h <- data.frame()
dccf_H <- data.frame()

#one-year window
#GARCH(1,1) for VAR(1) residuals, distribution for multispec can be mvnorm
for (i in seq(1, (nrow(data)-250), by=20))
{
  print(i)
  
  #VAR(1) model
  var_model <- VAR(data[i:(i+250), 2:5], p=1)
  var_coeffs <- rbind(var_coeffs, matrix(data=c(var_model$varresult$SPX$coefficients,
                                                var_model$varresult$DAX$coefficients,
                                                var_model$varresult$WIG20$coefficients,
                                                var_model$varresult$EURUSD$coefficients), nrow=1))
  
  #forecast from help VAR(1) models from a slightly different window
  if(i==2041)
    k_iter <- seq(0, 18, by=1)
  else
    k_iter <- seq(0, 19, by=1)
  
  for (k in k_iter)
  {
    var_model_1 <- VAR(data[(i+k):(i+k+250), 2:5], p=1)
    varf <- predict(var_model_1, n.ahead=1, ci=0.95)
    var_coeffs_ext <- rbind(var_coeffs_ext, matrix(data=c(var_model_1$varresult$SPX$coefficients,
                                                          var_model_1$varresult$DAX$coefficients,
                                                          var_model_1$varresult$WIG20$coefficients,
                                                          var_model_1$varresult$EURUSD$coefficients), nrow=1))
    varf_h <- rbind(varf_h, cbind(varf$fcst$SPX[,1],
                                  varf$fcst$DAX[,1],
                                  varf$fcst$WIG20[,1],
                                  varf$fcst$EURUSD[,1]))
  }
  
  #DCC-GARCH for VAR(1) residuals
  e <- resid(var_model)
  
  dcc <- dccfit(
    data=e,
    spec=dccspec(
      uspec=multispec(replicate(4, ugarchspec(distribution.model='sstd',
                                              mean.model=list(armaOrder=c(1,0)))
      )), dccOrder=c(1,1), distribution='mvnorm'),
    solver='nlminb', solver.control=list(tol=1e-24))
  dccf <- dccforecast(dcc, n.ahead=1)
  dccf_H <- rbind(dccf_H, array(unlist(dccf@mforecast$H)))
}

dccf_H <- data.frame(dccf_H)
colnames(dccf_H) <- c('H11', 'H12', 'H13', 'H14',
                      'H21', 'H22', 'H23', 'H24',
                      'H31', 'H32', 'H33', 'H34',
                      'H41', 'H42', 'H43', 'H44')

var_coeffs <- data.frame(var_coeffs)
colnames(var_coeffs) <- c('SPX_SPXl1', 'SPX_DAXl1', 'SPX_WIG20l1', 'SPX_EURUSD_l1', 'SPX_const',
                          'DAX_SPXl1', 'DAX_DAXl1', 'DAX_WIG20l1', 'DAX_EURUSD_l1', 'DAX_const',
                          'WIG20_SPXl1', 'WIG20_DAXl1', 'WIG20_WIG20l1', 'WIG20_EURUSD_l1', 'WIG20_const',
                          'EURUSD_SPXl1', 'EURUSD_DAXl1', 'EURUSD_WIG20l1', 'EURUSD_EURUSD_l1', 'EURUSD_const')

var_coeffs_ext <- data.frame(var_coeffs_ext)
colnames(var_coeffs_ext) <- c('SPX_SPXl1', 'SPX_DAXl1', 'SPX_WIG20l1', 'SPX_EURUSD_l1', 'SPX_const',
                          'DAX_SPXl1', 'DAX_DAXl1', 'DAX_WIG20l1', 'DAX_EURUSD_l1', 'DAX_const',
                          'WIG20_SPXl1', 'WIG20_DAXl1', 'WIG20_WIG20l1', 'WIG20_EURUSD_l1', 'WIG20_const',
                          'EURUSD_SPXl1', 'EURUSD_DAXl1', 'EURUSD_WIG20l1', 'EURUSD_EURUSD_l1', 'EURUSD_const')


varf_h <- data.frame(varf_h)
varf_h <- cbind(data$X[251:nrow(data)], varf_h)
colnames(varf_h) <- c('Date', 'SPX', 'DAX', 'WIG20', 'EURUSD')

write.csv(dccf_H, 'dccf_H_forecast_new.csv', row.names=FALSE)
write.csv(varf_h, 'varf_h_forecast_new.csv', row.names=FALSE)
write.csv(var_coeffs, 'var_coeffs_new.csv', row.names=FALSE)
write.csv(var_coeffs_ext, 'var_coeffs_ext_new.csv', row.names=FALSE)
