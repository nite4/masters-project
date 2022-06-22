library(bayesDccGarch)
library(BVAR)

setwd("~/Downloads/Praca magisterska/Kod")
set.seed(1238)

data <- read.csv('study_data.csv', header=TRUE, sep=',')

var_coeffs <- data.frame()
varf_h <- data.frame()
dccf_H <- data.frame()

#initial values for VAR prior distributions
eq_spx_coef <- c(0.041785, -0.205799, 0.141369, -0.031199, -0.072048)
eq_dax_coef <- c(-0.005848, 0.269408, -0.082196, -0.016119, -0.114188)
eq_wig20_coef <- c(-0.045685, 0.240677, -0.089691, -0.049037, 0.023186)
eq_eurusd_coef <- c(0.009460, -0.178725, 0.068866, -0.010172, 0.075879)

eq_spx_sd <- c(0.060872, 0.075459, 0.075366, 0.064829, 0.097287)
eq_dax_sd <- c(0.069953, 0.086716, 0.086609, 0.074500, 0.111800)
eq_wig20_sd <- c(0.065959, 0.081765, 0.081664, 0.070247, 0.105417)
eq_eurusd_sd <- c(0.034561, 0.042843, 0.042790, 0.036808, 0.055236)

#initial values for DCC-GARCH prior distributions
alpha_ini <- c(0.882605, 0.023122, 0.046696, 0.117395)
sigma_alpha_ini <- c(0.034267, 0.015472, 0.036756, 0.039619)
beta_ini <- c(0.027462, 0.114750, 0.067125, 0.869239)
sigma_beta_ini <- c(0.018817, 0.053159, 0.027637, 0.060134)
omega_ini <- c(0.00001, 0.919674, 0.978022, 0.012545)
sigma_omega_ini <- c(0.036146, 0.030620, 0.061809, 0.015016)
a_ini <- 0.012303
sigma_a_ini <- 0.1
b_ini <- 0.970017
sigma_b_ini <- 0.1

#control list of DCC-GARCH parameters
dcc_control <- list(mu_omega=omega_ini, mu_alpha=alpha_ini, mu_beta=beta_ini,
                    sigma_omega=sigma_omega_ini, sigma_alpha=sigma_alpha_ini, sigma_beta=sigma_beta_ini,
                    mu_a=a_ini, mu_b=b_ini, sigma_a=sigma_a_ini, sigma_b=sigma_b_ini,
                    print=FALSE)

#one-year window
for (i in seq(1, (nrow(data)-250), by=20))
{
  print(i)
  #VAR(1) model
  var_model <- bvar(data[i:(i+250), 2:5], lags=1, n_thin=5L)
  var_coeffs <- rbind(var_coeffs, matrix(data=coef(var_model), nrow=1))
  
  if(i==2041)
    k_iter <- seq(0, 18, by=1)
  else
    k_iter <- seq(0, 19, by=1)
  
  #VAR(1) forecast
  for (k in k_iter)
  {
    var_model_1 <- bvar(data[(i+k):(i+k+250), 2:5], lags=1,
                        n_thin=5L, n_draw=20000L, n_burn=1000L, verbose=FALSE)
    varf <- predict(var_model_1, horizon=1)
    varf_h <- rbind(varf_h, matrix(data=head(varf$fcast, 1), nrow=1))
  }
  
  e <- resid(var_model)
  
  #DCC-GARCH for BVAR residuals
  #Error distribution - Standardize Skewed Normal
  dcc <- bayesDccGarch(mY=e,
                       alpha_ini=alpha_ini, beta_ini=beta_ini, omega_ini=omega_ini,
                       a_ini=a_ini, b_ini=b_ini,
                       errorDist=1, control=dcc_control)  #tail_ini=nu
  dccf = predict(dcc, horizon=1)$H
  dccf_H <- rbind(dccf_H, dccf)
  
  #extracting next prior expected values
  new_ini <- head(matrix(tail(dcc$MC, 1), nrow=2), 1)
  alpha_ini <- c(new_ini[,4], new_ini[,8], new_ini[,12], new_ini[,16])
  beta_ini <- c(new_ini[,5], new_ini[,9], new_ini[,13], new_ini[,17])
  omega_ini <- c(new_ini[,3], new_ini[,7], new_ini[,11], new_ini[,15])
  a_ini <- new_ini[,18]
  b_ini <- new_ini[,19]
  #what about standard deviations?
  
  #updating the control list of arguments
  dcc_control <- list(mu_omega=omega_ini, mu_alpha=alpha_ini, mu_beta=beta_ini,
                     sigma_omega=sigma_omega_ini, sigma_alpha=sigma_alpha_ini, sigma_beta=sigma_beta_ini,
                      mu_a=a_ini, mu_b=b_ini, sigma_a=sigma_a_ini, sigma_b=sigma_b_ini,
                      print=FALSE)
}

dccf_H <- data.frame(dccf_H)
colnames(dccf_H) <- c('H11', 'H12', 'H13', 'H14',
                      'H21', 'H22', 'H23', 'H24',
                      'H31', 'H32', 'H33', 'H34',
                      'H41', 'H42', 'H43', 'H44')
rownames(dccf_H) <- NULL

var_coeffs <- data.frame(var_coeffs)
colnames(var_coeffs) <- c('SPX_const', 'SPX_SPXl1', 'SPX_DAXl1', 'SPX_WIG20l1', 'SPX_EURUSD_l1',
                          'DAX_const', 'DAX_SPXl1', 'DAX_DAXl1', 'DAX_WIG20l1', 'DAX_EURUSD_l1',
                          'WIG20_const', 'WIG20_SPXl1', 'WIG20_DAXl1', 'WIG20_WIG20l1', 'WIG20_EURUSD_l1',
                          'EURUSD_const', 'EURUSD_SPXl1', 'EURUSD_DAXl1', 'EURUSD_WIG20l1', 'EURUSD_EURUSD_l1')

varf_h <- data.frame(varf_h)
varf_h <- cbind(data$X[251:nrow(data)], varf_h)
colnames(varf_h) <- c('Date', 'SPX', 'DAX', 'WIG20', 'EURUSD')

write.csv(dccf_H, 'dccf_H_forecast_bayesian_new.csv', row.names=FALSE)
write.csv(varf_h, 'varf_h_forecast_bayesian_new.csv', row.names=FALSE)
write.csv(var_coeffs, 'var_coeffs_bayesian_new.csv', row.names=FALSE)
