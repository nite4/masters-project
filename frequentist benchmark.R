library(rugarch)
library(parallel)
library(rmgarch)

data <- read.csv('study_data.csv', header=TRUE, sep=',')

#for VAR benchmark model residuals - DCC-GARCH(1,1)
#to determinate prior distributions of unconditional alpha and beta

bench_resid <- read.csv('bench_resid.csv', header=TRUE, sep=',')
dcc_bench <- dccfit(
  data = bench_resid[,2:5], 
  spec = dccspec(
    uspec = multispec(replicate(4, ugarchspec(distribution.model='sstd')
    )), distribution='mvnorm'))
coeffs <- data.frame(dcc_bench@mfit$coef)
coeffs$dcc_bench.mfit.coef <- as.numeric(coeffs$dcc_bench.mfit.coef)
write.csv(coeffs, 'dcc_bench_coeffs.csv')