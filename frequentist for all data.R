library(rugarch)
library(parallel)
library(rmgarch)

data <- read.csv('study_data.csv', header=TRUE, sep=',')

head(data,5)

#VAR for all data
var_all <- VAR(data[,2:5], p=1)
e_all <- cbind(var_all$varresult$SPX$residuals,
               var_all$varresult$DAX$residuals,
               var_all$varresult$WIG20$residuals, 
               var_all$varresult$EURUSD$residuals)

#CCC-GARCH for all data
ccc <- cgarchfit(data=e_all,
                 spec=cgarchspec(uspec=multispec(replicate(4, ugarchspec(distribution.model='sstd')))))
print(ccc@mfit[["Rt"]])

write.csv(x=ccc@mfit[["Rt"]], file='ccc_garch.csv', row.names=FALSE)

#DCC-GARCH for all data
dcc <- dccfit(
  data = e_all, 
  spec = dccspec(
    uspec = multispec(replicate(4, ugarchspec(distribution.model='sstd')
    )), distribution = 'mvnorm'))

para12 <- c(); para13 <- c(); para14 <- c(); para23 <- c(); para24 <- c(); para34 <- c()
for (i in 1:length(dcc@mfit$R))
{
  para12 <- c(para12, dcc@mfit$R[[i]][1,2])
  para13 <- c(para13, dcc@mfit$R[[i]][1,3])
  para14 <- c(para14, dcc@mfit$R[[i]][1,4])
  para23 <- c(para23, dcc@mfit$R[[i]][2,3])
  para24 <- c(para24, dcc@mfit$R[[i]][2,4])
  para34 <- c(para34, dcc@mfit$R[[i]][3,4])
}

dcc_df <- data.frame(data[2:nrow(data),1], para12, para13, para14, para23, para24, para34)
colnames(dcc_df) <- c("Date", "SPX_DAX", "SPX_WIG20", "SPX_EURUSD", "DAX_WIG20", "DAX_EURUSD", 'WIG20_EURUSD')

write.csv(dcc_df, 'dcc_garch.csv', row.names=FALSE)