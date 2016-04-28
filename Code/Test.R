####################################################################
###
###            Test on out of sample dataset
### 
####################################################################
# selected params : trainwindow = 120, lag = 3, nbClust = 8, 
#                   fClust = 10, slide = 1, drift = FALSE

#exec files that contain the main functions for data reading and processing
source('Code/process_jodi.R')
source('Code/functions.R')

## load data jodi
Test_result <- loadJODI(rep = "Data/Test_datasets/") 
Test_dataoil <- Test_result$dataoil
rm(Test_result)

## load oil prices
Test_price <- LoadPrices("Data/Test_datasets/")

## Set the params 
Test_blag <- 3    #result[[j]]$nblag
Test_h <- 2     
Test_trainwindow <- 120   #result[[j]]$trainwindow

## selection counctries for covariates with recomputing
#colcountries <- clustboot(nbClust = 8, fClust = 10)
## selection counctries for covariates without recomputing
Test_colcountries <- colcountries <- c("KAZAKHSTAN", "PHILIPP", "AZERBAIJAN", "NIGERIA", "CHINA",
                                       "IRAN", "IRAQ", "NORWAY", "SARABIA",  "USA")

####################################################################
####################################################################
### Test dataset will start at March 2015
####################################################################
####################################################################
Train_end <-  2015+ 1/12

#Test_Cols <- colnames(CVdata$trainoil_lasso_covs)
Test_Cols <- c("AZERBAIJAN_Production_Thousand Metric Tons (kmt)_Lag1",  
    "IRAN_Refinery intake_Thousand Metric Tons (kmt)_Lag1",   
    "IRAQ_Closing stocks_Thousand Metric Tons (kmt)_Lag1",    
    "NIGERIA_Refinery intake_Thousand Metric Tons (kmt)_Lag1",
    "NORWAY_Refinery intake_Thousand Metric Tons (kmt)_Lag1", 
    "PHILIPP_Imports_Thousand Metric Tons (kmt)_Lag1",        
    "PHILIPP_Refinery intake_Thousand Metric Tons (kmt)_Lag1",
    "SARABIA_Exports_Thousand Metric Tons (kmt)_Lag1",        
    "AZERBAIJAN_Production_Thousand Metric Tons (kmt)_Lag2",  
    "IRAN_Refinery intake_Thousand Metric Tons (kmt)_Lag2",   
    "IRAQ_Closing stocks_Thousand Metric Tons (kmt)_Lag2",    
    "KAZAKHSTAN_Imports_Thousand Metric Tons (kmt)_Lag2",     
    "NIGERIA_Production_Thousand Metric Tons (kmt)_Lag2",     
    "NIGERIA_Refinery intake_Thousand Metric Tons (kmt)_Lag2",
    "NORWAY_Refinery intake_Thousand Metric Tons (kmt)_Lag2", 
    "PHILIPP_Exports_Thousand Metric Tons (kmt)_Lag2",       
    "SARABIA_Exports_Thousand Metric Tons (kmt)_Lag2")


#preparing Test data
Test_CVdata <- PreDataset(Test_price, Test_dataoil, Test_colcountries,
                          lag = 1:Test_blag, cleaning = FALSE)
Test_CVdata$prices_lag <- 
  window(Test_CVdata$prices_lag, start = Train_end - (Test_trainwindow-1)/12)
Test_CVdata$oil_lag_feats <- 
  window(Test_CVdata$oil_lag_feats, start = Train_end - (Test_trainwindow-1)/12)
Test_CVdata$prices <- 
  window(Test_CVdata$prices, start = Train_end - (Test_trainwindow-1)/12)
Test_CVdata$oil_feats <- 
  window(Test_CVdata$oil_feats, start = Train_end - (Test_trainwindow-1)/12)

#Keep only selected cols from the training
Test_CVdata$trainoil_lasso_covs <- Test_CVdata$oil_lag_feats[,Test_Cols]


#run the test
Test_AR2result <- DynamicModel(trainwindow = Test_trainwindow, 
                          slide = 0,
                          h = Test_h, 
                          Test_CVdata, order = c(2,0,0), drift = FALSE, method = "CSS")

Test_ARIMA111result <- DynamicModel(trainwindow = Test_trainwindow, 
                               slide = 0,
                               h = Test_h, 
                               Test_CVdata, order = c(1,1,1), drift = FALSE)

Test_Dynamicresult <- DynamicModel(trainwindow = Test_trainwindow, 
                                    slide = 0,
                                    h = Test_h, 
                                    Test_CVdata)

##Plot Multi-Step forecast
matplot(cbind(colMeans(Test_AR2result$MmaeDync,na.rm=TRUE),
              colMeans(Test_ARIMA111result$MmaeDync,na.rm=TRUE), 
              colMeans(Test_Dynamicresult$MmaeDync,na.rm=TRUE),
              colMeans(Test_ARIMA111result$MmaeArima,na.rm=TRUE),
              colMeans(Test_ARIMA111result$Mbench,na.rm=TRUE)
), type = 'b', xlab="horizon", ylab="MAE", main="Multi-step")
legend("topleft",legend=c("Dynamic AR2",  "Dynamic ARIMA(1,1,1)" , "Auto Dynamic",
                              "Auto Arima", "Bench"),
       col=1:6,lty=1)

##Plot one step forecast
boxplot(cbind("AR 2"=Test_AR2result$MmaeDync[,1],
              "ARIMA(1,1,1)" = Test_ARIMA111result$MmaeDync[,1], 
              "Auto Dynamic" = Test_Dynamicresult$MmaeDync[,1], 
              "Bench"=Test_ARIMA111result$Mbench[,1]),
        ylab="MAE", main = "One-step MAE")
abline(h =median(Test_ARIMA111result$Mbench[,1]), lty = 2)


#######################################################
### Confidence intervals and hypothesis tests
#######################################################
t.test(Test_ARIMA111result$MmaeDync[,1],Test_ARIMA111result$Mbench[,1])
t.test(Test_ARIMA111result$MmaeDync[,2],Test_ARIMA111result$Mbench[,2])
#No conclusion is possible
