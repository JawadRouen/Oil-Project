####################################################################
###
###           LASSO SELECTION
###
####################################################################

Lassofiltring <-   function(trainprices_lag,
                            trainoil_lag_feats,
                            graph = FALSE, centered = FALSE) {
  ##### without diff ######
  myRMSE <- function (data, lev = NULL, model = NULL) {
    cost <- sqrt(mean((data[, "pred"] - data[, "obs"]) ^ 2))
    names(cost) <- 'myRMSE'
    return(cost)
  }
  library(caret)
  library(glmnet)
  glmnetcv <-
    cv.glmnet(as.matrix(trainoil_lag_feats),
              as.matrix(trainprices_lag),
              alpha = 1)
  #define params for lasso cross-validation
  grid <- expand.grid(alpha = 1, lambda = glmnetcv$lambda)
  lmfit <-
    train(
      trainprices_lag ~ .,
      data = trainoil_lag_feats,
      method = 'glmnet',
      metric = 'myRMSE',
      # used to select the optimal model
      maximize = FALSE,
      trControl = trainControl(
        method = 'timeslice',
        initialWindow = round(3 * length(trainprices_lag) / 4),
        fixedWindow = TRUE,
        summaryFunction = myRMSE,
        horizon = 2
      ),
      preProcess = c("center", "scale", "nzv"),
      tuneGrid = grid
    )
  
  if (graph) {
    plot(lmfit)
    plot(mod)
    matplot(mod$lambda, t(mod$beta), type = "l")
  }
  
  mod <- lmfit$finalModel
  # find the index of the best model
  ind <-
    which(round(mod$tuneValue$lambda - mod$lambda, digits = 6) == 0)
  
  #retrun selected covariates
  if (centered){    
    covs <-
      ts(mod$call$x[, which(mod$beta[, ind] > 0)],
         frequency = 12,
         start = start(trainoil_lag_feats))}
  else {
    covs <-
      ts(trainoil_lag_feats[, which(mod$beta[, ind] > 0)],
         frequency = 12,
         start = start(trainoil_lag_feats))}
  
  return(covs)}


####################################################################
###
###            CROSS VALIDATION FUNCTION
###
####################################################################

DynamicModel <- function(trainwindow, slide, h, CVdata,
                         order = NULL, drift = TRUE, method = "CSS-ML") {
  
  trainprices  <- CVdata$prices
  trainprices_lag <- CVdata$prices_lag
  trainoil_lasso_covs <- CVdata$trainoil_lasso_covs
  
  str <- tsp(trainprices_lag)[1]
  datamin  <- str  + trainwindow / 12
  n <- length(window(trainprices_lag, start = datamin)) - h + 1

  MmaeArima <- OmaeArima <- maeLM <- matrix(NA, n, h)
  MmaeDync <- OmaeDync <- matrix(NA, n, h)
  Mbench <- Obench <- matrix(NA, n, h)
  AICc <- rep(0, n)
  arimaOrder <- DynamicOrder <- matrix(NA, n, 3)
  colnames(DynamicOrder) <- c("AR", "Diff", "MA")
  colnames(arimaOrder) <- c("AR", "Diff", "MA")
  
  #dynCoef <- matrix(NA, n, dim(trainoil_lasso_covs)[2] + 1)
  dynCoef<- list()
  
  for (i in 1:n) {
    istart <- str + slide * (i - 1) / 12
    imid  <-  datamin + (i - 1) / 12
    iend <-   datamin - 0.01 + (i + h - 1) / 12
    
    train <-  window(trainprices, start = istart, end = imid - 0.01)
    test <-   window(trainprices, start = imid, end = iend)
    xlong <-  window(trainprices, start = istart, end = iend)
    
    prices_lag <-  window(trainprices_lag, start = istart, end = imid - 0.01)
    prices_lag_test <- window(trainprices_lag, start = imid, end = iend)
    prices_lag_long <- window(trainprices_lag, start = istart, end = iend)
    #### lasso here ??####
    lasso_lag <- window(trainoil_lasso_covs, start = istart, end = imid - 0.01)
    lasso_lag_test <- window(trainoil_lasso_covs, start = imid, end = iend)
    lasso_lag_long <- window(trainoil_lasso_covs, start = istart, end = iend)
    
    ######################
    #arima auto
    ######################
    Mrefit <-
      auto.arima(train, seasonal = c(0, 0, 0))  # redifine the model
    arimaOrder[i, ] <- arimaorder(Mrefit)
    
    Monefit <- Arima(xlong, model = Mrefit)
    #forecast multistep
    Mfcast <- forecast(Mrefit, h = h)
    #forecast one-step
    Monefcast <- fitted(Monefit)[-(1:length(train))]
    #mean absolute value multi step
    MmaeArima[i, ] <- abs(Mfcast$mean - test)
    #mean absolute value one step
    OmaeArima[i, ] <- abs(Monefcast - test)
    
    ######################
    #arima Dynamic lasso
    ######################
    try({
      if (is.null(order)) {
        fitlagDynamic <-
          auto.arima(
            prices_lag,
            xreg = lasso_lag,
            seasonal = c(0, 0, 0),
            allowdrift = TRUE
          )
      } else {
        fitlagDynamic <-
          Arima(
            prices_lag,
            xreg = lasso_lag,
            order = order,
            seasonal = c(0, 0, 0),
            include.drift = drift,
            method = method
          )
      }
      DynamicOrder[i, ] <- arimaorder(fitlagDynamic)
      dynCoef[[i]] <- fitlagDynamic$coef
      
      if (colnames(fitlagDynamic$xreg)[1] == "drift"){
        outlagDynamic  <-
          Arima(prices_lag_long, 
                xreg = cbind(drift=1:length(prices_lag_long),lasso_lag_long), model = fitlagDynamic)
      } else outlagDynamic  <-
        Arima(prices_lag_long, xreg = lasso_lag_long, model = fitlagDynamic)
      
      #forecast multistep
      Mfcast <- forecast(fitlagDynamic, xreg = lasso_lag_test, h = h)
      #forecast one-step
      Monefcast <- fitted(outlagDynamic)[-(1:length(prices_lag))]
      #mean absolute value multi step
      MmaeDync[i, ] <- abs(Mfcast$mean - prices_lag_test)
      #mean absolute value one step
      OmaeDync[i, ] <- abs(Monefcast - prices_lag_test)
      #AIC out of simple
      AICc[i] <- outlagDynamic$aicc
    })
    
    #model linear
    LMfit <- tslm(train ~ trend + season)
    LMfcast <- forecast(LMfit, h = h)
    maeLM[i, ] <- abs(LMfcast$mean - test)
    
    #bench
    Mbench[i, ] <- abs(prices_lag_test - prices_lag[length(prices_lag)])
    Obench[i, ] <-
      abs(prices_lag_test - c(prices_lag[length(prices_lag)], prices_lag_test[1:(h-1)]))
  }
  
  list(
    MmaeArima = MmaeArima, OmaeArima = OmaeArima,
    maeLM = maeLM,
    MmaeDync = MmaeDync, OmaeDync = OmaeDync,
    Mbench = Mbench, Obench = Obench,
    arimaOrder = arimaOrder, DynamicOrder = DynamicOrder,
    dynCoef = dynCoef, AICc = AICc
  )
}


####################################################################
###
###            LAUNCH CROSS VALIDATION FOR ARIMA
###
####################################################################

#load data jodi
result <- loadJODI(rep = "Data/jodi/") 
countrydata <- result$countrydata
dataoil <- result$dataoil
rm(result)

#load oil prices
price <- LoadPrices(rep = "Data/opec/")
####################################################################
####################################################################
### Train dataset will stop at Feb 2015
####################################################################
####################################################################
price <- price[Date<"2015-03-01"]
####################################################################

rm(result)
result <- list()
#sliding window - 0 for no and 1 for yes
slide <- 1
#number of step test
h <- 2
for (nbClust in c(8,10)){  #clustering params
  for (fClust in 2*(4:(nbClust/2+1))){
    #get country to keep in the model
    colcountries <- clustboot(nbClust, fClust)
    for (nblag in c(1,2,3,4,5)) {
      #lag params
      lag <- 1:nblag
      #preparing data
      CVdata <- PreDataset(price, dataoil, colcountries, lag)
      #launch lasso filtring
      CVdata$trainoil_lasso_covs  <-
        Lassofiltring(CVdata$prices_lag, CVdata$oil_lag_feats)
      print(paste("Post lasso dim: ",dim(CVdata$trainoil_lasso_covs)))
      #number of months od training
      for (trainwindow in c(120)){ 
        print(paste("nbClust: ",nbClust,", fClust: ",fClust, 
                    ",nblag: ",nblag, ",trainwindow: ",trainwindow))
        result[[length(result)+1]] <- list(CVoutput = DynamicModel(trainwindow, slide, h, CVdata),
                                           trainwindow= trainwindow, slide = slide, nblag = nblag,
                                           h = h, nbClust = nbClust, fClust = fClust)
      }
    }
  }
}

#save(result, file="CVresult.rda")
load(file = "CVresult.rda")


######################
# Analysis by params 
######################
#mean analysis
MDY <- sapply(result, function(l)mean(l$CVoutput$MmaeDync[,1]))
rAICc <- sapply(result, function(l)l$CVoutput$AICc)
MAuto <- sapply(result, function(l)mean(l$CVoutput$MmaeArima[,1]))
Ben <- sapply(result, function(l)mean(l$CVoutput$Mbench[,1]))
#median analysis
mdMDY <- sapply(result, function(l)median(l$CVoutput$MmaeDync[,1]))
mdMAuto <- sapply(result, function(l)median(l$CVoutput$MmaeArima[,1]))
mdBen <- sapply(result, function(l)median(l$CVoutput$Mbench[,1]))
xlag <- sapply(result, function(l)l$nblag)
xtrainwindow <- sapply(result, function(l)l$trainwindow)
xfClust <- sapply(result, function(l)l$fClust)
xnbClust <- sapply(result, function(l)l$nbClust)
nbitemCV <- sapply(result, function(l)dim(l$CVoutput$MmaeDync)[1])
DynDf <-  1:length(nbitemCV)
DynDf <- sapply(result,function(l)length(l$CVoutput$dynCoef[[1]]))

#Gaussian Case
BIC = function(n, p, mae){
  # n number of observation and p number of covs
  return(n*log(mean(mae^2)) + log(n)*p)
}
# for small number of samples
AICc = function(n, p, mae){
  # n number of observation and p number of covs
  return(n*log(mean(mae^2)) + 2*p + 2*p*(p+1)/(n-p-2))
}
VarDY <- sapply(result, function(l)mean(l$CVoutput$MmaeDync[,1]^2))
BICDY <- 1:length(nbitemCV)
AICcDY <- 1:length(nbitemCV)
BICDY <-  BIC(nbitemCV,DynDf,VarDY)
AICcDY <-  AICc(nbitemCV,DynDf,VarDY)
params <- data.frame(cbind(1:length(xtrainwindow),xtrainwindow, xlag, xnbClust, 
                           xfClust, nbitemCV, DynDf, MDY,mdMDY,BICDY,AICcDY, 
                           MAuto, mdMAuto, Ben, mdBen, pre = (Ben-MDY)/Ben))


#ploting Mae od Dync model per params values
p <- ggplot(data = data.frame(params), aes(x = xlag, y = MDY))
p <- p + geom_smooth(aes(colour = factor(xtrainwindow)))
p <- p + facet_grid(. ~ xfClust)
p

# we note a signification change at lag 20 (when the samples are significant).

##means absolute error comparaison
boxplot(x = rAICc[[1]],xlim = c(0.5, length(rAICc)+ 0.5), ylim = c(760,835))
for(i in 2:length(rAICc)) boxplot(rAICc[[i]], at = i, add = TRUE)
legend("topright",legend=c("AICc", "mean AICc"),
       col=c('black','red'),lty=1)
text(x =1:length(rAICc), y = sapply(rAICc, mean), labels = params[,"nbitemCV"], col = 'red')

library(rgl)
#xlag  xnbClust xfClust xtrainwindow
par3d("windowRect"= c(0,0,1100,800))
plot3d(params$xfClust, params$xlag, params$MDY, col = params$xnbClust,
       xlab ="final Clust", ylab = "nb lag", zlab = "Mean Absolut Error")
legend3d("topleft", c("8","10"), col = c(8,10), pch = 19, cex = 1.5)

#historam of params that beats the bench
par(mfrow = c(2,2))
barplot(table(params[params$MDY<params$Ben,]$xfClust), main = "number of final Clusts")
barplot(table(params[params$MDY<params$Ben,]$xnbClust), main = "number of interm Clusts")
barplot(table(params[params$MDY<params$Ben,]$xlag), main = "Included lags")
par(mfrow = c(1,1))

#ploting Mae od Dync model per params values
p <- ggplot(data = params, aes(x = xlag, y = MDY))
#p <- p + geom_smooth(aes(colour = factor(xtrainwindow)))
p <- p + geom_point(aes(colour = factor(xtrainwindow),size = nbitemCV/20))
p <- p + facet_grid(. ~ xfClust)
p

##means absolute error comparaison
boxplot(x = result[[params[1,1]]]$CVoutput$MmaeDync[,1],xlim = c(0.5, nrow(params)+ 0.5))
for(i in 2:nrow(params)) boxplot(result[[params[i,1]]]$CVoutput$MmaeDync[,1], at = i, add = TRUE)
points(params[,"MDY"], type = 'b', col = 'red')
points(params[,"Ben"], type = 'b', col = 'green')
legend("topright",legend=c("Dynamic Dyn auto lasso", "Bench"),
       col=c('red','green'),lty=1)
text(x =1:nrow(params), y = params[,"MDY"]-1, labels = params[,"nbitemCV"])

## AICc/BIC comparaison
par(mfrow = c(2,1))
plot(params[,"BICDY"], type = 'b', col = 4,  main = "BIC")
text(x =1:nrow(params), y = params[,"BICDY"]-10, labels = params[,"nbitemCV"])

plot(params[,"AICcDY"], type = 'b', col = 2, main = "AICc")
text(x =1:nrow(params), y = params[,"AICcDY"]-10, labels = params[,"nbitemCV"])
par(mfrow = c(1,1))




#####################################################################################
#Plot the comparaison for selected values of params
#####################################################################################
params[which.min(params$AICcDY),]
j <- params[which.min(params$AICcDY),1]


MmaeArima <- result[[j]]$CVoutput$MmaeArima
OmaeArima <- result[[j]]$CVoutput$OmaeArima
maeLM <- result[[j]]$CVoutput$maeLM
MmaeDync <- result[[j]]$CVoutput$MmaeDync
OmaeDync <- result[[j]]$CVoutput$OmaeDync
Mbench <- result[[j]]$CVoutput$Mbench
Obench <- result[[j]]$CVoutput$Obench
arimaOrder <- result[[j]]$CVoutput$arimaOrder
DynamicOrder <- result[[j]]$CVoutput$DynamicOrder

matplot(cbind(colMeans(MmaeArima,na.rm=TRUE),colMeans(MmaeDync,na.rm=TRUE), 
              colMeans(Mbench,na.rm=TRUE)), type = 'b', xlab="horizon", ylab="MAE", main="Multi-step")
legend("bottomright",legend=c("ARIMA auto", "Dynamic Dyn auto lasso", "Bench"), col=1:4,lty=1)

boxplot(cbind("Arima auto"=MmaeArima[,1], "Dynamic Dyn auto lasso"=MmaeDync[,1],"Bench"=Mbench[,1]), 
        ylab="MAE", main = "One-step MAE moving model")
abline(h =mean(Mbench[,1]), col = "red", lty = 2)

#summary
summary(cbind("Arima auto"=MmaeArima[,1],"ARIMA Dyn Man Lasso"=MmaeDync[,1],
              "Bench"=Mbench[,1]))

##### One-step MAE on time axe
par(mfrow = c(2,1))
matplot(cbind("Dynamic Dyn auto lasso"=MmaeDync[,1],"Bench"=Mbench[,1]),
        lty=1,type = 'b', ylab="MAE", main = "One-step MAE on time axe")
abline(h =0, col = "green", lty = 2)
legend("topleft",legend=c("Dynamic Dyn auto lasso", "Bench"), col=1:2,lty=1)
      # One-step MAE moving average on time axe
matplot(cbind("Dynamic Dyn auto lasso"=ma(MmaeDync[,1],5),"Bench"=ma(Mbench[,1],3)),
        lty=1,type = 'b', ylab="MAE", main = "One-step MAE on time axe")
abline(h =0, col = "green", lty = 2)
legend("topleft",legend=c("Dynamic Dyn auto lasso", "Bench"), col=1:2,lty=1)

summary(cbind("Arima auto"=MmaeArima[,1],"ARIMA Dyn Man Lasso"=MmaeDync[,1],
              "Bench"=Mbench[,1]))
par(mfrow = c(1,1))



#####  RMSE
boxplot(cbind("Arima auto"=sqrt(rowMeans(OmaeArima^2)),
              "Dynamic Dyn auto lasso"=sqrt(rowMeans(OmaeDync^2)),"Bench"=sqrt(rowMeans(Obench^2)))
        , ylab="RMSE", main = "One-step RMSE")

boxplot(cbind("Arima auto"=rowMeans(OmaeArima),
              "Dynamic Dyn auto lasso"=rowMeans(OmaeDync), "Bench"=rowMeans(Obench),
              "LM"=rowMeans(maeLM)), ylab="MAE", main = "One-step MAE")

### analysis of arima orders  & coefs #######################
matplot(1:nrow(DynamicOrder),DynamicOrder,type="l", main="auto Dynamic Arima")
legend("bottomright",legend=c("AR","Diff","MA"),col=1:3,lty=1)


########################################################################################################    
########################################################################################################    
###
### comparing slide with not slide/ slide without clust
###
########################################################################################################   
########################################################################################################

#preparing data for feats without clustering 
CVdata <- PreDataset(price, dataoil, colcountries = NULL,  lag = 1:result[[j]]$nblag)
#launch lasso filtring
CVdata$trainoil_lasso_covs  <-
  Lassofiltring(CVdata$prices_lag, CVdata$oil_lag_feats)
#plot predictors
NoClustresult <- DynamicModel(trainwindow = result[[j]]$trainwindow, 
                              slide = 1,
                              h = result[[j]]$h, 
                              CVdata, method = "CSS")


colcountries <- clustboot(nbClust = result[[j]]$nbClust,
                          fClust = result[[j]]$fClust)
#preparing data
CVdata <- PreDataset(price, dataoil, colcountries,  lag = 1:result[[j]]$nblag)
#launch lasso filtring
CVdata$trainoil_lasso_covs  <-
  Lassofiltring(CVdata$prices_lag, CVdata$oil_lag_feats)
#plot predictors
plot(cbind(CVdata$prices_lag,CVdata$trainoil_lasso_covs))
Slideresult <- DynamicModel(trainwindow = result[[j]]$trainwindow, 
                            slide = 1,
                            h = result[[j]]$h, 
                            CVdata, method = "CSS")
NoSlideresult <- DynamicModel(trainwindow = result[[j]]$trainwindow, 
                              slide = 0,
                              h = result[[j]]$h, 
                              CVdata, method = "CSS")
AR2result <- DynamicModel(trainwindow = result[[j]]$trainwindow, 
                          slide = 0,
                          h = result[[j]]$h, 
                          CVdata, order = c(2,0,0), drift = FALSE, method = "CSS")
ARIMA111result <- DynamicModel(trainwindow = result[[j]]$trainwindow, 
                               slide = 0,
                               h = result[[j]]$h, 
                               CVdata, order = c(1,1,1), drift = FALSE)



##Plot Multi-Step forecast
matplot(cbind(colMeans(Slideresult$MmaeDync,na.rm=TRUE),
              colMeans(NoSlideresult$MmaeDync,na.rm=TRUE), 
              colMeans(Slideresult$Mbench,na.rm=TRUE),
              colMeans(NoClustresult$MmaeDync,na.rm=TRUE),
              colMeans(AR2result$MmaeDync,na.rm=TRUE),
              colMeans(ARIMA111result$MmaeDync,na.rm=TRUE)
), type = 'b', 
xlab="horizon", ylab="MAE", main="Multi-step")
legend("topleft",legend=c("Sliding window", "Not slinding window", 
                              "Bench",  "No Clust", "AR 2", "ARIMA (1,1,1)"),
       col=1:6,lty=1)

#comparing AR 2, ARIMA (1,1,1)
summary(AR2result$MmaeDync[,1])
summary(ARIMA111result$MmaeDync[,1]) 


##Plot one step forecast
boxplot(cbind("Sliding window"=Slideresult$MmaeDync[,1],
              "Not slinding window"=NoSlideresult$MmaeDync[,1], 
              "AR 2"=AR2result$MmaeDync[,1], 
              "ARIMA (1,1,1)"=ARIMA111result$MmaeDync[,1], 
              "Bench"=Slideresult$Mbench[,1]),
        ylab="MAE", main = "One-step MAE")
abline(h =median(NoSlideresult$Mbench[,1]), lty = 2)

##Plot all one step errors MA
plot(ma(Slideresult$MmaeDync[,1],5), type = 'b', col = 'red')
points(ma(NoSlideresult$MmaeDync[,1],5), type = 'b', col = 'green')
points(ma(Slideresult$Mbench[,1],5), type = 'b', col = 'black')
points(ma(NoClustresult$MmaeDync[,1],5), type = 'b', col = 'blue')
legend("bottomleft", c("Sliding window", "Not slinding window", "Bench", "No Clust"),
       lty=1,col =c('red', 'green','black','blue'))
abline(h =0, lty = 2)


#dynCoef minus the intercept
dynCoef <- ARIMA111result$dynCoef
dynCoef_tab <- do.call("rbind", dynCoef)

ARIMA_tab <-  dynCoef_tab[,1:2]
Cov_tab <- dynCoef_tab[,3:ncol(dynCoef_tab)]

boxplot(ARIMA_tab,cex.axis = .5,main = "Coefs values out of CV", horizontal = TRUE)
text(x= colMeans(ARIMA_tab), y = 1:dim(ARIMA_tab)[2],
     label = colnames(ARIMA_tab), cex = 0.8)
abline(v =0, col = "blue", lty = 3)
boxplot(Cov_tab,cex.axis = .5,main = "Coefs values out of CV", horizontal = TRUE)
text(x= colMeans(Cov_tab), y = 1:dim(Cov_tab)[2] -0.3,
     label = colnames(Cov_tab), cex = 0.8)
abline(v =0, col = "blue", lty = 3)


###
### trainwindow = 120, lag = 3, nbClust = 8, fClust = 10, slide = 0, drift = FALSE
### The best model is the AR2 with sliding window and without drift (others have been tested)
###

#####################################################################################################
#####################################################################################################
#####################################################################################################
###
###           Other models
###
#####################################################################################################
#####################################################################################################
#####################################################################################################

trainprices_lag <- CVdata$prices_lag
trainoil_lag_feats <- CVdata$trainoil_lasso_covs
trainoil_all_feats <- CVdata$oil_lag_feats

myRMSE <- function (data, lev = NULL, model = NULL) {
  cost <- sqrt(mean((data[, "pred"] - data[, "obs"]) ^ 2))
  names(cost) <- 'myRMSE'
  return(cost)
}
myMAE <- function (data, lev = NULL, model = NULL) {
  cost <- abs(data[, "pred"] - data[, "obs"])
  names(cost) <- 'myMAE'
  return(cost)
}
data <- data.frame(na.omit(cbind(trainoil_lag_feats, y = trainprices_lag,
                                 lag(trainprices_lag,-1),lag(trainprices_lag,-2))))

####################################################################
###
###           linear model /equivilant to ARIMA
###
####################################################################
# testing the equivalence for both models
coef(Arima(trainprices_lag, order = c(2,0,0), method = "CSS"))
coef(lm(y~. , data[18:20]))

lmfit <-
  train(
    y ~ .,
    data = data[18:20],
    method = 'lm',
    metric = 'myMAE',
    # used to select the optimal model
    maximize = FALSE,
    trControl = trainControl(
      method = 'timeslice',
      initialWindow = 120,
      fixedWindow = TRUE,
      summaryFunction = myMAE,
      horizon = 1
    )
  )
lmfit
summary(cbind("lm" = lmfit$resample[,1],
"Dyn" = AR2result$MmaeDync[,1],
"Arima" = AR2result$MmaeArima[,1]))

####################################################################
###
###           Random forest
###
####################################################################
rffit <-
  train(
    y ~ .,
    data = data,
    method = 'rf',
    metric = 'myMAE',
    # used to select the optimal model
    maximize = FALSE,
    trControl = trainControl(
      method = 'timeslice',
      initialWindow = 120,
      fixedWindow = TRUE,
      summaryFunction = myMAE,
      horizon = 1
    )
  )
rffit
####################################################################
###
###           glmnet with lasso and interection
###
####################################################################
glminterfit <-
  train(
    y ~ .*.,
    data = data,
    method = 'glmnet',
    metric = 'myMAE',
    # used to select the optimal model
    maximize = FALSE,
    trControl = trainControl(
      method = 'timeslice',
      initialWindow = 120,
      fixedWindow = TRUE,
      summaryFunction = myMAE,
      horizon = 1
    ),
    preProcess = c("center", "scale")
  )
glminterfit
####################################################################
###
###           Radial smoother
###
####################################################################
kernelfit <-
  train(
    y ~ .,
    data = data,
    method = 'krlsRadial',
    metric = 'myMAE',
    # used to select the optimal model
    maximize = FALSE,
    trControl = trainControl(
      method = 'timeslice',
      initialWindow = 120,
      fixedWindow = TRUE,
      summaryFunction = myMAE,
      horizon = 1
    ),
    preProcess = c("center", "scale")
  )
kernelfit
####################################################################
###
###           GAM model
###
####################################################################
library(mgcv)
gamfit <-
  train(
    y ~ .,
    data = data,
    method = 'gam',
    metric = 'myMAE',
    # used to select the optimal model
    maximize = FALSE,
    trControl = trainControl(
      method = 'timeslice',
      initialWindow = 120,
      fixedWindow = TRUE,
      summaryFunction = myMAE,
      horizon = 1
    ),
    preProcess = c("center", "scale")
  )
gamfit

####################################################################
###
###           neurol model
###
####################################################################
netfit <-
  train(
    y ~ .,
    data = data,
    method = 'avNNet',
    metric = 'myMAE',
    # used to select the optimal model
    maximize = FALSE,
    trControl = trainControl(
      method = 'timeslice',
      initialWindow = 120,
      fixedWindow = TRUE,
      summaryFunction = myMAE,
      horizon = 1
    ),
    repeats = 5,
    trace = FALSE,
    preProcess = c("center", "scale")
  )
netfit
####################################################################
###
###           SVM model
###
####################################################################
svmRadialfit <-
  train(
    y ~ .,
    data = data,
    method = 'svmRadial',
    metric = 'myMAE',
    # used to select the optimal model
    maximize = FALSE,
    trControl = trainControl(
      method = 'timeslice',
      initialWindow = 120,
      fixedWindow = TRUE,
      summaryFunction = myMAE,
      horizon = 1
    ),
    trace = FALSE,
    preProcess = c("center", "scale")
  )
svmRadialfit



# base on Ar covs, build oder type of models  (kalman)
# include new data. text analysis ? 
# deeplearning
# review the countries clust every n step
# review the lasso every t step

