
#############################################################
####
####          Preparing data for models building
####
#############################################################

PreDataset <- function(price, oildataset, colcountries = NULL, 
                       lag = 0:4, cleaning = TRUE){
  
  #Keep only data of countries that are in the colcountries list
  if (!is.null(colcountries)) 
    oildataset <-   oildataset[country %in% colcountries]
  
  #### get mean value and last per month
  y <- price[,c("lastprice","mean"):=.(max(day)==day, mean(RBRTE)), by=.(year,month)][lastprice==TRUE, .(year,month,RBRTE,mean)]
  melt <- melt(oildataset, id.vars = c("country","year","month"))
  x <- dcast(melt, year+month ~ country+variable, value.var = "value")
  
  #### remove cols with more than 100% NA
  ind <- apply(x,2,function(l)all(is.na(l)))
  x[,names(which(ind)) := NULL]
  
  if (cleaning){
    #### remove cols with more than 40% NA
    ind <- apply(x,2,function(l)mean(is.na(l))) > 0.4
    x[,names(which(ind)) := NULL]
    
    ##### remove nearzerovar columns
    x[,names(x)[nearZeroVar(x)] := NULL]
    
    ##### remove correlated columns
    coral <- cor(x[,-c(1:2), with = FALSE], use = 'pairwise.complete.obs')
    x[,colnames(coral)[findCorrelation(coral, cutoff = .8)] :=NULL]
  }
  
  ##### merging predictors and outputs
  data <- merge(y,x, by.x = c("year","month"), by.y = c("year","month"))
  
  #### create the time series
  prices <- ts(data[,RBRTE], frequency = 12, start = c(2002,1))
  oil_feats <- ts(data[,-c(1:4), with = FALSE], frequency = 12, start = c(2002,1))
  
  ######################################################################################################
  # Lagged predictors
  ######################################################################################################
  
  Advert <-lapply(lag, function(v)lag(oil_feats,-1*v))
  oil_lag_feats <- do.call("cbind",Advert)
  if (cleaning)oil_lag_feats <- na.omit(oil_lag_feats)
  colnames(oil_lag_feats) <- paste(rep(colnames(oil_feats),length(lag)),
                                   rep(lag, each = dim(oil_feats)[2]), sep = "_Lag")
  
  prices_lag <-
    window(prices,
           start = start(oil_lag_feats),
           end = end(oil_lag_feats))
  
  oil_lag_feats <-
    window(oil_lag_feats,
           start = start(prices_lag),
           end = end(prices_lag))
  
  list(prices = prices,oil_feats = oil_feats , prices_lag=prices_lag, oil_lag_feats=oil_lag_feats)
}


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

