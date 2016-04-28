
## Data type names
## The file contains data for all the categories below
pdts <- rbind(
  c("Crude oil", "CRUDEOIL"),
  c("Natural Gas Liquid", "NGL"),
  c("Other", "OTHERCRUDE"),
  c("Total", "TOTCRUDE"))
colnames(pdts)<- c("name","code")

## Time series name
## This time series are given for all contries of the file scope

flw <- rbind(
  c("Production", "PRODREFOUT"),
  c("From other sources", "OTHSOURCES" ),
  c("Receipts", "PRECEIPTS"),
  c("Imports", "TOTIMPSB" ),
  c("Exports", "TOTEXPSB" ),
  c("Products transferred/Backflows", "PTRANSFBF"), 
  c("Products transferred", "TRANSF"),
  c("Direct use", "DIRECTUSE" ),
  c("Interproduct transfers", "INTPRODTRANSF"),
  c("Stock change", "STCHANAT" ),
  c("Statistical difference", "STATDIFF"), 
  c("Refinery intake", "REFOBSDEM" ),
  c("Demand", "REFOBSDEM"),
  c("Closing stocks", "CSNATTER")) 
colnames(flw)<- c("name","code")

## Units of the time series.
unt <- rbind(c("Thousand Barrels per day (kb/d)", "KBD"),
             c("Thousand Barrels (kbbl)", "KBBL"),
             c("Thousand Kilolitres (kl)", "KL"),
             c("Thousand Metric Tons (kmt)", "TONS"),
             c("Conversion factor barrels/ktons", "CONVBBL"))
colnames(unt) <- c("name","code")


####################################################################
###
###           LOAD OIL PRICES
###
####################################################################
#Loading the oil prices using the data.table   
LoadPrices <- function(rep =".", file = "daily_prices.csv"){
  library(data.table)
  price <- fread(paste0(rep,file))
  price[,c("month", "year", "day") := .(month(Date), year(Date), mday(Date))]
  
  price
  }

####################################################################
###
###           LOAD JODI DATA
###
####################################################################

loadJODI <- function(rep = ".", file = "world_Primary_CSV.csv"){
  library(data.table)
  library(ggplot2)
  library(caret)
  
  #word table
  #country(109)*Product(7)*Flows(6)*Month(12)*Year(14)*Unit(5)
  #prim table
  #country(109)*Product(4)*Flows(10)*Month(12)*Year(14)*Unit(5)
  #secondery table
  #country(109)*Product(9)*Flows(10)*Month(12)*Year(14)*Unit(5)
  
  #Colour Codes
  #1 Results of the assessment show reasonable levels of comparability
  #2 Consult metadata/Use with caution
  #3 Data has not been assessed
  #4 Data under verification
  
  
  #read data from the primary file
  prim <- fread(paste0(rep,file))
  #second <- fread("world_Secondary_CSV.csv")
  prim <- data.table(prim, key = "country,product,flow")

  
  #change labels
  for(i in 1:dim(flw)[1])prim[flow == flw[i,2], flow := flw[i,1]]
  for(j in 1:dim(unt)[1])prim[unit == unt[j,2], unit := unt[j,1]]
  
  ############################
  ## CRUDEOIL
  ############################
  
  #keep only crudeoil data
  prim <- prim[product=="CRUDEOIL"][,product:=NULL]
  
  #creating new cols
  prim[, c("month","year", "date") := .(substr(date,1,3), as.numeric(substr(date, 4,8)),NULL)]
  
  #"Statistical difference" has very few values
  #"From other sources" is empty
  #Products transferred/Backflows
  prim <- prim[flow !="From other sources"]
  prim <- prim[flow !="Statistical difference"]
  prim <- prim[flow !="Products transferred/Backflows"]
  
  #"Stock change" measures started in 2009
  
  # Convertion factor changes from a country to another 
  # and also between two flows of the same contry. 
  # exact number depends on the grade of crude oil(quality).
  # Removing Convesrion factor for the dataset
  prim <- prim[unit != "Conversion factor barrels/ktons"]
  
  #Keep the column that are less ferequently empty or equal to 0
  select <- prim[,mean(is.na(quantity) | quantity==0),by=c("flow", "unit")][, .SD[which.min(V1)], by = .(flow)]
  #setting a key
  setkey(prim,flow, unit)
  prim <- prim[.(select[,flow],select[,unit])]
  
  #remove country with no data 
  select <- prim[,mean(is.na(quantity) | quantity==0),by=.(country)][V1<0.8, country]
  prim <- prim[country %in% select]
  
  #change the names of months by the numbers
  monthNum <- function(m){
    return <- 1:length(m)
    for (i in 1:length(m)) {return[i] <- which(toupper(month.abb)==m[i])}
    return}
  prim[,month:=monthNum(month)]
  
  #tidy data
  dataoil <- dcast(prim, country+year+month ~ flow+unit, value.var = "quantity")
  setkey(dataoil,country, year, month)
  
  #################################################
  ### cleaning
  ##################################################
  library(forecast)
  #function when value eq 0 or NA take previous value 
  takeprev <- function(l){
    for(i in 2:length(l)) if(is.na(l[i]) | l[i]==0)l[i] <- l[i-1]
    l
  }
  
  # Detect extremes
  # Used to correct extreme values
  # can be reviewed
  alpha = 1
  extrems <- function(l){
    exts <- abs((ma(l,order = 7)-l)/sd(diff(l), na.rm = TRUE))
    is.exts <- !is.na(exts) & (exts > alpha) 
    prev <- takeprev(l)
    l[ l==0 & is.exts] <- prev[l==0 & is.exts]
    l[ is.na(l)] <- prev[is.na(l)]
    l
  }
  
  Chosen.cols <- names(dataoil)[-c(1:3)]
  dataoil[,(Chosen.cols) := lapply(.SD, extrems), by = country, .SDcols = Chosen.cols]
  
  #####################################################
  #create time series
  #####################################################
  require (zoo)
  lct <- Sys.getlocale("LC_TIME"); Sys.setlocale("LC_TIME", "C")
  
  create_ts <- function(SD)
    zoo(SD[,-c(1:3), with = FALSE],
        order.by = as.Date(paste(SD$year,SD$month, "01", sep = "/"), "%Y/%m/%d"))
  
  #apply to data
  countr_l <- split(dataoil, dataoil$country)
  countrydata <- lapply(countr_l, create_ts)
  
  list(countrydata=countrydata,dataoil=dataoil)
}

####################################################################
###
###           CLUSTERING of the OIL DATA
###
####################################################################

clustboot <- function(nbClust, fClust, graph = FALSE) {
  ######################################################################################
  ##        dynamic time wrapping
  #####################################################################################
  ######################################################################################
  
  ###########################################
  #distance and clusts on all cols
  ###########################################
  library(dtw)
  m <- length(countrydata)
  result <- matrix(nrow = m, ncol = m)
  for (i in 1:m)
    for (j in 1:m) {
      if (is.na(result[j, i]))
        try(result[i, j] <-
              dtw(dist(countrydata[[i]], countrydata[[j]], method = "Euclidean"),
                  distance.only = T)$normalizedDistance)
      else
        result[i, j] <- result[j, i]
    }
  
  # hierarchical clustering
  naVals <- unique(as.vector(which(is.na(result), arr.ind = TRUE)))
  names(countrydata)[naVals]
  
  if (length(naVals) > 0) {
    distMatrix <- log(as.dist(result[-naVals, -naVals]) + 1)
  } else {
    distMatrix <- log(as.dist(result) + 1)
  }
  
  hc <- hclust(distMatrix)
  if (graph) {
    if (length(naVals) == 0)
    {
      plot(hc , labels = names(countrydata))
    } else
      plot(hc , labels = names(countrydata)[-naVals])
  }
  cdwtglob <- cutree(hc, nbClust)
  dwtglob <- split(names(countrydata), cdwtglob)
  
  ################################################
  #distance and clusts only production and Export
  ################################################
  
  #Stock change_Thousand Metric Tons (kmt) starts at 2009
  library(dtw)
  m <- length(countrydata)
  distProd <- matrix(nrow = m, ncol = m)
  distExport <- matrix(nrow = m, ncol = m)
  for (i in 1:m)
    for (j in 1:m) {
      if (is.na(distProd[j, i]))
        try(distProd[i, j] <-
              dtw(ts(countrydata[[i]][, 5]), ts(countrydata[[j]][, 5]))$normalizedDistance)
      else
        distProd[i, j] <- distProd[j, i]
      if (is.na(distExport[j, i]))
        try(distExport[i, j] <-
              dtw(ts(countrydata[[i]][, 3]), ts(countrydata[[j]][, 3]))$normalizedDistance)
      else
        distExport[i, j] <- distExport[j, i]
    }
  
  # hierarchical clustering
  #prod
  distMatrix <- log(as.dist(distProd) + 1)
  hcProd <- hclust(distMatrix)
  if (graph)
    plot(hcProd , labels = names(countrydata))
  
  if (graph)
    plot(hcProd$height, type = "h")
  hdelta <-
    c(hcProd$height) - c(0, hcProd$height[-length(hcProd$height)])
  if (graph)
    points(hdelta, col = "red", type = "l")
  
  #chose the distance to cut
  cdwtprod <- cutree(hcProd, nbClust)
  dwtprod <- split(names(countrydata), cdwtprod)
  
  #Export
  distMatrix <- log(as.dist(distExport) + 1)
  hcExp <- hclust(distMatrix)
  if (graph)
    plot(hcExp , labels = names(countrydata))
  
  if (graph)
    plot(hcExp$height, type = "h")
  hdelta <- c(hcExp$height) - c(0, hcExp$height[-length(hcExp$height)])
  if (graph)
    points(hdelta, col = "red", type = "l")
  #chose the distance to cut
  cdwtexp <- cutree(hcExp, nbClust)
  dwtexp <-  split(names(countrydata), cdwtexp)
  
  ######################################################################################
  ######################################################################################
  ####                  wavelets
  ######################################################################################
  ######################################################################################
  library(wavelets)
  
  genwt <- function (l) {
    wt <- wavelets::dwt(ts(l[complete.cases(l[, -2]),-2]))
    as.vector(rbind(do.call("rbind", wt@W), wt@V[[wt@level]]))
  }
  
  
  wtdata <- sapply(countrydata, genwt)
  
  #wav <- prcomp(t(wtdata), scale = TRUE)
  #biplot(wav)
  
  hcwave <- hclust(dist(t(wtdata)))
  if (graph)
    plot(hcwave)
  
  cdistwav <- cutree(hcwave, nbClust)
  distwav <-  split(names(countrydata), cdistwav)
  
  ###########################################
  # hclust based on correlation /dist p
  ###########################################
  
  country_cor <- cor(wtdata)
  
  distMatrix <- as.dist(country_cor)
  hccor <- hclust(distMatrix)
  if (graph)
    plot(hccor , labels = names(countrydata))
  
  cdistcor <- cutree(hccor, nbClust)
  distcor <-  split(names(countrydata), cdistcor)
  
  ######################################################################################
  ######################################################################################
  ####                  Functional principal component analysis
  ######################################################################################
  ######################################################################################
  require(fda)
  
  #cut the 6 last month
  #dataoil <- dataoil[!(year==2015 & month>6)]
  
  ##prod data
  melt <- melt(dataoil[], id.vars = c("country", "year", "month"))
  melt <- melt[variable == "Production_Thousand Metric Tons (kmt)"]
  if (graph)
    qplot(
      x = as.Date(paste(year, month, "01", sep = "/"), "%Y/%m/%d"),
      y = value,
      color = country,
      data = melt,
      geom = "line"
    )
  Prod <- dcast(melt, year + month ~ country, value.var = "value")
  Prod <- as.matrix(Prod[, -c(1, 2), with = FALSE])
  
  #Export data
  melt <- melt(dataoil[], id.vars = c("country", "year", "month"))
  melt <- melt[variable == "Exports_Thousand Metric Tons (kmt)"]
  if (graph)
    qplot(
      x = as.Date(paste(year, month, "01", sep = "/"), "%Y/%m/%d"),
      y = value,
      color = country,
      data = melt,
      geom = "line"
    )
  Exp <- dcast(melt, year + month ~ country, value.var = "value")
  Exp <- as.matrix(Exp[, -c(1, 2), with = FALSE])
  
  nb <- dim(Prod)[1]
  datatime <- 0:(nb - 1)
  datarange <- c(0, nb)
  Prod_Exp_data <- array(NA,
                         dim = c(nb, dim(Prod)[2], 2),
                         dimnames =  list(datatime, colnames(Prod), c("Production", "Exports")))
  
  Prod_Exp_data[, , 1] <- Prod
  Prod_Exp_data[, , 2] <- Exp
  
  #set up the harmonic acceleration operator
  harmaccelLfd <- vec2Lfd(c(0, (2 * pi / nb) ^ 2, 0), rangeval = datarange)
  #create fourier basis of range "datarange" with 31 basis functions.
  #Recall that a fourier basis has an odd number of basis functions.
  databasis <- create.fourier.basis(datarange, nbasis = 31)
  
  #Choose level of smoothing using
  #the generalized cross-validation criterion (GCV)
  dataLoglam <- seq(-1, 3, .2)
  nglam   <- length(dataLoglam)
  
  dataSmoothStats <- array(NA,
                           dim = c(nglam, 3),
                           dimnames = list(dataLoglam, c("log10.lambda", "df", "gcv")))
  dataSmoothStats[, 1] <- dataLoglam
  
  # smoothing each individual curve
  # loop through smoothing parameters
  for (ilam in 1:nglam) {
    dataSmooth <- smooth.basisPar(
      datatime,
      Prod_Exp_data,
      databasis,
      Lfdobj = harmaccelLfd,
      lambda = 10 ^ dataLoglam[ilam]
    )
    dataSmoothStats[ilam, "df"]  <- dataSmooth$df
    dataSmoothStats[ilam, "gcv"] <- sum(dataSmooth$gcv)
    # note: gcv is a matrix in this case
  }
  
  #  set up plotting arrangements for one and two panel displays
  #  allowing for larger fonts
  if (graph) {
    par(mfrow = c(2, 1))
    plot(
      dataLoglam,
      dataSmoothStats[, "gcv"],
      type = "b",
      xlab = "Log_10 lambda",
      ylab = "GCV Criterion",
      main = "data Smoothing",
      log = "y"
    )
    plot(
      dataLoglam,
      dataSmoothStats[, "df"],
      type = "b",
      xlab = "Log_10 lambda",
      ylab = "Degrees of freedom",
      main = "data Smoothing"
    )
    par(mfcol = c(1, 1))
  }
  
  #keep the best value of lambda
  lambda <- as.numeric(names(which.min(dataSmoothStats[, "gcv"])))
  datafd <- smooth.basisPar(
    datatime,
    Prod_Exp_data,
    databasis,
    Lfdobj = harmaccelLfd,
    lambda = 10 ^ (lambda)
  )$fd
  
  names(datafd$fdnames) <-
    c("Normalized time", "Country", "Mesures")
  datafd$fdnames[[3]] <-
    c("Production_Thousand Metric Tons (kmt)",
      "Exports_Thousand Metric Tons (kmt)")
  
  if (graph) {
    par(mfrow = c(2, 1))
    plot(datafd, cex = 1.2)
    par(mfcol = c(1, 1))
  }
  # For more plots look at fda::data Code demonstrations
  
  #compute the mean functions (in frequence)
  datameanfd <- mean.fd(datafd)
  
  #  we need the values of the two mean functions also
  datameanvec <- eval.fd(datatime, datameanfd)
  
  #  plot these functions and their first two derivatives
  
  if (graph) {
    par(mfcol = c(2, 3))
    plot(datameanfd)
    plot(datameanfd, Lfdobj = 1)
    plot(datameanfd, Lfdobj = 2)
    par(mfcol = c(1, 1))
  }
  
  ##  --- Principal components analysis ---
  #  do the PCA with varimax rotation
  datafdPar  <- fdPar(databasis, harmaccelLfd, lambda = 10 ^ (lambda))
  # nharm	the number of harmonics or principal components to compute
  nbprs = 4
  datapca.fd <- pca.fd(datafd, nharm = nbprs, datafdPar)
  #varimax rotation
  datapca.fd <- varmx.pca.fd(datapca.fd)
  
  if (graph) {
    plot(datapca.fd$values, type = 'b') #the complete set of eigenvalues
    plot(datapca.fd$varprop, type = 'b') #the proportion of variance explained by each eigenfunction
  }
  
  #  compute the values of the harmonics at time values for each angle
  dataharmmat <- eval.fd(datatime, datapca.fd$harmonics)
  lbpcs <- paste("PC", 1:nbprs, sep = "")
  dimnames(dataharmmat) <- list(dimnames(Prod_Exp_data)[[1]],
                                lbpcs, dimnames(Prod_Exp_data)[[3]])
  Prodharmmat  = dataharmmat[, , 1]
  Expharmmat = dataharmmat[, , 2]
  
  #  we need the values of the two mean functions also
  datameanvec = eval.fd(datatime, datameanfd)
  Prodmeanvec  = datameanvec[, , 1]
  Expemeanvec = datameanvec[, , 2]
  
  if (graph) {
    par(mfrow = c(1, 2))
    matplot(Prodharmmat, Expharmmat, type = "l")
    matplot(Prodmeanvec, Expemeanvec, type = "l")
    
    par(mfrow = c(2, 2))
    plot.pca.fd(datapca.fd, cycle = TRUE)# points are the mean and + and - are the PCs * +/-values
    
    par(mfrow = c(1, 2))
    plot(
      datapca.fd$scores[, , 1],
      xlab = "PCA1",
      ylab = "PCA2",
      main = "Product"
    )
    text(datapca.fd$scores[, , 1], datafd$fdnames[[2]])
    plot(
      datapca.fd$scores[, , 2],
      xlab = "PCA1",
      ylab = "PCA2",
      main = "Exports"
    )
    text(datapca.fd$scores[, , 2], datafd$fdnames[[2]])
    par(mfrow = c(1, 1))
    
    plot(datameanvec[, 1, 1],
         main = "Product",
         pch = 19,
         type = 'b')
    plot(datameanvec[, 1, 2],
         main = "Export",
         pch = 19,
         type = 'b')
    par(mfrow = c(1, 1))
  }
  
  hcProdfda <- hclust(dist(datapca.fd$scores[, , 1]))
  if (graph)
    plot(hcProdfda , labels = names(countrydata))
  cProdfda <- cutree(hcProdfda, nbClust)
  Prodfda <-  split(names(countrydata), cProdfda)
  
  hcExpfda <- hclust(dist(datapca.fd$scores[, , 2]))
  if (graph)
    plot(hcExpfda , labels = names(countrydata))
  cExpfda <- cutree(hcExpfda, nbClust)
  Expfda <-  split(names(countrydata), cutree(hcExpfda, nbClust))
  
  ###########################################################################################
  #######           kmeans/medoid
  ###########################################################################################
  #find the countries that will be included in the prediction model
  #weigthed clustering ?
  library(fpc)
  library(cluster)
  
  pamk.best <- list()
  pamk.best[[1]] <-
    pamk(
      datapca.fd$scores[, , 1],
      krange = 1:10,
      criterion = "asw",
      critout = TRUE
    )
  #number of clusters estimated by optimum average silhouette width:
  if (graph)
    plot(pam(datapca.fd$scores[, , 1], pamk.best[[1]]$nc))
  pamk.best[[2]] <-
    pamk(
      datapca.fd$scores[, , 2],
      krange = 1:10,
      criterion = "asw",
      critout = TRUE
    )
  #number of clusters estimated by optimum average silhouette width:
  if (graph)
    plot(pam(datapca.fd$scores[, , 2], pamk.best[[2]]$nc))
  
  if (graph) {
    plot(
      datapca.fd$scores[, , 1],
      xlab = "PCA1",
      ylab = "PCA2",
      main = "Product"
    )
    points(
      x = pamk.best[[1]]$pamobject$medoids[, 1],
      y = pamk.best[[1]]$pamobject$medoids[, 2],
      pch = 17,
      col = "red"
    )
    text(datapca.fd$scores[, , 1], datafd$fdnames[[2]])
    
    plot(
      datapca.fd$scores[, , 2],
      xlab = "PCA1",
      ylab = "PCA2",
      main = "Exports"
    )
    points(
      x = pamk.best[[2]]$pamobject$medoids[, 1],
      y = pamk.best[[2]]$pamobject$medoids[, 2],
      pch = 17,
      col = "red"
    )
    text(datapca.fd$scores[, , 2], datafd$fdnames[[2]])
  }
  #colDync <- c(paste(datafd$fdnames[[2]][pamk.best[[1]]$pamobject$id.med], "Imports", sep = "_"),
  #paste(datafd$fdnames[[2]][pamk.best[[2]]$pamobject$id.med], "Production", sep = "_"))
  
  kmProdfda <- datafd$fdnames[[2]][pamk.best[[1]]$pamobject$id.med]
  kmExpfda <- datafd$fdnames[[2]][pamk.best[[2]]$pamobject$id.med]
  
  #include functions with different frequence mesures(bayes hiarchical)
  
  ##############################################################################
  #
  #
  #         CLUSTERING SUM-UP
  #
  #
  #############################################################################
  #lets now concider clusters as categories
  Clust <-
    cbind(cdwtglob,
          cdwtprod,
          cdwtexp,
          cdistwav,
          cdistcor,
          cProdfda,
          cExpfda)
  
  factor_Clust <-
    data.frame(apply(Clust, 2, as.character), stringsAsFactors = TRUE)
  rownames(factor_Clust) <- names(cdistcor)
  
  #Multiple Correspondence Analysis of all clusters results
  require(FactoMineR)
  res.mca <- MCA(factor_Clust, ncp = 7, graph = TRUE)
  if (graph)
    plot(res.mca,
         invisible = c("var", "quali.sup", "quanti.sup"),
         cex = 0.7)
  
  #kmoid
  pamk.best  <-
    pamk(
      as.data.frame(res.mca$ind),
      krange = fClust ,
      criterion = "asw",
      critout = TRUE
    )
  
  if (graph) {
    plot(
      as.data.frame(res.mca$ind)[, 1:2],
      xlab = "PCA1",
      ylab = "PCA2",
      main = "Clust"
    )
    points(
      x = pamk.best$pamobject$medoids[, 1],
      y = pamk.best$pamobject$medoids[, 2],
      pch = 17,
      col = "red"
    )
    text(as.data.frame(res.mca$ind)[, 1:2], names(cdistcor))
  }
  groupMcakmoid <-
    split(names(countrydata), pamk.best$pamobject$clustering)
  
  #kmod - kmeans like for categorical data
  #library(klaR)
  ## run algorithm on x:
  #cl <- kmodes(factor_Clust, modes = fClust, iter.max = 100)
  #groupKmod <-  split(names(countrydata), cl$cluster)
  
  
  #hclust
  distclust <- dist(
    Clust,
    method = function(x, y)
      sum(x != y)
  )
  if (graph)
    plot(hclust(distclust))
  
  ######
  ##### keep the columns
  ######
  
  colcountries <- rownames(pamk.best$pamobject$medoids)
  
  colcountries
}
