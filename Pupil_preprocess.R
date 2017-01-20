#this will run a general linear model on the data with multiple linear regressions (similar
#to the de Gee paper)

Pupil_preprocess <- function(pupil,pup_times=TRUE,pupil_times,eyeRM,timings,startRM,endRM,sampleRate,ign=FALSE,NonInterest,Zscore=FALSE,pupilRTdiv=TRUE,deTrend=TRUE,lowPass=TRUE,low,highPass=FALSE,high,cEyes=TRUE,output,plotIT=TRUE,indivPIRF=FALSE){
  library("stats")
  library("zoo")
  library("signal")
  library("contrast")
  library("multcomp")
  library("plyr")
  
  tms <- read.table(timings, sep=" ", header=F, fill=T,stringsAsFactors=F)
  eys <- read.table(pupil, sep=" ", header=F, fill=T,stringsAsFactors=F)
  botched <- read.table(eyeRM, sep=" ",header=F,fill=T,stringsAsFactors=F)
  if(pup_times){
    pupil_tms <- read.table(pupil_times, sep=" ",header=F,fill=T,stringsAsFactors=F)
  }
  if(ign){
    ignores <- read.table(NonInterest, sep="\t",head=F,fill=T)
  }
  
  #add the time component to the eye data
  if(!pup_times){
    tm <- "time"
    eys[,tm] <- NA
    lngth <- nrow(eys)
    rate <- 1/sampleRate
    tmp <- seq(0,lngth*rate-rate,rate)
    eys$time <- tmp
  }else{tm <- "time"; eys[,tm] <- pupil_tms$V1}
  
  #remove all data not wanted at start and end of data
  last <- eys$time[nrow(eys)]
  eys <- subset(eys, time > startRM & time < last - endRM)
  
  ##add a linear detrend to each pupil data so it is centered at 0
  if(deTrend){
    m1 <- lm(eys$V1~index(eys$V1))
    detr1 <- zoo(resid(m1), index(eys$V1))
    eys$V1 <- detr1
    m2 <- lm(eys$V2~index(eys$V2))
    detr2 <- zoo(resid(m2), index(eys$V2))
    eys$V2 <- detr2
  }
  
  ############################Low Pass FILTER
  ###need to build a butterworth low pass filter of 4Hz, so nothing that happens faster
  ##than 4 times per second is included
  ##these numbers are based on the de Gee 2014 paper
  ##altering the ord (order) of the filter alters how steep the pass filter is, 3rd order
  ##based on the de Gee paper
  if(lowPass){
    Fs = sampleRate
    Fn = Fs/2
    Ws = c(low)/Fn
    ord = 3
    bf <- butter(ord, Ws,type="low")
    eye1 <- eys$V1
    eye2 <- eys$V2
    x <- length(eye1)
    filtered1 <- filtfilt(bf,c(rep(eye1[1],50),eye1,rep(eye1[x],50)))
    x <- length(eye2)
    filtered2 <- filtfilt(bf,c(rep(eye2[1],50),eye2,rep(eye2[x],50)))
    x <- length(filtered1)
    filtered1 <- filtered1[51:(length(filtered1)-50)]
    x <- length(filtered2)
    filtered2 <- filtered2[51:(length(filtered2)-50)]
    filt_length <- length(filtered1)
    eys$V1 <- filtered1
    eys$V2 <- filtered2
  }
  ##need to make a high-pass filter
  if(highPass){
    Fs = sampleRate
    Fn = Fs/2
    Ws = c(high)/Fn
    ord = 3
    bf <- butter(ord, Ws,type="high")
    eye1 <- eys$V1
    eye2 <- eys$V2
    x <- length(eye1)
    filtered1 <- filtfilt(bf,c(rep(eye1[1],50),eye1,rep(eye1[x],50)))
    x <- length(eye2)
    filtered2 <- filtfilt(bf,c(rep(eye2[1],50),eye2,rep(eye2[x],50)))
    x <- length(filtered1)
    filtered1 <- filtered1[51:(length(filtered1)-50)]
    x <- length(filtered2)
    filtered2 <- filtered2[51:(length(filtered2)-50)]
    filt_length <- length(filtered1)
    eys$V1 <- filtered1
    eys$V2 <- filtered2
  }
  
  ##use scale to z-score the entire time series
  if(Zscore){
    eys$V1 <- scale(eys$V1)
    eys$V2 <- scale(eys$V2)
  }
   
  #set unwanted timepoints to be NA (i.e. blinks or looking away from the screen)
  l <- nrow(botched)
  for(i in 1:l){
    if(botched$V1[i] == "L"){
      eys$V1[eys$time >= botched$V2[i] & eys$time <= botched$V2[i]+botched$V3[i]] <- NA
    }
    if(botched$V1[i] == "R"){
      eys$V2[eys$time >= botched$V2[i] & eys$time <= botched$V2[i]+botched$V3[i]] <- NA
    }
  }
  
  #combine the two eyes if desired, if only one eye has data (non-NA) then use only that eye
  if(cEyes){
    h<-"used_dilation"
    eys[,h] <- 0
    l <- nrow(eys)
    for(i in 1:l){
      if(is.na(eys$V1[i]) & !is.na(eys$V2[i])){
        eys$used_dilation[i] <- eys$V2[i]
      }
      if(!is.na(eys$V1[i]) & is.na(eys$V2[i])){
        eys$used_dilation[i] <- eys$V1[i]
      }
      if(!is.na(eys$V1[i]) & !is.na(eys$V2[i])){
        eys$used_dilation[i] <- (eys$V1[i] + eys$V2[i])/2
      }
      if(is.na(eys$V1[i]) & is.na(eys$V2[i])){
        eys$used_dilation[i] <- NA
      }
    }
  }
  #if they don't want to combine the eyes use the left eye
  if(!cEyes){
    h<-"used_dilation"
    eys[,h] <- 0
    eys$used_dilation <- eys$V1
  }

  ##convolve all regressors with canonical pupil response function
  ##modification of the de Gee 2014 pupil impulse response function, they used a tmax
  ##of 930ms which translates into 111.6 time points when sampled at 120Hz
  ##can vary tmax and w if you want to do the same as the paper which sampled 
  ##tmax of 500, 650, 800, 1100, 1250, 1400 and w of 4, 6, 8, 12, 14, 16
  if(!indivPIRF){
  tmax <- sampleRate*0.93
  w <- 10.1
  pirf <- matrix()
  for(t in 1:(4*sampleRate)){
    tmp <- (t^w)
    tmp2 <- (exp(-t*w/tmax))
    pirf[t] <- tmp*tmp2
  }
  ##figure out max of pirf to divide by so it can be scaled to 1
  mx <- max(pirf)
  ##sclae down the impulse response function so it peaks at 1
  pirf <- pirf/mx
  #build the regressors the individual wants and include the trials in the formula
  l <- nrow(tms)
  for(i in 1:l){
    tmp <- eys
    tmp2 <- "tmp2"
    tmp[,tmp2] <- 0
    tmp$tmp2[tmp$time >= tms$V2[i] & tmp$time <= (tms$V2[i] + tms$V3[i])] <- 1*tms$V4[i]
    if(pupilRTdiv){
      tmp$tmp2 <- tmp$tmp2/tms$V3[i]
    }
    len <- length(tmp$tmp2)
    tmp3 <- zapsmall(conv(tmp$tmp2,pirf))
    #uncomment the line below to scale the Betas up
    #tmp3 <- tmp3/100
    tmp4 <- tms$V1[i]
    eys[,tmp4] <- tmp3[1:len]
    if(!exists("model_formula")){model_formula <- tmp4}
    else model_formula <- paste(model_formula,"+",tmp4,sep=" ")
  }
  #build the regressors of non-interest
  if(ign){
  l <- nrow(ignores)
  for(i in 1:l){
    tmp <- "ignore"
    eys[,tmp] <- 0
    eys$ignore[eys$time >= ignores$V1[i] & eys$time <= (ignores$V1[i] + ignores$V2[i])] <- 1
    len <- length(eys$ignore)
    tmp2 <- zapsmall(conv(eys$ignore,pirf))
    eys$ignore <- tmp2[1:len]
    model_formula <- paste(model_formula,"+",tmp,sep=" ")
  }}
  #create data to output to txt file for import into a group multiple linear regression model (similar to de Gee)
  lm.formula <- paste("used_dilation ~ ", model_formula,sep=" ")
  #run the model and exclude the NA points (because interpolation sucks)
  lm1 <- lm(formula = lm.formula, data = eys, na.action = na.exclude)
  lm.output <- data.frame(coef(lm1))
  model.summary <- lm1$model
  residuals <- lm1$residuals
  fitted <- lm1$fitted.values
  }
  ##if an individualized model fit is desired for the boxcar convolution
  if(indivPIRF){
    rm(model_fit)
    for(t_value in seq(0.5,1.5,0.05)){
      for(dubya in seq(4,15,0.5)){
        tmax <- sampleRate*t_value
        w <- dubya
        pirf <- matrix()
        for(t in 1:(10*sampleRate)){
          tmp <- (t^w)
          tmp2 <- (exp(-t*w/tmax))
          pirf[t] <- tmp*tmp2
        }
        ##figure out max of pirf to divide by so it can be scaled to 1
        mx <- max(pirf)
        ##sclae down the impulse response function so it peaks at 1
        pirf <- pirf/mx
        #build the regressors the individual wants and include the trials in the formula
        l <- nrow(tms)
        for(i in 1:l){
          tmp <- eys
          tmp2 <- "tmp2"
          tmp[,tmp2] <- 0
          tmp$tmp2[tmp$time >= tms$V2[i] & tmp$time <= (tms$V2[i] + tms$V3[i])] <- 1*tms$V4[i]
          if(pupilRTdiv){
           tmp$tmp2 <- tmp$tmp2/tms$V3[i]
          }
          len <- length(tmp$tmp2)
          tmp3 <- zapsmall(conv(tmp$tmp2,pirf))
          tmp4 <- tms$V1[i]
          eys[,tmp4] <- tmp3[1:len]
          if(!exists("model_formula")){model_formula <- tmp4}
          else model_formula <- paste(model_formula,"+",tmp4,sep=" ")
        }
        #build the regressors of non-interest
        if(ign){
        l <- nrow(ignores)
        for(i in 1:l){
          tmp <- "ignore"
          eys[,tmp] <- 0
          eys$ignore[eys$time >= ignores$V1[i] & eys$time <= (ignores$V1[i] + ignores$V2[i])] <- 1
          len <- length(eys$ignore)
          tmp2 <- zapsmall(conv(eys$ignore,pirf))
          eys$ignore <- tmp2[1:len]
          model_formula <- paste(model_formula,"+",tmp,sep=" ")
        }}
        #create data to output to txt file for import into a group multiple linear regression model (similar to de Gee)
        lm.formula <- paste("used_dilation ~ ", model_formula,sep=" ")
        #run the model and exclude the NA points (because interpolation sucks)
        lm1 <- lm(formula = lm.formula, data = eys, na.action = na.exclude)
        r <- summary(lm1)$r.squared
        if(!exists("model_fit")){model_fit <- data.frame(t_value,w,r)}else{
          z <- c(t_value,w,r); model_fit <- rbind(model_fit,z)}
  }}
  write.table(model_fit, file = paste(output,"model_fits.txt",sep=""), sep = " ", quote = F, append = F, na = "NA", col.names = F, row.names = F)
  highest_r <- max(model_fit$r)
  t_value <- model_fit$t_value[model_fit$r == highest_r]
  w <- model_fit$w[model_fit$r == highest_r]
  tmax <- sampleRate*t_value
  pirf <- matrix()
  for(t in 1:(4*sampleRate)){
    tmp <- (t^w)
    tmp2 <- (exp(-t*w/tmax))
    pirf[t] <- tmp*tmp2
  }
  used_model <- model_fit[model_fit$r == highest_r,]
  write.table(used_model, file = paste(output,"used_PIRF.txt",sep=""), sep = " ", quote = F, append = F, na = "NA", col.names = F, row.names = F)
  ##figure out max of pirf to divide by so it can be scaled to 1
  mx <- max(pirf)
  ##sclae down the impulse response function so it peaks at 1
  pirf <- pirf/mx
  #build the regressors the individual wants and include the trials in the formula
  l <- nrow(tms)
  for(i in 1:l){
    tmp <- eys
    tmp2 <- "tmp2"
    tmp[,tmp2] <- 0
    tmp$tmp2[tmp$time >= tms$V2[i] & tmp$time <= (tms$V2[i] + tms$V3[i])] <- 1*tms$V4[i]
    if(pupilRTdiv){
      tmp$tmp2 <- tmp$tmp2/tms$V3[i]
    }
    len <- length(tmp$tmp2)
    tmp3 <- zapsmall(conv(tmp$tmp2,pirf))
    tmp4 <- tms$V1[i]
    eys[,tmp4] <- tmp3[1:len]
    if(!exists("model_formula")){model_formula <- tmp4}
    else model_formula <- paste(model_formula,"+",tmp4,sep=" ")
  }
  #build the regressors of non-interest
  if(ign){
  l <- nrow(ignores)
  for(i in 1:l){
    tmp <- "ignore"
    eys[,tmp] <- 0
    eys$ignore[eys$time >= ignores$V1[i] & eys$time <= (ignores$V1[i] + ignores$V2[i])] <- 1
    len <- length(eys$ignore)
    tmp2 <- zapsmall(conv(eys$ignore,pirf))
    eys$ignore <- tmp2[1:len]
    model_formula <- paste(model_formula,"+",tmp,sep=" ")
  }}
  #create data to output to txt file for import into a group multiple linear regression model (similar to de Gee)
  lm.formula <- paste("used_dilation ~ ", model_formula,sep=" ")
  #run the model and exclude the NA points (because interpolation sucks)
  lm1 <- lm(formula = lm.formula, data = eys, na.action = na.exclude)
  lm.output <- data.frame(coef(lm1))
  model.summary <- lm1$model
  residuals <- lm1$residuals
  fitted <- lm1$fitted.values
  }
  
  if(plotIT){
    l <- ncol(eys)
    rm(tmp)
    for(i in 4:l){
      if(!exists("tmp")){tmp <- plot(eys[,i]~eys$time,type="l",xlab="Time",ylab="Model")}else{
        lines((eys[,i]/max(eys[,i])~eys$time),col="red")}
    }
  }
  
  write.table(residuals, file = paste(output,"residuals.txt",sep=""), sep = " ", quote = F, append = F, na = "NA", col.names = F, row.names = F)
  write.table(fitted, file = paste(output,"fitted_values.txt",sep=""), sep = " ", quote = F, append = F, na = "NA", col.names = F, row.names = F)  
  write.table(model.summary, file = paste(output,"model_summary.txt",sep=""), sep = " ", quote = F, append = F, na = "NA", col.names = F, row.names = F)
  write.table(lm.output, file = paste(output,"lm_outputs.txt",sep=""), sep = " ", quote = F, append = F, na = "NA", col.names = F, row.names = T)
  write.table(lm.formula, file = paste(output,"model_formula.txt",sep=""), sep = " ", quote = F, append = F, na = "NA", col.names = F, row.names = F)  
  
  list(lm.output)
  
}


