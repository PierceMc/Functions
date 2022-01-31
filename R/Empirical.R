#' Variable Selection in a Cox Proportional Hazards Model
#'
#' @param formula Formula object with the global model
#' @param data Data to be passed to coxph
#' @param ties Treatment of tied events. Default="exact'
#' @param direction Direction of selection. Default='backwards', options 'backwards'
#' @param Alpha P-value threshold. Default = 0.05
#' @return coxph model object of final model selected.
#' @export
CoxModelSelector <- function(formula, data, ties='exact', direction='backwards', Alpha=0.05){
	model.formula <- formula
	model <- coxph(model.formula, data, ties=ties)
	model.Coef <- coef(model)
	model.Variables <- names(model.Coef)

	dropped <- data.frame()
	while(length(dropped)<length(model.Variables)){
		print(paste0("Dropping length is ", length(dropped)))
		if(length(dropped) > 0){
			dropping.formula <- as.formula(paste("~. -", paste(dropped, collapse=" -")))
			tmp.formula <- update(model.formula, dropping.formula)
		} else { 
			tmp.formula = model.formula
		}
		model.tmp <- coxph(tmp.formula, data, ties=ties)
		model.CovMat <- vcov(model.tmp)
		model.Coef <- coef(model.tmp)
		dyn.Variables <- names(model.Coef)
		scores <- data.frame()
		for(i in 1:length(dyn.Variables)){
			W.test <- wald.test(model.CovMat,model.Coef,c(i))
			W.score <- W.test$result$chi2[1]
			W.p <- W.test$result$chi2[3]
			scores<- rbind(scores, cbind(Variable=dyn.Variables[i],Wald=W.score, p.value=W.p))
		}
		zk <- scores[which(scores[,2] == min(scores[,2])),]
		if(zk[,3]<=Alpha){
			return(model.tmp)
		} else {
			print(paste("Dropping" , zk[,1], "."))
			print(zk)
			print("-------------")
			dropped <- c(dropped, zk[,1])
		}
	}
	print("Failed to find model")

}


#' Find the first year with defoliation equal to or greater than some threshold
#'
#' @param x vector including defoliation values.
#' @param nyears time span to the data. Default=13
#' @param Defol Threshold value. Default=2
#' @param maxyear Final year of data collection. Default=2020
#' @return Year of first defoliation
#' @export 
FirstDefol <- function(x, nyears=13, Defol=75, maxyear=2020){
	years <- (maxyear-nyears):maxyear
		names=character()
		for(i in years){
			names=c(names, paste0('mean_defoliation', i))
		}
		cols=x[names]#c((length(x)-nyears):length(x))
		year=min(which(x[names] >= Defol))+maxyear-nyears
		if(year == Inf){
			year=maxyear+1
		}
		return(year)
}

#' Find the first year with cumulative defoliation equal to or greater than some threshold
#'
#' @param x vector including defoliation values.
#' @param nyears time span to the data. Default=13
#' @param maxyear Final year of data collection. Default=2020
#' @param MortDefol Threshold value. Default=15
#' @return Year of mortality
#' @export 
Mortality <- function(x, nyears=13, maxyear=2020, MortDefol=600){
	years <- (maxyear-nyears):maxyear
	names=character()
	for(i in years){
		names=c(names, paste0('mean_defoliation', i))
	}
	data <- x[names]
	for(i in 1:length(data)){
		if(sum(data[1:i])>=MortDefol){
			return(years[i])
		}
		if(i == length(data)){
			return(maxyear+1)
		}
	}
}

#' Prepare output of main.py for cox analysis
#'
#' @param x Data from main.py
#' @param finalyear Final year of data collection. Default=2020
#' @param strata Size of bins for stratification in metres. Default=23000
#' @return Year of first defoliation
#' @export 
CoxDataPrep <- function(x, finalyear=2020, strata=23000){
	sampleddata <- x
	sampleddata[is.na(sampleddata)]=0
	sampleddata <- na.omit(sampleddata)
	sampleddata$time <- apply(sampleddata, 1, FirstDefol)
	sampleddata$event <- as.numeric(!(sampleddata$time ==finalyear))
	sampleddata$RoundedDistToOrigin <- mround(sampleddata$DistToOrigin, strata)
	sampleddata$First <- apply(sampleddata, 1, FirstDefol, Defol=1)
	sampleddata$CumulativeDefoliation <- apply(sampleddata, 1, Mortality)
	sampleddata$TimetilMortality <- sampleddata$CumulativeDefoliation - sampleddata$First
	sampleddata$Mortality <- as.numeric(!sampleddata$CumulativeDefoliation ==finalyear)
	return(sampleddata)
}

#' Stepwise selection for AIC/BIC (Wrapper for MASS package stepAIC())
#'
#' @param x Model object
#' @param direction Direction for searching
#' @param test AIC or BIC
#' @param n Number of observations (for BIC)
#' @return AICModel
#' @export 
stepwise <- function(x, direction='backwards', test,n=NULL){
	if(test=="AIC"){
		k=2
	} else if(test=="BIC"){
		k=log(n)
	}
	formula=x$formula
	result=stepAIC(x, formula, direction=direction, k=k)
	return(result)	
}


