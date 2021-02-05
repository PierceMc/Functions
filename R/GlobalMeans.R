#' Summary statistics for global temporal data (Output...)
#' 
#' Returns summary statistics for outputs of the model from Output... files 
#' 
#' 
#' @param workingdir Directory the data is in
#' @param vari Variable
#' @param i Value
#' @param r Total number of simulations
#' @param fun Function
#' @return The number of slices in the array
#' @export
Means <- function(workingdir, vari, i, r, fun="mean"){
	allvals=matrix(2,2,2)
	y=1
	while(nrow(allvals) <= 10){
		if(burnin > 0){
		allvals= read.csv(paste0(workingdir,"Output", vari, i ,  "run", y, ".csv"))[-c(1:burnin),]
			#if(allvals[nrow(allvals),"Parasitoids"] == 0){
			#	allvals=matrix(2,2,2)
			#	if(y > r){
			#		break
			#	}
			#	y=y+1
			#	next
			#}
		} else {
		allvals= read.csv(paste0(workingdir,"Output", vari, i ,  "run", y, ".csv"))
		}
		prims=allvals[,1]
		secs=allvals[,2]
		nons=allvals[,3]
		sbws=allvals[,4]
		paras=allvals[,5]
		cds=allvals[,6]
		y=y+1
		if(y > r){
			break
		}
	}
	periods <- numeric()
	peak <- numeric()
	if(y<=r){
		for(x in y:r){
			if(burnin > 0){
			newvals= read.csv(paste0(workingdir,"Output", vari, i ,  "run", x, ".csv"))[-c(1:burnin),]
			#if(newvals[1,"SBW"] == 0 | newvals[1,"Parasitoids"] == 0){
			#	next
			#}
			} else {
			newvals= read.csv(paste0(workingdir,"Output", vari, i ,  "run", x, ".csv"))
			}

			if(nrow(newvals) == reps){
				prims=cbind(prims, newvals[,1])
				secs=cbind(secs, newvals[,2])
				nons=cbind(nons, newvals[,3])
				sbws=cbind(sbws, newvals[,4])
				paras=cbind(paras, newvals[,5])
				cds=cbind(cds, newvals[,6])
			}
			if(fun=="PeriodMax" || fun=="PeriodMaxPeak"){
				max <- periodogrammax(newvals$SBW, 1, 2)
				periods <- c(periods, max)
				peak <- c(peak, periodogramdata(newvals$SBW, 1, 2)[max])
			} 
		}

		if(is.null(ncol(prims))) return(c(0,0,0,0,0,0))
		
		primsmean=apply(prims, 1, mean)
		primsse=se(prims, r)
		secsmean=apply(secs, 1, mean)
		secsse=se(secs, r)
		nonsmean=apply(nons, 1, mean)
		nonsse=se(nons, r)
		sbwsmean=apply(sbws, 1, mean)
		sbwsse=se(sbws, r)
		parasmean=apply(paras, 1, mean)
		parasse=se(paras, r)
		cdsmean=apply(cds, 1, mean)
		cdsse=se(cds, r)

		globmax=max(globmax, max(allvals))
		if(fun=='mean'){
			MeanPrim= signif(mean(primsmean, na.rm=T), 2)
			MeanSec= signif(mean(secsmean, na.rm=T), 2)
			MeanNon= signif(mean(nonsmean, na.rm=T), 2)
			MeanSBW = signif(mean(sbwsmean, na.rm=T), 2)
			MeanParasitoid = signif(mean(parasmean, na.omit=T), 2)
			MeanCD = signif(mean(cdsmean, na.rm=T), 2)
			MeanPrimSD= round(TotEr(primsmean, primsse), 3)
			MeanSecSD=  round(TotEr(secsmean, secsse), 3)
			MeanNonSD=  round(TotEr(nonsmean, nonsse), 3)
			MeanSBWSD =signif(TotEr(sbwsmean, sbwsse), 3)
			MeanParasitoidSD = signif(TotEr(parasmean, parasse), 3)
			MeanCDSD = round(TotEr(cdsmean, cdsse), 3)
			return(data.frame(c('Primary', 'Secondary', 'NonHost', 'Herbivore', 'Parasitoid','Defoliation'), c(MeanPrim, MeanSec, MeanNon, MeanSBW, MeanParasitoid, MeanCD), c(MeanPrimSD, MeanSecSD, MeanNonSD, MeanSBWSD, MeanParasitoidSD, MeanCDSD)))
		} else if(fun == "PeriodMax"){
			PeriodMean=round(mean(periods, na.rm=T), 2)
			PeriodSD=round(sd(periods), 2)
			return(c(PeriodMean, PeriodSD))
		} else if(fun == "PeriodMaxPeak"){
			return(c(signif(mean(peak), 3), signif(sd(peak), 2)))
		}
	} else {
		return(c(NA, NA, NA, NA, NA, NA)) 
	}
}
