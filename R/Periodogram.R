#' Calculate Standard Error
#'
#' @param x Vector to calculate se from
#' @return Standard Error
#' @export
se <- function(x) return(sd(x)/sqrt(length(x)))


#' Calculate Moving Average
#'
#' @param x Vector to calculate ma from
#' @param n window for moving average
#' @return Moving average
#' @export
ma <- function(x, n = 5){stats::filter(x, rep(1 / n, n), sides = 2)}

#' Calculate periodogram
#'
#' @param sbw Vector of values
#' @param smoothing window for moving average
#' @param spli should be 2
#' @param val Title for plot
#' @param averageline Line placement for an averageline 
#' @return periodogram plot (ggplot)
#' @export
periodogram<- function(sbw, smoothing=1, spli=2, val="", averageline=0){
	sbw = sbw-mean(sbw)
	sbw <- ma(sbw, smoothing)
	sbw <- sbw[-c(1:smoothing,(length(sbw)-smoothing):length(sbw))]
	Fs <- 1 
	N <- length(sbw)
	xdft <- (1/(N^(1/2)))*fft(sbw)
	p <- abs(xdft)^2
	p <- p[1:(N/spli+1)]
	p[2:(N/spli)] = spli*p[2:(N/spli)];
	freq<-seq(0,(Fs/spli),by=(Fs/N))
	P <- ggplot() + geom_line(aes(x=freq[2:(N/spli)], y=p[2:(N/spli)])) + xlab("Freq.") +ggtitle(val)
	if(data=F){
		return(P)
	} else {
		return(list(freq[2:(N/spli)], p[2:(N/spli)]))
	}
}

#' Find peak from periodogram
#'
#' @param sbw Vector of values
#' @param smoothing window for moving average
#' @param spli should be 2
#' @return Frequency of periodogram peak
#' @export
periodogrammax<- function(sbw, smoothing, spli){
	sbw = sbw-mean(sbw)
	sbw <- ma(sbw, smoothing)
	sbw <- sbw[-c(1:smoothing,(length(sbw)-smoothing):length(sbw))]
	Fs <- 1 #Sampling Rate
	N <- length(sbw)
	xdft <- (1/(N^(1/2)))*fft(sbw)
	p <- abs(xdft)^2
	p <- p[1:(N/spli+1)]
	p[2:(N/spli)] = spli*p[2:(N/spli)];
	freq<-seq(0,(Fs/spli),by=(Fs/N))

	opts=p[2:(N/spli)][-c(1:5)]
	res=freq[which(p[2:(N/spli)] == max(opts))]*N
	return(res[1])
}

#' Calculate periodogram 
#'
#' @param sbw Vector of values
#' @param smoothing window for moving average
#' @param spli should be 2
#' @param val Title for plot
#' @param averageline Line placement for an averageline 
#' @return periodogram plot (ggplot) wilth geom_line()
#' @export
periodogramlines<- function(sbw, smoothing, spli, val="", averageline=0){
	sbw = sbw-mean(sbw)
	sbw <- ma(sbw, smoothing)
	sbw <- sbw[-c(1:smoothing,(length(sbw)-smoothing):length(sbw))]
	Fs <- 1 #Sampling Rate
	N <- length(sbw)
	xdft <- (1/(N^(1/2)))*fft(sbw)
	p <- abs(xdft)^2
	p <- p[1:(N/spli+1)]
	p[2:(N/spli)] = spli*p[2:(N/spli)];
	freq<-seq(0,(Fs/spli),by=(Fs/N))
	#plot(freq[2:(N/spli)]*N,p[2:(N/spli)],type='l', xlab="Period (years)")
	P <- ggplot() + geom_line(aes(x=freq[2:(N/spli)], y=p[2:(N/spli)])) + xlab("Years") +ggtitle(val)
	return(P)
}

#' Calculate periodogram data
#'
#' @param sbw Vector of values
#' @param smoothing window for moving average
#' @param spli should be 2
#' @param val Title for plot
#' @param averageline Line placement for an averageline 
#' @return periodogram data
#' @export
periodogramdata<- function(sbw, smoothing, spli, val=""){
	sbw = sbw-mean(sbw)
	sbw <- ma(sbw, smoothing)
	sbw <- sbw[-c(1:smoothing,(length(sbw)-smoothing):length(sbw))]
	Fs <- 1 #Sampling Rate
	N <- length(sbw)
	xdft <- (1/(N^(1/2)))*fft(sbw)
	p <- abs(xdft)^2
	p <- p[1:(N/spli+1)]
	p[2:(N/spli)] = spli*p[2:(N/spli)];
	freq<-seq(0,(Fs/spli),by=(Fs/N))
	return(p[2:(N/spli)])
}

#' Calculate periodogram 
#'
#' @param sbw Vector of values
#' @param smoothing window for moving average
#' @param spli should be 2
#' @param val Title for plot
#' @param averageline Line placement for an averageline 
#' @return periodogram data structured for plotting
#' @export
periodogramplotdata<- function(sbw, smoothing, spli, val="", struc="list"){
	sbw = sbw-mean(sbw)
	sbw <- ma(sbw, smoothing)
	sbw <- sbw[-c(1:smoothing,(length(sbw)-smoothing):length(sbw))]
	Fs <- 1 #Sampling Rate
	N <- length(sbw)
	xdft <- (1/(N^(1/2)))*fft(sbw)
	p <- abs(xdft)^2
	p <- p[1:(N/spli+1)]
	p[2:(N/spli)] = spli*p[2:(N/spli)];
	freq<-seq(0,(Fs/spli),by=(Fs/N))
	if(struc == "list"){
		return(list(freq[2:(N/spli)], p[2:(N/spli)]))
	} else if(struc == "data.frame"){
		return(data.frame("Frequency"=freq[2:(N/spli)], "Power"=p[2:(N/spli)]))
	}

}

#' Calculate Average periodogram over whole landscape from aggregate data in 'Output files'
#'
#' @param DataDir Directory containing the data
#' @param Variable Name of the Variable in the filename
#' @param Values Vector of parameter values for the Moran's I to be calculated for
#' @param Species Species to use. Default=SBW
#' @param globr Number of Simulations to average
#' @param reps Length of simulations
#' @param burnin Burnin period
#' @param plotmax optional: Ymax height of plot
#' @param Names sometimes optional: Vection of names from data file
#' @param Title optional: Plot title
#' @param yl Remove y label
#' @param data - optional: Return plot data rather than grob
#' @param Density - optional: Value of total area for density calculations
#' @return Grob with the average peiodogram plotted if data = F, data.frame containing the data used to produce grob if data = T
#' @export
AveragePeriodogram <- function(DataDir, Variable, Value, Species, globr, reps, burnin, plotmax=0, Names=0, Title=T, yl=T, data=F, Density=F){
	if(globr < 1)stop('No data'); 
	tmpdata <- read.csv(paste0(DataDir, 'Output', Variable, Value, "run1.csv"), header=T)
	Hosts=F
	if(Species %in% c("Hosts", "Host", "host", "hosts")){Hosts = T}
	if(length(Names) > 1)(colnames(tmpdata)=Names)
	if(Hosts){
		tmp1=tmpdata[burnin:(reps-2),"Primary"]
		tmp2=tmpdata[burnin:(reps-2),"Secondary"]
		tmpdata <- tmp1+tmp2
	} else {
		tmpdata=tmpdata[burnin:(reps-2),Species]
	}
	if(is.numeric(Density)){
		print("Density")
		tmpdata <- tmpdata/Density
	}
	tmpPlotdata=periodogramplotdata(tmpdata, 1, 2)
	PlotData=tmpPlotdata[[1]]
	PlotData=cbind(PlotData, tmpPlotdata[[2]])
	for(i in 2:globr){
		tmpdata <- read.csv(paste0(DataDir, 'Output', Variable, Value, "run", i, ".csv"), header=T)
		if(length(Names) > 1)(colnames(tmpdata)=Names)
		if(Hosts){
			tmp1=tmpdata[burnin:(reps-2),"Primary"]
			tmp2=tmpdata[burnin:(reps-2),"Secondary"]
			tmpdata <- tmp1+tmp2
		} else {
			tmpdata=tmpdata[burnin:(reps-2),Species]
		}
		tmpPlotdata=periodogramdata(tmpdata, 1, 2)
		PlotData=cbind(PlotData, tmpPlotdata)
	}
	OutputPlotData=data.frame("Freq." = PlotData[,1], "Mean" = apply(PlotData[,-1], 1, mean, na.rm=T))
	sedat=data.frame(matrix(unlist(PlotData[,-1]), nrow=nrow(PlotData), byrow=T), stringsAsFactors=F)
	SE <-  apply(sedat, 1, se)
	SE[is.infinite(SE)] = 0
	OutputPlotData$SE <- SE
	if(plotmax == 0){
		OutputPlot <- ggplot(data=OutputPlotData, aes(y=Mean, x=Freq., ymin=Mean-SE, ymax=Mean+SE)) + geom_point() + geom_errorbar() + theme_bw() + ylab("Mean Power") + ggtitle(paste(Value))
	} else{
		OutputPlot <- ggplot(data=OutputPlotData, aes(y=Mean, x=Freq., ymin=Mean-SE, ymax=Mean+SE)) + geom_errorbar() + theme_bw() + ylab("Mean Power") + ggtitle(paste(Value)) + coord_cartesian(ylim=c(0,plotmax))
	}
	if(Title == F)(OutputPlot <- OutputPlot + ggtitle(""))
	if(yl == F)(OutputPlot <- OutputPlot + ylab(""))
	if(data==F){
		return(OutputPlot)
	} else {
		return(OutputPlotData)
	}

}





