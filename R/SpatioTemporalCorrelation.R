#' Calculate spatiotemporal cross correlation
#'
#' @param Array 3-D array of landscape snapshots (from OutbreakLandscape Files)
#' @param ms Number of smaples to take from the landscape
#' @param N Number of points in each sample
#' @param dclassess Distance classes to investigate
#' @param parallel Boolean. Calculate in parallel? Default=T
#' @param ncores Number of cores available for parallel computing. Default=4
#' @param verbose Print steps
#' @param timelag Increase timelag. Default=0
#' @param SecondArray 3-D array of landscape snapshots for another species (from OutbreakLandscape Files)
#' @return Data frame for correlogram of Rij
#' @export
rij.corr <- function(Array, ms, N, dclasses, parallel = T, ncores=4, verbose=F, timelag=0, SecondArray=NA){
	nd<-length(dclasses)
	tmpdata <- data.frame()
	starting=1
	if(!is.na(timelag)){
		starting=timelag+1
	}
	for(t in starting:nslice(Array)){
		if(verbose==T) print(paste("Time:", t))
		zt <- as.xyz(Array[,,t])
		xy <-zt[,1:2]
		df <- as.data.frame(zt[,3])
		sliceSP <- SpatialPointsDataFrame(coords=xy, data=df)
		names(sliceSP) <- 'z'
		if(!is.na(timelag)){
			ztminustl <- as.xyz(Array[,,(t-timelag)])
			xyminustl <-ztminustl[,1:2]
			dfminustl <- as.data.frame(ztminustl[,3])
			sliceSPminustl <- SpatialPointsDataFrame(coords=xyminustl, data=dfminustl)
			names(sliceSPminustl) <- 'z'
		}
		if(!is.na(SecondArray)){
			zt2 <- as.xyz(SecondArray[,,t])
			xy2 <-zt2[,1:2]
			df2 <- as.data.frame(zt2[,3])
			sliceSP2 <- SpatialPointsDataFrame(coords=xy2, data=df2)
			names(sliceSP2) <- 'z'
		}
		if(verbose==T) print(paste("Running in parallel with", ncores, "cores."))
		cl <- makeCluster(ncores)
		registerDoParallel(cl)
		if(is.na(SecondArray)){
			if(timelag==0){
				tmp <- foreach(s=c(1:ms), .combine=rbind, .packages=c("raster")) %dopar%{
					landscapesample=sliceSP[sample(c(1:nrow(sliceSP)), N),]
					internaltmp <- data.frame()
					for(d in 1:(length(dclasses)-1)){
						D <- dist(cbind(landscapesample$x, landscapesample$y))
						Dmat<-as.matrix(D)
						iHH<-(Dmat>dclasses[d]) & (Dmat<=dclasses[d+1])
						md<-mean(landscapesample$z)
						proddu<-(t(landscapesample$z-md) %*% iHH %*% (landscapesample$z-md))
						proddl<-(t((landscapesample$z-md)^2) %*% iHH %*% ((landscapesample$z-md)^2))
						internaltmp <- rbind(internaltmp, cbind("Dist"=dclasses[d], "Time"=t, "Sample"=s,  "upper"=proddu, "lower"=proddl))
						
					}
					return(internaltmp)
				}
				stopCluster(cl)
			} else  {
				tmp <- foreach(s=c(1:ms), .combine=rbind, .packages=c("raster")) %dopar%{
					sam <- sample(c(1:nrow(sliceSP)), N)
					landscapesampleitminuslag=sliceSPminustl[sam,]
					landscapesamplej=sliceSP[sam,]
					internaltmp <- data.frame()
					for(d in 1:(length(dclasses)-1)){
						D <- dist(cbind(landscapesampleitminuslag$x, landscapesampleitminuslag$y))
						Dmat<-as.matrix(D)
						iHH<-(Dmat>dclasses[d]) & (Dmat<=dclasses[d+1])
						mdi<-mean(landscapesampleitminuslag$z)
						mdj<-mean(landscapesamplej$z)
						proddu<-(t(landscapesampleitminuslag$z-mdi) %*% iHH %*% (landscapesamplej$z-mdj))
						proddl<-(t((landscapesampleitminuslag$z-mdi)^2) %*% iHH %*% ((landscapesamplej$z-mdj)^2))
						internaltmp <- rbind(internaltmp, cbind("Dist"=dclasses[d], "Time"=t, "Sample"=s,  "upper"=proddu, "lower"=proddl))
						
					}
					return(internaltmp)
				}
				stopCluster(cl)
			}
		tmpdata <- rbind(tmpdata, tmp)
		} else { 
			if(timelag==0){
				tmp <- foreach(s=c(1:ms), .combine=rbind, .packages=c("raster")) %dopar%{
					sam <- sample(c(1:nrow(sliceSP)), N)
					landscapesamplez=sliceSP[sam,]
					landscapesamplew=sliceSP2[sam,]
					internaltmp <- data.frame()
					for(d in 1:(length(dclasses)-1)){
						D <- dist(cbind(landscapesamplez$x, landscapesamplez$y))
						Dmat<-as.matrix(D)
						iHH<-(Dmat>dclasses[d]) & (Dmat<=dclasses[d+1])
						mdz<-mean(landscapesamplez$z)
						mdw<-mean(landscapesamplew$z)
						proddu<-(t(landscapesamplez$z-mdz) %*% iHH %*% (landscapesamplew$z-mdw))
						proddl<-(t((landscapesamplez$z-mdz)^2) %*% iHH %*% ((landscapesamplew$z-mdw)^2))
						internaltmp <- rbind(internaltmp, cbind("Dist"=dclasses[d], "Time"=t, "Sample"=s,  "upper"=proddu, "lower"=proddl))
						
					}
					return(internaltmp)
				}
			} else  {
				stop("Time lags on two variables is currently unsupported")
				#tmp <- foreach(s=c(1:ms), .combine=rbind, .packages=c("raster")) %dopar%{
				#	landscapesampleitminuslag=sliceSP[sample(c(1:(nrow(sliceSP)-timelag)), N),]
				#	landscapesamplej=sliceSP[sample(c(timelag:nrow(sliceSP)), N),]
				#	internaltmp <- data.frame()
				#	for(d in 1:(length(dclasses)-1)){
				#		D <- dist(cbind(landscapesample$x, landscapesample$y))
				#		Dmat<-as.matrix(D)
				#		iHH<-(Dmat>dclasses[d]) & (Dmat<=dclasses[d+1])
				#		md<-mean(landscapesample$z)
				#		proddu<-(t(landscapesampleitminuslag$z-md) %*% iHH %*% (landscapesamplej$z-md))
				#		proddl<-(t((landscapesampleitminuslag$z-md)^2) %*% iHH %*% ((landscapesamplej$z-md)^2))
				#		internaltmp <- rbind(internaltmp, cbind("Dist"=dclasses[d], "Time"=t, "Sample"=s,  "upper"=proddu, "lower"=proddl))
				#		
				#	}
				#	return(internaltmp)
				#}
			}
		tmpdata <- rbind(tmpdata, tmp)
		stopCluster(cl)
	}
	}
	names(tmpdata) <- c("Dist", "Time", "Sample", "upper", "lower")
	upper <- aggregate(upper~Dist+Time,tmpdata, mean)
	lower <- aggregate(lower~Dist+Time,tmpdata, mean)
	uppersd <- aggregate(upper~Dist+Time,tmpdata,sd)
	lowersd <- aggregate(lower~Dist+Time,tmpdata,sd)
	upper <- aggregate(upper~Dist,upper, sum)
	lower <- aggregate(lower~Dist,lower, function(x) sqrt(sum(x)))
	uppersd <- aggregate(upper~Dist,uppersd, sum)
	lowersd <- aggregate(lower~Dist,lowersd, function(x) sqrt(sum(x)))
	
	rij <- upper[,2]/lower[,2]
	SD <- (sqrt((uppersd[,2]/upper[,2])^2+(lowersd[,2]/lower[,2])^2))*rij

	return(data.frame("Dist"=dclasses[-length(dclasses)], "R"=rij, "SD"=SD))
	
}
