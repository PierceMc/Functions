#' Calculate Morans I with a least cost path weighted matrix
#'
#' @param Basedir Directory containing the data
#' @param Variable Name of the Variable in the filename
#' @param Values Vector of parameter values for the Moran's I to be calculated for
#' @param dclasses Distance classes for the correlogram
#' @param Species Species to use. Default=SBW
#' @param NSims How many model simulations. Default=32
#' @param Samples How many samples of the landscape to take. Default=32
#' @param N Points from the landscape to take each sample. Default=2000
#' @param Transition A function describing the cost of travel over non habitat patches. Default="transitionfunc"
#' @param ncores The number of cores available for calculation
#' @param verbose Boolean. Print status.
#' @return List of 2. Data frame with Moran's I correlogram data and List of example landscapes.
#' @export
Morans.LC <- function(BaseDir, Variable, Values, dclasses, Species='SBW', NSims=32, Samples=32, N=2000, Transition="transitionfunc", ncores=4, verbose=F){

	require(gdistance)
	require(doSNOW)
	require(parallel)
	require(doParallel)
	Output <- list()
	
	if(Transition != 'transitionfunc'){
		transitionfunc  <- Transition
	} else {
	}
	ms = Samples
	SampleLength <- 15

	MoransLeastCostWeighted <- data.frame()
	ExampleLandscapes <- list()
	ExamplePoints <- list()
	for(V in Values){
		for(s in 1:NSims){
			if(verbose==T) print(paste("Reading data for value", V, ", run", s))
			datafull <-read.csv(paste0(BaseDir, Variable, V, 'run', s, '.csv', Species, '.csv'))

			d <- datafull[1:SampleLength,]
			landscapeSize <- nrow(square(d[1,]))

			if(verbose==T) print('Creating 3D array')

			Array <- LandscapeDataToArray(d)

			if(verbose==T) print('Turn array into raster')

			flattenedArray <- apply(Array, c(1,2), sum)
			r1 <- raster(nrows=ncol(flattenedArray), ncols=ncol(flattenedArray), xmn=1, xmx=ncol(flattenedArray), ymn=1, ymx=landscapeSize)
			landscapeRaster <- raster(flattenedArray, xmn=1, xmx=ncol(flattenedArray), ymn=1, ymx=ncol(flattenedArray))
			landscapeRaster[landscapeRaster <1] = 0.1
			landscapeRaster[landscapeRaster >=1] = 1

			if(verbose==T) print("Raster to Transition layer")

			landscapeTransition <- transition(landscapeRaster, transitionFunction=get("transitionfunc"), directions=8, symm=F)


			landscapesquare <- Array[,,sample(c(1:nslice(Array)), 1)]
							  

			if(verbose==T) print("Converting square to xyz")

			landscapesmallxyz <- as.xyz(landscapesquare)
			landscapesmallxyz_no0 <- landscapesmallxyz[which(landscapesmallxyz[,3] != 0),]
			landscapesample_no0=landscapesmallxyz_no0[sample(c(1:nrow(landscapesmallxyz_no0)), N),]
			xy <- as.matrix(landscapesample_no0[,1:2])

			if(verbose==T) print('Running Parallel section...')
			cl <- makeCluster(ncores)
			registerDoParallel(cl)
			data.ms <- foreach(samples = c(1:ms), .combine=rbind, .packages=c('spdep', 'sp', 'gdistance')) %dopar% {

				source('~/Documents/Cellular/Analysis/Functions/MoranI.R', local=T)
				if(verbose==T) print(paste('Assessing sample', samples))
				landscapesample_no0=landscapesmallxyz_no0[sample(c(1:nrow(landscapesmallxyz_no0)), N),]
				xy <- as.matrix(landscapesample_no0[,1:2])
				df <- as.data.frame(landscapesample_no0[,3])
				names(df) <- "CD"
				landscapeSP <- SpatialPointsDataFrame(coords=xy, data=df)
				distanceasthewolfwalks <- costDistance(landscapeTransition, landscapeSP)
				datww.sym <- forceSymmetric(as.matrix(distanceasthewolfwalks), "L")
				w=as.matrix(distanceasthewolfwalks)
				Morans.LC.Res <- MoranI(w, df, dclasses, weight=T)
				Morans.LC.Res.DF <- cbind("distances"=Morans.LC.Res@distances, 'mI'=Morans.LC.Res@mI, 'number_pairs'=Morans.LC.Res@number_pairs, 'prob.adjust'=Morans.LC.Res@prob.adjust, 'varI'=Morans.LC.Res@varI); Morans.LC.Res.DF <- as.data.frame(Morans.LC.Res.DF); 
				Morans.LC.Res.DF$Simulation = s
				Morans.LC.Res.DF$Sample=samples
				Morans.LC.Res.DF$Value=V
				return(Morans.LC.Res.DF)
				names(Morans.LC.Res.DF) <- c('distances', 'mI', 'number_pairs', 'prob.adjust', 'varI', 'Simulation', 'Sample', 'Value')
			}
			stopCluster(cl)
			if(verbose==T) print('Saving')
			MoransLeastCostWeighted <- rbind(MoransLeastCostWeighted, data.ms)
			if(s == NSims){
				ExampleLandscapes[[as.character(V)]] <- landscapeRaster
				ExamplePoints[[as.character(V)]] <- xy
			}

		}
	}
		Output[["Morans.df"]] <- MoransLeastCostWeighted
		Output[["Landscapes"]] <- ExampleLandscapes
		Output[["Points"]] <- ExamplePoints
		return(Output)
}


transitionfunc <- function(j) ifelse(j[2]==1, 1, 0.4)

