#' @export
prepplotdata <- function(d, x, v, Temp=F){
	if(Temp==F){
		MoransData1 <- read.csv(paste0(d, "/MoransData", x, v[1],".csv"))
		MoransData3 <- read.csv(paste0(d, "/MoransData", x, v[2],".csv"))
		MoransData5 <- read.csv(paste0(d, "/MoransData", x, v[3], ".csv"))
	} else {
		MoransData1 <- read.csv(paste0(x, v[1],".csv"))
		MoransData3 <- read.csv(paste0(x, v[2],".csv"))
		MoransData5 <- read.csv(paste0(x, v[3], ".csv"))
	}

	MoransData1  <- cbind("V"=v[1], MoransData1[,-1])
	MoransData3  <- cbind("V"=v[2], MoransData3[,-1])
	MoransData5  <- cbind("V"=v[3], MoransData5[,-1])

	MoransData1$X=MoransData1$MoranDistance
	MoransData3$X=MoransData3$MoranDistance
	MoransData5$X=MoransData5$MoranDistance

	if(length(v) == 4){
		if(Temp==F){
			MoransData2 <- read.csv(paste0(d, "/MoransData", x, v[4], ".csv"))
		} else {
			MoransData2 <- read.csv(paste0(x, v[4], ".csv"))
		}
		MoransData2  <- cbind("V"=v[4], MoransData2[,-1])
		MoransData2$X=MoransData2$MoranDistance
		MoransData <- rbind(MoransData1, MoransData2, MoransData3, MoransData5)
	} else {
		MoransData <- rbind(MoransData1, MoransData3, MoransData5)
	}
	MoransData$V <- as.factor(MoransData$V)
	#names(MoransData) <- c("V","I","p","C","Cp","X")
	return(MoransData)
}

#' @export
prepplotdatalp <- function(d, x, v){
	out <- data.frame()
	for(i in v){
		tmp <- read.csv(paste0(d,'/MoransData', x, i, ".csv"))
		tmp  <- cbind("V"=i, tmp[,-1])
		tmp$X=tmp$MoranDistance
		out <- rbind(out, tmp)
	}
	out$V <- as.factor(out$V)
	return(out)
}
