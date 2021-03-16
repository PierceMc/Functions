#' Makes a line of data from OutbreakLandscape files into a landscape snapshot matrix
#'
#' @param d Single line from OutbreakLandscape file
#' @return Matrix of snapshot in tim from the model
#' @export
square=function(d, rows = 0, cols = 0){
	if(rows == cols){
		rows=sqrt(length(d))
		cols=sqrt(length(d))
	}

	output=matrix(d, rows, cols)
	#for(i in 1:rows){
		#nextline=d[1:cols]
		#names(nextline) <- c(1:rows)
		#output=rbind(output, nextline)
		#d=d[-(1:cols)]
	#}
	return(output)
}

#' Returns bins of sizes of patches.
#'
#' @param d Output from hosen()
#' @return Binned sizes of patches
#' @export
sizes <- function(d){
	output=cbind(1:max(d), 0)
	for(i in 1:max(d)){
		output[i,2] <- length(which(d == i))
	}
	return(output)
}

#' Identifies contiguous patches in the landscape
#'
#' @param squared Snapshot of landscape output
#' @return Matrix of labelled patches
#' @export
hosen=function(squared, mincd=0){
	squared <- as.matrix(squared)
	rows <- nrow(squared)
	cols <- ncol(squared)
	squared[which(squared < mincd)] = 0
	label=matrix(0, rows, cols)
	largest_label = 0;
	EqList=list()
	for(x in 1:cols){
	 for(y in 1:rows){
	   if(squared[y,x] >= 1){
	     left = squared[y,x-1]
	     above = squared[y-1,x]
	     if(x==1){left=0}
	     if(y==1){above=0}
	     if(left == 0 && above == 0){ #/* Neither a label above nor to the left. */
	       largest_label = largest_label + 1;        #/* Make a new, as-yet-unused cluster label. */
       		EqList[[largest_label]]=numeric()
	       label[y,x] = largest_label
	     } else if (left != 0 && above == 0) {#/* One neighbor, to the left. */
	       label[y,x] = label[y,x-1]
	     } else if (left == 0 && above != 0) {   #/* One neighbor, above. */
	       label[y,x] = label[y-1,x]
	     }  else{                               #         /* Neighbors BOTH to the left and above. */
	       label[y,x] = label[y-1,x]
       		leftlabel <- label[y,x-1]
	     	label[which(label == leftlabel)] = label[y,x]

       	       #if(!(label[y,x-1] %in% EqList[[label[y,x]]])){
       	       		#EqList[[label[y,x]]] <- c(EqList[[label[y,x]]], leftlabel)
	       #}
	     }
	   }
	 }
	}
	#if(length(EqList) > 1){
		#for(e in c(largest_label:1)){
			#label[which(label %in% EqList[[e]])] = e
		#}
	#}
	return(label)
}

#' Calculate periods over time of patch sizes
#'
#' @param OutbreakLandscape file
#' @return Periodogram data
#' @export
twodPeriodogram=function(d){
	rows=sqrt(length(d[1,]))
	cols=rows
	finaldata <- list()
	for(year in 1:(nrow(d)-1)){
		squareData <- square(d[year,], rows, cols)
		groupData <- hosen(squareData, rows)
		Sizes<- sizes(groupData)
		finaldata[[year]]<- periodogramdata(Sizes[,2], 1,2)
	}
	finaldata <- lapply(finaldata, FUN=function(x){x=c(x, rep(0, max(lengths(finaldata))-length(x)));return(x)})
	finaldata <- matrix(unlist(finaldata), nrow=length(finaldata), byrow=T)
	finaldata[which(is.na(finaldata))]  <- 0
	rownames(d)  <- c(1:nrow(d))
	return(finaldata)
}

#' Plot periods over time of patch sizes
#'
#' Don't really use anymore
#'
#' @param twodPeriodogram output
#' @return Plot
#' @export
plottwodPeriodogram <- function(d, type='image', theta = 320){
	if(type == 'image'){
		image(z=(d),x=c(1:nrow((d))),y=c(1:ncol(d)), xlab="Year", ylab="Outbreak Size")
	} else if( type == '3D') {
		persp3D(x=c(1:nrow(d)),y=c(1:ncol(d)), z=d, theta=theta, zlab="Frequency", ylab="Outbreak Size", xlab="Year", ticktype='detailed', xaxt='n')
	}
}

#' Histogram of patch sizes at one time step
#'
#' @param d OutbreakLandscape file
#' @param year row to generate histogram for
#' @param ret Return "plot" or "matrix"
#' @return Histogram of patch sizes
#' @export
SpatialHistogram <- function(d, year ,ret="plot"){
	rows=sqrt(length(d[1,]))
	cols=rows
	finaldata <- list()
	squareData <- square(d[year,], rows, cols)
	groupData <- hosen(squareData)
	Sizes<- sizes(groupData)
	if(ret == "plot") hist(Sizes[,2], breaks=seq(0,max(Sizes[,2])+2,1), xlab="Outbreak Size")
	if(ret == "matrix") return(Sizes)
}


#' Histogram of patch sizes at one time step from already squared data
#'
#' @param d Landscape snapshot
#' @return Histogram of patch sizes
#' @export
SpatialHistogramOnSquare <- function(d, ret="plot"){
	d <- hosen(d)
	Sizes<- sizes(d)

	output <- qplot(Sizes[,2], geom="histogram", fill=I("grey"), col=I("black")) + stat_bin(bins=10)#+ scale_x_continuous(trans="log10") + ylab("Count") + xlab("log Bin Size") + theme_bw() + stat_bin(bins=30)
	if(ret == "plot") hist(Sizes[,2], breaks=seq(0,max(Sizes[,2])+20,50), xlab="Outbreak Size")
	if(ret == "qplot") return(output)
}




