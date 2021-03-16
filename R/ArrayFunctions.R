#' Return number of slices in a 3D array
#'
#' @param x Three-Dimensional Array
#' @return The number of slices in the array
#' @export
nslice=function(x){
	return(length(x[1,1,]))
}

#' Convert Cellular model output to 3 dimensional array
#'
#' @param data Output file from Cellular Automaton (OutbreakLandscape)
#' @return 
#' @export
LandscapeDataToArray <- function(d,verbose=F){
	d=as.matrix(d)
	landscapeSize <- sqrt(length(d[1,]))
	testarray <- array(NaN, c(landscapeSize, landscapeSize, nrow(d)))
	for(i in 1:nrow(d)){
		if(verbose==T) print(i)
		testarray[,,i] <- as.matrix(square(d[i,]))
	}
	return(testarray)
}
