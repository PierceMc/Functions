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
LandscapeDataToArray <- function(data){
	landscapeSize <- sqrt(length(data[1,]))
	testarray <- array(NaN, c(landscapeSize, landscapeSize, nrow(data)))
	for(i in 1:nrow(data)){
		if(silent==F) print(i)
		testarray[,,i] <- as.matrix(square(data[i,]))
	}
	return(testarray)
}
