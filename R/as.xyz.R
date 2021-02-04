#' Converts a matrix of spatial values
#'
#' Takes a Matrix of spatial values and returns a 
#' data frame of x, y, and z (values).
#'
#' @param d Matrix of spatial values
#' @param raster Boolean. Inverts y values to account for raster standards
#' @return Data frame of x position, y position, and values (z).
#' @export
as.xyz=function(d, raster=T){
	rows=rep(c(1:nrow(d)), ncol(d))
	cols=rep(c(1:ncol(d)), each=nrow(d))
	rasterrows <- rep(c(nrow(d):1), ncol(d))
	values <- as.numeric(as.matrix(d))
	if(raster == T){
		output <- cbind(x=cols, y=rasterrows, z=values)
	} else {
		output <- cbind(x=cols, y=rows, z=values)
	}
	return(as.data.frame(output))
}
