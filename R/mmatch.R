#' Matches multiple columns in a dictionary to return a value
#'
#' Like base match but accepts multiple 'from' arguments
#'
#' @param from1 Vector of from data
#' @param from2 Vector of from data
#' @param to What the data will become
#' @param datafrom1 vector containing values from from1 from the dataframe
#' @param datafrom2 vector containing values from from2 from the dataframe
#' @return Vector of target values based on from values
#' @export
mmatch <- function(from1, from2, to, datafrom1, datafrom2){
	d <- as.data.frame(cbind(from1=from1, from2=from2, to=to))
	output <- rep(NA,length(datafrom1))
	for(one in unique(from1)){
		    for(two in unique(from2)){
			    tr=which(datafrom1==one & datafrom2==two)
			    output[tr] <- d[which(from1==one & from2 == two),"to"]
		    }
	}
	return(output)
}
