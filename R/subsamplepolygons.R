#' Return a rectangular sub-map of a polygon layer
#'
#' @param sizex X axis length
#' @param sizey Y axis length
#' @param polygons Polygon layer to subsample
#' @param rasterize [optional] Return a raster object insteard of a polygon
#' @param layer [for rasterize] Layer of the polygon to convert to raster
#' @return polygon or raster layer of size [siex, sizey] with the data from the input layer
#' @export
subsample <- function(sizex, sizey, polygons, rasterize=F, layer=""){
	mapsize <- extent(polygons)
	randomx <- sample(c(mapsize[1]:(mapsize[2]-sizex-1)),1)
	randomy <- sample(c(mapsize[3]:(mapsize[4]-sizey-1)),1)
	newextent <- extent(randomx, randomx+sizex, randomy, randomy+sizey)
	newpolygons <- crop(polygons, newextent)
	if(is.null(newpolygons)){
		newpolygons <- subsample(sizex, sizey, polygons)
	}
	if(rasterize==T){
		r <- raster(newpolygons, crs=proj4string(polygons), nrows=sizex,ncols=sizey)
		rs <- rasterize(newpolygons, r, newpolygons@data[layer])
		return(rs)
	} else {
		return(newpolygons)
	}
}
