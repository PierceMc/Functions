####' Round
####'
####' @param x numeric to round.
####' @param base Value to round to.
####' @return Rounded value(s) of x
####' @export
mround <- function(x,base){
        base*round(x/base)
}
