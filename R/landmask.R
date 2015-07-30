#' Initialize the land-masking function
#'
#' @return function
#'
#' @seealso \code{\link{surrogate_arsimulate}}
#'
#' @examples
#' lm=landmask_init() ## initialize land mask function
#' lm(c(100,-65)) ## test point lon,lat
#'
#' @export landmask_init

landmask_init=function() {
    land.mask=readPNG(system.file("extdata","land_mask-0.1-nosub.png",package="availability")) ## 0=land, 1=ocean
    land.lon=seq(from=-180,to=180,length.out=dim(land.mask)[2])
    land.lat=seq(from=0,to=-90,length.out=dim(land.mask)[1])
    test_point=function(pt) {
        lonidx=which.min(abs(land.lon-((pt[1]+180)%%360-180)))
        latidx=which.min(abs(land.lat-pt[2]))
        land.mask[latidx,lonidx]==1
    }
    test_point
}

## slightly faster, but has unexpected behaviour at edges, so don't use for now
##landmask_init=function() {
##    land.mask=readPNG(system.file("extdata","land_mask-0.1-nosub.png",package="availability")) ## 0=land, 1=ocean
##    land.mask=land.mask[nrow(land.mask):1,] ## reverse in latitude
##    land.lon=seq(from=-180,to=180,length.out=dim(land.mask)[2])
##    lon.bin=abs(land.lon[2]-land.lon[1])
##    land.lat=seq(from=-90,to=0,length.out=dim(land.mask)[1])
##    lat.bin=abs(land.lat[2]-land.lat[1])
##    test_point=function(pt) {
##        lonidx=.bincode((pt[1]+180)%%360-180,land.lon-lon.bin/2)
##        latidx=.bincode(pt[2],land.lat-lat.bin/2)
##        land.mask[latidx,lonidx]==1
##    }
##    test_point
##}

