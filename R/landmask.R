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
        lonidx=which.min(abs(land.lon-pt[1]))
        latidx=which.min(abs(land.lat-pt[2]))
        land.mask[latidx,lonidx]==1
    }
    test_point
}
