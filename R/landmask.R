##' Initialize the land-masking function
##'
##' @return function
##'
##' @seealso \code{\link{surrogate_arsimulate}}
##'
##' @examples
##' mask  <- landmask_init() ## initialize land mask function
##' mask(0,c(100,-65)) ## test point lon,lat
##'
##' @export landmask_init

landmask_init=function() {
  land.mask=readPNG(system.file("extdata","land_mask-0.1-nosub.png",package="availability")) ## 0=land, 1=ocean
  land.lon=seq(from=-180,to=180,length.out=dim(land.mask)[2])
  land.lat=seq(from=0,to=-90,length.out=dim(land.mask)[1])
  function(tm,pt) {
    lonidx=which.min(abs(land.lon-((pt[1]+180)%%360-180)))
    latidx=which.min(abs(land.lat-pt[2]))
    land.mask[latidx,lonidx]==1
  }

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

##' A land mask based on Gebco 08.
##'
##' Generate a land mask function based on Gebco 08 topography.  The
##' mask is constant, the \code{tm} argument to the mask is
##' ignored. The \code{land} argument determines whether the mask
##' function returns \code{TRUE} or \code{FALSE} for land.
##'
##' @title Land Mask
##' @param path the path to a folder containing gebco_08.tif
##' @param land logical - the value to return for land.
##' @return a logical indicating whether the point is land or sea.
##' @importFrom raster raster writeRaster extract
##' @export
gebcoMask <- function(path="c:/Gebco/",land=FALSE) {
  if(!file.exists(file.path(path,"gebco_08.grd")))
    writeRaster(raster(file.path(path,"gebco_08.tif")),file.path(path,"gebco_08.grd"))

  gebco <- raster(file.path(path,"gebco_08.grd"))

  if(land)
    function(tm,pt) pt[2] < 90 & pt[2] > -90 & extract(gebco,cbind((pt[1]+180)%%360-180,pt[2])) >= 0
  else
    function(tm,pt) pt[2] < 90 & pt[2] > -90 & extract(gebco,cbind((pt[1]+180)%%360-180,pt[2])) <= 0
}
