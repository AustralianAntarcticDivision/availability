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

##' A land mask based on ETOPO1
##'
##' Generate a land mask function based on ETOPO1 topography. The
##' etopo geotiff is not bundled with the package and must be
##' downloaded from
##' \url{https://www.ngdc.noaa.gov/mgg/global/global.html}.
##'
##' When the mask is intially created, a native raster (grd,gri)
##' version of the geotiff is created in the directory \code{tmp},
##' which must be writable. This file can be deleted when the
##' computation is finished.
##'
##' The \code{land} argument determines whether the mask function
##' returns \code{TRUE} or \code{FALSE} for land. The mask is constant
##' and the \code{tm} argument to the mask is ignored.
##'
##' @title Land Mask
##' @param basename the name of the etopo geotiff (without file extension).
##' @param path the path to a folder containing the etopo geotiff
##' @param tmp the path to a writeable folder
##' @param land the logical value to return for land.
##' @return a logical indicating whether the point is land or sea.
##' @importFrom raster raster writeRaster extract
##' @export
etopoMask <- function(basename="ETOPO1_Bed_c_geotiff",path=".",tmp=path,land=FALSE) {
  tif <- file.path(path,paste0(basename,".tif"))
  grd <- file.path(tmp,paste0(basename,".grd"))

  if(!file.exists(grd)) writeRaster(raster(tif),grd)
  etopo <- raster(grd)

  if(land)
    function(tm,pt) pt[2] < 90 & pt[2] > -90 & extract(etopo,cbind((pt[1]+180)%%360-180,pt[2])) >= 0
  else
    function(tm,pt) pt[2] < 90 & pt[2] > -90 & extract(etopo,cbind((pt[1]+180)%%360-180,pt[2])) <= 0
}
