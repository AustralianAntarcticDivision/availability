#' A land mask based on the GSHHS data set
#'
#' Generate a land mask function based on the Global Self-consistent, Hierarchical, High-resolution Geography Database.
#' The mask is provided at two (approximate) spatial resolutions: 0.1 degree and 0.05 degrees. The latter requires significantly more memory.
#' The mask is constant and the \code{tm} argument to the mask is ignored.
#'
#' @param res numeric: the spatial resolution of the mask, in degrees (either 0.1 or 0.05)
#' @return function that returns a logical indicating whether the point is at sea (TRUE) or on land (FALSE)
#' @seealso \code{\link{surrogateAR}}
#' @references Wessel, P., and W. H. F. Smith, A Global Self-consistent, Hierarchical, High-resolution Shoreline Database, J. Geophys. Res., 101, #B4, pp. 8741-8743, 1996. \url{https://www.ngdc.noaa.gov/mgg/shorelines/gshhs.html}
#' @examples
#' mask  <- gshhsMask() ## initialize land mask function
#' mask(0,c(100,-65)) ## test point lon,lat
#'
#' @export
gshhsMask <- function(res=0.1) {
    if (! res %in% c(0.1,0.05)) {
        res=0.1
    }
    land.mask=readPNG(system.file("extdata",paste0("land_mask_gshhs-",res,".png"),package="availability")) ## 0=land, 1=ocean
    if (length(dim(land.mask))>1) { land.mask=land.mask[,,1] }
    land.lon=seq(from=-180+res/2,to=180-res/2,length.out=dim(land.mask)[2])
    land.lat=seq(from=90-res/2,to=-90+res/2,length.out=dim(land.mask)[1])
    function(tm,pt) {
        lonidx=which.min(abs(land.lon-((pt[1]+180)%%360-180)))
        latidx=which.min(abs(land.lat-pt[2]))
        land.mask[latidx,lonidx]==1
    }
}

##' @rdname gshhsMask
##' @export
landmask_init <- gshhsMask

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


##' A land mask based on the full ETOPO1 dataset
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
