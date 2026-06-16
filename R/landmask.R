#' A land mask based on the GSHHS data set
#'
#' Generate a land mask function based on the Global Self-consistent, Hierarchical, High-resolution Geography Database.
#' The mask is provided at two (approximate) spatial resolutions: 0.1 degree and 0.05 degrees. The latter requires significantly more memory.
#' The mask is constant and the `tm` argument to the mask is ignored.
#'
#' @param res numeric: the spatial resolution of the mask, in degrees (either 0.1 or 0.05)
#' @param latmin numeric: southernmost latitude extent of the land mask
#' @param latmax numeric: northernmost latitude extent of the land mask
#' @return A function that returns a logical indicating whether the point is at sea (`TRUE`) or on land (`FALSE`)
#' @seealso [surrogateAR()], [surrogateAM()]
#' @references Wessel P, Smith WHF (1996) A Global Self-consistent, Hierarchical, High-resolution Shoreline Database. J. Geophys. Res. 101: 8741-8743. <https://www.ngdc.noaa.gov/mgg/shorelines/gshhs.html>
#' @examples
#' mask  <- gshhsMask() ## initialize land mask function
#' mask(0, c(100, -65)) ## test point lon,lat
#'
#' @export
gshhsMask <- function(res = 0.1, latmin = -90, latmax = 90) {
    if (!res %in% c(0.1, 0.05)) res <- 0.1
    land.mask <- readPNG(system.file("extdata", paste0("land_mask_gshhs-", res, ".png"), package = "availability")) ## 0 = land, 1 = ocean
    if (length(dim(land.mask)) > 1) land.mask <- land.mask[, , 1]
    land.lon <- seq(from = -180 + res / 2, to = 180 - res / 2, length.out = dim(land.mask)[2])
    land.lat <- seq(from = 90 - res / 2, to = -90 + res / 2, length.out = dim(land.mask)[1])
    function(tm, pt) {
        (pt[[2]] > latmin) & (pt[[2]] < latmax) & (land.mask[which.min(abs(land.lat-pt[[2]])), which.min(abs(land.lon-((pt[[1]]+180)%%360-180)))] > 0)
    }
}


##' @rdname gshhsMask
##' @export
landmask_init <- gshhsMask

## slightly faster, but has unexpected behaviour at edges, so don't use for now
##landmask_init <- function() {
##    land.mask <- readPNG(system.file("extdata", "land_mask-0.1-nosub.png", package = "availability")) ## 0 = land, 1 = ocean
##    land.mask <- land.mask[nrow(land.mask):1, ] ## reverse in latitude
##    land.lon <- seq(from = -180, to = 180, length.out = dim(land.mask)[2])
##    lon.bin <- abs(land.lon[2]-land.lon[1])
##    land.lat <- seq(from = -90, to = 0, length.out = dim(land.mask)[1])
##    lat.bin <- abs(land.lat[2]-land.lat[1])
##    test_point <- function(pt) {
##        lonidx <- .bincode((pt[1]+180)%%360-180, land.lon-lon.bin/2)
##        latidx <- .bincode(pt[2], land.lat-lat.bin/2)
##        land.mask[latidx, lonidx] > 0
##    }
##    test_point
##}


#' A land mask based on the ETOPO bathymetric dataset
#'
#' Generate a land mask function based on ETOPO topography. The ETOPO geotiff is not bundled with the package and must be
#' downloaded from <https://www.ncei.noaa.gov/products/etopo-global-relief-model>. You are free to choose the resolution appropriate to your usage (ETOPO 2022 is provided at 15-, 30-, and 60-arcsecond resolution). If you are working in the Southern Ocean you are advised to use the "Ice surface elevation geotiff" product rather than the bedrock elevation.
#'
#' The `land` argument determines whether the mask function returns `TRUE` or `FALSE` for land. The mask is constant and the `tm` argument to the mask is ignored.
#'
#' @title Land Mask
#' @param basename the name of the ETOPO geotiff (without file extension)
#' @param path the path to a folder containing the ETOPO geotiff
#' @param tif_path string: as an alternative to `basename` and `path`, provide the full path to the ETOPO geotiff (including file extension) in `tif_path`
#' @param land the logical value to return for land
#' @return A function that takes parameters `tm` and `pt` (locaton in longitude and latitude) and returns a logical value indicating whether the point is land or sea
#' @export
etopoMask <- function(basename = "ETOPO_2022_v1_30s_N90W180_surface", path = ".", tif_path, land = FALSE) {
    if (missing(tif_path)) tif_path <- file.path(path, basename, ".tif")
    etopo <- terra::rast(tif_path)
    if (land) {
        function(tm, pt) unname(pt[2] < 90 & pt[2] > -90 & terra::extract(etopo, cbind((pt[1] + 180) %% 360 - 180, pt[2])) >= 0)
    } else {
        function(tm, pt) unname(pt[2] < 90 & pt[2] > -90 & terra::extract(etopo, cbind((pt[1] + 180) %% 360 - 180, pt[2])) <= 0)
    }
}
