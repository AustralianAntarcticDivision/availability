#' \pkg{availability}
#'
#' Estimating geographic space available to animals based on telemetry data
#'
#' @name availability
#' @docType package
#' @import assertthat testthat
#' @importFrom png readPNG
#' @importFrom geosphere destPoint distVincentyEllipsoid finalBearing distVincentySphere
#' @importFrom mvtnorm rmvnorm
#' @importFrom sf sf_project st_as_sf st_coordinates st_crs `st_crs<-` st_is_longlat st_transform
#' @importFrom stats ar approxfun qnorm rnorm runif
#' @importFrom tmvtnorm rtmvnorm
"_PACKAGE"
