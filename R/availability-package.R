#' \pkg{availability}
#'
#' Estimating geographic space available to animals based on telemetry data
#'
#' @name availability
#' @author Ben Raymond \email{ben.raymond@@aad.gov.au}, Simon Wotherspoon, Ryan Reisinger
#' @docType package
#' @import assertthat testthat
#' @importFrom png readPNG
#' @importFrom geosphere destPoint distVincentyEllipsoid finalBearing distVincentySphere
#' @importFrom mvtnorm rmvnorm
#' @importFrom sf sf_project
#' @importFrom stats ar rnorm runif
#' @importFrom tmvtnorm rtmvnorm
"_PACKAGE"
