#' Fit first-order vector-autoregressive model to track
#'
#' @param x : a data.frame containing the observed track. It can optionally be an `sf` tibble (i.e. in projected coordinates). If not projected, it must contain columns named "lon" and "lat". In all cases it must contain a "date" column with timestamps in POSIX format
#' @param ... : additional arguments, currently not documented
#' @return An object of class "av_fit"
#' @seealso [surrogateAR()] for the original (differently-parameterized) version of this function, and [ar()] for the underlying function that fits the VAR model to the step speeds
#' @export
av_fit_ar <- function(x, ...) {
    dots <- list(...)
    if (inherits(x, "sf")) {
        stopifnot("x must be an sf-tibble with column 'date', or a data.frame with columns 'lon', 'lat', and 'date'" = "date" %in% names(x))
        xy <- st_coordinates(x)
        proj <- st_crs(x)
        ## distance increments dx, dy for each time step
        ## TODO deal with dateline wrapping
        dx <- diff(xy[, 1])
        dy <- diff(xy[, 2])
    } else if (is.data.frame(x)) {
        stopifnot("x must be an sf-tibble with column 'date', or a data.frame with columns 'lon', 'lat', and 'date'" = all(c("lon", "lat", "date") %in% names(x)))
        xy <- x[, c("lon", "lat")]
        proj <- NULL##st_crs("EPSG:4326")

        ## distance increments dx, dy for each time step
        nr <- nrow(xy)
        dx <- distVincentyEllipsoid(cbind(xy$lon[-1], xy$lat[-nr]), xy[-nr, ]) * sign(xy$lon[-1] - xy$lon[-nr])
        ## eastwards movement has positive sign

        ## the sign of dx will be wrong when the track crosses the date line
        ##  noting that we can have a track segment going from positive longitude (say 179E) to negative longitude (179W) crossing the date line, but also a segment crossing the zero line (1E to 1W) - the latter will be correct
        tempdx <- xy$lon[-nr] - xy$lon[-1] ## x difference on given coords
        tempdx2 <- xy$lon[-nr] - (xy$lon[-1] + 180) ## lon difference but with starting point shifted 180 degrees eastwards
        tempdx3 <- xy$lon[-nr] + 180 - xy$lon[-1] ## lon difference but with ending point shifted 180 degrees eastwards
        idx <- which((xy$lon[-nr] > 0 & xy$lon[-1] < 0 & abs(tempdx2) < abs(tempdx)) | (xy$lon[-nr] < 0 & xy$lon[-1] > 0 & abs(tempdx3) < abs(tempdx)))
        dx[idx] <- -dx[idx] ## change the sign of these

        dy <- distVincentyEllipsoid(cbind(xy$lon[-nr], xy$lat[-1]), xy[-nr, ]) * sign(xy$lat[-1] - xy$lat[-nr])
        ## northwards movement has positive sign

    } else {
        stop("x must be an sf-tibble with column 'date', or a data.frame with columns 'lon', 'lat', and 'date'")
    }

    ## convert to speed, units is m/s for lon-lat data or units/s for projected coordinates (in whatever units the projection uses)
    dt <- as.numeric(diff(x$date, units = "secs"))
    sx <- dx / dt
    sy <- dy / dt

    ## as-yet-undocumented transformation options
    if ("transform" %in% names(dots)) {
        if (identical(dots$transform, TRUE)) {
            ## use ecdf-quantile method to normalize x-speeds and y-speeds
            dots$transform <- list(
                structure(list(transform = approxfun(sx, ecdfq(sx), rule = 2), inverse = approxfun(ecdfq(sx), sx), rule = 2), class = "transform"),
                structure(list(transform = approxfun(sy, ecdfq(sy), rule = 2), inverse = approxfun(ecdfq(sy), sy), rule = 2), class = "transform")
            )
            dxdy <- data.frame(x = dots$transform[[1]]$transform(sx), y = dots$transform[[2]]$transform(sy))
        } else if (inherits(dots$transform, "transform")) {
            dxdy <- data.frame(x = dots$transform$transform(sx), y = dots$transform$transform(sy))
        } else if (is.list(dots$transform) && length(dots$transform) == 2 && all(vapply(dots$transform, inherits, "transform", FUN.VALUE = TRUE))) {
            dxdy <- data.frame(x = dots$transform[[1]]$transform(sx), y = dots$transform[[2]]$transform(sy))
        } else {
            warning("the `transform` parameter should be a transform object or a list of two such objects, ignoring")
            dxdy <- data.frame(x = sx, y = sy)
        }
    } else {
        dxdy <- data.frame(x = sx, y = sy)
    }
    ## Fit VAR1 to stepwise speeds
    structure(list(
        method = "var1",
        model = ar(dxdy, order.max = 1L, aic = FALSE),
        projection = proj,
        transform = dots$transform), ## will be NULL if not set
        class = "av_fit")
}


#' Simulate a track
#'
#' @param fit av_fit: an object of class `av_fit`, as returned by e.g. [av_fit_ar()]
#' @param x : a data.frame containing a "template" track. It can optionally be an `sf` tibble (i.e. in projected coordinates). If not projected, it must contain columns named "lon" and "lat". In all cases it must contain a "date" column with timestamps in POSIX format. This template is used to specify the projection and length of the simulated track, the location of fixed points (see `fixed`), and the time difference between successive points on the track
#' @param fixed logical: a vector (with length equal to the number of rows in `x`) indicating which locations in the template path are to be held fixed. By default the first and last points are fixed
#' @param point_check function: a function that accepts a parameters `tm` (date) and `pt` (location) and returns `TRUE` or `FALSE` indicating whether the location is acceptable. By default all points are accepted
#' @param partial logical: if `TRUE`, a partial track is returned if the sampling fails
#' @param random_rotation numeric: for the VAR-based method only. A two-element vector giving the upper and lower limits (radians) of the rotation applied to the model
#'
#' @return If `fit` was constructed with projected coordinates, an sf-tibble in the same projection with column `date`. Otherwise for unprojected (longitude-latitude) data, a data.frame with columns `lon`, `lat`, and `date`
#'
#' @export
av_sim <- function(fit, x, fixed, point_check, partial = FALSE, random_rotation) {
    stopifnot("`fit` should be an object of class \"av_fit\" as returned by e.g. `av_fit_ar()`, or an `ssm_df` object as returned by `aniMotum::fit_ssm()`" = inherits(fit, "av_fit") || inherits(fit, "ssm_df"))
    if (missing(fixed)) fixed <- rep(c(TRUE, FALSE, TRUE), c(1, nrow(x) - 2, 1))
    if (missing(point_check)) point_check <- function(...) TRUE
    if (inherits(fit, "ssm_df")) stop("av_sim_am not coded yet") ## return(av_sim_am(fit, x = x, fixed = fixed, point_check = point_check, partial = partial))
    ## if (!is.null(fit$projection) && !inherits(x, "sf")) stop("`x` should be an sf-tibble for simulating on projected coordinates")
    if (is.null(fit$projection) && inherits(x, "sf") && !st_is_longlat(x)) {
        ## we are asking for simulation on projected coords from a long-lat fit
        stop("cannot simulate on projected coordinates from a `fit` object that was fitted to longitude-latitude data")
    }
    if (fit$method == "var1") {
        av_sim_var1(fit, x = x, fixed = fixed, point_check = point_check, partial = partial, random_rotation = random_rotation)
    } else {
        stop("`fit` uses an unrecognized method: ", fit$method)
    }
}

av_sim_var1 <- function(fit, x, fixed, point_check, partial, random_rotation) {
    if (missing(random_rotation)) random_rotation <- c(-pi, pi)
    if (!is.null(fit$projection) && inherits(x, "sf") && !isTRUE(st_crs(x) == fit$projection)) warning("the projections of `fit` and `x` are different")
    ## ## lon-lat. Temporarily(?) call the old surrogateAR code
    ## mod <- fit$model
    ## attr(mod, "speed") <- TRUE ## fit$model is a speed-based one, this is just a temporary workaround to tell the sim code to scale by timestep to get from speed to distance
    ## attr(mod, "transform") <- fit$transform
    ## out <- surrogateAR(mod, xs = x[, c("lon", "lat")], ts = x$date, fixed = fixed, point.check = point_check, random.rotation = random_rotation, partial = partial)
    ## data.frame(lon = out$xs[, 1], lat = out$xs[, 2], date = x$date)
    av_sim_core(fit, x = x, fixed = fixed, point_check = point_check, random_rotation = random_rotation, partial = partial, speed = TRUE)
}

## av_sim_am <- function() {
## create fit object, create x if needed
        ## if (nrow(fit) > 1) {
        ##     warning("`fit` contains multiple models: only the first will be used")
        ##     fit <- fit[1, ]
        ## }

## @param fit
## @param ts the times at which the track is sampled
## @param xs the template sequence of states
## @param fixed a logical vector indicating which locations in the template path are to be held fixed.
## @param point.check function that accepts a state and returns boolean indicating whether the state is acceptable.
## @param random.rotation the upper and lower limits (radians) of the rotation applied to the VAR(1) model.
## @param partial if `TRUE`, a partial track is returned if the sampling fails.
## @param speed is `TRUE`, `fit` has been fitted to speeds not step lengths
## @return An array of states the define the simulated path.
av_sim_core <- function(fit, x, fixed, point_check, random_rotation, partial, speed) {
    proj <- NULL ## leave as NULL for simulation on long-lat coords
    if (inherits(x, "sf")) {
        xs <- st_coordinates(x)
        if (!st_is_longlat(x)) proj <- st_crs(x)
    } else if (is.data.frame(x)) {
        xs <- x[, c("lon", "lat")]
    } else {
        stop("`x` is an unexpected format")
    }
    ts <- x$date
    n <- nrow(xs)

    if (is.data.frame(xs)) xs <- as.matrix(xs)
    xs <- unname(xs[, 1:2, drop = FALSE])

    dt <- if (inherits(ts, "POSIXct")) c(0, as.numeric(diff(ts, units = "secs"))) else c(0, diff(ts))

    outf <- function(xs, ts) {
        if (is.null(proj)) {
            data.frame(lon = xs[, 1], lat = xs[, 2], date = ts)
        } else {
            out <- st_as_sf(data.frame(x = xs[, 1], y = xs[, 2], date = ts), coords = c("x", "y"))
            st_crs(out) <- proj
            out
        }
    }

    ## Simulate forward from k0.  Returns the index of the last fixed point reached if an acceptable next candidate cannot be found.
    ## parms will be a method-specific list of model parameters
    sample <- function(k0, parms) {
        if (fit$method == "var1") {
            ## initialize z
            z <- double(2)
            for (i in 1:100) z <- drop(parms$A %*% z) + drop(rnorm(2) %*% parms$U)
        }

        ## Simulate forward from k0
        pos <- xs[k0, ] ## starting position
        k <- k0 + 1

        ## Find remaining fixed points
        kfixed <- if(k <= n) (k:n)[fixed[k:n]] else integer(0)

        while (k <= n) {
            if (fixed[k]) {
                ## Skip fixed points
                while (k <= n && fixed[k]) {
                    pos <- xs[(k0 <- k), ]
                    k <- k + 1L
                }
                ## Find any remaining fixed points
                kfixed <- if (k <= n) (k:n)[fixed[k:n]] else integer(0)
                if (fit$method == "var1") {
                    ## re-initialize z
                    z <- double(2)
                    for (i in 1:100) z <- drop(parms$A %*% z) + drop(rnorm(2) %*% parms$U)
                }
            } else {
                ## Try at most 100 new candidates
                for (r in 1:100) {
                    if (fit$method == "var1") {
                        ## Trial z - simulate from centred VAR(1) model
                        z1 <- drop(parms$A %*% z) + drop(rnorm(2) %*% parms$U)
                        step <- z1 + parms$mu
                    } else {
                        stop("method not coded")
                    }
                    if (!is.null(fit$transform)) {
                        if (is.list(fit$transform)) {
                            step <- c(fit$transform[[1]]$inverse(step[1]),
                                      fit$transform[[2]]$inverse(step[2]))
                        } else {
                            step <- fit$transform$inverse(step)
                        }
                    }
                    if (speed) step <- step * dt[k]
                    if (is.null(proj)) {
                        ## take step in lon then lat
                        x1 <- as.vector(destPoint(destPoint(pos, 90, step[1]), 0, step[2]))
                    } else {
                        x1 <- as.vector(pos + step)
                    }
                    ## if we have a fixed point coming up, nudge towards it
                    ## TODO handle dateline if we are simulating in projected coords
                    if (length(kfixed)) {
                        this_nudge <- xs[kfixed[1], ] - x1
                        this_nudge[1] <- angle_normalise(this_nudge[1] / 180 * pi) / pi * 180
                        x1 <- x1 + (this_nudge) / (kfixed[1] - k + 1L)
                    }
                    ## Test current candidate
                    if (point_check(ts[k], x1)) {
                        ## Accept candidate
                        z <- z1
                        xs[k, ] <<- (pos <- x1)
                        k <- k + 1L
                        break
                    } else {
                        ## On failure return last fixed point
                        if (r == 100L) return(k0)
                    }
                }
            }
        }
        ## Return n+1 on success
        return(n + 1L)
    }

## `fit` has structure
##   list(method = "var1" or "am"
##        model = <AR model or ssm_df>,
##        projection = st_crs OR NULL,
##        transform = <transform object or NULL>)

    if (fit$method == "var1") {
        ## Try 100 rotations of the model
        for (roti in 1:100) {
            theta <- if (is.null(random_rotation)) 0 else runif(1, random_rotation[1], random_rotation[2])
            model0 <- rotateVAR1(fit$model, theta)
            parms <- list(A = unname(model0$ar[1, , ]), U = chol(unname(model0$var.pred)), mu = as.vector(model0$x.mean))
            k <- 1L
            fails <- 0
            for (i in 1:50) {
                knew <- if (i < 25) sample(k, parms) else sample(1, parms)
                if (knew == k) fails <- fails + 1L ## failed to find valid new point
                k <- knew
                if (k > n) return(outf(xs = xs, ts = ts))
            }
            if (fails == 50) warning("Failed to find acceptable point at step ", k, ". If surrogateAR fails to return a track, this might indicate a location from which it is not possible to step to another, valid location.")
        }
    } else {
        stop("method '", fit$method, "' not coded yet")
    }
    ## Return partial track or NULL
    if (partial) {
        xs[k:n, ] <- NA
        outf(xs = xs, ts = ts)
    } else {
        NULL
    }
}
