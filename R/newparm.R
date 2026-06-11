#' Fit first-order vector-autoregressive model to track
#'
#' @param x : a data.frame containing the observed track. It can optionally be an `sf`-tibble (i.e. in projected coordinates). If not projected, it must contain columns named "lon" and "lat". In both cases it must contain a "date" column with timestamps in POSIX format
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
        ## TODO deal with dateline wrapping, pole crossing
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
                structure(list(transform = approxfun(sx, ecdfq(sx), rule = 2), inverse = approxfun(ecdfq(sx), sx, rule = 2)), class = "transform"),
                structure(list(transform = approxfun(sy, ecdfq(sy), rule = 2), inverse = approxfun(ecdfq(sy), sy, rule = 2)), class = "transform")
            )
        }
        dxdy <- apply_transform(sx, sy, transform = dots$transform)
    } else {
        dxdy <- data.frame(x = sx, y = sy)
    }
    ## Fit VAR1 to stepwise speeds
    structure(list(
        method = "var1",
        model = ar(dxdy, order.max = 1L, aic = FALSE),
        projection = proj,
        transform = dots$transform), ## will be NULL if not set
        data = x,
        class = "av_fit")
}


#' Simulate a track
#'
#' @param fit av_fit: an object of class `av_fit`, as returned by e.g. [av_fit_ar()]
#' @param x : a data.frame containing a "template" track. It can optionally be an `sf`-tibble (i.e. in projected coordinates). If not projected, it must contain columns named "lon" and "lat". In all cases it must contain a "date" column with timestamps in POSIX format. This template is used to specify the projection and length of the simulated track, the location of fixed points (see `fixed`), and the time difference between successive points on the track
#' @param fixed logical: a vector (with length equal to the number of rows in `x`) indicating which locations in the template path are to be held fixed. By default the first and last points are fixed
#' @param point_check function: a function that accepts a parameters `tm` (date) and `pt` (location) and returns `TRUE` or `FALSE` indicating whether the location is acceptable. By default all points are accepted. The `pt` parameter passed to the `point_check` function will be in the same projection as `x` (or will be longitude and latitude if `x` is not projected)
#' @param partial logical: if `TRUE`, a partial track is returned if the sampling fails
#' @param random_rotation numeric: for the VAR-based method only. A two-element vector giving the upper and lower limits (radians) of the rotation applied to the model
#' @param force logical: if `TRUE`, force the simulation to go ahead even if internal checks fail
#'
#' @return If `x` was an `sf`-tibble in projected coordinates, the returned value will also be an `sf`-tibble in the same projection with column `date`. Otherwise for unprojected (longitude-latitude) data, the returned value will be a data.frame with columns `lon`, `lat`, and `date`
#'
#' @export
av_sim <- function(fit, x, fixed, point_check, partial = FALSE, random_rotation, force = FALSE) {
    stopifnot("`fit` should be an object of class \"av_fit\" as returned by e.g. `av_fit_ar()`, or an `ssm_df` object as returned by `aniMotum::fit_ssm()`" = inherits(fit, "av_fit") || inherits(fit, "ssm_df"))
    if (missing(point_check)) point_check <- function(...) TRUE

    ## fit can be:
    ## - an ssm_df object, in which case it has been fitted to projected data, or
    ## - an av_fit object e.g. VAR model, which could have been fitted to long-lat or projected coords

    ## if our fit was made on projected coordinates, test that those coordinates are E/N aligned
    fit_alignment_check <- NA
    fit_proj <- if (inherits(fit, "ssm_df")) st_crs(aniMotum::grab(fit, what = "predicted", as_sf = TRUE)) else fit$projection
    if (!is.null(fit_proj)) {
        temp <- if (inherits(fit, "ssm_df")) aniMotum::grab(fit, what = "predicted", as_sf = TRUE) else fit$data
        fit_alignment_check <- st_is_longlat(temp) || check_proj_is_lonlat_aligned(temp)
    }

    ## if we are simulating on projected coordinates, test that those coordinates are E/N aligned
    sim_alignment_check <- NA
    if (inherits(x, "sf")) {
        sim_alignment_check <- st_is_longlat(x) || check_proj_is_lonlat_aligned(x)
    }

    ## check simulation on long-lat from a projected fit
    if (!is.null(fit_proj) && (!inherits(x, "sf") || st_is_longlat(x))) {
        if (!isTRUE(force) && !fit_alignment_check) stop("the simulation is to be done in longitude-latitude coordinates, but the model was fitted to projected coordinates that are not cardinal and rectilinear. If you are confident that this is OK, you can force the simulation to run with `av_sim(..., force = TRUE)`")
    }

    ## check simulation on projected coords from a projected fit (do the projections match?)
    if (inherits(x, "sf") && !st_is_longlat(x) && !is.null(fit_proj) && !isTRUE(st_crs(x) == fit_proj)) {
        if (!isTRUE(force)) {
            stop("the projections of `fit` and `x` are different. If you are confident that this is OK, you can force the simulation to run with `av_sim(..., force = TRUE)`")
        }
    }

    ## check simulation on projected coordinates from a fit made on long-lat data
    if (is.null(fit_proj) && inherits(x, "sf") && !st_is_longlat(x)) {
        if (!isTRUE(force) && !sim_alignment_check) stop("the model was fitted to longitude-latitude coordinates, but the simulation is to be done on projected coordinates that are not cardinal and rectilinear. If you are confident that this is OK, you can force the simulation to run with `av_sim(..., force = TRUE)`")
    }

    if (inherits(fit, "ssm_df")) return(av_sim_am(fit, x = x, fixed = fixed, point_check = point_check, partial = partial, alignment_checks = list(fit = fit_alignment_check, sim = sim_alignment_check)))
    if (missing(fixed)) fixed <- rep(c(TRUE, FALSE, TRUE), c(1, nrow(x) - 2, 1))
    if (fit$method == "var1") {
        av_sim_var1(fit, x = x, fixed = fixed, point_check = point_check, partial = partial, random_rotation = random_rotation, alignment_checks = list(fit = fit_alignment_check, sim = sim_alignment_check))
    } else {
        stop("`fit` uses an unrecognized method: ", fit$method)
    }
}

av_sim_var1 <- function(fit, x, fixed, point_check, partial, random_rotation, alignment_checks) {
    if (missing(random_rotation)) random_rotation <- c(-pi, pi)
    ## ## lon-lat. Temporarily(?) call the old surrogateAR code
    ## mod <- fit$model
    ## attr(mod, "speed") <- TRUE ## fit$model is a speed-based one, this is just a temporary workaround to tell the sim code to scale by timestep to get from speed to distance
    ## attr(mod, "transform") <- fit$transform
    ## out <- surrogateAR(mod, xs = x[, c("lon", "lat")], ts = x$date, fixed = fixed, point.check = point_check, random.rotation = random_rotation, partial = partial)
    ## data.frame(lon = out$xs[, 1], lat = out$xs[, 2], date = x$date)
    av_sim_core(fit, x = x, fixed = fixed, point_check = point_check, random_rotation = random_rotation, partial = partial, alignment_checks = alignment_checks)
}

av_sim_am <- function(fit, x, fixed, point_check, partial, random_rotation, alignment_checks) {
    if (nrow(fit) > 1) {
        warning("`fit` contains multiple models: only the first will be used")
        fit <- fit[1, ]
    }
    if (missing(x)) x <- aniMotum::grab(fit, what = "predicted")
    if (missing(fixed)) fixed <- rep(c(TRUE, FALSE, TRUE), c(1, nrow(x) - 2, 1))
    av_sim_core(list(method = "aniMotum", model = fit, projection = st_crs(aniMotum::grab(fit, what = "predicted", as_sf = TRUE)), transform = NULL), x = x, fixed = fixed, point_check = point_check, random_rotation = random_rotation, partial = partial, alignment_checks = alignment_checks)
}

## @param fit
## @param ts the times at which the track is sampled
## @param xs the template sequence of states
## @param fixed a logical vector indicating which locations in the template path are to be held fixed.
## @param point.check function that accepts a state and returns boolean indicating whether the state is acceptable.
## @param random.rotation the upper and lower limits (radians) of the rotation applied to the VAR(1) model.
## @param partial if `TRUE`, a partial track is returned if the sampling fails.
## @return An array of states the define the simulated path.
av_sim_core <- function(fit, x, fixed, point_check, random_rotation, partial, alignment_checks) {
    proj <- NULL ## leave as NULL for simulation on long-lat coords
    fit_unit_scale <- sim_unit_scale <- 1 ## scalar to apply to units of m in fit coords and sim coords. For long-lat the scalar is 1 because calculations are done in m. If projected coords are in units of km then then the unit_scale value will be 1000
    if (!is.null(fit$projection)) {
        crd_units <- convert_to_base(fit$projection$ud_unit)
        if (isTRUE(deparse_unit(crd_units) == "m")) {
            fit_unit_scale <- as.numeric(crd_units)
        } else {
            warning("the coordinate system of `fit` has an unexpected base unit ('", deparse_unit(crd_units), "'), check the scale of the simulated track")
        }
    }
    if (inherits(x, "sf")) {
        xs <- st_coordinates(x)
        if (!st_is_longlat(x)) proj <- st_crs(x)
        crd_units <- convert_to_base(st_crs(x, parameters = TRUE)$ud_unit)
        if (isTRUE(deparse_unit(crd_units) == "m")) {
            sim_unit_scale <- as.numeric(crd_units)
        } else {
            warning("the coordinate system of `x` has an unexpected base unit ('", deparse_unit(crd_units), "'), check the scale of the simulated track")
        }
    } else if (is.data.frame(x)) {
        xs <- x[, c("lon", "lat")]
    } else {
        stop("`x` is an unexpected format")
    }
    ts <- x$date
    n <- nrow(xs)

    if (is.data.frame(xs)) xs <- as.matrix(xs)
    xs <- unname(xs[, 1:2, drop = FALSE])

    if (!is.null(proj)) {
        ## helper functions to transform from projected coordinates (i.e. the projection we are simulating on, which might be different to the projection that was used to fit the model) to long-lat and vice-versa
        ## ll2xy <- function(lonlat) {
        ##     sf_project(from = "EPSG:4326", to = proj, pts = matrix(lonlat, ncol = 2, byrow = TRUE), authority_compliant = FALSE)
        ## }
        xy2ll <- function(xy) {
            sf_project(from = proj, to = "EPSG:4326", pts = matrix(xy, ncol = 2, byrow = TRUE), authority_compliant = FALSE)
        }
    }

    dt <- if (inherits(ts, "POSIXct")) c(0, as.numeric(diff(ts, units = "secs"))) else c(0, diff(ts))

    outf <- function(xs, ts) {
        if (is.null(proj)) {
            data.frame(lon = xs[, 1], lat = xs[, 2], date = ts)
        } else {
            ## return as an sf-tibble in projected coords
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
        pos <- xs[k0, ] ## starting position, either lon-lat or x-y
        k <- k0 + 1L

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
                    step <- rep(NA_real_, 2)
                    if (fit$method == "var1") {
                        ## Trial z - simulate from centred VAR(1) model
                        z1 <- drop(parms$A %*% z) + drop(rnorm(2) %*% parms$U)
                        thisv <- z1 + parms$mu
                    } else if (fit$method == "aniMotum") {
                        if (parms$am_model == "crw") {
                            ## following https://github.com/ianjonsen/aniMotum/blob/cae6bb0c69669fc8c427362d9308facde0bb4ac0/R/sim_fit.R#L217
                            ## note: vmin, vmax, and sigma are km/h. Our v is also km/h, but dt is secs not hours
                            thisv <- rtmvnorm(1, v[k - 1, ], sigma = parms$Sigma * dt[k] / 3600, lower = parms$vmin, upper = parms$vmax)
                        } else {
                            stop("aniMotum rw models are not supported yet")
                            ## rw, see https://github.com/ianjonsen/aniMotum/blob/cae6bb0c69669fc8c427362d9308facde0bb4ac0/R/sim_fit.R#L268
                        }
                    } else {
                        stop("method not coded")
                    }
                    if (!is.null(fit$transform)) {
                        ## v needs to be transformed (and saved for next iter), because the model was fitted on transformed v's
                        thisv <- apply_transform(thisv[1], thisv[2], transform = fit$transform, inverse = TRUE, df = FALSE)
                    }
                    v[k, ] <<- thisv
                    step <- if (fit$method == "var1") thisv * fit_unit_scale else if (fit$method == "aniMotum") thisv / 3.6 else stop("method not coded")
                    ## for VAR, scale the units if needed (will only matter if the VAR model was fitted to projected data not in units of m); for aniMotum, convert km/h to m/s
                    step <- step * dt[k] ## from m/s to m
                    ## `step` is the x, y step we are taking, in projected units or m for long-lat
                    if (is.null(proj)) {
                        ## take step in lon then lat
                        ## if the fit was projected, we've already tested that it was E/N aligned, and the step values are in m, so this should be OK
                        x1 <- as.vector(destPoint(destPoint(pos, 90, step[1]), 0, step[2]))
                        ## destPoint already handles pole-crossing but we need to adjust v if we did so, done below
                    } else {
                        ## in projected coords
                        ## adjust step length for simulation units: our step values are in m, but the projected units might be e.g. km
                        x1 <- as.vector(pos + step / sim_unit_scale)
                        ## TODO handle pole crossing, dateline crossing, out-of-bounds
                    }
                    ## if we have a fixed point coming up, nudge towards it
                    ## TODO handle dateline if we are simulating in projected coords
                    if (length(kfixed)) {
                        this_nudge <- xs[kfixed[1], ] - x1
                        if (is.null(fit$projection)) this_nudge[1] <- angle_normalise(this_nudge[1] / 180 * pi) / pi * 180
                        x1 <- x1 + (this_nudge) / (kfixed[1] - k + 1L)
                        if (is.null(fit$projection)) x1[1] <- angle_normalise(x1[1] / 180 * pi) / pi * 180
                    }
                    ## test current candidate
                    if (point_check(ts[k], x1)) {
                        ## Accept candidate
                        ## did we change from southerly to northerly heading (e.g. crossed the pole)? If we did, and if the velocity is parameterized as easterly/northerly components, then we need to reverse the sign of the northerly component for the next simulation step
                        ## we can only do this if the model was fitted to long-lat (unprojected) data, or its projection had easterly/northerly alignment
                        if (is.null(fit$projection) || isTRUE(alignment_checks$fit)) {
                            went_from <- if (is.null(proj)) pos else xy2ll(pos)
                            went_to <- if (is.null(proj)) x1 else xy2ll(x1)
                            ## those are the start and end locations of this step, either the longlat point, or transformed from xy to ll if we are simulating on projected coords
                            if (isTRUE(changed_ns(went_from, went_to))) {
                                ## cat("changed N/S from:", went_from, "to", went_to, "with v:", if (fit$method == "var1") z1 else v[k, ], "\n")
                                if (fit$method == "var1") {
                                    z1 <- (z1 + parms$mu) * c(1, -1) - parms$mu
                                } else if (fit$method == "aniMotum") {
                                    if (parms$am_model == "crw") {
                                        thisv <- v[k, ]
                                        thisv[2] <- -thisv[2]
                                        v[k, ] <<- thisv
                                    } else {
                                        stop("method not coded")
                                    }
                                }
                            }
                        }
                        if (fit$method == "var1") z <- z1
                        pos <- x1 ## new point
                        xs[k, ] <<- x1
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
##   list(method = "var1" or "aniMotum"
##        model = <AR model or ssm_df>,
##        projection = st_crs OR NULL,
##        data (only for method "var1") = original data <data.frame or sf-tibble> used to fit the model,
##        transform = <transform object or NULL>)

    for (roti in seq_len(if (fit$method == "var1") 100 else 1)) { ## Try 100 rotations of the model for VAR method
        if (fit$method == "var1") {
            theta <- if (is.null(random_rotation)) 0 else runif(1, random_rotation[1], random_rotation[2])
            model0 <- rotateVAR1(fit$model, theta)
            parms <- list(A = unname(model0$ar[1, , ]), U = chol(unname(model0$var.pred)), mu = as.vector(model0$x.mean))
            v <- matrix(NA_real_, nrow = n, ncol = 2) ## velocities, km/h
        } else if (fit$method == "aniMotum") {
            parms <- list(am_model = fit$model$ssm[[1]]$pm) ## aniMotum model type as a string, "crw", "rw", etc
            if (parms$am_model == "crw") {
                Sigma <- diag(2 * fit$model$ssm[[1]]$par[c("D_x", "D_y"), 1])
                Sigma[1, 2] <- Sigma[2, 1] <- fit$model$ssm[[1]]$par["rho_p", 1] * sqrt(Sigma[1, 1]) * sqrt(Sigma[2, 2])
                parms$Sigma <- Sigma
                uv <- aniMotum::grab(fit$model, what = "predicted")[, c("u", "v")]
                parms$vmin <- c(min(uv$u, na.rm = TRUE), min(uv$v, na.rm = TRUE)) ## km/h
                parms$vmax <- c(max(uv$u, na.rm = TRUE), max(uv$v, na.rm = TRUE)) ## km/h
                v <- matrix(NA_real_, nrow = n, ncol = 2) ## velocities, km/h
                v[1, ] <- as.numeric(uv[1, ])
                if (any(is.na(v[1, ]))) v[1, ] <- c(0, 0) ## fallback
            } else if (parms$am_model == "rw") {
                Sigma <- diag(2) * c(fit$model$ssm[[1]]$par["sigma_x", 1], fit$model$ssm[[1]]$par["sigma_y", 1]) ^ 2
                Sigma[!Sigma] <- prod(Sigma[1, 1]^0.5, Sigma[2, 2]^0.5) * fit$model$ssm[[1]]$par["rho_p", 1]
                parms$Sigma <- Sigma
                xy <- aniMotum::grab(fit$model, what = "predicted")[, c("x", "y")]
                parms$vmin <- c(min(diff(xy$x), na.rm = TRUE), min(diff(xy$y), na.rm = TRUE))
                parms$vmax <- c(max(diff(xy$x), na.rm = TRUE), max(diff(xy$y), na.rm = TRUE))
            } else {
                stop("unsupported aniMotum model type, must be 'crw' or 'rw'")
            }
        } else {
            stop("method '", fit$method, "' not coded yet")
        }
        k <- 1L
        fails <- 0
        for (i in 1:50) {
            knew <- if (i < 25) sample(k, parms = parms) else sample(1, parms = parms)
            if (knew == k) fails <- fails + 1L ## failed to find valid new point
            k <- knew
            if (k > n) return(outf(xs = xs, ts = ts))
        }
    }
    if (fails == 50) warning("Failed to find acceptable point at step ", k, ". If av_sim fails to return a track, this might indicate a location from which it is not possible to step to another, valid location.")
    ## Return partial track or NULL
    if (partial) {
        xs[k:n, ] <- NA
        outf(xs = xs, ts = ts)
    } else {
        NULL
    }
}
