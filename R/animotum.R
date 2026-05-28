#' Generate new tracks from a fitted aniMotum model
#'
#' Given a fitted aniMotum crw or rw model, and optionally a template track, this function generates a new track of the same length that coincides with the template at the start point and optionally other specified points along the track. This function works similarly to [aniMotum::sim_fit()] but with different handling of land masking and fixed points.
#'
#' Additional constraints can be placed on the path by rejection sampling through the function `point.check`.  This function must accept a state and return a boolean indicating whether the point is acceptable.  For example, the track can be constrained to the ocean by supplying a `point.check` function that compares the state to a land mask and returns `FALSE` for states corresponding to locations that fall on land.
#'
#' @title aniMotum bridge sampler
#' @param fit a fitted model generated with [aniMotum::fit_ssm()]
#' @param xs matrix: a two-column matrix giving the the template sequence of locations (longitude and latitude). If not provided, the "predicted" (i.e. filtered and interpolated) track from the `fit` object will be used
#' @param ts POSIXct: the times at which the track is sampled
#' @param fixed logical: a vector indicating which locations in the template path are to be held fixed
#' @param point.check logical: a function that accepts a state and returns `TRUE` or `FALSE` indicating whether the state is acceptable
#' @param partial if `TRUE`, a partial track is returned if the sampling fails
#' @return An array of states the define the simulated path
#' @export
surrogateAM <- function(fit, xs, ts, fixed = rep(c(TRUE, FALSE, TRUE), times = c(1, nrow(xs) - 2, 1)), point.check = function(tm, pt) TRUE, partial = FALSE) {

    stopifnot("fit must be an ssm_df object as returned by aniMotum::fit_ssm()" = inherits(x, "ssm_df"))
    if (nrow(fit) > 1) {
        warning("`fit` contains multiple models: only the first will be used")
        fit <- fit[1, ]
    }
    if (missing(xs)) {
        xs <- aniMotum::grab(fit, what = "predicted")
        ts <- xs$date
        xs <- xs[, c("lon", "lat")]
        fixed <- rep(c(TRUE, FALSE, TRUE), times = c(1, nrow(xs) - 2, 1))
    } else if (missing(ts)) {
        ts <- seq_len(nrow(xs))
    }

    if (is.data.frame(xs)) xs <- as.matrix(xs)

    xs <- unname(xs[, 1:2, drop = FALSE])
    n <- nrow(xs)
    if (inherits(ts, "POSIXct")) {
        dt <- as.numeric(difftime(ts, c(as.POSIXct(NA), ts[-length(ts)]), units = "hours"))
        dt[1] <- 0
    } else {
        dt <- c(0, diff(ts))
    }
    K <- 1L ## first model only
    model <- fit$ssm[[K]]$pm ## model type as a string, "crw", "rw", etc

    ## we'll do the simulation in projected xy space
    prj.sim <- "+proj=merc +lon_0=0 +datum=WGS84 +units=km +no_defs"

    ll2xy <- function(lon, lat) {
        ## st_coordinates(st_transform(st_as_sf(data.frame(lon = lon, lat = lat), coords = c("lon","lat"), crs = 4326), crs = prj.sim))
        sf_project(from = "EPSG:4326", to = prj.sim, pts = matrix(c(lon, lat), ncol = 2, byrow = TRUE), authority_compliant = FALSE)
    }
    xy2ll <- function(x, y) {
        ## st_coordinates(st_transform(st_as_sf(data.frame(x = x, y = y), coords = c("x","y"), crs = prj.sim), crs = 4326))
        sf_project(from = prj.sim, to = "EPSG:4326", pts = matrix(c(x, y), ncol = 2, byrow = TRUE), authority_compliant = FALSE)
    }

    ## get parameters from model fit object
    if (model == "crw") {
        Sigma <- diag(2) * 2 * fit$ssm[[K]]$par[c("D_x","D_y"), 1]
        Sigma[1,2] <- Sigma[2,1] <- fit$ssm[[K]]$par["rho_p",1] * sqrt(Sigma[1,1]) * sqrt(Sigma[2,2])
        vmin <- c(min(loc$u, na.rm = TRUE), min(loc$v, na.rm = TRUE))
        vmax <- c(max(loc$u, na.rm = TRUE), max(loc$v, na.rm = TRUE))
        v <- matrix(NA_real_, nrow = n, ncol = 2) ## velocities
        v[1, ] <- c(0, 0) ## initially zero? almost zero?
    } else if (model == "rw") {
        Sigma <- diag(2) * c(fit$ssm[[K]]$par["sigma_x", 1], fit$ssm[[K]]$par["sigma_y", 1]) ^ 2
        Sigma[!Sigma] <- prod(Sigma[1, 1]^0.5, Sigma[2, 2]^0.5) * fit$ssm[[K]]$par["rho_p", 1]
        vmin <- c(min(diff(loc$x), na.rm = TRUE), min(diff(loc$y), na.rm = TRUE))
        vmax <- c(max(diff(loc$x), na.rm = TRUE), max(diff(loc$y), na.rm = TRUE))
    } else {
        stop("unsupported model type, must be 'crw' or 'rw'")
    }

    ## Simulate forward from k0.  Returns the index of the last fixed point reached if an acceptable next candidate cannot be found.
    sample <- function(k0, Sigma, vmin, max) {
        ## Simulate forward from k0
        ll <- xs[k0, ] ## position at k0, lon lat
        k <- k0 + 1L

        ## Find remaining fixed points
        kfixed <- if(k <= n) (k:n)[fixed[k:n]] else integer(0)

        while (k <= n) {
            if (fixed[k]) {
                ## Skip fixed points
                while (k <= n && fixed[k]) {
                    ll <- xs[(k0 <- k), ]
                    k <- k + 1L
                }
                ## Find any remaining fixed points
                kfixed <- if (k <= n) (k:n)[fixed[k:n]] else integer(0)
            } else {
                x <- ll2xy(ll[[1]], ll[[2]]) ## position in projected coords
                ## Try at most 100 new candidate points
                for (r in 1:100) {
                    if (model == "crw") {
                        ## following https://github.com/ianjonsen/aniMotum/blob/cae6bb0c69669fc8c427362d9308facde0bb4ac0/R/sim_fit.R#L217
                        v[k, ] <<- tmvtnorm::rtmvnorm(1, v[k - 1, ], sigma = Sigma * dt[k], lower = vmin, upper = vmax)
                        ## trial new position
                        x1 <- x + v[k, ] * dt[k]
                    } else {
                        ## rw, see https://github.com/ianjonsen/aniMotum/blob/cae6bb0c69669fc8c427362d9308facde0bb4ac0/R/sim_fit.R#L268
                        if (k == 2) {
                            x1 <- rmvnorm(1, x, sigma = Sigma * dt[k]^2)
                        } else {
                            dxy <- rtmvnorm(100, x[k - 1, ] - x[k - 2,], sigma = Sigma * dt[i]^2,
                                            lower = vmin, upper = vmax, algorithm = "gibbs", burn.in.samples = 100)
                            dxy <- dxy[which(!is.na(dxy))[1],]
                            if(all(is.na(dxy))) dxy <- c(1,1)
                            x1 <- x + dxy
                        }
                    }
                    ll1 <- xy2ll(x1[1], x1[2]) ## convert new location to lon/lat
                    if (length(kfixed)) {
                        this_nudge <- xs[kfixed[1], ] - ll1
                        this_nudge[1] <- angle_normalise(this_nudge[1] / 180 * pi) / pi * 180
                        ll1 <- ll1 + (this_nudge) / (kfixed[1] - k + 1L)
                    }
                    ## Test current candidate
                    if (point.check(ts[k], ll1)) {
                        ## Accept candidate
                        ll <- ll1
                        xs[k, ] <<- ll1
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

    k <- 1L
    fails <- 0
    for (i in 1:50) {
        knew <- if (i < 25) sample(k, Sigma, vmin, vmax) else sample(1, Sigma, vmin, vmax)
        if (knew == k) fails <- fails + 1L ## failed to find valid new point
        k <- knew
        if(k > n) return(list(xs = xs, ts = ts))
    }
    if (fails == 50) warning("Failed to find acceptable point at step ", k, ". If surrogateAM fails to return a track, this might indicate a location from which it is not possible to step to another, valid location.")
    ## Return partial track or NULL
    if (partial) {
        xs[k:n, ] <- NA
        list(xs = xs, ts = ts)
    }
}
