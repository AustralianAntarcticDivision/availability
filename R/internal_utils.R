ecdfq <- function(z) qnorm((rank(z) - 0.5) / length(z))

## For VAR model, construct fit that would have been obtained had the data been rotated by angle theta (radians)
rotateVAR1 <- function(model, theta) {
    if (abs(theta) > 1e-09) {
        nms <- names(model$x.mean)
        R <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
        model$ar[1, , ] <- R %*% model$ar[1, , ] %*% t(R)
        model$var.pred <- R %*% model$var.pred %*% t(R)
        model$x.mean <- model$x.mean %*% t(R)
        names(model$x.mean) <- nms
    }
    model
}

## detect if we've changed from a southerly to northerly (or vice-versa) heading, along a geodesic path
## p0, p1 are long-lat locations
changed_ns <- function(p0, p1) {
    ## calculate the direction of departure and arrival for a geodesic path between p0 and p1
    chk <- geosphere::geodesic_inverse(p0, p1)
    ## if the sign of the northerly component of a1 differs from a2, then we started off heading south and ended up heading north (or vice-versa)
    chk <- cos(chk[2:3] / 180 * pi)
    ## note that for almost-easterly or almost-westerly paths, a great circle will switch from slightly northwards to slightly southwards or vice-versa. Does this matter? Consider increasing the abs(chk) threshold?
    all(abs(chk) > 1e-08) && sign(chk[1]) != sign(chk[2])
}

## helper function to apply a simulation model's transform function (or its inverse), which is used to make the velocity data more normal before fitting the model
apply_transform <- function(x, y, transform, inverse = FALSE, df = TRUE) {
    if (inherits(transform, "transform")) {
        outx <- if (inverse) transform$inverse(x) else transform$transform(x)
        outy <- if (inverse) transform$inverse(y) else transform$transform(y)
    } else if (is.list(transform) && length(transform) == 2 && all(vapply(transform, inherits, "transform", FUN.VALUE = TRUE))) {
        outx <- if (inverse) transform[[1]]$inverse(x) else transform[[1]]$transform(x)
        outy <- if (inverse) transform[[2]]$inverse(y) else transform[[2]]$transform(y)
    } else {
        warning("the `transform` parameter should be a transform object or a list of two such objects, ignoring")
        outx <- x
        outy <- y
    }
    if (df) {
        data.frame(x = outx, y = outy)
    } else {
        stopifnot(length(outx) == 1, length(outy) == 1)
        c(outx, outy)
    }
}

## given projected coordinates, do the projected axes strictly align along easterly/northerly directions?
## `proj` can be an sf object, or a projection string (in this case also need to provide x_range and y_range)
check_proj_is_lonlat_aligned <- function(proj, x_range, y_range) {
    ## construct an xy grid in projected coords
    if (inherits(proj, "sf")) {
        bb <- sf::st_bbox(proj)
        proj <- st_crs(proj)
        xy <- expand.grid(x = seq(bb$xmin, bb$xmax, length.out = 101), y = seq(bb$ymin, bb$ymax, length.out = 100))
    } else {
        xy <- expand.grid(x = seq(x_range[1], x_range[2], length.out = 101), y = seq(y_range[1], y_range[2], length.out = 100))
    }
    tryCatch({
        ## project it to long-lat
        ll <- sf_project(from = proj, to = "EPSG:4326", pts = xy, authority_compliant = FALSE)
        x <- matrix(ll[, 1], ncol = 101, byrow = TRUE)
        y <- matrix(ll[, 2], ncol = 101, byrow = TRUE)
        all(abs(apply(x, 2, diff)) < 1e-08) && ## columns of constant x are projected to the same longitude value
            all(abs(apply(y, 1, diff)) < 1e-08) ## rows of constant y are projected to the same latitude value
    }, error = function(e) FALSE)
}
