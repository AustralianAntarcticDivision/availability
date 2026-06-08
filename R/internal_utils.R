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
