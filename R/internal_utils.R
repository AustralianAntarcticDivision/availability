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
