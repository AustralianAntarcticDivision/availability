


##' Compute surrogate track by randomizing the steps of an observed
##' track
##'
##' Converts a track to a distance/bearing representation, and then
##' reconstructs a new track by randomly perturbing the bearings of
##' each increment, and otpionally, randomly reordering the
##' increments.
##'
##' @title Randomized track
##' @param lonlat a 2-column matrix or dataframe with longitude and
##' latitude of each point
##' @param rotate a 2-element numeric vector giving the lower and
##' upper limits of the random rotation to apply to the randomized
##' track
##' @param reorder should the track steps be randomly reordered.
##' @return A dataframe with columns
##' \item{\code{lon}}{the longitude of the randomized track}
##' \item{\code{lat}}{the latitude of the randomized track}
##' @author Ben Raymond
##' @export
randomize_track <- function(lonlat,rotate=c(-pi,pi),reorder=FALSE) {
  ## distance and final bearing for each step
  db <- calc_distbearing(lonlat)
  ## Randomize the step order
  if(reorder)
    db <- db[sample.int(nrow(db)),]
  if (!is.null(rotate)) {
    db[,2] <- db[,2]+runif(1,min=min(rotate)/pi*180,max=max(rotate)/pi*180)
  }
  for (k in 2:nrow(lonlat)) {
    lonlat[k,] <- destPoint(lonlat[k-1,],db[k-1,2],db[k-1,1])
  }
  data.frame(lon=lonlat[,1],lat=lonlat[,2])
}


##' Calculate distances in metres and bearing between successive points along track.
##'
##' This is an internal function used by \code{\link{randomize_track}}.
##' @title Distance and Bearing
##' @param lonlat a 2-column matrix or dataframe with longitude and
##' latitude of each point
##' @return A dataframe with columns
##' \item{\code{lon}}{the longitude of the randomized track}
##' \item{\code{lat}}{the latitude of the randomized track}
##' @author Ben Raymond
calc_distbearing <- function(lonlat) {
  nr <- nrow(lonlat)
  dst <- distVincentySphere(lonlat[-nr,],lonlat[-1,])
  ## Why is this the final bearing??
  brg <- finalBearing(lonlat[-nr,],lonlat[-1,])
  data.frame(distance=dst,bearing=brg)
}



##' Fit first-order vector-autoregressive model to track
##'
##' This function fits the vector AR(1) model used to model the track
##' incements in \code{\link{surrogate_arsimulate}}.
##' @title VAR(1) track model
##' @param lonlat a 2-column matrix or dataframe with longitude and
##' latitude of each point
##' @return An object of class "ar"
##' @seealso \code{\link{ar}}, \code{\link{surrogate_arsimulate}}
##' @author Ben Raymond
##' @export
surrogate_arfit <- function(lonlat) {

  ## fixed at 1st order model for now. Might allow this as a param,
  ## but needs surrogate_arsimulate code updated to handle higher
  ## model orders first
  model.order=1

  ## Calculate distance increments dx, dy for each time step
  nr <- nrow(lonlat)
  dx <- distVincentyEllipsoid(lonlat[-nr,],cbind(lonlat[-1,1],lonlat[-nr,2]))*sign(lonlat[-nr,1]-lonlat[-1,1])
  dy <- distVincentyEllipsoid(lonlat[-nr,],cbind(lonlat[-nr,1],lonlat[-1,2]))*sign(lonlat[-nr,2]-lonlat[-1,2])

  ## Fit AR1 to distance increments
  dxdy <- data.frame(dx=dx,dy=dy)
  ar(dxdy,order.max=model.order,aic=FALSE)
}


##' Simulate track from fitted vector autoregressive model
##'
##' Note that land masking uses a built-in land mask image, and it
##' only covers the southern hemisphere. A future version will do
##' something about this.
##'
##' @title Simulated VAR(1) tracks
##' @param arfit fitted object of class "ar" as returned by
##' \code{\link{surrogate_arfit}}
##' @param n number of points to simulate
##' @param startlonlat 2-element vector of starting longitude and latitude
##' @param fixed a dataframe or matrix in which the first column is
##' the index (from 1:n) of each fixed point, and the second and third
##' columns give the associated longitude and latitude
##' @param endlonlat a 2-element vector with ending longitude and
##' latitude. If NULL, no end constraint is imposed except for land
##' masking (if land masking is used). This is a simple way of
##' imposing a return-to-starting-location constraint; for more
##' complex constraints use the \code{fixed} argument
##' @param do.test.land a logical or function. If TRUE, use the
##' included land mask to avoid land. Alternatively, a function can be
##' passed that returns TRUE (point is okay, not on land) or FALSE
##' (point is on land) for a given lon,lat. Note that land masking is
##' ignored for fixed points. Note also that it is possible to create
##' a sitation where tracks are difficult or impossible to simulate,
##' because a fixed point is sufficiently far onto land that the track
##' cannot reach it.
##' @param random.rotation a 2-element vector giving the range of
##' the rotation to apply to the randomized track (values in
##' radians). use \code{random.rotation=NULL} for no such rotation. The
##' angle can be restricted using
##' \code{random.rotation=c(min.angle,max.angle)} - this may speed up
##' computation by avoiding impossible angles (e.g. tracks over a land
##' mass)
##' @param verbose an integer 0-3, if >0 spit out extra information
##' which may be helpful if things don't work as expected. Larger
##' numbers mean more output
##' @param return.all.points if TRUE, return points that were proposed
##' but rejected due to land masking (may be helpful for debugging). If
##' TRUE, the returned data.frame will have an extra column named
##' "valid"
##' @param intermediate.tries when land-masking, try how many times to
##' find a valid point at each step before giving up and starting
##' again? Higher values may improve overall run-time, but too-high
##' values may yield tracks that aren't a good representation of the
##' fitted model
##' @return 2 or 3 column dataframe with the longitude and latitude
##' of simulated track points (and point validity, if
##' return.all.points is TRUE)
##' @author Ben Raymond
##' @importFrom geosphere destPoint
##' @export
surrogate_arsimulate <- function(arfit,n,startlonlat,fixed=NULL,endlonlat=NULL,
                                 do.test.land=TRUE,random.rotation=c(-pi,pi),
                                 verbose=0,return.all.points=FALSE,intermediate.tries=10) {

  if(!is.null(endlonlat) && !is.null(fixed))
    stop("only one of fixed or endlonlat can be supplied")

  ## change to inbuilt land-mask
  if(is.logical(do.test.land) && do.test.land)
    do.test.land <- landmask_init()

  if (!is.null(random.rotation)) {
    this.rotation <- 0
    ## apply rotation to arfit parms
    for (ntries in 1:100) {
      rotate.by <- angle_normalise(runif(1)*diff(range(random.rotation))+min(random.rotation))
      if(verbose>0) cat(sprintf("Rotating track by %.1f degrees\n",rotate.by/pi*180))
      Rm <- matrix(c(cos(rotate.by),-sin(rotate.by),sin(rotate.by),cos(rotate.by)),nrow=2,byrow=TRUE)
      rotated.arfit <- arfit
      rotated.arfit$ar <- Rm%*%matrix(arfit$ar,2,2)%*%t(Rm)
      rotated.arfit$var.pred <- Rm%*%as.matrix(arfit$var.pred)%*%t(Rm)
      rotated.arfit$x.mean <- as.vector(arfit$x.mean)%*%t(Rm)
      ## call simulate on rotated parms
      simtrack <- Recall(arfit=rotated.arfit,n=n,startlonlat=startlonlat,fixed=fixed,endlonlat=endlonlat,
                         do.test.land=do.test.land,random.rotation=NULL,
                         verbose=verbose,return.all.points=return.all.points,intermediate.tries=intermediate.tries)
      if (dim(simtrack)[1]>0) {
        return(simtrack)
      }
    }
    return(simtrack)
  }
  ## convert to "fixed" format
  if(!is.null(endlonlat)) {
    endlonlat <- as.numeric(endlonlat)
    fixed <- data.frame(index=n,lon=endlonlat[1],lat=endlonlat[2])
  }

  ## add starting point as a fixed point
  if(is.null(fixed)) {
    fixed <- data.frame(index=1,lon=startlonlat[1],lat=startlonlat[2])
  } else {
    if(is.matrix(fixed)) {
      fixed <- data.frame(index=fixed[,1],lon=fixed[,2],lat=fixed[,3])
    }
    fixed <- rbind(data.frame(index=1,lon=startlonlat[1],lat=startlonlat[2]),fixed)
  }
  ## ensure ascending order by index
  fixed <- fixed[order(fixed$index),]

  A <- matrix(arfit$ar,2,2,byrow=FALSE)
  fitted.var <- as.matrix(arfit$var.pred)
  fitted.mean <- as.vector(arfit$x.mean)
  tempchol <- chol(fitted.var) ## calculate chol decomposition once
  xsim <- matrix(0,n,2)
  ## burnin for 100 steps
  for (k in 1:100) {
    thisrand <- drop(rnorm(2) %*% tempchol)
    xsim[1,] <- t(A %*% t(xsim[1,]-fitted.mean))+fitted.mean+thisrand ## simulated dx,dy for this time step
  }

  simtrack <- matrix(0,1,3)
  simtrack[1,1:2] <- as.numeric(startlonlat)
  simtrack[1,3] <- 1 ## valid
  sidx <- 1 ## pointer into simtrack matrix of last valid point
  for (k in 2:n) {
    ## k is index into xsim, the model of x- and y- step lengths/speeds
    point.okay <- TRUE
    if (verbose>0) cat(sprintf("Step %d, current location is %.3f, %.3f\n",k,simtrack[sidx,1],simtrack[sidx,2]))
    for (land.tries in 1:intermediate.tries) {
      thisrand <- matrix(rnorm(2) %*% tempchol,nrow=1)
      xsim[k,]  <- t(A %*% t(xsim[k-1,]-fitted.mean))+fitted.mean+thisrand ## simulated dx,dy for this time step
      ## calculate track point from steps
      tempx <- destPoint(simtrack[sidx,1:2],90,xsim[k,1]) # x step
      tempx[1] <- angle_normalise(tempx[1]/180*pi)/pi*180 ## ensure longitude is in range -180 to 180
      tempy <- destPoint(simtrack[sidx,1:2],0,xsim[k,2]) # y step
      ## note that destPoint will handle the case where a step crosses 90S or 90N
      simtrack <- rbind(simtrack,c(tempx[1],tempy[2],NA)) ## 3rd entry (valid) is NA for now
      ## as we get closer to the next fixed point, increasingly nudge the random point towards the designated fixed location
      next_fixed <- if(any(fixed$index>=k)) which.max(fixed$index>=k) else NA
      if(verbose>1) cat(sprintf("  proposed point %d is at %.3f, %.3f\n",k,simtrack[nrow(simtrack),1],simtrack[nrow(simtrack),2]))
      if(!is.na(next_fixed)) {
        if(verbose>1) cat(sprintf("    the next fixed point is at %.3f, %.3f in %d steps time\n",fixed$lon[next_fixed],fixed$lat[next_fixed],fixed$index[next_fixed]-k))
        ## next_fixed is row index into fixed
        a <- diag(1/(fixed$index[next_fixed]-k+1),2)
        simtrack[nrow(simtrack),1:2] <- simtrack[nrow(simtrack),1:2]+a%*%(c(fixed$lon[next_fixed],fixed$lat[next_fixed])-simtrack[nrow(simtrack),1:2])
        if (verbose>1) cat(sprintf("    the proposed point has been nudged to %.3f, %.3f because of the next fixed point\n",simtrack[nrow(simtrack),1],simtrack[nrow(simtrack),2]))
      }
      if(is.function(do.test.land)) {
        if(is.na(next_fixed) || fixed$index[next_fixed]!=k) {
          point.okay=do.test.land(simtrack[nrow(simtrack),])
        } else {
          if(verbose>1) cat("    not checking land-mask for this point, because it is a fixed point.\n")
        }
      }
      if(point.okay) {
        if (verbose>1 & is.function(do.test.land)) cat(sprintf("    this proposed point does not lie on land, accepting\n"))
        simtrack[nrow(simtrack),3] <- 1
        sidx=nrow(simtrack) ## update pointer to valid location
        break
      } else {
        if (verbose>1) cat(sprintf("    this proposed point lies on land\n"))
        simtrack[nrow(simtrack),3] <- 0
      }
    }
    if(!point.okay) {
      ## could not find a valid point at this step: give up
      if (verbose>0) cat(sprintf("  could not find valid point, abandoning this track and starting again\n"))
      simtrack <- NULL
      break
    }
  }
  if(return.all.points) {
    data.frame(lon=simtrack[,1],lat=simtrack[,2],valid=as.logical(simtrack[,3]))
  } else {
    data.frame(lon=simtrack[simtrack[,3]>0,1],lat=simtrack[simtrack[,3]>0,2])
  }
}




##' Fold angles (radians) into [-pi,pi)
##'
##' This is an internal function used by
##' \code{\link{surrogate_arsimulate}}.
##' @title Fold angles into [-pi,pi)
##' @param x a vector of angles
##' @return a vector of folded angles
##' @author Ben Raymond
angle_normalise <- function(x) {
  ## normalize angle in radians to range [-pi,pi)
  (x+pi)%%(2*pi)-pi
}


