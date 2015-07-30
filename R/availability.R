#' Compute surrogate track by randomizing the steps of an observed track
#'
#' Assumes that steps are equally spaced in time. Probably not a particular good way of doing it since it destroys the autocorrelation
#' structure of the track.
#'
#' @author Ben Raymond
#'
#' @param lonlat array: 2-column matrix or data.frame with longitude, latitude of each point
#' @param rotate array: 2-element numeric giving the lower and upper limits of the random rotation to apply to the randomized track
#'
#' @return 2-column data.frame with longitude,latitude of randomized track points
#'
#' @export randomize_track

randomize_track=function(lonlat,rotate=c(-pi,pi)) {
    db=calc_distbearing(lonlat) ## distance and bearing at each step
    db=db[sample.int(nrow(db)),] ## randomize the step order
    if (!is.null(rotate)) {
        db[,2]=db[,2]+runif(1,min=min(rotate)/pi*180,max=max(rotate)/pi*180)
    }
    temp=lonlat
    for (k in 2:nrow(temp)) {
        temp[k,]=destPoint(temp[k-1,],db[k-1,2],db[k-1,1])
    }
    temp
}


#' Fit first-order vector-autoregressive model to track
#'
#' @author Ben Raymond
#'
#' @param lonlat array: 2-column matrix or data.frame with longitude, latitude of each point
#'
#' @return object of class "ar"
#'
#' @seealso \code{\link{ar}} \code{\link{surrogate_arsimulate}}
#'
#' @export surrogate_arfit

# @param model.order numeric: the order of the AR model to fit, default=1

surrogate_arfit=function(lonlat) {
    ## calculate dx, dy separately at each time step
    model.order=1 ## fixed at 1st order model for now. Might allow this as a param, but needs surrogate_arsimulate code updated to handle higher model orders first
    nr=nrow(lonlat)
    dx=distVincentyEllipsoid(lonlat[1:(nr-1),],cbind(lonlat[-1,1],lonlat[1:(nr-1),2]))*sign(lonlat[1:(nr-1),1]-lonlat[-1,1])
    dy=distVincentyEllipsoid(lonlat[1:(nr-1),],cbind(lonlat[1:(nr-1),1],lonlat[-1,2]))*sign(lonlat[1:(nr-1),2]-lonlat[-1,2])

    dxdy=data.frame(dx=dx,dy=dy)
    ar(dxdy,order.max=model.order,aic=FALSE)
}


#' Simulate track from fitted vector autoregressive model
#'
#' Note that land masking uses a built-in land mask image, and it only covers the southern hemisphere. A future version will do something about this.
#'
#' @author Ben Raymond
#'
#' @param arfit : fitted object of class "ar" as returned by \code{\link{surrogate_arfit}}
#' @param n numeric: number of points to simulate
#' @param startlonlat numeric: 2-element array with starting longitude and latitude
#' @param fixed data.frame or matrix: first column is the index (from 1:n) of each fixed point, and the second and third columns give the associated longitude and latitude
#' @param endlonlat numeric: 2-element array with ending longitude and latitude. If NULL, no end constraint is imposed except for land masking (if land masking is used). This is a simple way of imposing a return-to-starting-location constraint; for more complex constraints use the \code{fixed} argument
#' @param do.test.land logical or function: if TRUE, use the included land mask to avoid land. Alternatively, a function can be passed that returns TRUE (point is okay, not on land) or FALSE (point is on land) for a given lon,lat. Note that land masking is ignored for fixed points. Note also that it is possible to create a sitation where tracks are difficult or impossible to simulate, because a fixed point is sufficiently far onto land that the track cannot reach it.
#' @param random.rotation numeric: 2-element array giving the range of the rotation to apply to the randomized track (values in radians). use \code{random.rotation=NULL} for no such rotation. The angle can be restricted using \code{random.rotation=c(min.angle,max.angle)} - this may speed up computation by avoiding impossible angles (e.g. tracks over a land mass)
#' @param verbose numeric: 0-3, if >0 spit out extra information which may be helpful if things don't work as expected. Larger numbers mean more output
#' @param return.all.points logical: if TRUE, return points that were proposed but rejected due to land masking (may be helpful for debugging). If TRUE, the returned data.frame will have an extra column named "valid"
#' @param intermediate.tries numeric: when land-masking, try how many times to find a valid point at each step before giving up and starting again? Higher values may improve overall run-time, but too-high values may yield tracks that aren't a good representation of the fitted model
#'
#' @return 2- or 3-column data.frame with longitude,latitude of simulated track points (and point validity, if return.all.points is TRUE)
#'
#' @seealso \code{\link{surrogate_arfit}}
#'
#' @export surrogate_arsimulate

surrogate_arsimulate=function(arfit,n,startlonlat,fixed=NULL,endlonlat=NULL,do.test.land=TRUE,random.rotation=c(-pi,pi),verbose=0,return.all.points=FALSE,intermediate.tries=10) {
    if (!is.null(endlonlat) && !is.null(fixed)) {
        stop("only one of fixed or endlonlat can be supplied")
    }
    if (is.logical(do.test.land) && do.test.land) {
        ## change to function for inbuilt land-mask
        do.test.land=landmask_init()
    }
    if (! is.null(random.rotation)) {
        this.rotation=0
        ## apply rotation to arfit parms
        for (ntries in 1:100) {
            rotate.by=angle.normalise(runif(1)*diff(range(random.rotation))+min(random.rotation))
            if (verbose>0) cat(sprintf("Rotating track by %.1f degrees\n",rotate.by/pi*180))
            Rm=matrix(c(cos(rotate.by),-sin(rotate.by),sin(rotate.by),cos(rotate.by)),nrow=2,byrow=TRUE)
            rotated.arfit=arfit
            rotated.arfit$ar=Rm %*% matrix(arfit$ar,nrow=2) %*% t(Rm)
            rotated.arfit$var.pred=Rm %*% as.matrix(arfit$var.pred) %*% t(Rm)
#            print(matrix(arfit$x.mean,nrow=1))
            rotated.arfit$x.mean=matrix(arfit$x.mean,nrow=1) %*% t(Rm)
            ## call simulate on rotated parms
            simtrack=surrogate_arsimulate(arfit=rotated.arfit,n=n,startlonlat=startlonlat,fixed=fixed,endlonlat=endlonlat,do.test.land=do.test.land,random.rotation=NULL,verbose=verbose,return.all.points=return.all.points,intermediate.tries=intermediate.tries)
            if (dim(simtrack)[1]>0) {
                break
            }
        }
        return(simtrack)
    }
    if (!is.null(endlonlat)) {
        endlonlat=as.numeric(endlonlat)
        fixed=data.frame(index=n,lon=endlonlat[1],lat=endlonlat[2]) ## convert to "fixed" format
    }
    ## add starting point as a fixed point
    if (is.null(fixed)) {
        fixed=data.frame(index=1,lon=startlonlat[1],lat=startlonlat[2])
    } else {
        if (is.matrix(fixed)) {
            fixed=data.frame(index=fixed[,1],lon=fixed[,2],lat=fixed[,3])
        }
        fixed=rbind(data.frame(index=1,lon=startlonlat[1],lat=startlonlat[2]),fixed)
    }
    fixed=fixed[order(fixed$index),] ## ensure ascending order by index
    #if (do.test.land) {
    #    land.mask=readPNG(system.file("extdata","land_mask-0.1-nosub.png",package="availability")) ## 0=land, 1=ocean
    #    land.lon=seq(from=-180,to=180,length.out=dim(land.mask)[2])
    #    land.lat=seq(from=0,to=-90,length.out=dim(land.mask)[1])
    #}
    A=matrix(arfit$ar,ncol=2,byrow=FALSE)
    fitted.var=as.matrix(arfit$var.pred)
    fitted.mean=matrix(arfit$x.mean,nrow=1)
    tempchol=chol(fitted.var) ## calculate chol decomposition once
    xsim=matrix(0,n,2)
    for (k in 1:100) { ## burnin period
        thisrand=matrix(rnorm(2) %*% tempchol,nrow=1)
#        print(thisrand)
#        print(A)
#        print(xsim[,1]-fitted.mean)
        xsim[1,]=t(A %*% t(xsim[1,]-fitted.mean))+fitted.mean+thisrand ## simulated dx,dy for this time step
#        stop()
    }

    simtrack=matrix(0,1,3)
    simtrack[1,1:2]=as.numeric(startlonlat)
    simtrack[1,3]=TRUE ## valid
    sidx=1 ## pointer into simtrack matrix of last valid point
    for (k in 2:n) {
        ## k is index into xsim, the model of x- and y- step lengths/speeds
        point.okay=TRUE
        if (verbose>0) cat(sprintf("Step %d, current location is %.3f, %.3f\n",k,simtrack[sidx,1],simtrack[sidx,2]))
        for (land.tries in 1:intermediate.tries) {
            thisrand=matrix(rnorm(2) %*% tempchol,nrow=1)
            xsim[k,]=t(A %*% t(xsim[k-1,]-fitted.mean))+fitted.mean+thisrand ## simulated dx,dy for this time step
            ## calculate track point from steps
            tempx=destPoint(simtrack[sidx,1:2],90,xsim[k,1]) # x step
            tempx[1]=angle.normalise(tempx[1]/180*pi)/pi*180 ## ensure longitude is in range -180 to 180
            tempy=destPoint(simtrack[sidx,1:2],0,xsim[k,2]) # y step
            ## note that destPoint will handle the case where a step crosses 90S or 90N
            simtrack=rbind(simtrack,c(tempx[1],tempy[2],NA)) ## 3rd entry (valid) is NA for now
            ## as we get closer to the next fixed point, increasingly nudge the random point towards the designated fixed location
            next_fixed=if (any(fixed$index>=k)) which.max(fixed$index>=k) else NA
            ##    cat(sprintf("%d: next_fixed=%d (%d) [diff=%.3f]\n",k,next_fixed,fixed$index[next_fixed],1/(fixed$index[next_fixed]-k+1)))
            if (verbose>1) cat(sprintf("  proposed point %d is at %.3f, %.3f\n",k,simtrack[nrow(simtrack),1],simtrack[nrow(simtrack),2]))
            if (!is.na(next_fixed)) {
                if (verbose>1) cat(sprintf("    the next fixed point is at %.3f, %.3f in %d steps time\n",fixed$lon[next_fixed],fixed$lat[next_fixed],fixed$index[next_fixed]-k))
                ## next_fixed is row index into fixed
                a=diag(1/(fixed$index[next_fixed]-k+1),2)
                simtrack[nrow(simtrack),1:2]=simtrack[nrow(simtrack),1:2]+a%*%(c(fixed$lon[next_fixed],fixed$lat[next_fixed])-simtrack[nrow(simtrack),1:2])
                if (verbose>1) cat(sprintf("    the proposed point has been nudged to %.3f, %.3f because of the next fixed point\n",simtrack[nrow(simtrack),1],simtrack[nrow(simtrack),2]))
            }
            if (is.function(do.test.land)) {
                if (is.na(next_fixed) || fixed$index[next_fixed]!=k) {
                    ## test if point over land, but not if this is a fixed point
                    #lonidx=which.min(abs(land.lon-simtrack[nrow(simtrack),1]))
                    #latidx=which.min(abs(land.lat-simtrack[nrow(simtrack),2]))
                                        #point.okay=land.mask[latidx,lonidx]==1
                    point.okay=do.test.land(simtrack[nrow(simtrack),])
                } else {
                    if (verbose>1) cat("    not checking land-mask for this point, because it is a fixed point.\n")
                }
            }
            if (point.okay) {
                if (verbose>1 & is.function(do.test.land)) cat(sprintf("    this proposed point does not lie on land, accepting\n"))
                simtrack[nrow(simtrack),3]=TRUE
                sidx=nrow(simtrack) ## update pointer to valid location
                break
            } else {
                if (verbose>1) cat(sprintf("    this proposed point lies on land\n"))
                simtrack[nrow(simtrack),3]=FALSE
            }
        }
        if (! point.okay) {
            ## could not find a valid point at this step: give up
            if (verbose>0) cat(sprintf("  could not find valid point, abandoning this track and starting again\n"))
            simtrack=NULL
            break
        }
    }
    if (return.all.points) {
        data.frame(lon=simtrack[,1],lat=simtrack[,2],valid=simtrack[,3])
    } else {
        data.frame(lon=simtrack[simtrack[,3],1],lat=simtrack[simtrack[,3],2])
    }
}



## internal helper functions

# Calculate x- and y-step distances in metres between successive points
#
# @author Ben Raymond
#
# @param lonlat array: 2-column matrix or data.frame with longitude, latitude of each point
#
# @return data.frame with dx and dy in metres

calc_dxdy=function(lonlat) {
    ## calculate dx, dy separately at each time step
    nr=nrow(lonlat)
    dx=distVincentyEllipsoid(lonlat[1:(nr-1),],cbind(lonlat[-1,1],lonlat[1:(nr-1),2]))*sign(lonlat[1:(nr-1),1]-lonlat[-1,1])
    dy=distVincentyEllipsoid(lonlat[1:(nr-1),],cbind(lonlat[1:(nr-1),1],lonlat[-1,2]))*sign(lonlat[1:(nr-1),2]-lonlat[-1,2])
    data.frame(dx=dx,dy=dy)
}

# Calculate distances in metres and bearing between successive points
#
# @author Ben Raymond
#
# @param lonlat array: 2-column matrix or data.frame with longitude, latitude of each point
#
# @return data.frame with distance in metres and bearing in degrees

calc_distbearing=function(lonlat) {
    nr=nrow(lonlat)
    dst=distVincentySphere(lonlat[1:(nr-1),],lonlat[-1,])
    brg=finalBearing(lonlat[1:(nr-1),],lonlat[-1,])
    data.frame(distance=dst,bearing=brg)
}


angle.normalise=function(x) {
    ## normalize angle to range [-pi,pi)
    (x+pi)%%(2*pi)-pi
}


