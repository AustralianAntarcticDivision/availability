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
#' @param endlonlat numeric: 2-element array with ending longitude and latitude. If NULL, no end constraint is imposed except for land masking (if \code{do.test.land} is TRUE). This is a simple way of imposing a return-to-starting-location constraint; for more complex constraints use the \code{fixed} argument
#' @param do.test.land logical: use the included land mask to avoid land? Note that land masking is ignored for fixed points. Note also that it is possible to create a sitation where tracks are difficult or impossible to simulate, because a fixed point is sufficiently far onto land that the track cannot reach it.
#' @param random.rotation numeric: 2-element array giving the range of the rotation to apply to the randomized track (values in radians). use \code{random.rotation=NULL} for no such rotation. The angle can be restricted using \code{random.rotation=c(min.angle,max.angle)} - this may speed up computation by avoiding impossible angles (e.g. tracks over a land mass)
#' @param verbose logical: if TRUE, spit out extra information which may be helpful if things don't work as expected
#'
#' @return 2-column data.frame with longitude,latitude of simulated track points
#'
#' @seealso \code{\link{surrogate_arfit}}
#'
#' @export surrogate_arsimulate

surrogate_arsimulate=function(arfit,n,startlonlat,fixed=NULL,endlonlat=NULL,do.test.land=TRUE,random.rotation=c(-pi,pi),verbose=FALSE) {
    if (!is.null(endlonlat) && !is.null(fixed)) {
        stop("only one of fixed or endlonlat can be supplied")
    }
    if (! is.null(random.rotation)) {
        this.rotation=0
        ## apply rotation to arfit parms
        for (ntries in 1:100) {
            rotate.by=angle.normalise(runif(1)*diff(range(random.rotation))+min(random.rotation))
            Rm=matrix(c(cos(rotate.by),-sin(rotate.by),sin(rotate.by),cos(rotate.by)),nrow=2,byrow=TRUE)
            rotated.arfit=arfit
            rotated.arfit$ar=Rm %*% matrix(arfit$ar,nrow=2) %*% t(Rm)
            rotated.arfit$var.pred=Rm %*% as.matrix(arfit$var.pred) %*% t(Rm)
#            print(matrix(arfit$x.mean,nrow=1))
            rotated.arfit$x.mean=matrix(arfit$x.mean,nrow=1) %*% t(Rm)
            ## call simulate on rotated parms
            simtrack=surrogate_arsimulate(arfit=rotated.arfit,n=n,startlonlat=startlonlat,fixed=fixed,endlonlat=endlonlat,do.test.land=do.test.land,random.rotation=NULL,verbose=verbose)
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
    if (do.test.land) {
        land.mask=readPNG(system.file("extdata","land_mask-0.1-nosub.png",package="availability")) ## 0=land, 1=ocean
        land.lon=seq(from=-180,to=180,length.out=dim(land.mask)[2])
        land.lat=seq(from=0,to=-90,length.out=dim(land.mask)[1])
    }
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

    simtrack=matrix(0,n,2)
    simtrack[1,]=as.numeric(startlonlat)
    for (k in 2:n) {
        point.okay=TRUE
        if (verbose) cat(sprintf("Step %d, current location is %.3f, %.3f\n",k,simtrack[k-1,1],simtrack[k-1,2]))
        for (land.tries in 1:10) {
            thisrand=matrix(rnorm(2) %*% tempchol,nrow=1)
            xsim[k,]=t(A %*% t(xsim[k-1,]-fitted.mean))+fitted.mean+thisrand ## simulated dx,dy for this time step
            ## calculate track point from steps
            temp=destPoint(simtrack[k-1,],90,xsim[k,1]) # x step
            simtrack[k,1]=temp[1]
            simtrack[k,1]=angle.normalise(simtrack[k,1]/180*pi)/pi*180 ## ensure is in range -180 to 180
            temp=destPoint(simtrack[k-1,],0,xsim[k,2]) # y step
            simtrack[k,2]=temp[2]
            ##if (!is.null(endlonlat)) {
            ##    ## as we get closer to the end of our track, increasingly nudge the random point towards the designated ending location
            ##    a=diag(1/(n-k+1),2)
            ##    simtrack[k,]=simtrack[k,]+a%*%(endlonlat-simtrack[k,])
            ##}

            ## as we get closer to the next fixed point, increasingly nudge the random point towards the designated fixed location
            next_fixed=if (any(fixed$index>=k)) which.max(fixed$index>=k) else NA
            ##    cat(sprintf("%d: next_fixed=%d (%d) [diff=%.3f]\n",k,next_fixed,fixed$index[next_fixed],1/(fixed$index[next_fixed]-k+1)))
            if (verbose) cat(sprintf("  proposed point %d is at %.3f, %.3f\n",k,simtrack[k,1],simtrack[k,2]))
            if (!is.na(next_fixed)) {
                if (verbose) cat(sprintf("    the next fixed point is at %.3f, %.3f in %d steps time\n",fixed$lon[next_fixed],fixed$lat[next_fixed],fixed$index[next_fixed]-k))
                ## next_fixed is row index into fixed
                a=diag(1/(fixed$index[next_fixed]-k+1),2)
                simtrack[k,]=simtrack[k,]+a%*%(c(fixed$lon[next_fixed],fixed$lat[next_fixed])-simtrack[k,])
                if (verbose) cat(sprintf("    the proposed point has been nudged to %.3f, %.3f because of the next fixed point\n",simtrack[k,1],simtrack[k,2]))
            }
            if (do.test.land) {
                if (fixed$index[next_fixed]!=k) {
                    ## test if point over land, but not if this is a fixed point
                    lonidx=which.min(abs(land.lon-simtrack[k,1]))
                    latidx=which.min(abs(land.lat-simtrack[k,2]))
                    point.okay=land.mask[latidx,lonidx]==1
                } else {
                    if (verbose) cat("    not checking land-mask for this point, because it is a fixed point.\n")
                }
            }
            if (point.okay) {
                if (verbose & do.test.land) cat(sprintf("    this proposed point does not lie on land, accepting\n"))
                break
            } else {
                if (verbose) cat(sprintf("    this proposed point lies on land\n"))
            }
        }
        if (! point.okay) {
            ## could not find a valid point at this step: give up
            if (verbose) cat(sprintf("  could not find valid point, abandoning this track and starting again\n"))
            simtrack=NULL
            break
        }
    }
    data.frame(lon=simtrack[,1],lat=simtrack[,2])
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


