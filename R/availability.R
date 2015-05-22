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


#' Fit vector-autoregressive model to track
#'
#' @author Ben Raymond
#'
#' @param lonlat array: 2-column matrix or data.frame with longitude, latitude of each point
#' @param model.order numeric: the order of the AR model to fit, default=1
#'
#' @return object of class "ar"
#'
#' @seealso \code{\link{ar}} \code{\link{surrogate_arsimulate}}
#'
#' @export surrogate_arfit

surrogate_arfit=function(lonlat,model.order=1) {
    ## calculate dx, dy separately at each time step
    nr=nrow(lonlat)
    dx=distVincentyEllipsoid(lonlat[1:(nr-1),],cbind(lonlat[-1,1],lonlat[1:(nr-1),2]))*sign(lonlat[1:(nr-1),1]-lonlat[-1,1])
    dy=distVincentyEllipsoid(lonlat[1:(nr-1),],cbind(lonlat[1:(nr-1),1],lonlat[-1,2]))*sign(lonlat[1:(nr-1),2]-lonlat[-1,2])

    dxdy=data.frame(dx=dx,dy=dy)
    ar(dxdy,order.max=model.order,aic=FALSE)
}


#' Simulate track from fitted vector autoregressive model
#'
#' @author Ben Raymond
#'
#' @param arfit : fitted object of class "ar" as returned by \code{\link{surrogate_arfit}}
#' @param n numeric: number of points to simulate
#' @param startlonlat numeric: 2-element array with starting longitude and latitude
#' @param endlonlat numeric: 2-element array with ending longitude and latitude. If NULL, no end constraint is imposed except for land masking (if do.test.land is TRUE)
#' @param do.test.land logical: use the included land mask to avoid land?
#' @param random.rotation numeric: 2-element array giving the range of the rotation to apply to the randomized track (values in radians). use random.rotation=NULL for no such rotation. The angle can be restricted using random.rotation=c(min.angle,max.angle) - this may speed up computation by avoiding impossible angles (e.g. tracks over a land mass)
#'
#' @return 2-column data.frame with longitude,latitude of simulated track points
#'
#' @seealso \code{\link{surrogate_arfit}}
#'
#' @export surrogate_arsimulate

surrogate_arsimulate=function(arfit,n,startlonlat,endlonlat=NULL,do.test.land=TRUE,random.rotation=c(-pi,pi)) {
    ## random.rotation=c(-pi,pi) will apply random rotation to parms before simulating
    ## use random.rotation=NULL for no such rotation
    ## angle can be restricted using random.rotation=c(min.angle,max.angle) - this may speed up
    ## computation by avoiding impossible angles (e.g. tracks over the continent)

    if (!is.null(endlonlat)) {
        endlonlat=as.numeric(endlonlat)
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
            simtrack=surrogate_arsimulate(arfit=rotated.arfit,n=n,startlonlat=startlonlat,endlonlat=endlonlat,do.test.land=do.test.land,random.rotation=NULL)
            if (dim(simtrack)[1]>0) {
                break
            }
        }
        return(simtrack)
    }
    if (do.test.land) {
        land.mask=readPNG(system.file("extdata","land_mask-0.1-nosub.png",package="availability")) ## 0=land, 1=ocean
        land.lon=seq(from=-180,to=180,length.out=dim(land.mask)[2])
        land.lat=seq(from=0,to=-90,length.out=dim(land.mask)[1])
    }
    A=matrix(arfit$ar,nrow=2) ## fill down columns
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
        for (land.tries in 1:10) {
            thisrand=matrix(rnorm(2) %*% tempchol,nrow=1)
            xsim[k,]=t(A %*% t(xsim[k-1,]-fitted.mean))+fitted.mean+thisrand ## simulated dx,dy for this time step
            ## calculate track point from steps
            temp=destPoint(simtrack[k-1,],90,xsim[k,1]) # x step
            simtrack[k,1]=temp[1]
            simtrack[k,1]=angle.normalise(simtrack[k,1]/180*pi)/pi*180 ## ensure is in range -180 to 180
            temp=destPoint(simtrack[k-1,],0,xsim[k,2]) # y step
            simtrack[k,2]=temp[2]
            if (!is.null(endlonlat)) {
                ## as we get closer to the end of our track, increasingly nudge the random point towards the designated ending location
                a=diag(1/(n-k+1),2)
                simtrack[k,]=simtrack[k,]+a%*%(endlonlat-simtrack[k,])
            }
            if (do.test.land) {
                ## test if point over land
                lonidx=which.min(abs(land.lon-simtrack[k,1]))
                latidx=which.min(abs(land.lat-simtrack[k,2]))
                point.okay=land.mask[latidx,lonidx]==1
            }
            if (point.okay) {
                break
            }
        }
        if (! point.okay) {
            ## could not find a valid point at this step: give up
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


