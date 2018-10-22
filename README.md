# Availability

Estimating geographic space available to animals based on telemetry data.

## Installation
```{r}
library(devtools)
install_github("AustralianAntarcticDataCentre/availability", build_vignettes = TRUE)
```

## Example usage
Using the vector-AR method (see [1]):
```{r}
library(availability)
## load observed track, 2-column matrix of longitude and latitude
## the track points should be equally sampled in time
realtrack <- ...
arf <- surrogateARModel(realtrack) ## fit AR model to track
st <- surrogateAR(arf, realtrack) ## simulate new track
```

Or using the crawl-based track simulator:
```{r}
library(availability)
library(crawl)
## fit crawl to your raw data
fit <- crwMLE(...)
## regularly-spaced times for which you want positions
predTime <- ...
## extract predicted positions at those times
predObj <- crwPredict(fit, predTime = predTime, speedEst = TRUE, flat = TRUE)
## keep only regularly-interpolated locations
pr <- data.frame(date = predTime, predObj[predObj$locType == "p", ])
## construct the corresponding transition and covariance matrices of the state space model for time increment `dt`
dt <- ...
model <- surrogateCrawlModel(fit, dt)
## and finally simulate the track
stcrw <- surrogateCrawl(model, as.matrix(pr[, c("mu.x", "mu.y", "nu.x", "nu.y")]), pr$date)
```

## Crawl note
This package was developed with version 1 of the crawl package. It hasn't yet been tested with crawl v2.

## Vignette
More detailed usage examples are in the package vignette (do `vignette("availability")` in R, provided that you did `install_github(..., build_vignettes = TRUE)`). The vignette is also [available as a PDF here](./vignettes/availability.pdf?raw=true).

## References
[1] Raymond B *et al.* (2014) Important marine habitat off East Antarctica revealed by two decades of multi-species predator tracking. *Ecography* **38**:121–129. [doi:10.1111/ecog.01021](http://dx.doi.org/10.1111/ecog.01021)
