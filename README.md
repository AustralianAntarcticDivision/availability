# availability

An R package for estimating geographic space available to animals based on telemetry data.

## Installation
```{r}
library(remotes)
install_github("AustralianAntarcticDivision/availability")
```

## Minimal usage examples

```{r}
library(availability)
```

Using the vector-AR method (see [1]):

```{r}
## load observed track, 2-column matrix of longitude and latitude
## the track points should be equally sampled in time
realtrack <- ... ## your data here

arf <- surrogateARModel(realtrack) ## fit AR model to track

st <- surrogateAR(arf, realtrack) ## simulate new track
```

Or using the crawl-based track simulator:

```{r}
library(crawl)

## fit a crawl model to your raw track data
fit <- crwMLE(...)

## regularly-spaced times for which you want positions
time_step <- 3 ## e.g. using a time step of 3 hours
predTime <- seq(my_starting_date, my_ending_date, by = time_step*3600)

## extract predicted positions at those times
predObj <- crwPredict(fit, predTime = predTime, speedEst = TRUE, flat = TRUE)

## keep only regularly-interpolated locations
pr <- data.frame(date = predTime, predObj[predObj$locType == "p", ])

## construct the corresponding transition and covariance matrices of the
##  state space model
model <- surrogateCrawlModel(fit, time_step)

## and finally simulate the track
stcrw <- surrogateCrawl(model, as.matrix(pr[, c("mu.x", "mu.y", "nu.x", "nu.y")]), pr$date)
```

## Crawl note

Please note: this package was developed with version 1 of the `crawl` package. It should also work with `crawl` v2, but note that v2 only works with projected coordinates (not longitude and latitude).


## Vignette
More detailed usage examples are in the [package vignette](https://australianantarcticdatacentre.github.io/availability/articles/availability.html).

## References
[1] Raymond B *et al.* (2015) Important marine habitat off East Antarctica revealed by two decades of multi-species predator tracking. *Ecography*. [doi:10.1111/ecog.01021](https://doi.org/10.1111/ecog.01021)

[2] Reisinger RR *et al* (2018) Habitat modelling of tracking data from multiple marine top predators reveals important habitat in the Southern Indian Ocean. *Diversity and Distributions*. [doi:10.1111/ddi.12702](https://doi.org/10.1111/ddi.12702)
