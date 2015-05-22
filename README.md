# Availability

Estimating geographic space available to animals based on telemetry data. This is entirely experimental code, use at your own risk. This code is likely to change in the future and in particular may be reimplemented to use existing R package functionality.

**Installation:**
```
require(devtools)
devtools::install_github("AustralianAntarcticDivision/availability")
```

**Example usage:**
```
library(availability)
realtrack=... ## load observed track, 2-column matrix of longitude and latitude. The track points should be equally sampled in time
arf=surrogate_arfit(realtrack)
st=surrogate_arsimulate(arf,nrow(realtrack),startlonlat=realtrack[1,],do.test.land=TRUE,random.rotation=c(-pi,pi))
```

*Note*: `do.test.land` currently will not work in the northern hemisphere, because the included land mask only covers south of the equator.
