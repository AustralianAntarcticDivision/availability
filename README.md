# Availability

Estimating geographic space available to animals based on telemetry data.

**Installation:**
```
require(devtools)
devtools::install_github("AustralianAntarcticDataCentre/availability")
```

**Example usage:**
```
library(availability)
realtrack=... ## load observed track, 2-column matrix of longitude and latitude. The track points should be equally sampled in time
arf=surrogateARModel(realtrack) ## fit AR model to track
(arfit,tp$xs,tp$ts)
## simulate new track
st=surrogateAR(arf,realtrack)
```
