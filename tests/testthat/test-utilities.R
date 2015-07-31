context("Test various utility functions")

test_that("calc_distbearing works as expected", {
    expect_equal(names(calc_distbearing(data.frame(lon=100:101,lat=-65:-64))),c("distance","bearing"))
})

                                        #> calc_distbearing(data.frame(lon=100:110,lat=-65:-55))
#   distance  bearing
#1  121193.6 22.86787
#2  121894.7 23.62824
#3  122610.3 24.37369
#4  123339.4 25.10412
#5  124080.7 25.81948
#6  124833.1 26.51971
#7  125595.7 27.20480
#8  126367.2 27.87476
#9  127146.6 28.52962
#10 127932.7 29.16943
