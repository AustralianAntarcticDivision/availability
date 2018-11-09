data <- read.csv("data-raw/MAR2015_DMS_PTT_7344_7344.csv", header = TRUE)
data$date <- as.POSIXct(strptime(paste0(data$year, "-", data$month, "-", data$day, " ", data$time), format = "%Y-%m-%d %H:%M:%S", tz = as.character(data$time_zone[1])))
data <- data[!duplicated(data$date), ]
ptt_data <- list(DMS = data)

data <- read.csv("data-raw/MAR2015_SES_PTT_BB151_24647.csv",header = TRUE)
data$date <- as.POSIXct(strptime(paste0(data$year, "-", data$month, "-", data$day, " ", data$time), format = "%Y-%m-%d %H:%M:%S", tz = as.character(data$time_zone[1])))
data <- data[!duplicated(data$date), ]
ptt_data$SES <- data

data <- read.csv("data-raw/MAR2015_AFS_PTT_PP048_38622.csv",header = TRUE)
data$date <- as.POSIXct(strptime(paste0(data$year, "-", data$month, "-", data$day, " ", data$time), format = "%Y-%m-%d %H:%M:%S", tz = as.character(data$time_zone[1])))
data <- data[!duplicated(data$date), ]
ptt_data$AFS <- data

usethis::use_data(ptt_data, overwrite = TRUE)
