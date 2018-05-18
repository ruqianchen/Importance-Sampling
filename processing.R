dat <- read.csv("user_network_data.csv")
datdate <- dat$date_followed
# datdate <- as.Date(datdate)
dat$year <- format(as.Date(datdate,  "%m/%d/%Y"), "20%y")
dat$month <- format(as.Date(datdate,  "%m/%d/%Y"), "%m")
dat$day <- format(as.Date(datdate,  "%m/%d/%Y"), "%d")
# reference https://stackoverflow.com/questions/9749598/extract-month-and-year-from-a-zooyearmon-object
# reference https://stackoverflow.com/questions/9508747/add-correct-century-to-dates-with-year-provided-as-year-without-century-y
dat
write.csv(dat, "user_network_data_year.csv")
