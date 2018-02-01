library(devtools)
WD <- "~/Dropbox/research/Caleb/imaging/SSKmeans/package"
setwd(WD)

create("SSKmeans")


setwd(WD)
setwd("SSKmeans")
document()

install()
