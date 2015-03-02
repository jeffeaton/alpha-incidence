####  Notes: chains 10:15 correspond to prior variance IG(1, 0.005), run on 28 January 2015
####         chains 20:25 correspond to prior variance IG(1, 0.00005), run on 3 February 2015

setwd(Sys.getenv("WORK"))
library(rstan)

lapply(paste("karongam-fit", 10:15, ".RData", sep=""), load, globalenv())
lapply(paste("karongaf-fit", 10:15, ".RData", sep=""), load, globalenv())
lapply(paste("kisesam-fit", 10:15, ".RData", sep=""), load, globalenv())
lapply(paste("kisesaf-fit", 10:15, ".RData", sep=""), load, globalenv())
lapply(paste("manicm-fit", 10:15, ".RData", sep=""), load, globalenv())
lapply(paste("manicf-fit", 10:15, ".RData", sep=""), load, globalenv())
lapply(paste("masakam-fit", 10:15, ".RData", sep=""), load, globalenv())
lapply(paste("masakaf-fit", 10:15, ".RData", sep=""), load, globalenv())
lapply(paste("rakaim-fit", 10:15, ".RData", sep=""), load, globalenv())
lapply(paste("rakaif-fit", 10:15, ".RData", sep=""), load, globalenv())
lapply(paste("umkhanm-fit", 10:15, ".RData", sep=""), load, globalenv())
lapply(paste("umkhanf-fit", 10:15, ".RData", sep=""), load, globalenv())

karongam.fit <- sflist2stanfit(lapply(paste("karongam.fit", 10:15, sep=""), get))
karongaf.fit <- sflist2stanfit(lapply(paste("karongaf.fit", 10:15, sep=""), get))
kisesam.fit <- sflist2stanfit(lapply(paste("kisesam.fit", 10:15, sep=""), get))
kisesaf.fit <- sflist2stanfit(lapply(paste("kisesaf.fit", 10:15, sep=""), get))
manicm.fit <- sflist2stanfit(lapply(paste("manicm.fit", 10:15, sep=""), get))
manicf.fit <- sflist2stanfit(lapply(paste("manicf.fit", 10:15, sep=""), get))
masakam.fit <- sflist2stanfit(lapply(paste("masakam.fit", 10:15, sep=""), get))
masakaf.fit <- sflist2stanfit(lapply(paste("masakaf.fit", 10:15, sep=""), get))
rakaim.fit <- sflist2stanfit(lapply(paste("rakaim.fit", 10:15, sep=""), get))
rakaif.fit <- sflist2stanfit(lapply(paste("rakaif.fit", 10:15, sep=""), get))
umkhanm.fit <- sflist2stanfit(lapply(paste("umkhanm.fit", 10:15, sep=""), get))
umkhanf.fit <- sflist2stanfit(lapply(paste("umkhanf.fit", 10:15, sep=""), get))


save(karongam.fit, karongaf.fit, kisesam.fit, kisesaf.fit,
     manicm.fit, manicf.fit, masakam.fit, masakaf.fit,
     rakaim.fit, rakaif.fit, umkhanm.fit, umkhanf.fit,
     file="alpha-incidence-fits_2015-01-29.RData")

############################################################
############################################################

lapply(paste("karongam-fit", 20:25, ".RData", sep=""), load, globalenv())
lapply(paste("karongaf-fit", 20:25, ".RData", sep=""), load, globalenv())
lapply(paste("kisesam-fit", 20:25, ".RData", sep=""), load, globalenv())
lapply(paste("kisesaf-fit", 20:25, ".RData", sep=""), load, globalenv())
lapply(paste("manicm-fit", 20:25, ".RData", sep=""), load, globalenv())
lapply(paste("manicf-fit", 20:25, ".RData", sep=""), load, globalenv())
lapply(paste("masakam-fit", 20:25, ".RData", sep=""), load, globalenv())
lapply(paste("masakaf-fit", 20:25, ".RData", sep=""), load, globalenv())
lapply(paste("rakaim-fit", 20:25, ".RData", sep=""), load, globalenv())
lapply(paste("rakaif-fit", 20:25, ".RData", sep=""), load, globalenv())
lapply(paste("umkhanm-fit", 20:25, ".RData", sep=""), load, globalenv())
lapply(paste("umkhanf-fit", 20:25, ".RData", sep=""), load, globalenv())

karongam.fit <- sflist2stanfit(lapply(paste("karongam.fit", 20:25, sep=""), get))
karongaf.fit <- sflist2stanfit(lapply(paste("karongaf.fit", 20:25, sep=""), get))
kisesam.fit <- sflist2stanfit(lapply(paste("kisesam.fit", 20:25, sep=""), get))
kisesaf.fit <- sflist2stanfit(lapply(paste("kisesaf.fit", 20:25, sep=""), get))
manicm.fit <- sflist2stanfit(lapply(paste("manicm.fit", 20:25, sep=""), get))
manicf.fit <- sflist2stanfit(lapply(paste("manicf.fit", 20:25, sep=""), get))
masakam.fit <- sflist2stanfit(lapply(paste("masakam.fit", 20:25, sep=""), get))
masakaf.fit <- sflist2stanfit(lapply(paste("masakaf.fit", 20:25, sep=""), get))
rakaim.fit <- sflist2stanfit(lapply(paste("rakaim.fit", 20:25, sep=""), get))
rakaif.fit <- sflist2stanfit(lapply(paste("rakaif.fit", 20:25, sep=""), get))
umkhanm.fit <- sflist2stanfit(lapply(paste("umkhanm.fit", 20:25, sep=""), get))
umkhanf.fit <- sflist2stanfit(lapply(paste("umkhanf.fit", 20:25, sep=""), get))

save(karongam.fit, karongaf.fit,
     kisesam.fit, kisesaf.fit,
     manicm.fit, manicf.fit,
     masakam.fit, masakaf.fit,
     rakaim.fit, rakaif.fit,
     umkhanm.fit, umkhanf.fit,
     file="alpha-incidence-fits_2015-02-04.RData")
