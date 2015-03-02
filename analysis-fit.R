#####################
####  Load data  ####
#####################

source('incidence-functions.R')

karongam.dat <- prepare.data('Karonga', 'Male')
karongaf.dat <- prepare.data('Karonga', 'Female')
kisesam.dat <- prepare.data('Kisesa', 'Male')
kisesaf.dat <- prepare.data('Kisesa', 'Female')
manicm.dat <- prepare.data('Manicaland', 'Male')
manicf.dat <- prepare.data('Manicaland', 'Female')
masakam.dat <- prepare.data('Masaka', 'Male')
masakaf.dat <- prepare.data('Masaka', 'Female')
rakaim.dat <- prepare.data('Rakai', 'Male')
rakaif.dat <- prepare.data('Rakai', 'Female')
umkhanm.dat <- prepare.data('uMkhanyakude', 'Male')
umkhanf.dat <- prepare.data('uMkhanyakude', 'Female')

##############################
####  Analysis functions  ####
##############################

library(adegenet)
library(RColorBrewer)

cred.region <- function(x, y, ...)
  polygon(c(x, rev(x)), c(y[1,], rev(y[2,])), ...)


plot.fit <- function(fit, dat, name, max.inc=0.07, ages=c(20, 30, 40, 50), times = c(1999, 2003, 2007, 2011)){
  ##
  fit.sample <- extract(fit)
  ##
  par(tcl=-0.25, mgp=c(2, 0.5, 0), mar=c(3.1, 3.1, 3.1, 1.1))
  ##
  plot(NA, NA, xlim=range(dat$x_time), ylim=c(0, max.inc), xlab="Year", ylab="HIV incidence rate", main=name)
  for(idx in 1:4){
    post.incrate <- exp(dat$X_time %*% apply(fit.sample$coef_time_age, 1, "%*%", dat$X_age[dat$x_age[-1]==ages[idx],]))
    col <- brewer.pal(4, "Dark2")[idx]
    cred.region(dat$x_time[-1]-dat$dt/2, apply(post.incrate, 1, quantile, c(0.1, 0.9)),
                col=transp(col, 0.3), border=NA)
    lines(dat$x_time[-1]-dat$dt/2, apply(post.incrate, 1, mean), lty=1, lwd=2, col=col)
    lines(dat$x_time[-1]-dat$dt/2, apply(post.incrate, 1, median), lty=2, lwd=2, col=col)
  }
  legend("topright", paste("age", ages), lwd=2, col= brewer.pal(4, "Dark2"))
  ##
  plot(NA, NA, xlim=range(dat$x_age), ylim=c(0, max.inc), xlab="Age", ylab="HIV incidence rate", main=name)
  for(idx in 1:4){
    post.incrate <- exp(dat$X_age %*% apply(fit.sample$coef_time_age, 1, function(Bmat) dat$X_time[min(which(dat$x_time[-1]>=times[idx])),] %*% Bmat))
    col <- brewer.pal(4, "Dark2")[idx]
    cred.region(dat$x_age[-1]-dat$dt/2, apply(post.incrate, 1, quantile, c(0.1, 0.9)),
                col=transp(col, 0.3), border=NA)
    lines(dat$x_age[-1]-dat$dt/2, apply(post.incrate, 1, mean), lty=1, lwd=2, col=col)
    lines(dat$x_age[-1]-dat$dt/2, apply(post.incrate, 1, median), lty=2, lwd=2, col=col)
  }
  legend("topright", legend=times, lwd=2, col= brewer.pal(4, "Dark2"))
  ##
  return(NULL)
}


#####################
####  Load fits  ####
#####################

load("alpha-incidence-fits_2015-01-29.RData")

quartz(h=6, w=6, pointsize=8)

pdf("ALPHA-incidence-fits_2015-01-29.pdf", h=6, w=6, pointsize=8)

par(mfrow=c(2,2), cex=1)

plot.fit(karongam.fit, karongam.dat, "Karonga men", 0.02, times=c(2006, 2008, 2010, 2011.5))
plot.fit(karongaf.fit, karongaf.dat, "Karonga women", 0.02, times=c(2006, 2008, 2010, 2011.5))

plot.fit(kisesam.fit, kisesam.dat, "Kisesa men", 0.03, times=c(1996, 2001, 2006, 2011))
plot.fit(kisesaf.fit, kisesaf.dat, "Kisesa women", 0.03, times=c(1996, 2001, 2006, 2011))

plot.fit(manicm.fit, manicm.dat, "Manicaland men", 0.05)
plot.fit(manicf.fit, manicf.dat, "Manicaland women", 0.05)

plot.fit(masakam.fit, masakam.dat, "Masaka men", 0.02, times=c(1992, 1998, 2004, 2010))
plot.fit(masakaf.fit, masakaf.dat, "Masaka women", 0.02, times=c(1992, 1998, 2004, 2010))

plot.fit(rakaim.fit, rakaim.dat, "Rakai men", 0.03, times = c(2000, 2003.5, 2007, 2010.5))
plot.fit(rakaif.fit, rakaif.dat, "Rakai women", 0.03, times = c(2000, 2003.5, 2007, 2010.5))


plot.fit(umkhanm.fit, umkhanm.dat, "uMkhanyakude men", 0.2, times = c(2003, 2006, 2009, 2012))
plot.fit(umkhanf.fit, umkhanf.dat, "uMkhanyakude women", 0.2, times = c(2003, 2006, 2009, 2012))

dev.off()



##########################################################################
####  Outputs with gamma(1.0, 0.00005) prior, censored Africa Centre  ####
##########################################################################


load("alpha-incidence-fits_2015-02-04.RData")

umkhanm.dat <- prepare.data('uMkhanyakude', 'Male', max.time=2013.0)
umkhanf.dat <- prepare.data('uMkhanyakude', 'Female', max.time=2013.0)

quartz(h=6, w=6, pointsize=8)

pdf("ALPHA-incidence-fits_2015-02-04.pdf", h=6, w=6, pointsize=8)

par(mfrow=c(2,2), cex=1)

plot.fit(karongam.fit, karongam.dat, "Karonga men", 0.02, times=c(2006, 2008, 2010, 2011.5))
plot.fit(karongaf.fit, karongaf.dat, "Karonga women", 0.02, times=c(2006, 2008, 2010, 2011.5))

plot.fit(kisesam.fit, kisesam.dat, "Kisesa men", 0.03, times=c(1996, 2001, 2006, 2011))
plot.fit(kisesaf.fit, kisesaf.dat, "Kisesa women", 0.03, times=c(1996, 2001, 2006, 2011))

plot.fit(manicm.fit, manicm.dat, "Manicaland men", 0.05)
plot.fit(manicf.fit, manicf.dat, "Manicaland women", 0.05)

plot.fit(masakam.fit, masakam.dat, "Masaka men", 0.02, times=c(1992, 1998, 2004, 2010))
plot.fit(masakaf.fit, masakaf.dat, "Masaka women", 0.02, times=c(1992, 1998, 2004, 2010))

plot.fit(rakaim.fit, rakaim.dat, "Rakai men", 0.03, times = c(2000, 2003.5, 2007, 2010.5))
plot.fit(rakaif.fit, rakaif.dat, "Rakai women", 0.03, times = c(2000, 2003.5, 2007, 2010.5))


plot.fit(umkhanm.fit, umkhanm.dat, "uMkhanyakude men", 0.2, times = c(2003, 2006, 2009, 2012))
plot.fit(umkhanf.fit, umkhanf.dat, "uMkhanyakude women", 0.2, times = c(2003, 2006, 2009, 2012))

dev.off()
