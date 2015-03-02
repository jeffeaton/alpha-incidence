library(foreign)
library(lubridate)

####  Stata stset information  ####

## . st
## -> stset exit, id(idno) failure(failure) time0(entry) origin(time dob)
##                scale(365.25)

##                 id:  idno
##      failure event:  failure != 0 & failure < .
## obs. time interval:  (entry, exit]
##  exit on or before:  failure
##     t for analysis:  (time-origin)/365.25
##             origin:  time dob


## alpha <- read.dta("../ALPHA -jeff-140618.dta")
## alpha <- read.dta("~/Documents/Data/ALPHA/Gates/ALPHA11Gates_ongoing_ready_yearsplit_agesplit_all_pooled - 140619_vStata12.dta")
## alpha <- read.dta("~/Documents/Data/ALPHA/Gates/ALPHA11Gates_ongoing_ready_yearsplit_agesplit_all_pooled - 140915_vStata12.dta")
alpha <- read.dta("~/Documents/Data/ALPHA/Gates/ALPHA11Gates_ongoing_ready_yearsplit_agesplit_all_pooled - 141216_vStata12.dta")
alpha$dob <- round(alpha$dob)
alpha$entry <- round(alpha$entry)
alpha$exit <- round(alpha$exit)
alpha$last_neg_date <- as.Date(alpha$last_neg_date, origin="1960-01-01")
alpha$frst_pos_date <- as.Date(alpha$frst_pos_date, origin="1960-01-01")

alpha$site <- factor(alpha$site)

dat <- alpha[,c("idno", "sex", "dob", "site", "entry", "exit", "failure", "_st", "_d", "_origin", "_t", "_t0",
                "last_neg_date", "frst_pos_date", "first_seen_date", "last_seen_date", "hivstatus_broad", "hivstatus_detail")]
names(dat)[names(dat) == "_st"] <- "st"
names(dat)[names(dat) == "_d"] <- "d"
names(dat)[names(dat) == "_origin"] <- "origin"
names(dat)[names(dat) == "_t"] <- "t"
names(dat)[names(dat) == "_t0"] <- "t0"

dat$dur <- dat$exit - dat$entry

dat.entry <- as.Date(tapply(dat$entry, dat$idno, min), origin="1970-01-01")
dat.entry <- data.frame(idno=as.integer(names(dat.entry)), first.entry=dat.entry)
dat.exit <- as.Date(tapply(dat$exit, dat$idno, max), origin="1970-01-01")
dat.exit <- data.frame(idno=as.integer(names(dat.exit)), last.exit=dat.exit)

dat.failure <- with(subset(dat, failure == 1, c("idno", "exit")), tapply(exit, idno, max))
dat.failure <- data.frame(idno=as.integer(names(dat.failure)), exit.failure = dat.failure)

dat.firsthivepis <- as.Date(with(subset(dat, hivstatus_detail %in% c("Negative", "Positive")),
                                 tapply(entry, idno, min)), origin="1970-01-01")
dat.firsthivepis <- data.frame(idno=as.integer(names(dat.firsthivepis)), firsthivepis=dat.firsthivepis)

datm <- dat; nrow(datm)
datm <- merge(datm, dat.entry); nrow(datm)
datm <- merge(datm, dat.exit); nrow(datm)
datm <- merge(datm, dat.failure, by="idno", all.x=TRUE, suffixes=c("", ".failure")); nrow(datm)
datm <- merge(datm, dat.firsthivepis, all.x=TRUE); nrow(datm)

datw <- subset(datm, exit==last.exit); nrow(datw)
table(table(datw$idno))
datw <- datw[order(datw$failure, decreasing=TRUE),]
datw <- subset(datw, !duplicated(datw$idno)); nrow(datw)

with(datw, table(site, last_neg_date < first.entry))
with(datw, table(site, frst_pos_date < first.entry))
with(datw, table(site, firsthivepis < first.entry))

with(datw, table(site, last_neg_date > last.exit))
with(datw, table(site, frst_pos_date > last.exit))
with(datw, table(site, firsthivepis < last.exit))

datw$last_neg_date <- with(datw, replace(last_neg_date, last_neg_date < first.entry, NA))
datw$frst_pos_date <- with(datw, pmax(frst_pos_date, first.entry))

datw$frst_pos_date <- with(datw, replace(frst_pos_date, frst_pos_date > last.exit, NA))
datw$last_neg_date <- with(datw, pmin(last_neg_date, last.exit))

datw$firsttest <- with(datw, pmin(firsthivepis, last_neg_date, frst_pos_date, na.rm=TRUE))

with(datw, table(site, firsttest < first.entry))
with(datw, table(site, firsttest <= last_neg_date))
with(datw, table(site, firsttest <= frst_pos_date))
with(datw, table(site, last_neg_date <= frst_pos_date))

datw$failure <- with(datw, replace(failure, is.na(failure), 0))

datw <- setNames(datw[,c("idno", "site", "sex", "dob", "first.entry", "last.exit", "failure", "firsttest", "last_neg_date", "frst_pos_date")],
                 c("idno", "site", "sex", "dob", "entry", "exit", "death", "firsttest", "lastneg", "firstpos"))
datw$seroconvdate <- as.Date((as.integer(datw$lastneg) + as.integer(datw$firstpos))/2, origin="1970-01-01")
datw$seroconv <- !is.na(datw$seroconvdate)
datw$seroconv[is.na(datw$lastneg)] <- NA

datw <- subset(datw, !is.na(dob) & sex %in% 1:2); nrow(datw)
datw$sex <- factor(datw$sex, 1:2, c("Male", "Female"))


## convert dates to decimals
decimal_date_na <- function(x){ y <- rep(NA, length(x));  y[!is.na(x)] <- decimal_date(x[!is.na(x)]); return(y) }

datw$dob <- decimal_date(datw$dob)
datw$entry <- decimal_date(datw$entry)
datw$exit <- decimal_date(datw$exit)
datw$firsttest <- decimal_date_na(datw$firsttest)
datw$lastneg <- decimal_date_na(datw$lastneg)
datw$firstpos <- decimal_date_na(datw$firstpos)
datw$seroconvdate <- decimal_date_na(datw$seroconvdate)

## save(datw, file="episode-data-formatted_20141031.RData")
save(datw, file="episode-data-formatted_20150204.RData")
