library(splines)
library(rstan)
set_cppo("fast")  # for best running speed

## Load dataset
load("episode-data-formatted_20141031.RData")
incdat <- subset(datw, !is.na(lastneg), select=c("idno", "site", "sex", "dob", "firsttest", "lastneg", "firstpos"))
incdat <- subset(incdat, lastneg > firsttest | firstpos > firsttest)
incdat <- subset(incdat, is.na(firstpos) | firstpos > lastneg) ## !!! Investigate this more


## Declare functions

prepare.data <- function(sites, sexes, dt=0.1,
                         min.time=NULL, max.time=NULL, nk.time=NULL, fixcoef.time=NULL,
                         min.age=15.0, max.age=55.0, nk.age=20L, fixcoef.age=8L){

  ## ################## ##
  ##  Prepare the data  ##
  ## ################## ##

  dat <- subset(incdat, site %in% sites & sex %in% sexes)

  ## round dates to time steps [NOTE: think more about this]
  dat$dobTS <- round(dat$dob / dt)
  dat$firsttestTS <- round(dat$firsttest / dt)
  dat$lastnegTS <- round(dat$lastneg / dt)
  dat$firstposTS <- round(dat$firstpos / dt)

  min.timeTS <- ifelse(is.null(min.time), min(dat$firsttestTS))
  max.timeTS <- ifelse(is.null(max.time), max(dat$lastnegTS, dat$firstposTS, na.rm=TRUE), round(max.time / dt))
  min.ageTS <- round(min.age / dt)
  max.ageTS <- round(max.age / dt)

  ## truncate data age range [NOTE: slight bias in censoring, should fix]
  dat <- subset(dat, lastnegTS - dobTS >= min.ageTS)
  dat <- subset(dat, firsttestTS - dobTS <= max.ageTS)

  startageTS <- dat$firsttestTS - dat$dobTS
  startageTS[startageTS < min.ageTS] <- min.ageTS
  dat$firsttestTS <- dat$dobTS + startageTS
  
  firstposageTS <- dat$firstposTS - dat$dobTS
  dat$firstposTS[firstposageTS > max.ageTS] <- NA
  
  lastnegageTS <- dat$lastnegTS - dat$dobTS
  lastnegageTS[lastnegageTS > max.ageTS] <- max.ageTS
  dat$lastnegTS <- dat$dobTS + lastnegageTS

  rm(startageTS, firstposageTS, lastnegageTS)

  ## truncate data time range
  dat <- subset(dat, lastnegTS >= min.timeTS)  # NOTE: bias in this censoring
  dat <- subset(dat, firsttestTS <= max.timeTS)
  
  dat$firsttestTS[dat$firsttestTS < min.timeTS] <- min.timeTS
  dat$firstposTS[dat$firstposTS > max.timeTS] <- NA ## BIAS!!!!
  dat$lastnegTS[dat$lastnegTS > max.timeTS] <- max.timeTS

  ## calculate array indices
  dat$firsttest.aIDX <- as.integer(dat$firsttestTS - dat$dobTS - min.ageTS) + 1
  dat$lastneg.aIDX <- as.integer(dat$lastnegTS - dat$dobTS - min.ageTS) + 1
  dat$firstpos.aIDX <- as.integer(dat$firstposTS - dat$dobTS - min.ageTS) + 1
  dat$firsttest.tIDX <- as.integer(dat$firsttestTS - min.timeTS) + 1
  dat$lastneg.tIDX <- as.integer(dat$lastnegTS - min.timeTS) + 1
  dat$firstpos.tIDX <- as.integer(dat$firstposTS - min.timeTS) + 1

  ## ###################### ##
  ##  Prepare spline model  ##
  ## ###################### ##

  x.time <- seq(min.timeTS*dt, max.timeTS*dt, dt)
  x.age <- seq(min.ageTS*dt, max.ageTS*dt, dt)

  if(is.null(nk.time)){ # place knots every half year for span of data
    k.time <- seq(0.5 * (min.timeTS * dt) %/% 0.5 - 1.5, -0.5 * (-max.timeTS * dt) %/% 0.5 + 1.5, 0.5)
    nk.time <- length(k.time) - 4L
  } else {
    time.dur <- max.timeTS * dt - min.timeTS * dt
    k.time <- seq(min.time*dt - 3*time.dur/(nk.time-3), max.time*dt + 3*time.dur/(nk.time-3), time.dur/(nk.time-3))
  }

  age.dur <- max.age - min.age
  k.age <- seq(min.ageTS*dt - 3*age.dur/(nk.age-3), max.ageTS*dt + 3*age.dur/(nk.age-3), age.dur/(nk.age-3))

  X.time <- splineDesign(k.time, x.time[-1]-dt/2)
  P.time <- diff(diag(nk.time), diff=1)

  X.age <- splineDesign(k.age, x.age[-1]-dt/2)
  P.age <- diff(diag(nk.age), diff=1)

  if(is.null(fixcoef.time))
    fixcoef.time <- as.integer(nk.time / 2)


  ## ################### ##
  ##  List of Stan data  ##
  ## ################### ##

  stan.data <- list(N=nrow(dat),
                    Nseroconv=sum(!is.na(dat$firstposTS)),
                    STEPS_time=length(x.time),
                    STEPS_age=length(x.age),
                    dt=dt,
                    nk_time=nk.time,
                    x_time=x.time,
                    X_time=X.time,
                    P_time=P.time,
                    fixcoef_time_idx=fixcoef.time,
                    nk_age=nk.age,
                    x_age=x.age,
                    X_age=X.age,
                    P_age=P.age,
                    fixcoef_age_idx=fixcoef.age,
                    firsttest_tIDX=dat$firsttest.tIDX,
                    lastneg_tIDX=dat$lastneg.tIDX,
                    lastneg_seroconv_tIDX=dat$lastneg.tIDX[!is.na(dat$firstpos.tIDX)],
                    firstpos_tIDX=dat$firstpos.tIDX[!is.na(dat$firstpos.tIDX)],
                    firsttest_aIDX=dat$firsttest.aIDX,
                    lastneg_aIDX=dat$lastneg.aIDX,
                    lastneg_seroconv_aIDX=dat$lastneg.aIDX[!is.na(dat$firstpos.aIDX)],
                    firstpos_aIDX=dat$firstpos.aIDX[!is.na(dat$firstpos.aIDX)])

  return(stan.data)
}


stan.code <- "
functions{

  vector cumsum(vector x){

     vector[rows(x)+1] val;
     val[1] <- 0;
     for(i in 2:rows(val))
       val[i] <- val[i-1] + x[i-1];

     return val;
  }

  matrix diagCumSum(matrix x){

    matrix[rows(x)+1, cols(x)+1] val;

    for(i in 1:rows(val))
      val[i,1] <- 0;
    for(j in 2:cols(val)){
      val[1,j] <- 0;
      for(i in 2:rows(val))
        val[i,j] <- val[i-1,j-1] + x[i-1,j-1];
    }

    return val;
  }
}
data {
        // spline model setup
        int<lower=1> STEPS_time;
        int<lower=1> STEPS_age;
        real<lower=0> dt;
        int<lower=5> nk_time;
        matrix[STEPS_time-1, nk_time] X_time;
        matrix[nk_time-1, nk_time] P_time;
        int fixcoef_time_idx;

        int<lower=5> nk_age;
        matrix[STEPS_age-1, nk_age] X_age;
        matrix[nk_age-1, nk_age] P_age;
        int fixcoef_age_idx;

        // data
        int<lower=1> N;
        int<lower=0> Nseroconv;

        int<lower=1, upper=STEPS_time> firsttest_tIDX[N];
        int<lower=1, upper=STEPS_time> lastneg_tIDX[N];
        int<lower=1, upper=STEPS_time> lastneg_seroconv_tIDX[Nseroconv];
        int<lower=1, upper=STEPS_time> firstpos_tIDX[Nseroconv];

        int<lower=1, upper=STEPS_age> firsttest_aIDX[N];
        int<lower=1, upper=STEPS_age> lastneg_aIDX[N];
        int<lower=1, upper=STEPS_age> lastneg_seroconv_aIDX[Nseroconv];
        int<lower=1, upper=STEPS_age> firstpos_aIDX[Nseroconv];
}
parameters {
        matrix[nk_time, nk_age] coef_time_age;
        real<lower=0> sigma2_time;
        real<lower=0> sigma2_age;
        real<lower=0> sigma2_time_age;
}
transformed parameters{
}
model {
        vector[nk_time] coef_time;
        row_vector[nk_age] coef_age;

        matrix[STEPS_time-1, STEPS_age-1] incrate_time_age;
        matrix[STEPS_time, STEPS_age] logcumavoid;

        1/sigma2_time ~ gamma(1.0, 0.00005);
        1/sigma2_age ~ gamma(1.0, 0.00005);
        1/sigma2_time_age ~ gamma(1.0, 0.00005);

        coef_time <- col(coef_time_age, fixcoef_age_idx);
        coef_age <- row(coef_time_age, fixcoef_time_idx) - coef_time_age[fixcoef_time_idx, fixcoef_age_idx] ;

        P_time * coef_time ~ normal(0, sqrt(sigma2_time));
        P_age * coef_age' ~ normal(0, sqrt(sigma2_age));
        to_vector(P_time * (coef_time_age - rep_matrix(coef_time, nk_age) - rep_matrix(coef_age, nk_time)) * P_age') ~ normal(0, sqrt(sigma2_time_age));

        incrate_time_age <- exp(X_time * coef_time_age * X_age');
        logcumavoid <- -diagCumSum(incrate_time_age) * dt;

        for(i in 1:N)
          increment_log_prob(logcumavoid[lastneg_tIDX[i], lastneg_aIDX[i]] - logcumavoid[firsttest_tIDX[i], firsttest_aIDX[i]]);

        for(i in 1:Nseroconv)
          increment_log_prob(log1m_exp(logcumavoid[firstpos_tIDX[i], firstpos_aIDX[i]] - logcumavoid[lastneg_seroconv_tIDX[i], lastneg_seroconv_aIDX[i]]));
}
"
