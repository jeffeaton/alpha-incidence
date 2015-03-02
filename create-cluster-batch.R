
create.cluster.scripts <- function(label, site, sex, chains, iter=1000, max.time="NULL"){

  for(i in chains){
    fileConn <- file(paste("cluster/", label, "-fit", i, ".R", sep=""))
    writeLines(c("setwd(paste(Sys.getenv('HOME'), '/ALPHA-incidence/', sep=''))",
                 "source('incidence-functions.R')",
                 paste(label, ".dat <- prepare.data('", site, "', '", sex, "', max.time=", max.time, ")", sep=""),
                 paste(label, ".fit", i, " <- stan(model_code = stan.code, data = ", label, ".dat, iter = ", iter, ", chains=1, chain_id=", i, ", refresh=10)", sep=""),
                 paste("save(", label, ".fit", i, ", file=paste(Sys.getenv('WORK'), '/", label, "-fit", i, ".RData', sep=''))", sep="")), fileConn)
    close(fileConn)

    fileConn <- file(paste("cluster/", label, "-fit", i, ".pbs", sep=""))
    writeLines(c("#!/bin/sh", 
                 "#PBS -l walltime=48:00:00",
                 "#PBS -l select=01:ncpus=1:mem=4gb",
                 "#PBS -j oe",
                 "",
                 "module load R/3.1.0",
                 "module load intel-suite",
                 "",
                 paste("R CMD BATCH --no-restore --no-save $HOME/ALPHA-incidence/", label, "-fit", i, ".R $WORK/", label, "-fit", i, ".Rout", sep=""),
                 "",
                 "qstat -f $PBS_JOBID"), fileConn)
    close(fileConn)
  }

  return(NULL)
}


create.cluster.scripts("karongam", "Karonga", "Male", 20:25)
create.cluster.scripts("karongaf", "Karonga", "Female", 20:25)

create.cluster.scripts("kisesam", "Kisesa", "Male", 20:25)
create.cluster.scripts("kisesaf", "Kisesa", "Female", 20:25)

create.cluster.scripts("manicm", "Manicaland", "Male", 20:25)
create.cluster.scripts("manicf", "Manicaland", "Female", 20:25)

create.cluster.scripts("masakam", "Masaka", "Male", 20:25)
create.cluster.scripts("masakaf", "Masaka", "Female", 20:25)

create.cluster.scripts("rakaim", "Rakai", "Male", 20:25)
create.cluster.scripts("rakaif", "Rakai", "Female", 20:25)

create.cluster.scripts("umkhanm", "uMkhanyakude", "Male", 20:25, max.time=2013.0)
create.cluster.scripts("umkhanf", "uMkhanyakude", "Female", 20:25, max.time=2013.0)
