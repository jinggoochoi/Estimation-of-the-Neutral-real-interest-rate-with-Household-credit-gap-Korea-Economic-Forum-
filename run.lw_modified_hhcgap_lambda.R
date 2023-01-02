
setwd("d:/LW_replication");

##----------------------------------------------------------------------------------##
## File:        run.lw_modified_hhcgap_lambda.R
##
## Description: calculate lambda_g and lambda_z for estimating NRIR by LW model(2003),
##              which includes hhcgap term in IS curve
##
## original file : http://www.frbsf.org/economic-research/economists/LW_replication.zip
## Modified by Choi
##----------------------------------------------------------------------------------##

rm(list=ls())

##------------------------------------------------------------------------------##
## Load required packages and source all programs to be used in HLW estimation.
##------------------------------------------------------------------------------##

if (!require("tis")) {install.packages("tis"); library("tis")} ## Time series package
if (!require("nloptr")) {install.packages("nloptr"); library("nloptr")} ## Optimization
if (!require("mFilter")) {install.packages("mFilter"); library("mFilter")} ## HP filter

## Source all R programs

source("kalman.log.likelihood.R")
source("utilities.R")
source("kalman.states.R")
source("median.unbiased.estimator.stage1.R")
source("median.unbiased.estimator.stage2.R")

source("log.likelihood.wrapper_hhc.R")
source("kalman.states.wrapper_hhc.R")
source("calculate.covariance_hhc.R")
source("kalman.standard.errors_hhc.R")
source("format.output_hhc.R")
source("rstar.stage1_modified_hhc.R")
source("rstar.stage2_modified_hhc.R")
source("unpack.parameters.stage1_modified_hhc.R")
source("unpack.parameters.stage2_modified_hhc.R")

##------------------------------------------------------------------------------##
## Define variables
##------------------------------------------------------------------------------##

## Upper bound on a_3 parameter (slope of the IS curve)

a.r.constraint.stage2 <- -0.0025
a.r.constraint.stage3 <- -0.0025 
c.constraint.lb<-0
c.constraint.ub<-1

## Lower bound on b_2 parameter (slope of the Phillips curve)

b.y.constraint.stage1 <- NA
b.y.constraint.stage2 <- 0.025
b.y.constraint.stage3 <- 0.025
rho.z.constraint.lb<-0
rho.z.constraint.ub<-1

## Set the start and end dates of the estimation sample as well as the data start date (format is c(year,quarter))
## due to using output gap estimated from the research department, data.start = c(1991,1)

data.start   <- c(1987,1)
sample.start <- c(1995,1) 
sample.end   <- c(2022,4)

## The estimation process uses data beginning 8 quarters prior to the sample start
est.data.start    <- shiftQuarter(sample.start,-8)

## Set start index for y for initialization of state vector
g.pot.start.index <- 1 + ti(shiftQuarter(sample.start,-3),'quarterly')-ti(est.data.start,'quarterly')

## Set column names for CSV output
output.col.names <- c("Date","rstar (filtered)","g (filtered)","z (filtered)","output gap (filtered)","","All results are output from the Stage 3 model.",rep("",12),"Standard Errors","Date","y*","r*","g","","rrgap","","Date","rstar (smoothed)","g (smoothed)","z (smoothed)","output gap (smoothed)")

## Set number of iterations for Monte Carlo standard error procedure
niter <- 5000

## Because the MC standard error procedure is time consuming, we include a run switch
## Set run.se to TRUE to run the procedure
run.se <- TRUE


##------------------------------------------------------------------------------##
## Read in data and compute inflation expectations series
##------------------------------------------------------------------------------##
data <- read.table("LW_input_data_kr_rp.csv",
                   sep=',',header=TRUE,stringsAsFactors=FALSE, na.strings=".")

## Get series beginning in data.start
log.output.all                      <- tis(data$gdp.log, start=data.start, tif='quarterly')
inflation.all                       <- tis(data$inflation, start=data.start, tif='quarterly')
hhc.all    <- tis(data$hhc, start=data.start, tif='quarterly') 
relative.import.price.inflation.all <- tis(data$import, start=data.start,
                                           tif='quarterly') - inflation.all
nominal.interest.rate.all           <- tis(data$interest, start=data.start, tif='quarterly')
output.gap <- tis(data$output.gap, start=data.start, tif='quarterly')
pvc.all    <- tis(data$pvc, start=data.start, tif='quarterly') 

## Get data in vector form beginning at est.data.start (set above)
log.output                      <- as.numeric(window(log.output.all, start=est.data.start))
inflation                       <- as.numeric(window(inflation.all, start=est.data.start))
hhc                      <- as.numeric(window(hhc.all, start=est.data.start))
pvc                      <- as.numeric(window(pvc.all, start=est.data.start))
relative.import.price.inflation <- as.numeric(window(relative.import.price.inflation.all,
                                                     start=est.data.start))
nominal.interest.rate           <- as.numeric(window(nominal.interest.rate.all, start=est.data.start))
output.gap <- as.numeric(window(output.gap, start=est.data.start)) 
real.interest.rate              <- nominal.interest.rate - inflation


##------------------------------------------------------------------------------##
## Run estimation
##------------------------------------------------------------------------------##

## Running the stage 1 model
out.stage1 <- rstar.stage1_modified_hhc(log.output,
                                    output.gap,
                                    inflation,
                                    hhc,
                                    relative.import.price.inflation,
                                    b.y.constraint=b.y.constraint.stage1)


## Median unbiased estimate of lambda_g
lambda.g.series <- median.unbiased.estimator.stage1(out.stage1$potential.smoothed)
lambda.g <- (lambda.g.series[1] + lambda.g.series[2] + lambda.g.series[3]) / 3

## Running the stage 2 model
out.stage2 <- rstar.stage2_modified_hhc(log.output,
                           inflation,
                           relative.import.price.inflation,
                           real.interest.rate,
                           hhc,
                           output.gap,
                           lambda.g,
                           a.r.constraint=a.r.constraint.stage2,
                           b.y.constraint=b.y.constraint.stage2)


## Median unbiased estimate of lambda_z
lambda.z.series <- median.unbiased.estimator.stage2(out.stage2$y, out.stage2$x)
lambda.z <- (lambda.z.series[1] + lambda.z.series[2] + lambda.z.series[3]) / 3

#####################################################################################################
