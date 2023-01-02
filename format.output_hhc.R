##------------------------------------------------------------------------------##
## File:        format.output.R
##
## Description: This generates a dataframe to be written to a CSV containing
##              one- and two-sided estimates, parameter values, standard errors,
##              and other statistics of interest.
##------------------------------------------------------------------------------##
format.output_hhc <- function(estimation, one.sided.est, two.sided.est, real.rate, start, end, run.se = TRUE) {
    output <- data.frame(matrix(NA,dim(one.sided.est)[1],32))
    
    output[,1]   <- seq(from = (as.Date(ti(shiftQuarter(start,-1),'quarterly'))+1), to = (as.Date(ti(shiftQuarter(end,-1),tif='quarterly'))+1), by = 'quarter')
    output[,2:5] <- one.sided.est

    output[1,7]    <- "Parameter Point Estimates"
    output[2,7:24] <- c("a_1","a_2","a_3","a_4","b_1","b_2","b_3","b_4","r_1","r_2","r_3","sigma_y","sigma_pi","sigma_hhc","sigma_y*","e1","rho_1","a_1 + a_2")
    output[3,7:23] <- estimation$theta
    output[3,24]   <- estimation$theta[1] + estimation$theta[2]
    ## Include standard errors in output only if run.se switch is TRUE
    if (run.se) {
        output[4,7]    <- "T Statistics"
        output[5,7:23] <- estimation$se$t.stats
        
        output[8,7]    <- "Average Standard Errors"
        output[9,7:9]  <- c("y*","r*","g")
        output[10,7:9] <- estimation$se$se.mean        
    }
    
    output[14,7] <- "Signal-to-noise Ratios"
    output[15,7] <- "lambda_g"; output[15,8] <- lambda.g
    output[16,7] <- "lambda_z"; output[16,8] <- lambda.z
    output[14,11] <- "Log Likelihood"; output[15,11] <- estimation$log.likelihood

    output[19,7] <- "State vector: [y_{t}* y_{t-1}* y_{t-2}* g_{t-1} g_{t-2} z_{t-1} z_{t-2}]"
    output[20,7] <- "Initial State Vector"
    output[21,7:13] <- estimation$xi.00
    output[22,7] <- "Initial Covariance Matrix"
    output[23:29,7:13] <- estimation$P.00

    if (run.se) {
        output[,25]    <- seq(from = (as.Date(ti(shiftQuarter(start,-1),'quarterly'))+1), to = (as.Date(ti(shiftQuarter(end,-1),tif='quarterly'))+1), by = 'quarter')
        output[,26:27] <- estimation$se$se
    }
    
    output[,29] <- real.rate[9:length(real.rate)] - estimation$rstar.filtered

    output[,31]    <- seq(from = (as.Date(ti(shiftQuarter(start,-1),'quarterly'))+1), to = (as.Date(ti(shiftQuarter(end,-1),tif='quarterly'))+1), by = 'quarter')
    output[,32:35] <- two.sided.est

    return(output)
}
