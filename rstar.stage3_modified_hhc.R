##------------------------------------------------------------------------------##
## File:        rstar.stage3.R
##
## Description: This file runs the model in the third stage of the LW estimation.
##------------------------------------------------------------------------------##
rstar.stage3_modified_hhc <- function(log.output,
                                      inflation,
                                      output.gap,
                                      relative.import.price.inflation,
                                      real.interest.rate,
                                      hhc,
                                      lambda.g,
                                      lambda.z,
                                      a.r.constraint=NA,
                                      b.y.constraint=NA,
                                      c.constraint.lb=NA,
                                      c.constraint.ub=NA,
                                      rho.z.constraint.lb=NA,
                                      rho.z.constraint.ub=NA,
                                      run.se=FALSE,
                                      xi.00=NA, P.00=NA) {
  
  stage <- 3
  
  ## Data must start 8 quarters before estimation period
  
  T<-length(log.output)-8
  ## Original output gap estimate
  #x.og <- cbind(rep(1,T+4), 1:(T+4), c(rep(0,56),1:(T+4-56)), c(rep(0,142),1:(T+4-142)))
  #y.og <- log.output[5:(T+8)]
  #output.gap <- (y.og - x.og %*% solve(t(x.og) %*% x.og, t(x.og) %*% y.og)) * 100
  
  
  ## Initialization of state vector for Kalman filter using HP trend of log output
  ## Pulled into unpack.parameters.stage3.R
  ## exported from HLW(2016)
  log.output.hp.trend <- hpfilter(log.output,freq=1600,type="lambda",drift=FALSE)$trend
  g.pot <- log.output.hp.trend[(g.pot.start.index):length(log.output.hp.trend)]
  g.pot.diff <- diff(g.pot)
  z <- real.interest.rate[-1:-6]-400*g.pot.diff
  z.t <- sum(z[1:20])/20
  g.t <- sum(g.pot[1:20])/20
  g.diff <- sum(g.pot.diff[1:20])/20
  #xi.00.gpot <- c(100*g.pot[3:1],100*g.pot.diff[2:1],z.t,z.t)
  xi.00.gpot <- c(g.t,g.t,g.t,g.diff,g.diff,z.t,z.t)
  ## xi.00.gpot: Initial vector for state variables(y(t)*,y(t-1)*,y(t-2)*,g(t),g(t-1),z(t),z(t-1))
  
  #b.pot <- solve(t(x.og) %*% x.og, t(x.og) %*% y.og)
  #g.pot <- x.og %*% solve(t(x.og) %*% x.og, t(x.og) %*% y.og)
  #xi.00.gpot <- c(100*g.pot[5:3],100*b.pot[2],100*b.pot[2],0,0)
  
  ## IS curve
  y.is <- output.gap[9:(T+8)]
  x.is <- cbind(output.gap[8:(T+7)], output.gap[7:(T+6)],
                (real.interest.rate[9:(T+8)] + real.interest.rate[8:(T+7)])/2,
                hhc[9:(T+8)])
  b.is <- solve(t(x.is) %*% x.is, t(x.is) %*% y.is)
  r.is <- y.is - x.is %*% b.is
  s.is <- sqrt(sum(r.is^2) / (T-dim(x.is)[2]))
  
  ## Phillips curve
  y.ph <- inflation[9:(T+8)]
  x.ph <- cbind(inflation[8:(T+7)],
                (inflation[7:(T+6)]+inflation[6:(T+5)]+inflation[5:(T+4)])/3,
                (inflation[4:(T+3)]+inflation[3:(T+2)]+inflation[2:(T+1)]+inflation[1:T])/4,
                output.gap[8:(T+7)],
                relative.import.price.inflation[9:(T+8)])
  b.ph <- solve(t(x.ph) %*% x.ph, t(x.ph) %*% y.ph)
  r.ph <- y.ph - x.ph %*% b.ph
  s.ph <- sqrt(sum(r.ph^2) / (T-dim(x.ph)[2]))
  
  ## hhc curve
  y.hc <- hhc[9:(T+8)]
  x.hc <- cbind(hhc[8:(T+7)],
                hhc[7:(T+6)],
                (real.interest.rate[8:(T+7)] + real.interest.rate[7:(T+6)])/2)
  b.hc <- solve(t(x.hc) %*% x.hc, t(x.hc) %*% y.hc)
  r.hc <- y.hc - x.hc %*% b.hc
  s.hc <- sqrt(sum(r.hc^2) / (T-dim(x.hc)[2]))
  
  y.data <- cbind(100 * log.output[9:(T+8)],
                  inflation[9:(T+8)],
                  hhc[9:(T+8)])
  x.data <- cbind(100 * log.output[8:(T+7)],
                  100 * log.output[7:(T+6)],
                  real.interest.rate[9:(T+8)],
                  real.interest.rate[8:(T+7)],                  
                  inflation[8:(T+7)],
                  (inflation[7:(T+6)]+inflation[6:(T+5)]+inflation[5:(T+4)])/3,
                  (inflation[4:(T+3)]+inflation[3:(T+2)]+inflation[2:(T+1)]+inflation[1:T])/4,
                  hhc[8:(T+7)],
                  hhc[7:(T+6)],
                  relative.import.price.inflation[9:(T+8)])
  
  ###################a_1,a_2,a_3,a_4,b_1,b_2,b_3,b_4,r_1,r_2,r_3,sigma_y,sigma_pi,sigma_hhc,sigma_y*,e_1,rho_z###############
  initial.parameters <- c(b.is, b.ph[1:2], b.ph[4:5], b.hc, s.is, s.ph, s.hc, 0.5, 0.9, 1)
  ###################a_1,a_2,a_3,a_4,b_1,b_2,b_3,b_4,r_1,r_2,r_3,sigma_y,sigma_pi,sigma_hhc,sigma_y*,e_1,rho_z###############
  
  theta.lb <- c(rep(-Inf,length(initial.parameters)))
  theta.ub <- c(rep(Inf,length(initial.parameters)))
  
  ## Set a lower bound for the Phillips curve slope (b_3) of b.y.constraint, if not NA
  if (!is.na(b.y.constraint)) {
    print(paste("Setting a lower bound of b_y >",as.character(b.y.constraint),"in Stage 3"))
    if (initial.parameters[7] < b.y.constraint) {
      initial.parameters[7] <- b.y.constraint
    }
    theta.lb[7] <- b.y.constraint
  }
  
  ## Set an upper bound for the IS curve slope (a_3) of a.r.constraint, if not NA
  if (!is.na(a.r.constraint)) {
    print(paste("Setting an upper bound of a_r <",as.character(a.r.constraint),"in Stage 3"))
    if (initial.parameters[3] > a.r.constraint) {
      initial.parameters[3] <- a.r.constraint
    }
    theta.ub[3] <- a.r.constraint      
  }
  
  if (!is.na(c.constraint.ub)) {
    print(paste("Setting an upper bound of c <=",as.character(c.constraint.ub),"in Stage 3"))
    if (initial.parameters[16] >= c.constraint.ub) {
      initial.parameters[16] <- c.constraint.ub
    }
    theta.ub[16] <- c.constraint.ub      
  }
  
  if (!is.na(c.constraint.lb)) {
    print(paste("Setting an lower bound of c >=",as.character(c.constraint.lb),"in Stage 3"))
    if (initial.parameters[16] <= c.constraint.lb) {
      initial.parameters[16] <- c.constraint.lb
    }
    theta.lb[16] <- c.constraint.lb      
  }
  
  if (!is.na(rho.z.constraint.ub)) {
    print(paste("Setting an upper bound of rho.z <=",as.character(rho.z.constraint.ub),"in Stage 3"))
    if (initial.parameters[17] >= rho.z.constraint.ub) {
      initial.parameters[17] <- rho.z.constraint.ub
    }
    theta.ub[17] <- rho.z.constraint.ub      
  }
  
  if (!is.na(rho.z.constraint.lb)) {
    print(paste("Setting an lower bound of rho.z >=",as.character(rho.z.constraint.lb),"in Stage 3"))
    if (initial.parameters[17] <= rho.z.constraint.lb) {
      initial.parameters[17] <- rho.z.constraint.lb
    }
    theta.lb[17] <- rho.z.constraint.lb      
  }
  
  
  P.00.gpot <- calculate.covariance_hhc(initial.parameters, theta.lb, theta.ub, y.data, x.data, stage, lambda.g, lambda.z, xi.00=xi.00.gpot)
  
  f <- function(theta) {return(-log.likelihood.wrapper_hhc(theta, y.data, x.data, stage,
                                                           lambda.g, lambda.z,
                                                           xi.00=NA, P.00=NA,
                                                           xi.00.gpot, P.00.gpot)$ll.cum)}
  nloptr.out <- nloptr(initial.parameters, f, eval_grad_f=function(x) {gradient(f, x)},
                       lb=theta.lb,ub=theta.ub,
                       opts=list("algorithm"="NLOPT_LD_MMA","xtol_rel"=1.0e-8,"maxeval"=5000))
  theta <- nloptr.out$solution
  
  if (nloptr.out$status==-1 | nloptr.out$status==5) {
    print("Look at the termination conditions for nloptr in Stage 3")
    stop(nloptr.out$message)
  }
  
  ## Get xi.00 and P.00
  init.vals <- unpack.parameters.stage3_modified_hhc(theta, y.data, x.data, lambda.g, lambda.z,
                                                     xi.00=NA, P.00=NA, xi.00.gpot, P.00.gpot)
  xi.00 <- init.vals$xi.00
  P.00  <- init.vals$P.00
  print(paste("GLS switch value:",init.vals$gls))
  
  log.likelihood <- log.likelihood.wrapper_hhc(theta, y.data, x.data, stage, lambda.g, lambda.z, xi.00, P.00)$ll.cum
  states <- kalman.states.wrapper_hhc(theta, y.data, x.data, stage, lambda.g, lambda.z, xi.00, P.00)
  
  ## If run.se = TRUE, compute standard errors for estimates of the states and report run time
  if (run.se) {
    ptm <- proc.time()
    se <- kalman.standard.errors_hhc(T, states, theta, y.data, x.data, stage, lambda.g, lambda.z, xi.00, P.00, niter, a.r.constraint, b.y.constraint)
    print("Standard error procedure run time")
    print(proc.time() - ptm)
  }
  
  ## One-sided (filtered) estimates
  trend.filtered      <- states$filtered$xi.tt[,4] * 4
  z.filtered          <- states$filtered$xi.tt[,6]
  rstar.filtered      <- trend.filtered*theta[16] + z.filtered*theta[17]
  potential.filtered  <- states$filtered$xi.tt[,1]/100
  output.gap.filtered <- y.data[,1] - (potential.filtered * 100)
  
  ## Two-sided (smoothed) estimates
  trend.smoothed      <- states$smoothed$xi.tT[,4] * 4
  z.smoothed          <- states$smoothed$xi.tT[,6]
  rstar.smoothed      <- trend.smoothed*theta[16] + z.smoothed*theta[17]
  potential.smoothed  <- states$smoothed$xi.tT[,1]/100
  output.gap.smoothed <- y.data[,1] - (potential.smoothed * 100)
  
  ## Save variables to return
  return.list <- list()
  return.list$rstar.filtered      <- rstar.filtered
  return.list$trend.filtered      <- trend.filtered
  return.list$z.filtered          <- z.filtered
  return.list$potential.filtered  <- potential.filtered
  return.list$output.gap.filtered <- output.gap.filtered
  return.list$rstar.smoothed      <- rstar.smoothed
  return.list$trend.smoothed      <- trend.smoothed
  return.list$z.smoothed          <- z.smoothed
  return.list$potential.smoothed  <- potential.smoothed
  return.list$output.gap.smoothed <- output.gap.smoothed
  return.list$theta               <- theta
  return.list$log.likelihood      <- log.likelihood
  return.list$states              <- states
  return.list$xi.00               <- xi.00
  return.list$P.00                <- P.00
  return.list$y.data              <- y.data
  return.list$initial.parameters  <- initial.parameters
  return.list$init.vals           <- init.vals
  if (run.se) { return.list$se    <- se }
  return(return.list)
}
