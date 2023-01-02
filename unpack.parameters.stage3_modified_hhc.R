##------------------------------------------------------------------------------##
## File:        unpack.parameters.stage3.R
##
## Description: This file generates coefficient matrices for the stage 3
##              state-space model for the given parameter vector.
##
## Stage 3 parameter vector: [a_1(ygap(-1)), a_2(ygap(-2)), a_3(r_t-1,2), a_4(hhc), 
##                            b_1(pi(-1)), b_2(pi(-2~-4)), b_3(ygap(-1)), b_5(imp-inf),  
##                            r_1(hhc_t-1), r_2(hhc_t-2), r_3(r_t-1,2),
##                            sigma_y~, sigma_pi, sigma_hc, sigma_y*, e_1, rho_z] 
##------------------------------------------------------------------------------##
unpack.parameters.stage3_modified_hhc <- function(parameters, y.data, x.data, lambda.g, lambda.z, xi.00=NA, P.00=NA, xi.00.gpot=NA, P.00.gpot=NA) {
  n.state.vars <- 7
    
  A         <- matrix(0, 3, 10)
  A[1, 1:2] <- parameters[1:2] ## a_1, a_2
  A[1, 3:4] <- parameters[3]/2 ## a_3
  A[1, 9  ] <- parameters[4]   ## a_4
  A[2, 1  ] <- parameters[7]   ## b_3
  A[2, 5:6] <- parameters[5:6] ## b_1, b_2
  A[2, 7  ] <- 1 - parameters[5] - parameters[6] ## b_1,b_2
  A[2, 8] <- parameters[8] ##b_5
  A[3, 9:10] <- parameters[9:10]   ## r_1,r_2
  A[3, 3:4] <- parameters[11]/2   ## r_3
  A         <- t(A)
  
  H         <- matrix(0, 3, 7)
  H[1, 1  ] <- 1
  H[1, 2:3] <- -parameters[1:2]                  ## a_1,a_2
  H[1, 4:5] <- -2*parameters[3]*parameters[16]   ## a_3 (annualized)
  H[1, 6:7] <- -parameters[3]/2*parameters[17]   ## a_3*rho_z
  H[2, 2]   <- -parameters[7]                    ## b_3
  H[3, 4:5] <- -2*parameters[11]*parameters[16]   ## r_3 (annualized)
  H[3, 6:7] <- -parameters[11]/2*parameters[17]   ## r_3*rho_z
  H         <- t(H)

  R         <- diag(c(parameters[12]^2, parameters[13]^2, parameters[14]^2)) ## sigma_y,sigma_pi,sigma_hc

  Q         <- matrix(0, 7, 7)
  Q[1, 1]   <- (1+lambda.g^2)*parameters[15]^2 ## sigma_4
  Q[1, 4]   <- Q[4, 1] <- Q[4, 4] <- (lambda.g*parameters[15])^2 ## sigma_4
  Q[6, 6]   <- (lambda.z*parameters[12]/parameters[3])^2 ##sigma_1, a_3
  
  F <- matrix(0, 7, 7)
  F[1, 1] <- F[1, 4] <- F[2, 1] <- F[3, 2]  <- F[5,4]<- 1
  F[4,4] <- parameters[16]
  F[6,6] <-  parameters[17]
  F[7,6] <- 1
  cons <- matrix(0, n.state.vars, 1)
  
  ##  Starting values for xi.00 and P.00
  if (any(is.na(xi.00))) {
      x <- rbind(t(H), t(H)%*%F, t(H)%*%F%*%F, t(H)%*%F%*%F%*%F, t(H)%*%F%*%F%*%F%*%F)
      om <- matrix(0,15,15)
      
      om[7:9, 4:6] <- t(H) %*% F %*% Q %*% H
      om[10:12, 4:6] <- t(H) %*% F %*% F %*% Q %*% H
      om[13:15,4:6] <- t(H) %*% F %*% F %*% F %*% Q %*% H
      om[10:12, 7:9] <- t(H) %*% (F %*% F %*% Q %*% t(F) + F %*% Q) %*% H
      om[13:15,7:9] <- t(H) %*% F %*% (F %*% F %*% Q %*% t(F) + F %*% Q) %*% H
      om[13:15,10:12] <- t(H) %*% F %*% (F %*% F %*% Q %*% t(F) %*% t(F) + F %*% Q %*% t(F) + Q) %*% H
      om           <- om + t(om)
      om[1:3, 1:3] <- R
      om[4:6, 4:6] <- t(H) %*% Q %*% H + R
      om[7:9, 7:9] <- t(H) %*% (F %*% Q %*% t(F) + Q) %*% H + R
      om[10:12, 10:12] <- t(H) %*% (F %*% F %*% Q %*% t(F) %*% t(F) + F %*% Q %*% t(F) + Q) %*% H + R
      om[13:15, 13:15] <- t(H) %*% (F %*% F %*% F %*% Q %*% t(F) %*% t(F) %*% t(F) + F %*% F %*% Q %*% t(F) %*% t(F)
                                    + F %*% Q %*% t(F) + Q) %*% H + R
      if ( (det(om) <= 10^(-3)) | (det(t(x) %*% solve(om, x)) <= 10^(-3)) ) {
          gls   <- FALSE
          xi.00 <- xi.00.gpot
          P.00  <- P.00.gpot
      }
      else {
          gls <- TRUE
          p1 <- t(x) %*% solve(om, x)
          yy <- c(y.data[1,], y.data[2,], y.data[3,], y.data[4,], y.data[5,])
          tmp <- c(t(A) %*% x.data[1,],
                   t(A) %*% x.data[2,],
                   t(A) %*% x.data[3,],
                   t(A) %*% x.data[4,],
                   t(A) %*% x.data[5,])
          xi.00 <- solve(p1, t(x)) %*% solve(om, yy - tmp)
          tmp <- yy - tmp - x %*% xi.00
          P.00 <- solve(p1, diag(nrow=4) * (sum(tmp^2) / (length(yy) - n.state.vars) ))
      }
  }
  else { gls <- FALSE }
  return(list("xi.00"=xi.00, "P.00"=P.00, "F"=F, "Q"=Q, "A"=A, "H"=H, "R"=R, "cons"=cons, "x.data"=x.data, "y.data"=y.data,"gls"=gls))
}



