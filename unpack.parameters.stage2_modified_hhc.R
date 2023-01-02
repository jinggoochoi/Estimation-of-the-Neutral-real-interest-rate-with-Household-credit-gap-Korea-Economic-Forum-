##------------------------------------------------------------------------------##
## File:        unpack.parameters.stage2.R
##
## Description: This file generates coefficient matrices for the stage 2
##              state-space model for the given parameter vector.
##
## Stage 2 parameter vector: [a_1(ygap(-1)), a_2(ygap(-2)), a_3(r_t-1,2), a_4(hhc), a_5(1 column), a_6(g)
##                            b_1(pi(-1)), b_2(pi(-2~-4)), b_4(ygap(-1)), b_5(imp-inf),  r_1(hhc_t-1), r_2(hhc_t-2), r_3(r_t-1,2),
##                            sigma_y~, sigma_pi, sigma_hc, sigma_y*]  
##------------------------------------------------------------------------------##
unpack.parameters.stage2_modified_hhc <- function(parameters, y.data, x.data, lambda.g, xi.00=NA, P.00=NA) {
  n.state.vars <- 4
  
  A         <- matrix(0, 3, 11)
  A[1, 1:2] <- parameters[1:2] ## a_1, a_2
  A[1, 3:4] <- parameters[3]/2 ## a_3
  A[1, 8]   <- parameters[4]   ##a_4
  A[1, 11 ] <- parameters[5]   ## a_5
  A[2, 1  ] <- parameters[9]   ## b_4
  A[2, 5:6] <- parameters[7:8] ## b_1, b_2
  A[2, 7  ] <- 1 - parameters[7] - parameters[8] ## b_1,b_2
  A[2, 10] <- parameters[10] ##b_5
  A[3, 3:4] <- parameters[13]/2 ## r_3
  A[3, 8]  <- parameters[11] ## r_1  
  A[3, 9]  <- parameters[12] ## r_2  
  A         <- t(A)
  
  H         <- matrix(0, 3, 4)
  H[1, 1  ] <- 1
  H[1, 2:3] <- -parameters[1:2] ## a_1,a_2
  H[1, 4  ] <- parameters[6]    ## a_6
  H[2, 2]   <- -parameters[9]   ## b_4
  H         <- t(H)
  
  R         <- diag(c(parameters[14]^2, parameters[15]^2, parameters[16]^2)) ## sigma_1,sigma_2,sigma_6
  Q         <- matrix(0, 4, 4)
  Q[1, 1]   <- parameters[17]^2              ## sigma_4
  Q[4, 4]   <- (lambda.g * parameters[17])^2 ## sigma_4
  
  F <- matrix(0, 4, 4)
  F[1, 1] <- F[1, 4] <- F[2, 1] <- F[3, 2] <- F[4,4] <- 1
  
  cons <- matrix(0, n.state.vars, 1)  
  
  ## Starting values for xi.00 and P.00
  if (any(is.na(xi.00))) {
    x  <- rbind(t(H), t(H) %*% F, t(H) %*% F %*% F, t(H) %*% F %*% F %*% F, t(H) %*% F %*% F %*% F %*% F)
    om <- matrix(0, 15, 15)
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
    p1 <- t(x) %*% solve(om, x)
    yy <- c(y.data[1,], y.data[2,], y.data[3,], y.data[4,], y.data[5,])
    tmp <- c(t(A) %*% x.data[1,],
             t(A) %*% x.data[2,],
             t(A) %*% x.data[3,],
             t(A) %*% x.data[4,],
             t(A) %*% x.data[5,])
    xi.00 <- solve(p1, t(x)) %*% solve(om, yy - tmp)
    tmp <- yy - tmp - x %*% xi.00
    P.00 <- solve(p1, (diag(nrow=4) * sum(tmp^2) / (length(yy) - n.state.vars)))
  }
  return(list("xi.00"=xi.00, "P.00"=P.00, "F"=F, "Q"=Q, "A"=A, "H"=H, "R"=R, "cons"=cons, "x.data"=x.data, "y.data"=y.data))
}
