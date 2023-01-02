##------------------------------------------------------------------------------##
## File:        unpack.parameters.stage1.R
##
## Description: This file generates coefficient matrices for the stage 1
##              state-space model for the given parameter vector.
##
## Stage 1 parameter vector: [a_1(ygap(-1)), a_2(ygap(-2)), a_3(hhc), b_1(pi(-1)), b_2(pi(-2~-4)), 
##                            b_4(output gap), b_5(imp-inf), r_1(hhc(-1)), r_2(hhc(-2)), g, sigma_y, sigma_pi, sigma_hhc, sigma_4]  
##------------------------------------------------------------------------------##
unpack.parameters.stage1_modified_hhc <- function(parameters, y.data, x.data, xi.00=NA, P.00=NA) {
  n.state.vars <- 3
    
  A         <- matrix(0, 8, 3)
  A[1:3, 1] <- parameters[1:3]
  A[1, 2]   <- parameters[6]
  A[5:6, 2] <- parameters[4:5]
  A[7, 2]   <- 1-sum(A[4:5, 2])
  A[8, 2] <- parameters[7]
  A[3, 3] <- parameters[8]
  A[4, 3] <- parameters[9]
  
  H         <- matrix(0, 3, 3)
  H[1, 1]   <- 1
  H[2:3, 1] <- -parameters[1:2]
  H[2, 2]   <- -parameters[6]

  R         <- diag(c(parameters[11]^2, parameters[12]^2, parameters[13]^2))
  Q         <- matrix(0, 3, 3)
  Q[1, 1]   <- parameters[14]^2

  F <- matrix(0, 3, 3)
  F[1, 1] <- F[2, 1] <- F[3, 2] <- 1

  cons <- matrix(0, 3, 1)
  cons[1, 1] <- parameters[10]
  
  if (any(is.na(xi.00))) {
      xi.00 <- rep(0, n.state.vars)
      
      ##  Starting values for xi.00 and P.00
      x  <- rbind(t(H), t(H) %*% F, t(H) %*% F %*% F, t(H) %*% F %*% F %*% F)
      om <- matrix(0, 12, 12)
      om[7:9, 4:6] <- t(H) %*% F %*% Q %*% H
      om[10:12, 4:6] <- t(H) %*% F %*% F %*% Q %*% H
      om[10:12, 7:9] <- t(H) %*% (F %*% F %*% Q %*% t(F) + F %*% Q) %*% H
      om           <- om + t(om)
      om[1:3, 1:3] <- R
      om[4:6, 4:6] <- t(H) %*% Q %*% H + R
      om[7:9, 7:9] <- t(H) %*% (F %*% Q %*% t(F) + Q) %*% H + R
      om[10:12, 10:12] <- t(H) %*% (F %*% F %*% Q %*% t(F) %*% t(F) + F %*% Q %*% t(F) + Q) %*% H + R
      
      p1 <- t(x) %*% solve(om, x) ## x' * inv(om) * x
      yy <- c(y.data[1,], y.data[2,], y.data[3,], y.data[4,])
      tmp <- c(t(A) %*% x.data[1,],
               (t(A) %*% x.data[2,] + t(H) %*% cons),
               (t(A) %*% x.data[3,] + t(H) %*% cons + t(H) %*% F %*% cons),
               (t(A) %*% x.data[4,] + t(H) %*% (diag(n.state.vars)+F+F%*%F) %*% cons))
      xi.00 <- solve(p1, t(x)) %*% solve(om, yy - tmp)
      tmp <- yy - tmp - x %*% xi.00
      P.00 <- solve(p1, (diag(nrow=3) * sum(tmp^2) / 3))
  }
  return(list("xi.00"=xi.00, "P.00"=P.00, "F"=F, "Q"=Q, "A"=A, "H"=H, "R"=R, "cons"=cons, "x.data"=x.data, "y.data"=y.data))
}
