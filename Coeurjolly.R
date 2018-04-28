library("dvfBm")

whittST <- function(fBm, Hprel) {
  if(missing(fBm)) fBm <- circFBM(n = 500, H = 0.7)
  Hprel <- 0.72
  # if(missing(Hprel)) {
  #   Db4 <- c(0.4829629, -0.8365163, 0.22414386, 0.12940952)
  #   Hprel <- VaPkolST(fBm=fBm, k=2, a=Db4, M=5, 0)$Hols
  # }
  fGn <- c(fBm[1], diff(fBm))
  n <- length(fGn)
  nstar <- trunc((n - 1)/2)
  perio <- (Mod(fft(fGn, inverse = F))) ** 2/(2 * pi * n)
  perio <- perio[1:nstar]
  global <- vector("list", 2)
  global[[1]] <- perio
  global[[2]] <- n
  # NOTE: command does not work in R
  assign("global", global, frame = 1)  ##  ## ---------------------
  ## criterion to minimize
  ## ---------------------
  Q <- function(Htry) {
    perio <- global[[1]]
    n <- global[[2]]
    spd <- spdFGN(Htry, n)
    result <- (4 * pi)/n * sum(perio/spd)
    drop(result)
  }
  epsilon <- 1.0000000000000001e-05
  Hest <- nlminb(objective = Q, start = Hprel, lower = epsilon, upper = 1 -
  epsilon)$parameters

  return(Hest)
}


spdFGN <- function(Htry, n) {
  if(missing(Htry)) Htry <- 0.7
  if(missing(n)) n <- 500
  alpha <- 2 * Htry + 1
  nstar <- trunc((n - 1)/2)
  clambda <- (sin(pi * Htry) * gamma(alpha))/pi ##
  ## -----------------------------------------
  ## computation of spd at Fourier frequencies
  ## -----------------------------------------
  j <- 2 * pi * ((-300):300)
  spd <- rep(0, nstar)
  for(k in (1:nstar)) {
    lambda <- (2 * pi * k)/n
    stocksum <- sum(abs(j + lambda) ** ( -alpha))
    spd[k] <- clambda * (1 - cos(lambda)) * stocksum
  }
  ## ----------------------
  ## renormalization of spd
  ## ----------------------
  theta <- exp(2/n * sum(log(spd)))
  spd <- spd/theta
  drop(spd)
}

# print(spdFGN())



VaPkolST <- function(fBm, k, a, M, llplot) {
     N <- length(fBm)   ##
## ------------------------------
## dilatation m times of filter a
## ------------------------------
     dilatation <- function(a, m)
     {
           la <- length(a)
           am <- rep(0, m * la - 1)
           am[seq(1, m * la - 1, by = m)] <- a
           drop(am)
     }
## ---------------------------------------------
## estimation of H by a simple linear regression
## ---------------------------------------------
     la <- length(a)
     Vam <- rep(0, M)
SNkam <- rep(0, M)
Vam <- filter(fBm, a, sides = 1)
Vam <- Vam[ - (1:(la - 1))]
SNkam[1] <- mean(abs(Vam)**k)
for(m in (2:M)) {
      am <- dilatation(a, m)
      Vam <- filter(fBm, am, sides = 1)
      lam <- m * la - 1
      Vam <- Vam[ - (1:(lam - 1))]
      SNkam[m] <- mean(abs(Vam)**k)
}
LN <- log(SNkam)
m <- 1:M
Reg <- lsfit(log(m), LN, intercept = T)
Hols <- Reg$coef[2]/k ##
mean.eps <- rep(0, M)
     for(m in (1:M)) {  ## approximation by expansion of psi(z)-log(z)
           z <- 0.5 * (N - m * (la - 1))
           mean.eps[m] <- -1/2/z - 1/12/z**2 + 1/120/z**4 - 1/252/z**6
     }
     m <- 1:M
     A <- log(m) - mean(log(m))
     norm.A <- as.vector(t(A) %*% A)
     bias <-  - t(A) %*% mean.eps/2/norm.A ##
     Ek <- (2**(k/2) * gamma(0.5 * k + 0.5))/gamma(0.5)
      pia0 <- piaH(a, Hols, 0)
      thetaols <- Reg$coef[1]
      Cols <- N**(Hols)/Ek**(1/k)/sqrt(pia0) * exp(thetaols/k) ##
      if(llplot == 1) {
            par(mfrow = c(1, 1))
            Hchar <- as.character(round(Hols, 4))
            Hleg <- paste(c("Hest", Hchar), collapse = "  =  ")
            plot(log(m), log(SNkam), main = "Regression of
                        log( SN(k,aˆm) ) on log(m)")
                        abline(Reg)
                        legend(min(log(m)), max(log(SNkam)), c("", Hleg, ""))
                      }
drop(list(Hols = Hols, Cols = Cols, LN = LN, bias = bias))
}

piaH <- function(a, H, i) {
  l <- length(a) - 1
  d <- l + 1
  mat <- matrix(rep(0, d * d), ncol = d)
  for(q in (0:l)) {
            for(r in (0:l)) {
                       z <- a[q + 1] * a[r + 1] * abs(q - r + i)**(2 * H)
                       mat[q + 1, r + 1] <- -0.5 * z
}
  }
  piaH.i <- sum(apply(mat, 1, sum))
  drop(piaH.i)
}

print(whittST())
