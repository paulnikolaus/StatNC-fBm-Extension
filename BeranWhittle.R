## Purpose: Calculation of the spectral density f of
## normalized fractional Gaussian noise with self-similarity parameter H=h
## at the Fourier frequencies 2*pi*j/m (j=1,...,(m-1)).
##
## Remarks:
## -------
## 1. cov(X(t),X(t+k)) = integral[ exp(iuk)f(u)du ]
## 2. f=theta1*fspec and integral[log(fspec)]=0.
## -------------------------------------------------------------------------
## INPUT: m = sample size
##        h = self-similarity parameter
##
## OUTPUT: list(fspec=fspec,theta1=theta1)
## -------------------------------------------------------------------------
## Author: Jan Beran;  modified: Martin Maechler, Date: Sep 95.

fspecFGN <- function(eta, m) {
  # ---------parameters for the calculation of f--------
  h <- eta[1]
  nsum <- 200
  # nsum = m
  hh <- -2 * h - 1
  const <- 1 / pi * sin(pi * h) * gamma(-hh)
  j <- 2 * pi * c(0:nsum)
  #------   x = 2*pi*(j-1)/m (j=1,2,...,(n-1)/2)   -----
  #------   Fourier frequencies   ----------------------
  mhalfm <- trunc((m - 1) / 2)
  x <- 1:mhalfm
  x <- 2 * pi / m * x
  #-----   calculation of f at Fourier frequencies   -------
  fspec <- matrix(0, mhalfm)
  for (i in seq(1:mhalfm)) {
    lambda <- x[i]
    fi <- matrix(lambda, (nsum + 1))
    fi <- abs(j + fi) ** hh + abs(j - fi) ** hh
    fi[1] <- fi[1] / 2
    fi <- (1 - cos(lambda)) * const * fi
    fspec[i] <- sum(fi)
  }
  #----    adjusted spectrum (such that int(log(fspec))=0 ---
  logfspec <- log(fspec)
  fint <- 2 / (m) * sum(logfspec)
  theta1 <- exp(fint)
  fspec <- fspec / theta1
  drop(list(fspec = fspec, theta1 = theta1))
}


###################
#
# Splus-functions and program for the calculation of Whittle's
# estimator and the goodness of fit statistic defined in Beran
# (1992). The models are fractional Gaussian noise or
# fractional ARIMA. The data series may be divided
# into subseries for which the parameters are fitted separately.
#
###################
#FFFFFFFFFFFFFFFFF
# Functions
#FFFFFFFFFFFFFFFFF
###################
# CetaFGN
#
# Covariance matrix of hat {eta}
# for fGn.
###################

CetaFGN <- function(eta) {
  M <- length(eta)
  # size of steps in Riemann sum: 2*pi/m
  m <- 10000
  # trunc(): integer part of (m-1)/2
  mhalfm <- trunc((m - 1) / 2)
  # size of delta for numerical calculation of derivative
  delta <- 0.000000001
  # partial derivatives of log f (at each Fourier frequency)
  lf <- matrix(1, ncol = M, nrow = mhalfm)
  f_O <- fspecFGN(eta, m)$fspec
  for (j in (1:M)) {
    etaj <- eta
    etaj[j] <- etaj[j] + delta
    fj <- fspecFGN(etaj, m)$fspec
    lf[, j] <- log(fj / f_O) / delta
  }
  # Calculate D
  Djl <- matrix(1, ncol = M, nrow = M)
  for (j in (1:M)) {
    for (l in (1:M)) {
      Djl[j, l] <- 2 * 2 * pi / m * sum(lf[, j] * lf[, l])
    }
  }
  # Result
  drop(matrix(4 * pi * solve(Djl), ncol = M, nrow = M, byrow = T))
}

###################
# Qeta
# -------------------------------
# Function for the calculation of A, B and
# Tn = A/B ** 2
# where A = 2pi / n sum 2 * [I(lambda <- j) / f(lambda <- j)],
# B = 2pi/n sum 2*[I(lambda) <- j]/f(lambda <- j)]**2 and
# the sum is taken over all Fourier frequencies
# lambda <- j = 2pi*j/n (j=1,...,(n-1)/2)
# f is the spectral density of fractional Gaussion noise or fractional
# ARIMA(p, d, q) with self similarity parameter H=h.
# cov(X(t), X(t+k)) = integral(exp(iuk)f(u)du)
#
# NOTE: yper[1] ,ist be the periodogram I(lambda<-1) at
# the frequency 2pi/n (i.e. not the frequency zero!)
# INPUT: h
# (n, nhalfm = trunc[(n - 1) / 2]) and
# nhalfm-dimensional vector yper must be defined.)
# Tn is the goodness of fit test statistic
# Tn = A / B ** 2 defined in Beran (1992),
# z is the standardized test statistic,
# pval the corresponding p-value P(w>z).
# theta1 is the scale parameter such that
# f = theta1 * fspec and integral(log[fspec]) = 0.

Qeta <- function(eta) {
  cat("in function Qeta", fill = T)
  h <- eta[1]
  # spectrum at Fourier frequencies

  fspec <- fspecFGN(eta, n)
  theta1 <- fspec$theta1
  fspec <- fspec$fspec

  # if (imodel == 1) {
  #   fspec <- fspecFGN(eta, n)
  #   theta1 <- fspec$theta1
  #   fspec <- fspec$fspec
  # } else {
  #   fspec <- fspecARIMA(eta, p, q, n)
  #   theta1 <- fspec$theta1
  #   fspec <- fspec$fspec
  # }

  # Tn = A / (B ** 2)
  yf <- yper / fspec
  yfyf <- yf ** 2
  A <- 2 * (2 * pi / n) * sum(yfyf)
  B <- 2 * (2 * pi / n) * sum(yf)
  Tn <- A / (B ** 2)
  z <- sqrt(n) * (pi * Tn - 1) / sqrt(2)
  pval <- 1 - pnorm(z)
  theta1 <- B / (2 * pi)
  fspec <- fspec

  Qresult <- list(n = n, h = h, eta = eta, A = A, B = B, Tn = Tn, z = z,
    pval = pval, theta1 = theta1, fspec = fspec)
drop(Qresult)
}

###################
# definition of the periodogram
#################
per <- function(z) {
  n <- length(z)
  (Mod(fft(z)) ** 2 / (2 * pi * n)) [1:(n %/% 2 + 1)]
}


###################
# definition of function to be minimized
###################
Qmin <- function (etatry) {
  result <- Qeta(etatry)$B
  cat("etatry=", etatry, "B=", result, sep = " " , fill = T)
  drop(result)
}

# #MMMMMMMMMMMMM # Main program #MMMMMMMMMMMMM
# # read data
# cat(" in which file are the data ?")
# filedata <- readline()
# cat(" total number of observations ?")
# nmax <- scan(n = 1)
# cat("first and last observation to be considered (istart, iend) ?")
#
# startend <- c(scan(n = 2)) # > we only look at
# istart <- startend[1] # observations
# iend <- startend[2] # istart,istart+1,...,end
# cat("into how many subseries do you divide the data ?")
# nloop <- scan(n = 1)
# n <- trunc((iend - istart + 1) / nloop)
# nhalfm <- trunc((n - 1) / 2)
# # choose model
# cat("model: fr.G.noise (1) or fractional ARIMA(p,d,q) (2)?")
#
#
# p <- 0
# q <- 0
#
# # imodel <- scan(n = 1)
# # p <- 0
# # q <- 0
# # if (imodel == 2) {
# #   cat(" order of AR ?") #
# #   p <- scan(n = 1) #
# #   cat("order of MA ?")
# #   q <- scan(n = 1)
# # }
#
# # initialize h
# cat(" initial estimate of h=?")
# h <- scan(n = 1)
# eta <- c(h)
# # initialize AR parameter
# if (p > 0) {
#   cat(" initial estimates of AR parameters=?")
#   eta[2:(p + 1)] <- scan(n = -p)
# }
# # initialize MA parameter
# if (q > 0) {
#   cat(" initial estimates of MA parameters=?")
#   eta[(p + 2):(p + q + 1)] <- scan(n = q)
# }
# M <- length(eta)
# # loop
# thetavector <- c()
# i0 <- istart
# for (iloop in (1:nloop)) {
#   h <- max(0.2, min(h, 0.9)) # avoid extreme initial values
#   eta[1] <- h
#   i1 <- i0 + n - 1
#   y <- c(scan(filedata, n = nmax))[i0:i1] # read only y[i0:i1]
#   # standardize data
#   vary <- var(y)
#   y <- (y - mean(y)) / sqrt(var(y))
#   # periodogram of data
#   yper <- per(y)[2:(nhalfm + 1)]
#   # find estimate
#   s <- 2 * (1.-h)
#   etatry <- eta
#   # nlmin not part of R. Use nlm() or nlminb() instead
#   result <- nlmin(Qmin, etatry, xc.to1 = 0.0000001, init.step = s)
#   eta <- result$x
#   thetal <- Qeta(eta)$thetal
#   theta <- c(thetal, eta)
#   thetavector <- c(thetavector, theta)
#   # calculate goodness of fit statistic
#   Qresult <- Qeta(eta)
#   # output
#   M <- length(eta)
#
#   SD <- CetaFGN(eta)
#   SD <- matrix(SD, nco1 = M, nrow = M, byrow = T) / n
#
#   # if (imodel == 1) {
#   #   SD <- CetaFGN(eta)
#   #   SD <- matrix(SD, ncol = M, nrow = M, byrow = T) / n
#   # } else {
#   #   SD <- CetaARIMA(eta, p, q)
#   #   SD <- matrix(SD, ncol = M, nrow, M, byrow = T) / n
#   # }
#   Hlow <- eta[1] - 1.96 * sqrt(SD[1, 1])
#   Hup <- eta[1] + 1.96 * sqrt(SD[1, 1])
#   cat("theta=", theta, fill = T)
#   cat("H=", eta[1], fill = T)
#   cat("95%-C.I. for H: [", Hlow, ",", Hup, "]", fill = T)
#   etalow <- c()
#   etaup <- c()
#   for (i in (1:length(eta))){
#     etalow <- c(etalow, eta[i] - 1.96 * sqrt(SD[i, i]))
#     etaup <- c(etaup, eta[i] + 1.96 * sqrt(SD[i, i]))
#   }
#   cat("95%-C.I.:", fill = T)
#   print(cbind(etalow, etaup), fill = T)
#   cat("periodogram is in yper", fill = T)
#   fest <- QresultStheta1 * Qresult$fspec
#   cat(" spectral density is in fest", fill = T)
#   # next subseries
#   i0 <- i0 + n
# }
