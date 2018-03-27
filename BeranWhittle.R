fspecFGN <- function(eta,m) {
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
  ##---------parameters for the calculation of f--------
  h <- eta[1]
  nsum <- 200
  #nsum = m
  hh <- -2*h-1
  const <- 1/pi*sin(pi*h)*gamma(-hh)
  j <- 2*pi*c(0:nsum)
  ##------   x = 2*pi*(j-1)/m (j=1,2,...,(n-1)/2)   -----
  ##------   Fourier frequencies   ----------------------
  mhalfm <- trunc((m-1)/2)
  x <- 1:mhalfm
  x <- 2*pi/m*x
  ##-----   calculation of f at Fourier frequencies   -------
  fspec <- matrix(0,mhalfm)
  for(i in seq(1:mhalfm)) {
    lambda <- x[i]
    fi <- matrix(lambda,(nsum+1))
    fi <- abs(j+fi)**hh+abs(j-fi)**hh
    fi[1] <- fi[1]/2
    fi <- (1-cos(lambda))*const*fi
    fspec[i] <- sum(fi)
  }
  ##----    adjusted spectrum (such that int(log(fspec))=0 ---
  logfspec <- log(fspec)
  fint <- 2/(m)*sum(logfspec)
  theta1 <- exp(fint)
  fspec <- fspec/theta1
  drop(list(fspec=fspec,theta1=theta1))
}


###################
#
#Splus-functions and program for the calculation of Whittle's
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
#Covariance matrix of hat {eta}
# for fGn.
###################

CetaFGN <- function(eta)
{
M <- length(eta)
# size of steps in Riemann sum: 2*pi/m
m <- 10000
mhalfm <- trunc((m-1)/2)
# size of delta for numerical calculation of derivative
delta <- 0.000000001
# partial derivatives of log f (at each Fourier frequency)
lf<- matrix(1, ncol = M, nrow = mhalfm)
fO <- fspecFGN(eta,m)$fspec
for(j in (1:M))
{
etaj <- eta
etaj[j] <- etaj[j] + delta
fj <- fspecFGN(etaj,m)$fspec
lf[,j] <- log(fj/fO)/delta
}
# CalculateD
Djl <- matrix(1, ncol = M, nrow = M)
for(j in (1:M))
{for(l in (1:M))
{Djl[j,l]<-2*2*pi/m*sum(lf[,j]*lf[,l])
}
}
# Result
drop(matrix(4*pi*solve(Djl), ncol = M, nrow = M, byrow = T))
}
