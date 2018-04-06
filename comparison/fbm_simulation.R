##### fbm_simulation.R #####

library("fArma")
library("foreach")

# Comparison of the fBm simulations in the fArma package

method_to_description <- function(method) {
  if (method == "mvn") {
      return("Numerical approximation of the stochastic integral")
    } else if (method == "chol") {
      return("Choleski’s decomposition of the covariance matrix")
    } else if (method == "lev") {
      return ("Levinson method")
    } else if (method == "circ") {
      return("Wood and Chan method")
    } else if (method == "wave") {
      return("Wavelet synthesis")
    } else {
      stop("This method is not in the documentation")
    }
}

get_execution_time <- function(method, hurst, sample_length,
                               iterations = 100) {
  print(method_to_description(method))
  data <- foreach(i = 1:iterations, .combine = 'c') %do% {
    system.time(fbmSim(n = sample_length, H = hurst, method = method,
                       doplot = FALSE))[3]
  }

  return(paste0("overall time: ", sum(data)))
}

h <- 0.7
sample_length <- 2 ** 10
iterations <- 100

# mvn is way to slow
print(get_execution_time(method = "mvn", hurst = h,
                         sample_length = sample_length,
                         iterations = iterations))
# chol is also too slow
print(get_execution_time(method = "chol", hurst = h,
                        sample_length = sample_length,
                        iterations = iterations))
print(get_execution_time(method = "lev", hurst = h,
                        sample_length = sample_length,
                        iterations = iterations))
# circ is by far the fastest method
print(get_execution_time(method = "circ", hurst = h,
                         sample_length = sample_length,
                         iterations = iterations))
print(get_execution_time(method = "wave", hurst = h,
                        sample_length = sample_length,
                        iterations = iterations))

# results:

#"Numerical approximation of the stochastic integral"
#"overall time: 108.168"
#"Choleski’s decomposition of the covariance matrix"
#"overall time: 31.163"
#"Levinson method"
#"overall time: 6.77200000000002"
#"Wood and Chan method"
#"overall time: 0.164999999999907"
#"Wavelet synthesis"
#"overall time: 6.31800000000021"
