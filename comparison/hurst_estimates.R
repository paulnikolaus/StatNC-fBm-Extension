##### hurst_estimates.R #####

library("fArma")
library("foreach")

source("../simulation.R")

# Comparison of the hurst estimations in the fArma package
# Functions are based on
# Taqqu, Teverovsky, and Willinger,
# Estimators for Long-Range Dependence: An Empirical Study,
# 1995

flow_example <- build_flow(arrival_rate = 1.0, hurst = 0.7,
                           sample_length = 2 ** 12, std_dev = 1.0)

fgn_sample <- (flow_example - 1.0) / 1.0

estimate_random_data <- function(estimator_fun) {
  flow_example <- build_flow(arrival_rate = 1.0, hurst = 0.7,
                             sample_length = 2 ** 12, std_dev = 1.0)

  fgn_sample <- (flow_example - 1.0) / 1.0

  return(estimator_fun(x = fgn_sample)@hurst$"H")
}

get_execution_time <- function(fgn_data, estimator_fun, iterations = 100) {
  print(estimator_fun(x = fgn_data)@method)
  data <- foreach(i = 1:iterations, .combine = "c") %do% {
    system.time(estimate_random_data(estimator_fun))[3]
  }

  return(paste0("overall time: ", sum(data)))
}

print(get_execution_time(fgn_data = fgn_sample, estimator_fun = aggvarFit))
print(get_execution_time(fgn_data = fgn_sample, estimator_fun = diffvarFit))
print(get_execution_time(fgn_data = fgn_sample, estimator_fun = absvalFit))
# higuchiFit is quite slow
# print(get_execution_time(fgn_data = fgn_sample, estimator_fun = higuchiFit))
# pengFit is very slow
# print(get_execution_time(fgn_data = fgn_sample, estimator_fun = pengFit))
# rsFit is really fast
# print(get_execution_time(fgn_data = fgn_sample, estimator_fun = rsFit))
print(get_execution_time(fgn_data = fgn_sample, estimator_fun = perFit))
# print(get_execution_time(fgn_data = fgn_sample, estimator_fun = boxperFit))

mean_and_var_estimate <- function(fgn_data, estimator_fun, iterations = 100) {
  print(estimator_fun(x = fgn_data)@method)
  data <- foreach(i = 1:iterations, .combine = "c") %do% {
    estimate_random_data(estimator_fun)
  }

  return(paste0("mean: ", mean(data), " , standard deviation: ", sd(data)))
}

print(mean_and_var_estimate(fgn_data = fgn_sample, estimator_fun = aggvarFit))
print(mean_and_var_estimate(fgn_data = fgn_sample, estimator_fun = diffvarFit))
print(mean_and_var_estimate(fgn_data = fgn_sample, estimator_fun = absvalFit))
#' print(mean_and_var_estimate(fgn_data = fgn_sample,
#'                             estimator_fun = higuchiFit))
#' print(mean_and_var_estimate(fgn_data = fgn_sample, estimator_fun = pengFit))
#' rsFit has a high standard deviations
#' print(mean_and_var_estimate(fgn_data = fgn_sample, estimator_fun = rsFit))
print(mean_and_var_estimate(fgn_data = fgn_sample, estimator_fun = perFit))
#' boxperFit tends to strongly underestimate
#' print(mean_and_var_estimate(fgn_data = fgn_sample,
#'                             estimator_fun = boxperFit))

# results:

#"Aggregated Variance Method"
#"overall time: 1.195"
#"Differenced Aggregated Variance"
#"overall time: 1.03299999999999"
#"Absolute Moment - No. 1"
#"overall time: 0.970000000000011"
#"Higuchi Method"
#"overall time: 7.198"
#"Peng Method"
#"overall time: 71.414"
#"R/S Method"
#
#"overall time: 0.504000000000019"
#"Periodogram Method"
#"overall time: 0.522000000000048"
#"Boxed Periodogram"
#"overall time: 0.545000000000016"
#"Aggregated Variance Method"
#"mean: 0.677286469854631 , standard deviation: 0.0426196817324083"
#"Differenced Aggregated Variance"
#"mean: 0.759534413200037 , standard deviation: 0.0675571365827665"
#"Absolute Moment - No. 1"
#"mean: 0.688786925869913 , standard deviation: 0.0508917313229061"
#"Higuchi Method"
#"mean: 0.666372132581795 , standard deviation: 0.0354403572715817"
#"Peng Method"
#"mean: 0.685477964846515 , standard deviation: 0.0260295596672761"
#"R/S Method"
#"mean: 0.695710477045381 , standard deviation: 0.0873768452040556"
#"Periodogram Method"
#"mean: 0.710900428971295 , standard deviation: 0.0317333692886241"
#"Boxed Periodogram"
#"mean: 0.64353694400691 , standard deviation: 0.0265709820541408"
