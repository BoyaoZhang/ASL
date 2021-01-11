library(parallel)
library(ggplot2)
library(reshape2)
# simulate data
set.seed(123)

n <- 1000

for (i in 1:3) {
  assign(paste0("x", i), runif(n, -1, 1))
}

mu <- 1 + x1 + 2 * x2 - x3
sigma <- exp(.1 + .1 * x1 - .2 * x2 + .1 * x3)

y <- rnorm(n, mu, sigma)
data <- as.matrix(as.data.frame(mget(paste0("x", 1:3))))

folds <- mboost::cv(rep(1, length(y)), type = "kfold", B = 10)

# 10-folds CV
cv_saasl <- cv_risk(y, data, folds = folds, method = "SAASL")
# plot CV
cvrisk_plot(cv_saasl$cvlike)

# fit GaussianLSS
b_saasl <- boost_gaussianLSS(y, data, m_stop = 1000, center_x = TRUE, method = "SAASL")
# get the coefficients at stopping iteration
b_saasl$mu_mat[, cv_saasl$mstop]
b_saasl$si_mat[, cv_saasl$mstop]

matplot(t(b_saasl$mu_mat), type = "l")
matplot(t(b_saasl$si_mat), type = "l")
