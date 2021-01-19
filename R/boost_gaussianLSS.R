#' get residuals and estimated coefficients
#' @param y response variable
#' @param X feature matrix
#' @keywords internal
#' @export
getRes <- function(y, X) {
  blong <- NULL
  for (i in 1:ncol(X)) {
    x1 <- cbind(1, X[, i])
    tmp <- solve(t(x1) %*% x1) %*% (t(x1) %*% y)
    blong <- c(blong, tmp)
  }
  pred <- t(blong[seq(1, length(blong), 2)] + blong[seq(2, length(blong), 2)] * t(X))
  res <- colSums((y - pred)^2)
  epsilon <- colMeans(y - pred)
  epsilon_vec <- y - pred
  return(list(blong = blong,
              res = res,
              epsilon = epsilon,
              epsilon_vec = epsilon_vec))
}

#' @title Boosting GaussianLSS Model
#' @description Fit boosting GaussianLSS model with different kinds of step lengths.
#' @param y response variable
#' @param data feature matrix
#' @param m_stop stopping iteration, default is 1000
#' @param center_x center the feature matrix? default is TRUE
#' @param weights a vector of weights indicating the training (1) and testing (0) data,
#'       default is NULL, i.e. all observations as training data
#' @param method step length method, should be one of {"FSL", "ASL", "SAASL", "SAASL05"},
#'       default is "SAASL"
#' @return \code{mu} the estiamted coefficients of mean \cr
#'         \code{sigma} the estimated coefficients of the predictor of standard deviation \cr
#'         \code{mu_mat} matrix of the estimated coefficients of mu in each iteartion \cr
#'         \code{si_mat} matrix of the estimated coefficients of sigma in each iteration \cr
#'         \code{v_mu} step length of mu in each iteartion \cr
#'         \code{v_mu} step length of sigma in each iteration \cr
#'         \code{muvsi} which distribution paramter is updated in each iteration \cr
#'         \code{v_mu_var} which variable is used for updating mu in each iteration \cr
#'         \code{v_si_var} which variable is used for updating sigma in each iteartion \cr
#'         \code{like_train} positive likelihood of the training data in each iteartion \cr
#'         \code{like_test} positive likelihood of the test data in each iteration
#' @importFrom stats sd dnorm optimize
#' @export
boost_gaussianLSS <- function(y, data, m_stop = 1000, center_x = TRUE, weights = NULL,
                              method = "SAASL") {
  do_m <- do_s <- TRUE

  muvsi <- rep(0, m_stop)

  # Split data into train and test sets. If weights=NULL, then data = train_set = test_set
  if (!is.null(weights)) {
    weights <- as.logical(weights)
    y_train <- y[weights]
    data_train <- data[weights, ]
    y_test <- y[!weights]
    data_test <- data[!weights, ]
  }
  else {
    y_train <- y
    data_train <- data
    y_test <- y
    data_test <- data
  }

  if (center_x) {
    data_train <- apply(data_train, FUN = function(x) scale(x, scale = FALSE, center = TRUE), MARGIN = 2)
    data_test <- apply(data_test, FUN = function(x) scale(x, scale = FALSE, center = TRUE), MARGIN = 2)
  }
  else {
    data_train <- as.matrix(data_train)
    data_test <- as.matrix(data_test)
  }

  # add intercept
  X_train <- cbind(1, data_train)
  X_test <- cbind(1, data_test)
  p <- ncol(data_train)

  # initialize beta_mu and beta_si
  beta_mu <- matrix(c(mean(y), rep(0, p)), ncol = 1)
  beta_si <- matrix(c(log(sd(y)), rep(0, p)), ncol = 1)

  # initialize other outputs
  mu_train <- mu_test <- list()
  sigma_train <- sigma_test <- list()
  ngr_mu <- list()
  ngr_si <- list()

  mu_mat <- si_mat <- matrix(nrow = ncol(X_train), ncol = m_stop)

  v_mu <- v_mu_var <- vector(length = m_stop)
  v_si <- v_si_var <- vector(length = m_stop)

  Xbs_train <- X_train %*% beta_si
  Xbs_test <- X_test %*% beta_si
  Xbm_train <- rep(mean(y_train), length(y_train))
  Xbm_test <- rep(mean(y_test), length(y_test))
  step.length.mu <- step.length.si <- 0
  update.si <- update.mu <- 0

  like_train <- like_test <- vector(length = m_stop)
  h_mu <- hs2 <- vector(length = m_stop)

  # boosting start
  for (m in 1:m_stop) {
    mu_train[[m]] <- Xbm_train
    mu_test[[m]] <- Xbm_test
    sigma_train[[m]] <- exp(Xbs_train)
    sigma_test[[m]] <- exp(Xbs_test)

    # potential update of mu
    if (do_m) {
      # calculate negative gradient
      ngr_mu[[m]] <- (1/sigma_train[[m]]^2) * (y_train - mu_train[[m]])
      RS <- getRes(as.vector(ngr_mu[[m]]), data_train)

      # select the best performing variable
      best_mu <- which.min(RS$res)
      v_mu_var[m] <- best_mu
      beta_best_mu <- RS$blong[c(best_mu*2-1,best_mu*2)]
      update.mu_train <- X_train[,best_mu+1]*beta_best_mu[2] + beta_best_mu[1]
      update.mu_test <- X_test[,best_mu+1]*beta_best_mu[2] + beta_best_mu[1]

      # find the optimal step length by line search
      if (method == "ASL") {
        phi_mu <- function(v) {
          ret <- -sum(dnorm(y_train, mu_train[[m]] + v * update.mu_train, sigma_train[[m]], log = TRUE))
          return(ret)
        }
        v <- optimize(phi_mu, interval = c(-1, 10))$min
      }
      else if (method %in% c("SAASL", "SAASL05")) {
        h_mu[m] <- sum(update.mu_train^2)
        hs2[m] <- sum(update.mu_train^2/sigma_train[[m]]^2)

        v <- h_mu[m] / hs2[m]
      }
      else if (method == "FSL") {
        v <- 1
      }

      # the adaptive step length is 10% of the optimal step length
      v_mu[m] <- v
      step.length.mu <- .1 * v_mu[m]

      # update the predictors
      beta_mu_st <- beta_mu
      beta_mu_st[1] <- beta_mu[1]+step.length.mu*beta_best_mu[1]
      beta_mu_st[1+best_mu] <- beta_mu[1+best_mu]+step.length.mu*beta_best_mu[2]
    }

    # potential update of sigma
    if (do_s) {
      ngr_si[[m]] <- (-1 + exp(-2 * Xbs_train) * ((y_train - mu_train[[m]])^2))
      RS <- getRes(as.vector(ngr_si[[m]]), data_train)

      best_si <- which.min(RS$res)
      v_si_var[m] <- best_si
      beta_best_si <- RS$blong[c(best_si*2-1,best_si*2)]
      update.si_train <- X_train[,best_si+1]*beta_best_si[2] + beta_best_si[1]
      update.si_test <- X_test[,best_si+1]*beta_best_si[2] + beta_best_si[1]

      if (method %in% c("ASL", "SAASL", "SAASL05")) {
        phi_si <- function(v) {
          ret <- -sum(dnorm(y_train, mu_train[[m]], exp(Xbs_train + v * update.si_train), log = TRUE))
          return(ret)
        }
        v <- optimize(phi_si, interval = c(-1, 1))$min
      }
      else if (method == "FSL") {
        v <- 1
      }

      v_si[m] <- v
      step.length.si <- .1 * v_si[m]

      beta_si_st <- beta_si
      beta_si_st[1] <- beta_si[1]+step.length.si*beta_best_si[1]
      beta_si_st[1+best_si] <- beta_si[1+best_si]+step.length.si*beta_best_si[2]
    }

    # select the best performing distribution parameter
    if (do_m & do_s) {
      # compare the positive likelihood
      like_mu <- sum(dnorm(y_test, mu_test[[m]] + step.length.mu * update.mu_test, sigma_test[[m]], log = TRUE))
      like_si <- sum(dnorm(y_test, mu_test[[m]], exp(Xbs_test + step.length.si * update.si_test), log = TRUE))

      # if like_mu > like_si, then update mu
      if (like_mu > like_si) {
        beta_mu <- beta_mu_st
        Xbm_train <- Xbm_train + step.length.mu * update.mu_train
        Xbm_test <- Xbm_test + step.length.mu * update.mu_test
        muvsi[m] <- 1
      }
      # else update sigma
      else {
        beta_si <- beta_si_st
        Xbs_train <- Xbs_train + step.length.si * update.si_train
        Xbs_test <- Xbs_test + step.length.si * update.si_test
        muvsi[m] <- 2
      }
    }

    # # if either paramter stops updating, update the others
    # else {
    #   if (do_m) {
    #     beta_mu <- beta_mu_st
    #     Xbm_train <- Xbm_train + step.length.mu * update.mu_train
    #     Xbm_test <- Xbm_test + step.length.mu * update.mu_test
    #     muvsi[m] <- 1
    #   }
    #   if (do_s) {
    #     beta_si <- beta_si_st
    #     Xbs_train <- Xbs_train + step.length.si * update.si_train
    #     Xbs_test <- Xbs_test + step.length.si * update.si_test
    #     muvsi[m] <- 2
    #   }
    # }

    # update the predictors of the best performing paramter
    mu_mat[, m] <- beta_mu
    si_mat[, m] <- beta_si
    like_train[m] <- sum(dnorm(y_train, mu_train[[m]], sigma_train[[m]], log = TRUE))
    like_test[m] <- sum(dnorm(y_test, mu_test[[m]], sigma_test[[m]], log = TRUE))
  }

  # format the output
  rownames(beta_mu) <- c("(Intercept)", colnames(data))
  rownames(beta_si) <- c("(Intercept)", colnames(data))

  return(list(mu = t(beta_mu),
              sigma = t(beta_si),
              mu_mat = mu_mat,
              si_mat = si_mat,
              v_si = v_si,
              v_mu = v_mu,
              muvsi = muvsi,
              v_mu_var = v_mu_var,
              v_si_var = v_si_var,
              like_train = like_train,
              like_test = like_test,
              call = sys.call()))
}
