#' @title Cross-Validation for boosting GaussianLSS
#' @description Cross-validation for boosting GaussianLSS
#' @param y response variable
#' @param data feature matrix
#' @param m_stop stopping iteration, default is 1000
#' @param center_x center the feature matrix? default is TRUE
#' @param folds a matrix indicating the folds of CV, specified by mboost::cv,
#'       default is 10-folds CV.
#' @param method step length method, should be one of {"FSL", "ASL", "SAASL", "SAASL05"},
#'       default is "SAASL"
#' @param papply parallel computation, default is mclapply
#' @importFrom mboost cv
#' @importFrom parallel mclapply
#' @export
cv_risk = function(y, data, m_stop = 1000, center_x = T,
                   folds = cv(rep(1, length(y)), type = "kfold", B = 10),
                   method = "SAASL",
                   papply = mclapply) {
  force(m_stop)
  force(folds)

  mod = papply(1:ncol(folds), function(i) {
    if (is.vector(data)) {
      data = as.matrix(data)
      colnames(data) = colnames(data)
    }
    boost_gaussianLSS(y, data, m_stop = m_stop, center_x = center_x, weights = as.logical(folds[, i]),
                      method = method)
  })

  likelihood = -sapply(mod, `[[`, "like_test")

  output = list()
  output$cvlike = t(likelihood)
  output$mstop = which.min(apply(output$cvlike, 2, mean, na.rm = T))
  return(output)
}

#' @title cv_risk plot
#' @description Plot the negative log-likelihood of the cv_risk output.
#' @param cv_risk likelihood of the cv_risk output
#' @importFrom reshape2 melt
#' @importFrom ggplot2 ggplot aes geom_line geom_vline ggtitle ylab theme_light
#' @export
cvrisk_plot = function(cv_risk) {
  melt_cr = melt(cv_risk)
  colnames(melt_cr) = c("folds", "Iteration", "Likelihood")
  mean_like = apply(cv_risk, 2, mean, na.rm = T)

  ggplot() +
    geom_line(data = melt_cr, aes(x = get("Iteration"), y = get("Likelihood"), group = get("folds")),
              color = "gray") +
    geom_line(aes(x = 1:ncol(cv_risk), y = mean_like)) +
    geom_vline(xintercept = which.min(mean_like)) +
    ggtitle(paste("CV: mstop =", which.min(mean_like))) +
    ylab("Negative log-Likelihoood") +
    theme_light()
}
