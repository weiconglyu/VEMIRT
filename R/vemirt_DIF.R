# Written by Weicong Lyu

best_vemirt_DIF <- function(all, criterion) {
  ic <- sapply(all, `[[`, criterion)
  best <- which.min(ic)
  if (length(ic) > 1)
    if (best == 1)
      warning('The optimal lambda0 may be smaller than ', all[[best]]$lambda0, '.')
    else if (best == length(ic))
      warning('The optimal lambda0 may be larger than ', all[[best]]$lambda0, '.')
  best
}

new_vemirt_DIF <- function(all, criterion) {
  best <- best_vemirt_DIF(all, criterion)
  structure(list(fit = all[[best]], best = best, all = all), class = 'vemirt_DIF')
}

#' Extract Parameter Estimates from DIF Analysis
#'
#' @param object An object of class \code{vemirt_DIF}
#' @param criterion Information criterion for model selection, one of \code{'AIC'}, \code{'BIC'} and \code{'GIC'}, otherwise use the criterion specified when fitting the model(s)
#'
#' @usage coef(object, criterion = NULL)
#'
#' @seealso \code{\link{em_DIF}}, \code{\link{gvemm_DIF}}, \code{\link{lrt_DIF}}, \code{\link{print.vemirt_DIF}}
#' @export
coef.vemirt_DIF <- function(object, criterion = NULL, ...) {
  if (length(criterion) == 1 && criterion %in% c('AIC', 'BIC', 'GIC'))
    object$all[[best_vemirt_DIF(object$all, criterion)]]
  else
    object$fit
}

#' Print DIF Items
#'
#' @param x An object of class \code{vemirt_DIF}
#' @param criterion Information criterion for model selection, one of \code{'AIC'}, \code{'BIC'} and \code{'GIC'}, otherwise use the criterion specified when fitting the model(s)
#'
#' @usage print(x, criterion = NULL)
#'
#' @seealso \code{\link{em_DIF}}, \code{\link{gvemm_DIF}}, \code{\link{lrt_DIF}}, \code{\link{coef.vemirt_DIF}}
#' @export
print.vemirt_DIF <- function(x, criterion = NULL, ...) {
  fit <- coef.vemirt_DIF(x, criterion)
  cat(sep = '', 'lambda0 = ', fit$lambda0, '\n')
  for (k in 2:nrow(fit$beta)) {
    cat(sep = '', '\nGroup ', k, ':\n')
    dif <- cbind(fit$gamma[k, , ], fit$beta[k, ]) != 0
    cat(sep = '', '\t\t', paste0('a', 1:(ncol(dif) - 1), collapse = '\t'), '\tb\n')
    for (i in 1:nrow(dif))
      cat(sep = '', 'Item ', i, '\t', paste0(c(' ', 'X')[dif[i, ] + 1], collapse = '\t'), '\n')
  }
  invisible(fit)
}
