###################
#####  Print  #####
###################

#' @title Print APML1 Model Summary
#' @description
#' Prints a formatted summary of APML1 (Adaptive Penalized Model L1) objects.
#' Displays model type, penalty information, and regularization path details.
#'
#' @param x An APML1 model object
#' @param digits Number of significant digits to display (default: 4)
#' @param ... Additional arguments passed to print methods
#'
#' @return No return value, prints model summary to console
#'
#' @examples
#' \dontrun{
#' # After fitting an APML1 model
#' model <- APML1(x, y, family = "gaussian", penalty = "Lasso")
#' print(model)
#' }
#'
#' @method print APML1
#' @export
print.APML1 = function(x, digits = 4, ...) {
  #cat("\nCall: ", deparse(x$call))
  tem = switch(
    x$penalty,
    "Lasso" = "Lasso (L1)",
    "Enet" = "Enet (L1 + L2)",
    "Net" = "Net (L1 + Laplacian)"
  )

  if (x$family == "gaussian") {
    cat("\nRegularized linear regression: ", tem)
  } else if (x$family == "binomial") {
    cat("\nRegularized logistic regression: ", tem)
  } else if (x$family == "cox") {
    cat("\nRegularized Cox model: ", tem)
  }

  if (x$penalty %in% c("Lasso", "Enet") & x$adaptive[1]) {
    cat("\nAdaptive: Beta (L1)")
  } else if (x$penalty == "Net" & sum(x$adaptive) == 1) {
    cat("\nAdaptive: ", c("Beta (L1)", "sign (Laplacian)")[x$adaptive])
  } else if (x$penalty == "Net" & sum(x$adaptive) == 2) {
    cat("\nAdaptive: Beta (L1), sign (Laplacian)")
  }

  cat("\n\n\nThe path of lambda:\n\n")
  if (!any(colnames(x$fit) %in% "cvm")) {
    print(signif(x$fit, digits))
  } else if (any(colnames(x$fit) %in% "cvm") & is.null(x$lambda.opt)) {
    #print(signif(x$fit,digits))
    switch(
      x$family,
      "gaussian" = print(cbind(signif(x$fit[-6], digits), x$fit[6])),
      "binomial" = print(cbind(signif(x$fit[-6], digits), x$fit[6])),
      "cox" = print(cbind(signif(x$fit[-5], digits), x$fit[5]))
    )
  } else if (any(colnames(x$fit) %in% "cvm") & !is.null(x$lambda.opt)) {
    #print(signif(x$fit,digits))
    switch(
      x$family,
      "gaussian" = print(cbind(signif(x$fit[-6], digits), x$fit[6])),
      "binomial" = print(cbind(signif(x$fit[-6], digits), x$fit[6])),
      "cox" = print(cbind(signif(x$fit[-5], digits), x$fit[5]))
    )
    cat("\n\nTuning the number of non-zeros with lambda:\n\n")
    print(signif(x$fit0, digits))
  }
  cat("\n")
}
