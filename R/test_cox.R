#' @keywords internal
test_cox <- function(X, Y, network, alpha, cell_num, n = 100, nfold = 10, ...) {
  # * SigBridgeR Config
  dots <- rlang::list2(...)
  seed <- dots$seed %||% SigBridgeRUtils::getFuncOption("seed")
  verbose <- dots$verbose %||% SigBridgeRUtils::getFuncOption("verbose")

  set.seed(seed)
  m1 <- sum(Y[, 2] == 1)
  m2 <- sum(Y[, 2] == 0)
  index0 <- sample(cut(seq(m1 + m2), breaks = nfold, labels = FALSE))
  #index1 <- sample(cut(seq(m1), breaks = nfold, labels = F))
  #index2 <- sample(cut(seq(m2), breaks = nfold, labels = F))

  background <- numeric(n)
  if (verbose) {
    ts_cli$cli_alert_info("Perform cross-validation on X with true label")
    cli::cli_progress_bar("CV with true labels", total = nfold)
  }
  for (j in seq_len(nfold)) {
    c_index_test_real[j] <- RunCVFold(
      j = j,
      X = X,
      Y = Y,
      index0 = index0,
      network = network,
      alpha = alpha,
      cell_num = cell_num,
      seed = seed
    )

    if (verbose) {
      cli::cli_progress_update()
    }
  }

  if (verbose) {
    cli::cli_progress_done()
    ts_cli$cli_alert_info("Perform cross-validation on X with permutated label")
    cli::cli_progress_bar("CV with permutated labels", total = n)
  }

  background <- numeric(n)
  for (i in seq_len(n)) {
    background[i] <- RunPermutation(
      i = i,
      X = X,
      Y = Y,
      index0 = index0,
      network = network,
      alpha = alpha,
      cell_num = cell_num,
      nfold = nfold,
      seed = seed
    )

    if (verbose) {
      cli::cli_progress_update()
    }
  }
  if (verbose) {
    cli::cli_progress_done()
  }

  statistic <- Matrix::mean(c_index_test_real)
  p <- sum(background > statistic) / n

  if (verbose) {
    cli::cli_text(sprintf(
      "Test statistic = %s",
      formatC(statistic, format = "f", digits = 3)
    ))
    cli::cli_text(sprintf(
      "Reliability significance test p = %s",
      formatC(p, format = "f", digits = 3)
    ))
  }

  c_index_test_back <- matrix(
    background,
    ncol = 1,
    dimnames = list(paste0("Permutation_", 1:n), "Mean_Concordance")
  )

  list(
    statistic = statistic,
    p = p,
    c_index_test_real = c_index_test_real,
    c_index_test_back = c_index_test_back
  )
}


#' @title Fit a Cox proportional hazards model with network regularization
#'
#' @param X_train A matrix of predictor variables for training
#' @param Y_train A matrix of response variables for training (survival data)
#' @param network A network/adjacency matrix used for regularization
#' @param alpha The elastic net mixing parameter (0 <= alpha <= 1)
#' @param cell_num The target number of non-zero coefficients
#' @param seed Random seed for reproducibility
#'
#' @export
#' @return A numeric vector of fitted coefficients
#'
#' @details This function fits a Cox proportional hazards model using network-regularized
#'   regression (APML1 algorithm). It automatically selects the model with the number
#'   of non-zero coefficients closest to the specified cell_num. The fitting process
#'   is repeated until successful convergence.
FitModel <- function(X_train, Y_train, network, alpha, cell_num, seed) {
  fit <- NULL

  while (is.null(fit$fit)) {
    set.seed(seed)
    fit <- APML1(
      X_train,
      Y_train,
      family = "cox",
      penalty = "Net",
      alpha = alpha,
      Omega = network,
      nlambda = 100
    )
  }

  index <- which.min(abs(fit$fit$nzero - cell_num))
  as.numeric(fit$Beta[, index])
}


#' Calculate Concordance Index (C-index) for survival predictions
#'
#' Computes the C-index (concordance statistic) which evaluates the predictive accuracy
#' of survival models by comparing predicted risk scores with actual survival times.
#'
#' @param X_test A numeric matrix of test data features/observations
#' @param Y_test A numeric matrix with two columns: survival times (column 1) and event status (column 2, 1=event, 0=censored)
#' @param Coefs A numeric vector of coefficients/weights for the features
#'
#' @return A numeric value between 0 and 1 representing the concordance index, where:
#' \itemize{
#'   \item 1 indicates perfect prediction
#'   \item 0.5 indicates random prediction
#'   \item <0.5 indicates worse than random
#' }
#' @export
#' @details The function fits a Cox proportional hazards model to the test data using
#' the linear predictor (X_test %*% Coefs) as the sole covariate, then computes the
#' concordance statistic from the fitted model.
CIndexCalc <- function(X_test, Y_test, Coefs) {
  prediction <- as.vector(X_test %*% Coefs)
  test_dt <- data.table::data.table(
    OS_time = Y_test[, 1],
    Status = Y_test[, 2],
    Prediction = prediction
  )

  res.cox <- survival::coxph(
    survival::Surv(OS_time, Status) ~ Prediction,
    data = test_dt
  )

  survival::concordance(res.cox)$concordance
}


#' @title Run cross-validation fold for model evaluation
#'
#' @param j Integer indicating the current fold index
#' @param X Matrix of training data predictors
#' @param Y Matrix of training data responses
#' @param index0 Vector indicating fold assignments for each observation
#' @param network Network structure used for model fitting
#' @param alpha Regularization parameter for model fitting
#' @param cell_num Number of cells/units in the model
#' @param seed Random seed for reproducibility
#'
#' @return Computed concordance index (C-index) for the test fold
#' @export
#'
#'
RunCVFold <- function(j, X, Y, index0, network, alpha, cell_num, seed) {
  c_index <- which(index0 == j)

  X_train <- X[-c_index, , drop = FALSE]
  Y_train <- Y[-c_index, , drop = FALSE]
  X_test <- X[c_index, , drop = FALSE]
  Y_test <- Y[c_index, , drop = FALSE]

  Coefs <- FitModel(X_train, Y_train, network, alpha, cell_num, seed)
  CIndexCalc(X_test, Y_test, Coefs)
}

#' @title Run permutation test for cross-validation results
#' @description Performs permutation test by shuffling response matrix Y and running cross-validation,
#' returning the mean performance metric across all folds. Used to assess significance of model results.
#'
#' @param i Integer seed offset for random number generation
#' @param X Input feature matrix
#' @param Y Response matrix (will be permuted)
#' @param index0 Vector of initial feature indices
#' @param network Network/graph structure for the model
#' @param alpha Regularization parameter
#' @param cell_num Number of cells/samples
#' @param nfold Number of cross-validation folds
#' @param seed  Random seed for reproducibility
#'
#' @return Mean cross-validation index across all folds for permuted data
#' @export
#'
RunPermutation <- function(
  i,
  X,
  Y,
  index0,
  network,
  alpha,
  cell_num,
  nfold,
  seed
) {
  set.seed(i + 100)
  Y_perm <- Y[sample(nrow(Y)), , drop = FALSE]

  c_indices <- numeric(nfold)
  for (j in seq_len(nfold)) {
    c_indices[j] <- RunCVFold(
      j = j,
      X = X,
      Y = Y_perm,
      index0 = index0,
      network = network,
      alpha = alpha,
      cell_num = cell_num,
      seed = seed
    )
  }

  mean(c_indices)
}
