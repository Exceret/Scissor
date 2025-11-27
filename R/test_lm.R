#' @title Perform linear model testing with cross-validation and permutation tests
#'
#' @param X A matrix or Matrix object containing the input features
#' @param Y A vector containing the response variable
#' @param network A network object used in the analysis
#' @param alpha A numeric value for the significance level
#' @param cell_num Number of cells (observations)
#' @param n Number of permutations to perform (default: 100)
#' @param nfold Number of folds for cross-validation (default: 10)
#' @param ... Additional parameters including:
#'   \itemize{
#'     \item{seed}{Random seed (default: from SigBridgeRUtils options)}
#'     \item{verbose}{Whether to show progress messages (default: from SigBridgeRUtils options)}
#'   }
#'
#' @return A list containing:
#'   \itemize{
#'     \item{statistic}{Mean MSE test statistic}
#'     \item{p}{P-value from reliability test}
#'     \item{MSE_test_real}{Vector of MSE values from real label cross-validation}
#'     \item{MSE_test_back}{List of MSE values from permutation tests}
#'   }
#'
#' @export
#'
test_lm <- function(
    X,
    Y,
    network,
    alpha,
    cell_num,
    n = 100,
    nfold = 10,
    ...
) {
    dots <- rlang::list2(...)
    seed <- dots$seed %||% SigBridgeRUtils::getFuncOption("seed")
    verbose <- dots$verbose %||% SigBridgeRUtils::getFuncOption("verbose")

    set.seed(seed)
    m <- nrow(X)
    index0 <- sample(cut(seq_len(m), breaks = nfold, labels = FALSE))

    if (verbose) {
        ts_cli$cli_alert_info(
            "Performing {nfold}-fold cross-validation on X with true labels"
        )
    }
    X <- Matrix::Matrix(X)
    MSE_test_real <- numeric(nfold)

    if (verbose) {
        cli::cli_progress_bar("CV with true labels", total = nfold)
    }
    for (j in seq_len(nfold)) {
        MSE_test_real[j] <- ComputeFold(
            j,
            X,
            Y,
            index0,
            network,
            alpha,
            cell_num
        )
        if (verbose) {
            cli::cli_progress_update()
        }
        Sys.sleep(1 / 100)
    }
    if (verbose) {
        cli::cli_progress_done()
    }

    if (verbose) {
        ts_cli$cli_alert_info(
            "Perform cross-validation on X with permutated label"
        )
    }

    MSE_test_back <- purrr::map(
        seq_len(n),
        ~ PermutationValidate(
            i = .x, # i
            X = X,
            Y = Y,
            index0 = index0,
            network = network,
            alpha = alpha,
            cell_num = cell_num,
            nfold = nfold,
            m = m,
            seed = seed
        ),
        .progress = if (verbose) 'Permutation Validation' else FALSE
    )

    names(MSE_test_back) <- paste0("n_", seq_len(n))
    background <- purrr::map_dbl(MSE_test_back, mean)

    statistic <- Matrix::mean(background, na.rm = TRUE)
    p <- sum(background < statistic, na.rm = TRUE) / sum(!is.na(background))

    if (verbose) {
        cli::cli_h3("Results")
        cli::cli_alert_success(
            "Test statistic (mean MSE) = {.val {round(statistic, 3)}}"
        )
        cli::cli_alert_success(
            "Reliability test p-value = {.val {round(p, 3)}}"
        )
    }

    if (p >= 0.05) {
        cli::cli_warn(
            "Model is NOT significantly better than random (p >= 0.05)"
        )
    }

    list(
        statistic = statistic,
        p = p,
        MSE_test_real = MSE_test_real,
        MSE_test_back = MSE_test_back
    )
}


#' @title Compute cross-validation fold performance for network-regularized regression
#'
#' @description
#' Performs network-regularized regression (APML1) on training data and evaluates
#' performance on test fold. Uses elastic net penalty with network regularization.
#' Automatically selects lambda based on target sparsity (cell_num)
#'
#' @param j Integer indicating the current fold index
#' @param X Matrix of predictor variables for all samples
#' @param Y Vector of response values for all samples
#' @param index0 Vector indicating fold assignments for samples
#' @param network Network/penalty matrix for regularization
#' @param alpha Elastic net mixing parameter (0 = ridge, 1 = lasso)
#' @param cell_num Target number of non-zero coefficients
#' @param seed Random seed for reproducibility
#'
#' @export
#' @return Mean squared error for the test fold
#'
ComputeFold <- function(j, X, Y, index0, network, alpha, cell_num, seed) {
    c_index <- which(index0 == j)
    X_train <- X[-c_index, , drop = FALSE]
    Y_train <- Y[-c_index]
    X_test <- X[c_index, , drop = FALSE]
    Y_test <- Y[c_index]

    # 拟合模型
    fit <- NULL
    while (is.null(fit$fit)) {
        set.seed(seed)
        fit <- APML1(
            X_train,
            Y_train,
            family = "gaussian",
            penalty = "Net",
            alpha = alpha,
            Omega = network,
            nlambda = 100
        )
    }

    index <- which.min(abs(fit$fit$nzero - cell_num))
    Coefs <- as.numeric(fit$Beta[, index])

    mean((Y_test - X_test %*% Coefs)^2)
}


#' @title Perform permutation validation for network-based analysis
#'
#' @description This function performs permutation validation by randomly shuffling
#' the response vector and computing cross-validated MSE for each fold. Used for
#' assessing significance in network-based predictive modeling.
#'
#' @param i Integer iteration counter for permutation testing
#' @param X Input feature matrix
#' @param Y Response vector
#' @param index0 Vector of initial indices
#' @param network Network/graph object used for analysis
#' @param alpha Regularization parameter
#' @param cell_num Number of cells/samples
#' @param nfold Number of cross-validation folds
#' @param m Total number of permutations
#' @param seed Random seed for reproducibility
#'
#' @return Numeric vector containing MSE (Mean Squared Error) values for each fold
#' @export
#'
PermutationValidate <- function(
    i,
    X,
    Y,
    index0,
    network,
    alpha,
    cell_num,
    nfold,
    m,
    seed
) {
    set.seed(i + seed)
    Y2 <- Y[sample(m)]

    mse_vec <- numeric(nfold)
    names(mse_vec) <- paste0("nfold_", seq_len(nfold))
    for (j in seq_len(nfold)) {
        mse_vec[j] <- ComputeFold(
            j,
            X,
            Y2,
            index0,
            network,
            alpha,
            cell_num,
            seed
        )
    }

    mse_vec
}
