#' @keywords internal
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
    parallel <- dots$parallel %||% SigBridgeRUtils::getFuncOption("parallel")
    workers <- dots$workers %||% SigBridgeRUtils::getFuncOption("workers")
    parallel_type <- dots$parallel.type %||%
        SigBridgeRUtils::getFuncOption("parallel.type")

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

    if (verbose) {
        ts_cli$cli_alert_info("Permutation test")
    }

    MSE_test_back <- if (parallel) {
        SigBridgeRUtils::plan(parallel_type, workers = workers)
        on.exit(SigBridgeRUtils::plan("sequential"), add = TRUE)

        SigBridgeRUtils::future_map(
            seq_len(n),
            ~ PermutationValidate(
                .x, # i
                X,
                Y,
                index0,
                network,
                alpha,
                cell_num,
                nfold,
                m,
                seed
            ),
            .options = SigBridgeRUtils::furrr_options(seed = TRUE),
            .progress = verbose
        )
    } else {
        purrr::map(
            seq_len(n),
            ~ PermutationValidate(
                .x, # i
                X,
                Y,
                index0,
                network,
                alpha,
                cell_num,
                nfold,
                m,
                seed
            ),
            .progress = verbose
        )
    }
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

#' @keywords internal
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
