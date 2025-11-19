#' @keywords internal
test_cox <- function(X, Y, network, alpha, cell_num, n = 100, nfold = 10, ...) {
    # * SigBridgeR Config
    dots <- rlang::list2(...)
    seed <- dots$seed %||% SigBridgeRUtils::getFuncOption("seed")
    verbose <- dots$verbose %||% SigBridgeRUtils::getFuncOption("verbose")

    set.seed(seed)
    m1 <- sum(Y[, 2] == 1)
    m2 <- sum(Y[, 2] == 0)
    index0 <- sample(cut(seq(m1 + m2), breaks = nfold, labels = F))
    #index1 <- sample(cut(seq(m1), breaks = nfold, labels = F))
    #index2 <- sample(cut(seq(m2), breaks = nfold, labels = F))

    Matrix::print("|**************************************************|")
    ts_cli$cli_alert_info("Perform cross-validation on X with true label")
    c_index_test_real <- NULL
    cli::cli_progress_bar("CV with true labels", total = nfold)
    for (j in 1:nfold) {
        c_index <- Matrix::which(index0 == j)
        #c_index <- c(Matrix::which(Y[,2] == 1)[Matrix::which(index1 == j)], Matrix::which(Y[,2] == 0)[Matrix::which(index2 == j)])
        X_train <- X[-c_index, ]
        Y_train <- Y[-c_index, ]
        fit <- NULL
        while (is.null(fit$fit)) {
            set.seed(123)
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
        Coefs <- as.numeric(fit$Beta[, index])
        Cell1 <- Coefs[Matrix::which(Coefs > 0)]
        Cell2 <- Coefs[Matrix::which(Coefs < 0)]

        X_test <- X[c_index, ]
        Y_test <- Y[c_index, ]
        test_data <- data.frame(cbind(Y_test, X_test %*% Coefs))
        colnames(test_data) <- c("OS_time", "Status", "Prediction")
        res.cox <- survival::coxph(
            survival::Surv(OS_time, Status) ~ Prediction,
            data = test_data
        )
        c_index_test_real[j] <- survival::concordance(res.cox)$concordance

        #pb1$tick()
        Sys.sleep(1 / 100)
        cli::cli_progress_update()
        if (j == nfold) {
            cli::cli_progress_done()
            cat("Finished!\n")
        }
    }

    Matrix::print("|**************************************************|")
    ts_cli$cli_alert_info("Perform cross-validation on X with permutated label")
    c_index_test_back <- list()
    cli::cli_progress_bar("CV with permutated labels", total = n)
    for (i in 1:n) {
        set.seed(i + 100)
        c_index_test_back[[i]] <- matrix(
            0,
            nfold,
            1,
            dimnames = list(paste0("Testing_", 1:nfold), "Concordance")
        )
        Y2 <- Y[sample(nrow(Y)), ]
        for (j in 1:nfold) {
            c_index <- Matrix::which(index0 == j)
            #c_index <- c(Matrix::which(Y2[,2] == 1)[Matrix::which(index1 == j)], Matrix::which(Y2[,2] == 0)[Matrix::which(index2 == j)])
            X_train <- X[-c_index, ]
            Y_train <- Y2[-c_index, ]
            fit <- NULL
            while (is.null(fit$fit)) {
                set.seed(123)
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
            Coefs <- as.numeric(fit$Beta[, index])
            Cell1 <- Coefs[Matrix::which(Coefs > 0)]
            Cell2 <- Coefs[Matrix::which(Coefs < 0)]

            X_test <- X[c_index, ]
            Y_test <- Y2[c_index, ]
            test_data <- data.frame(cbind(Y_test, X_test %*% Coefs))
            colnames(test_data) <- c("OS_time", "Status", "Prediction")
            res.cox <- survival::coxph(
                survival::Surv(OS_time, Status) ~ Prediction,
                data = test_data
            )
            c_index_test_back[[i]][j] <- survival::concordance(
                res.cox
            )$concordance
        }
        #pb2$tick()
        Sys.sleep(1 / 100)
        cli::cli_progress_update()
        if (i == n) {
            cli::cli_process_done()
            cat("Finished!\n")
        }
    }
    statistic <- Matrix::mean(c_index_test_real)
    background <- NULL
    for (i in 1:n) {
        background[i] <- Matrix::mean(c_index_test_back[[i]][, 1])
    }
    p <- sum(background > statistic) / n

    Matrix::print(sprintf(
        "Test statistic = %s",
        formatC(statistic, format = "f", digits = 3)
    ))
    Matrix::print(sprintf(
        "Reliability significance test p = %s",
        formatC(p, format = "f", digits = 3)
    ))

    return(list(
        statistic = statistic,
        p = p,
        c_index_test_real = c_index_test_real,
        c_index_test_back = c_index_test_back
    ))
}
