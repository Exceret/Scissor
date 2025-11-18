#' @keywords internal
test_logit <- function(X, Y, network, alpha, cell_num, n = 100, nfold = 10) {
    set.seed(2)
    m1 <- sum(Y == 1)
    m2 <- sum(Y == 0)
    index1 <- sample(cut(seq(m1), breaks = nfold, labels = F))
    index2 <- sample(cut(seq(m2), breaks = nfold, labels = F))

    Matrix::print("|**************************************************|")
    Matrix::print("Perform cross-validation on X with true label")
    AUC_test_real <- NULL
    cli::cli_progress_bar("CV with true labels", total = nfold)
    for (j in 1:nfold) {
        c_index <- c(
            Matrix::which(Y == 1)[Matrix::which(index1 == j)],
            Matrix::which(Y == 0)[Matrix::which(index2 == j)]
        )
        X_train <- X[-c_index, ]
        Y_train <- Y[-c_index]
        fit <- NULL
        while (is.null(fit$fit)) {
            set.seed(123)
            fit <- APML1(
                X_train,
                Y_train,
                family = "binomial",
                penalty = "Net",
                alpha = alpha,
                Omega = network,
                nlambda = 100
            )
        }
        index <- which.min(abs(fit$fit$nzero - cell_num))
        Coefs <- as.numeric(fit$Beta[2:(ncol(X_train) + 1), index])
        Cell1 <- Coefs[Matrix::which(Coefs > 0)]
        Cell2 <- Coefs[Matrix::which(Coefs < 0)]

        X_test <- X[c_index, ]
        Y_test <- Y[c_index]
        score_test <- 1 / (1 + exp(-X_test %*% Coefs - fit$Beta[1, index]))[, 1]
        AUC_test_real[j] <- pROC::roc(
            Y_test,
            score_test,
            direction = "<",
            quiet = T
        )$auc

        #pb1$tick()
        Sys.sleep(1 / 100)
        cli::cli_progress_update()
        if (j == nfold) {
            cli::cli_progress_done()
            cat("Finished!\n")
        }
    }

    Matrix::print("|**************************************************|")
    Matrix::print("Perform cross-validation on X with permutated label")
    AUC_test_back <- list()
    cli::cli_progress_bar("CV with permutated labels", total = n)

    for (i in 1:n) {
        set.seed(i + 100)
        AUC_test_back[[i]] <- matrix(
            0,
            nfold,
            1,
            dimnames = list(paste0("Testing_", 1:nfold), "AUC")
        )
        Y2 <- sample(Y)
        names(Y2) <- rownames(X)
        for (j in 1:nfold) {
            c_index <- c(
                Matrix::which(Y2 == 1)[Matrix::which(index1 == j)],
                Matrix::which(Y2 == 0)[Matrix::which(index2 == j)]
            )
            X_train <- X[-c_index, ]
            Y_train <- Y2[-c_index]
            fit <- NULL
            while (is.null(fit$fit)) {
                set.seed(123)
                fit <- APML1(
                    X_train,
                    Y_train,
                    family = "binomial",
                    penalty = "Net",
                    alpha = alpha,
                    Omega = network,
                    nlambda = 100
                )
            }
            index <- which.min(abs(fit$fit$nzero - cell_num))
            Coefs <- as.numeric(fit$Beta[2:(ncol(X_train) + 1), index])
            Cell1 <- Coefs[Matrix::which(Coefs > 0)]
            Cell2 <- Coefs[Matrix::which(Coefs < 0)]

            X_test <- X[c_index, ]
            Y_test <- Y2[c_index]
            score_test <- 1 /
                (1 + exp(-X_test %*% Coefs - fit$Beta[1, index]))[, 1]
            AUC_test_back[[i]][j] <- pROC::roc(
                Y_test,
                score_test,
                direction = "<",
                quiet = T
            )$auc
        }
        #pb2$tick()
        Sys.sleep(1 / 100)
        cli::cli_progress_update()
        if (i == n) {
            cli::cli_progress_done()
            cat("Finished!\n")
        }
    }
    statistic <- Matrix::mean(AUC_test_real)
    background <- NULL
    for (i in 1:n) {
        background[i] <- Matrix::mean(AUC_test_back[[i]][, 1])
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
        AUC_test_real = AUC_test_real,
        AUC_test_back = AUC_test_back
    ))
}
