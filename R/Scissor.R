#' @title Optimized Scissor Algorithm for Seurat ver5
#' @description
#' Scissor.v5 from `https://doi.org/10.1038/s41587-021-01091-3`and `https://github.com/sunduanchen/Scissor/issues/59`
#' Another version of Scissor.v5() to optimize memory usage and execution speed in preprocess.
#'
#' @references
#' Sun D, Guan X, Moran AE, Wu LY, Qian DZ, Schedin P, et al. Identifying phenotype-associated subpopulations by integrating bulk and single-cell sequencing data. Nat Biotechnol. 2022 Apr;40(4):527–38.
#'
#' @section LICENSE:
#' Licensed under the GNU General Public License version 3 (GPL-3.0).
#' A copy of the license is available at <https://www.gnu.org/licenses/gpl-3.0.en.html>.
#'
#' @family scissor
#' @export
#'
Scissor.v5.optimized <- function(
    bulk_dataset,
    sc_dataset,
    phenotype,
    tag = NULL,
    alpha = NULL,
    cutoff = 0.2,
    family = c("gaussian", "binomial", "cox"),
    Save_file = "Scissor_inputs.RData",
    Load_file = NULL,
    verbose = SigBridgeRUtils::getFuncOption("verbose"),
    seed = SigBridgeRUtils::getFuncOption("seed"),
    ...
) {
    if (verbose) {
        ts_cli$cli_alert_info(
            cli::col_green("Scissor start...")
        )
    }

    if (is.null(Load_file)) {
        if (verbose) {
            ts_cli$cli_alert_info("Start from raw data...")
        }
        common <- intersect(
            rownames(bulk_dataset),
            rownames(sc_dataset)
        )
        if (length(common) == 0) {
            cli::cli_abort(c(
                "x" = "There is no common genes between the given single-cell and bulk samples. Please check Scissor inputs."
            ))
        }

        if (inherits(sc_dataset, "Seurat")) {
            # sc_exprs <- as.matrix(sc_dataset@assays$RNA$data)
            sc_exprs <- SeuratObject::LayerData(sc_dataset, layer = "data")

            if ("RNA_snn" %chin% names(sc_dataset@graphs)) {
                # network <- as.matrix(sc_dataset@graphs$RNA_snn)
                network <- SeuratObject::Graphs(
                    sc_dataset,
                    slot = "RNA_snn"
                )

                if (verbose) {
                    cli::cli_alert_info(
                        "Using {.val RNA_snn} graph for network."
                    )
                }
            } else if ("integrated_snn" %chin% names(sc_dataset@graphs)) {
                # network <- as.matrix(sc_dataset@graphs$integrated_snn)
                network <- SeuratObject::Graphs(
                    sc_dataset,
                    slot = "integrated_snn"
                )

                if (verbose) {
                    cli::cli_alert_info(
                        "Using {.val integrated_snn} graph for network."
                    )
                }
            } else {
                cli::cli_abort(c(
                    "x" = "No `RNA_snn` or `integrated_snn` graph in the given Seurat object. Please check Scissor inputs."
                ))
            }
        } else {
            sc_exprs <- Matrix::Matrix(as.matrix(sc_dataset))
            Seurat_tmp <- SCPreProcess(
                sc_dataset,
                quality_control.pattern = c("^MT-"),
                verbose = FALSE
            )
            network <- SeuratObject::Graphs(
                Seurat_tmp,
                slot = "RNA_snn"
            )
            rm(Seurat_tmp)
        }
        Matrix::diag(network) <- 0
        network <- as.matrix((network != 0) * 1)

        # bulk_mat <- as.matrix(bulk_dataset[common, ])
        bulk_mat <- Matrix::Matrix(as.matrix(bulk_dataset[common, ]))

        # sc_mat <- as.matrix(sc_exprs[common, ])
        sc_mat <- sc_exprs[common, ] # A dgCMatrix object

        dataset0 <- Matrix::cbind2(bulk_mat, sc_mat) # much smaller
        if (verbose) {
            ts_cli$cli_alert_info(
                "Normalizing quantiles of data"
            )
        }

        dataset1 <- SigBridgeRUtils::normalize.quantiles(
            as.matrix(dataset0),
            keep.names = TRUE
        )
        # rownames(dataset1) <- common
        # colnames(dataset1) <- c(colnames(bulk_mat), colnames(sc_mat))
        if (verbose) {
            ts_cli$cli_alert_info(
                "Subsetting data"
            )
        }

        n_bulk <- ncol(bulk_mat)
        # gene-sample
        Expression_bulk <- dataset1[, seq_len(n_bulk), drop = FALSE]
        # gene-cell
        Expression_cell <- dataset1[, (n_bulk + 1):ncol(dataset1), drop = FALSE]

        gc(verbose = FALSE)
        if (verbose) {
            ts_cli$cli_alert_info(
                "Calculating correlation"
            )
        }

        X <- stats::cor(Expression_bulk, Expression_cell)

        quality_check <- colQuantiles(X, probs = seq(0, 1, 0.25))
        if (verbose) {
            cli::cli_text(
                strrep("-", floor(getOption("width") / 2)),
                "\n",
                sep = ""
            )
            cli::cli_text("Five-number summary of correlations:")
            quality_check %>%
                asplit(2) %>%
                purrr::map_dbl(mean) %>%
                round(digits = 6) %>%
                paste(sep = " ", collapse = " ") %>%
                cli::cli_text()
            cli::cli_text(
                strrep("-", floor(getOption("width") / 2)),
                "\n",
                sep = ""
            )
        }
        # median
        if (quality_check[3] < 0.01) {
            cli::cli_warn(
                "The median correlation between the single-cell and bulk samples is relatively low."
            )
        }

        FamilyProcessor <- list(
            binomial = function() {
                Y <- as.numeric(phenotype)
                z <- table(Y)
                if (length(z) != length(tag)) {
                    cli::cli_abort(
                        "x" = "The length differs between tags and phenotypes. Please check Scissor inputs and selected regression type."
                    )
                }
                if (verbose) {
                    cli::cli_alert_info(
                        "Current phenotype contains {.val {z[1]}} {tag[1]} and {.val {z[2]}} {tag[2]} samples."
                    )
                    ts_cli$cli_alert_info(
                        "Perform logistic regression on the given phenotypes..."
                    )
                }
                Y
            },
            gaussian = function() {
                Y <- as.numeric(phenotype)
                z <- table(Y)
                if (length(z) != length(tag)) {
                    cli::cli_abort(
                        "x" = "The length differs between tags and phenotypes. Please check Scissor inputs and selected regression type.",
                        "i" = "length of tags: {.val {length(tag)}}",
                        "i" = "length of phenotypes: {.val {length(z)}}"
                    )
                }
                if (verbose) {
                    tmp <- paste(z, tag)

                    cli::cli_alert_info(
                        "Current phenotype contains: {.val {length(tmp)}} samples."
                    )
                    cli::cli_text("Sample examples:")
                    cli::cli_bullets(c(
                        " " = "{.val {head(tmp, 5)}}",
                        " " = "... ({length(tmp)-6} more samples)",
                        " " = "{.val {tail(tmp, 1)}}"
                    ))
                    ts_cli$cli_alert_info(
                        "Perform linear regression on the given phenotypes..."
                    )
                }
                Y
            },
            cox = function() {
                Y <- as.matrix(phenotype)
                if (ncol(Y) != 2) {
                    cli::cli_abort(
                        "x" = "The size of survival data is wrong. Please check Scissor inputs and selected regression type."
                    )
                }
                if (verbose) {
                    ts_cli$cli_alert_info(
                        "Perform cox regression on the given clinical outcomes..."
                    )
                }
                Y
            }
        )

        Y <- FamilyProcessor[[family]]()

        if (!is.null(Save_file)) {
            save(
                X,
                Y,
                network,
                Expression_bulk,
                Expression_cell,
                file = Save_file
            )
            if (verbose) {
                ts_cli$cli_alert_success(
                    "Statistics data saved to {.file {Save_file}}."
                )
            }
        }
    } else {
        # Load data from previous work
        if (verbose) {
            ts_cli$cli_alert_info(
                "Loading data from {.file {Load_file}}"
            )
        }
        load(Load_file)
    }

    # garbage collection
    rm(
        Expression_bulk,
        Expression_cell,
        sc_dataset,
        bulk_dataset,
        phenotype
    )
    gc(verbose = FALSE)
    if (verbose) {
        ts_cli$cli_alert_info("Screening...")
    }

    alpha <- alpha %||%
        c(0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

    results <- list()

    for (i in seq_along(alpha)) {
        set.seed(seed)

        fit0 <- APML1(
            X,
            Y,
            family = family,
            penalty = "Net",
            alpha = alpha[i],
            Omega = network,
            nlambda = 100,
            nfolds = min(10, nrow(X))
        )

        fit1 <- APML1(
            X,
            Y,
            family = family,
            penalty = "Net",
            alpha = alpha[i],
            Omega = network,
            lambda = fit0$lambda.min
        )

        if (family == "binomial") {
            Coefs <- as.numeric(fit1$Beta[2:(ncol(X) + 1)])
        } else {
            Coefs <- as.numeric(fit1$Beta)
        }

        pos_mask <- Coefs > 0
        neg_mask <- Coefs < 0
        cells <- colnames(X)
        Cell1 <- cells[pos_mask]
        Cell2 <- cells[neg_mask]
        percentage <- (length(Cell1) + length(Cell2)) / length(cells)

        if (verbose) {
            cli::cli_h2("At alpha = {.val {alpha[i]}}")
            cli::cli_text(sprintf(
                "Scissor identified {.val {%d}} Scissor+ cells and {.val {%d}} Scissor- cells.",
                length(Cell1),
                length(Cell2)
            ))
            cli::cli_text(sprintf(
                "The percentage of selected cell is: {.val {%s}}%%",
                round(percentage * 100, digits = 3)
            ))
        }

        if (percentage < cutoff) {
            ts_cli$cli_alert_info(cli::col_green("Scissor Ended."))
            break
        }
    }

    list(
        para = list(
            alpha = alpha,
            lambda = fit0$lambda.min,
            family = family,
            Coefs = Coefs # for miscellaneous informationss
        ),
        Coefs = Coefs, # for cell evaluation
        Scissor_pos = Cell1,
        Scissor_neg = Cell2,
        X = X,
        Y = Y,
        network = network
    )
}
