# ? Package startup messages
.onAttach <- function(libname, pkgname) {
  pkg_version <- utils::packageVersion(pkgname)

  msg <- cli::cli_fmt(cli::cli_alert_success(
    "{.pkg {pkgname}} v{pkg_version} loaded"
  ))
  packageStartupMessage(msg)
  invisible()
}

#' @keywords internal
ts_cli <- SigBridgeRUtils::CreateTimeStampCliEnv()

#' @importFrom data.table `%chin%`
NULL
