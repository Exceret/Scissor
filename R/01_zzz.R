# ? Package startup messages
.onAttach <- function(libname, pkgname) {
    pkg_version <- utils::packageVersion(pkgname)

    msg <- cli::cli_fmt(cli::cli_alert_success(
        "{.pkg {pkgname}} v{pkg_version} loaded"
    ))
    packageStartupMessage(msg)
    invisible()
}

.onLoad <- function(libname, pkgname) {
    # Add timestamp to cli functions
    assign(
        "ts_cli",
        SigBridgeRUtils::CreateTimeStampCliEnv(),
        envir = asNamespace(pkgname)
    )

    invisible()
}
