#' sdpMig2: stochastic dynamic programming model of optimal bird migration
#'
#' @keywords internal
#' @useDynLib sdpMig2, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @import methods
"_PACKAGE"

# computeOutputs() and pltNetworkGG() use data.table syntax (:=); this tells data.table that the package is data.table-aware
.datatable.aware <- TRUE
