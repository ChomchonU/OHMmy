# Internal helpers shared across OHMmy -----------------------------------------

#' Make a string safe for use in file names
#' @keywords internal
#' @noRd
sanitize <- function(x) gsub("[^A-Za-z0-9_.-]+", "_", x)

#' Stop with an informative message when suggested packages are missing
#' @keywords internal
#' @noRd
.check_suggested <- function(pkgs, fun) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(fun, " requires the following package(s): ", paste(missing, collapse = ", "),
         ". Install them, e.g. with BiocManager::install(c(",
         paste0('"', missing, '"', collapse = ", "), ")).", call. = FALSE)
  }
  invisible(TRUE)
}
