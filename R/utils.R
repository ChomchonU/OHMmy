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

#' Would FindMarkers() accept this SCT assay with UMI re-correction?
#'
#' Mirrors the check in Seurat's FindMarkers() for SCT assays: every model's stored
#' median UMI must equal the smallest observed median UMI.
#' @keywords internal
#' @noRd
.sct_models_consistent <- function(seurat_obj, assay = "SCT") {
  sct <- seurat_obj[[assay]]
  if (length(levels(sct)) <= 1) return(TRUE)
  observed <- vapply(SCTResults(sct, slot = "cell.attributes"),
                     function(x) stats::median(x[, "umi"]), numeric(1))
  stored <- tryCatch(unlist(SCTResults(sct, slot = "median_umi")), error = function(e) NULL)
  if (is.null(stored) || anyNA(observed)) return(FALSE)
  all(stored == min(observed))
}
