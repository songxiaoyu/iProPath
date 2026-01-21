#' Extract rows matching a given key from a list of data frames
#'
#' @param df_list list of data frames
#' @param match_key_g value of matching_key to extract
#' @return data frame
extract_data_g <- function(df_list, match_key_g) {

  # ---- guard: df_list can be NULL / empty ----
  if (is.null(df_list) || length(df_list) == 0L) {
    return(NULL)
  }

  # keep only non-null data frames
  df_list <- Filter(Negate(is.null), df_list)
  if (length(df_list) == 0L) return(NULL)

  # subset each df by matching_key
  out_list <- lapply(df_list, function(df) {
    df <- as.data.frame(df)  # handle tibble-like inputs safely

    if (!("matching_key" %in% names(df))) {
      stop("extract_data_g: column 'matching_key' not found in one element of df_list.")
    }
    if (nrow(df) == 0L) return(NULL)

    sub <- df[df$matching_key == match_key_g, , drop = FALSE]
    if (nrow(sub) == 0L) return(NULL)
    sub
  })

  out_list <- Filter(function(x) !is.null(x) && nrow(x) > 0L, out_list)
  if (length(out_list) == 0L) return(NULL)

  # bind rows (base)
  do.call(rbind, out_list)
}


#' Drop index columns and transpose data frame
#'
#' @param df A data frame containing covariates or mediator values.
#'   Can be \code{NULL} or have zero rows.
#'
#' @param drop_cols Character vector of column names to be removed
#'   before transposition. Defaults to \code{c("unique_index", "matching_key")}.
#'
#' @param prefix Optional character string used as a prefix for column
#'   names when \code{name_style = "prefix_col"}.
#'
#' @param name_style Character string specifying how column names are generated
#'   after transposition. One of:
#'   \describe{
#'     \item{\code{"keep"}}{Keep default column identifiers (row indices).}
#'     \item{\code{"prefix_col"}}{Prefix column identifiers with \code{prefix}.}
#'     \item{\code{"mediator_row_key"}}{Construct names using
#'       \code{mediator_tag}, row names, and \code{matching_key}.}
#'   }
#'
#' @param mediator_tag Character string used when
#'   \code{name_style = "mediator_row_key"} to label mediator columns.
#'
#' @param sep Character string used to separate components in constructed
#'   column names.
#'
#' @return A numeric matrix with transposed values and defensively
#'   constructed column names. If input is \code{NULL} or empty, a
#'   \code{0 x 0} matrix is returned.
#'
#' @export
dropindex_and_transpose <- function(
    df,
    drop_cols = c("unique_index", "matching_key"),
    prefix = NULL,
    name_style = c("keep", "prefix_col", "mediator_row_key"),
    mediator_tag = "Mediator",
    sep = "_"
) {
  name_style <- match.arg(name_style)

  # ---- guard: NULL / empty -> return NULL (not 0x0 matrix) ----
  if (is.null(df) || NROW(df) == 0L) {
    return(NULL)
  }

  df <- as.data.frame(df)

  keep <- setdiff(colnames(df), drop_cols)
  if (length(keep) == 0L) {
    return(NULL)
  }

  m <- as.matrix(df[, keep, drop = FALSE])
  out <- t(m)

  if (NROW(out) == 0L || NCOL(out) == 0L) {
    return(NULL)
  }

  row_id <- rownames(df)
  if (is.null(row_id) || all(row_id == "")) {
    row_id <- as.character(seq_len(nrow(df)))
  }

  key <- if ("matching_key" %in% names(df)) as.character(df$matching_key) else rep("", nrow(df))

  # base-R NULL-coalesce: if prefix is NULL, treat as ""
  pref <- if (is.null(prefix)) "" else as.character(prefix)

  cn <- switch(
    name_style,
    keep = row_id,
    prefix_col = paste0(pref, row_id),
    mediator_row_key = paste(mediator_tag, row_id, key, sep = sep)
  )

  # fallback
  if (length(cn) != ncol(out)) {
    cn <- as.character(seq_len(ncol(out)))
    if (!is.null(prefix)) cn <- paste0(pref, cn)
  }

  colnames(out) <- cn
  out
}



#' Time a code block and print elapsed time
#'
#' A lightweight helper to measure and report the running time of an
#' expression using \code{message()}. Mainly intended for internal
#' progress/diagnostic messages in package functions.
#'
#' @param label A short character string used as the message label.
#' @param expr An expression to be evaluated.
#'
#' @return The result of evaluating \code{expr}.
#'
#' @keywords internal
.time_block <- function(label, expr) {
  t0 <- Sys.time()
  out <- force(expr)
  message(sprintf(
    "[%s] %.2f min",
    label,
    as.numeric(difftime(Sys.time(), t0, units = "mins"))
  ))
  out
}


