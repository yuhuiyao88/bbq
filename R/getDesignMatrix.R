#' Sustained step-design matrix (persistent shifts)
#'
#' Builds a lower-triangular, blocky design matrix where each column turns on
#' at a specific index and then stays on (1) through the end. Useful for modeling
#' *sustained* level shifts that persist once they occur (e.g., interventions,
#' regime changes, or structural breaks).
#'
#' @param n Integer \eqn{\ge 1}. Number of rows (time points or observations).
#' @param l Integer \eqn{\ge 1}. Nominal block length. Columns are created every
#'   \code{l} positions after the warm-up; the last block may be shorter to reach \code{n}.
#' @param w Integer \eqn{\ge 0}. Warm-up (initial zeros). No columns start before
#'   row \code{w + 1}. Defaults to \code{0}.
#'
#' @details
#' The number of columns is
#' \deqn{n_{\mathrm{col}} = \left\lfloor\frac{n - w}{l}\right\rfloor + \mathbf{1}\{(n-w) \bmod l > 0\}.}
#' For column \code{i}, rows \code{start_i:n} are set to 1, where
#' \code{start_1 = w + 1}, \code{start_{i+1} = start_i + l}.
#'
#' @return A numeric \code{n x ncol} matrix with 0/1 entries. Each column
#'   represents a sustained shift starting at its block’s first row.
#'
#' @examples
#' # 20 points, blocks of length 5, no warm-up:
#' S <- getSustainedShift(n = 20, l = 5)
#' dim(S)     # 20 x 4
#' S[1:8, ]   # first rows show the staggered "on" times
#'
#' # With warm-up of 3 rows:
#' S2 <- getSustainedShift(n = 20, l = 5, w = 3)
#'
#' @seealso \code{\link{getIsolatedShift}}
#' @export
getSustainedShift <- function(n, l, w = 0) {
  nrw <- n
  ncl <- as.integer((n - w) / l) + as.numeric((((n - w) %% l) > 0))# - 1
  out <- matrix(0, nrow = nrw, ncol = ncl)
  sta <- w
  end <- w
  for (i in 1:ncl) {
    sta <- end + 1
    end <- ifelse(i == ncl, nrw, sta + (l - 1))
    out[sta:nrw, i] <- 1 }
  out
}


#' Isolated step-design matrix (finite-length shifts)
#'
#' Builds a block-diagonal-like design matrix where each column is 1 for a
#' *finite* run of length \code{l} (except possibly the last, which may be shorter),
#' and 0 elsewhere. Useful for modeling *transient* or *windowed* effects
#' (e.g., short interventions, windowed basis functions).
#'
#' @param n Integer \eqn{\ge 1}. Number of rows (time points or observations).
#' @param l Integer \eqn{\ge 1}. Block length for each column. The last block
#'   may be shorter if it would exceed \code{n}.
#' @param w Integer \eqn{\ge 0}. Warm-up (initial zeros). No columns start before
#'   row \code{w + 1}. Defaults to \code{0}.
#'
#' @details
#' The number of columns matches \code{getSustainedShift()}:
#' \deqn{n_{\mathrm{col}} = \left\lfloor\frac{n - w}{l}\right\rfloor + \mathbf{1}\{(n-w) \bmod l > 0\}.}
#' For column \code{i}, rows \code{start_i:end_i} are set to 1, where
#' \code{start_1 = w + 1}, \code{start_{i+1} = start_i + l}, and
#' \code{end_i = min(start_i + l - 1, n)}.
#'
#' @return A numeric \code{n x ncol} matrix with 0/1 entries. Each column
#'   represents an isolated shift of length \code{l} (last column possibly shorter).
#'
#' @examples
#' # 20 points, blocks of length 5:
#' H <- getIsolatedShift(n = 20, l = 5)
#' dim(H)     # 20 x 4
#' colSums(H) # each column has about 5 ones
#'
#' # With warm-up of 3 rows:
#' H2 <- getIsolatedShift(n = 20, l = 5, w = 3)
#'
#' @seealso \code{\link{getSustainedShift}}
#' @export
getIsolatedShift <- function(n, l, w = 0) {
  nrw <- n
  ncl <- as.integer((n - w) / l) + as.numeric((((n - w) %% l) > 0))# - 1
  out <- matrix(0, nrow = nrw, ncol = ncl)
  sta <- w
  end <- w
  for (i in 1:ncl) {
    sta <- end + 1
    end <- ifelse(i == ncl, nrw, sta + (l - 1))
    out[sta:end, i] <- 1 }
  out
}


#' Gradual-drift design matrix (unbounded linear trends)
#'
#' Builds a design matrix where each column ramps up linearly from its onset
#' and continues to grow without bound. Useful for modeling *gradual drifts*
#' or *trending* effects that accumulate over time (e.g., wear, degradation,
#' or progressive regime changes).
#'
#' @param n Integer \eqn{\ge 1}. Number of rows (time points or observations).
#' @param l Integer \eqn{\ge 1}. Nominal block length. Columns are created every
#'   \code{l} positions after the warm-up; the last block may be shorter to reach \code{n}.
#' @param w Integer \eqn{\ge 0}. Warm-up (initial zeros). No columns start before
#'   row \code{w + 1}. Defaults to \code{0}.
#'
#' @details
#' The number of columns matches \code{getSustainedShift()} and
#' \code{getIsolatedShift()}:
#' \deqn{n_{\mathrm{col}} = \left\lfloor\frac{n - w}{l}\right\rfloor + \mathbf{1}\{(n-w) \bmod l > 0\}.}
#' For column \code{i}, rows \code{r >= start_i} are set to \code{r - start_i + 1}
#' (i.e. 1, 2, 3, \ldots), where \code{start_1 = w + 1},
#' \code{start_{i+1} = start_i + l}. Rows before \code{start_i} are zero.
#'
#' @return A numeric \code{n x ncol} matrix with non-negative entries. Each column
#'   represents a linear ramp starting at its block's first row.
#'
#' @examples
#' # 20 points, blocks of length 5, no warm-up:
#' D <- getGradualDrift(n = 20, l = 5)
#' dim(D)     # 20 x 4
#' D[1:8, ]   # first rows show staggered linear ramps
#'
#' # With warm-up of 3 rows:
#' D2 <- getGradualDrift(n = 20, l = 5, w = 3)
#'
#' @seealso \code{\link{getSustainedShift}}, \code{\link{getIsolatedShift}}
#' @noRd
getGradualDrift <- function(n, l, w = 0) {
  nrw <- n
  ncl <- as.integer((n - w) / l) + as.numeric((((n - w) %% l) > 0))
  out <- matrix(0, nrow = nrw, ncol = ncl)
  sta <- w
  end <- w
  for (i in 1:ncl) {
    sta <- end + 1
    end <- ifelse(i == ncl, nrw, sta + (l - 1))
    out[sta:nrw, i] <- seq_len(nrw - sta + 1)
  }
  out
}
