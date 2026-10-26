# =============================================================================
# BBQ visualization (ggplot2)
# =============================================================================
# Graph types:
#   (1) plotQuantileProcess() - data with fitted quantile bands over time; the
#       localized change-points as circles, colored by source when a comparator
#       method's change-points are supplied; no block-onset rules by default
#   (2) plotQSSProcess() / plotLmomProcess() - shape profiles over time, bands only
#   (3) plotGammaHeatmap()    - block-shift coefficient diagnosis (heatmap)
#   (4) plotBBQSummary()      - (1) left, (2) over (3) right: the talk's figure
# A Date `time` vector gets yearly breaks with rotated labels.
# ggplot2 is an optional (Suggests) dependency; each function checks for it at
# call time so the package still loads and fits without ggplot2 installed.

# UA / crimson palette shared across the three graphs
.bbq_pal <- list(ink = "#2A2123", steel = "#5B666D", crimson = "#9E1B32",
                 brick = "#76232F", rose = "#C46A78", gray = "#9AA5B1",
                 # Muted, low-chroma Morandi-style diverging pair for the heatmap:
                 # positive shifts read red, negative shifts read blue. There is no
                 # canonical hex for "Morandi red/blue"; these are representative
                 # tones, matched in lightness, and overridable per call.
                 morandi_red = "#AD6A6C", morandi_blue = "#6E8CA0")

.bbq_need_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("ggplot2 is required for BBQ plots. Install it with install.packages(\"ggplot2\").",
         call. = FALSE)
  }
}

.bbq_theme <- function() {
  ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      strip.placement = "outside",
      strip.text.y.left = ggplot2::element_text(angle = 0)
    )
}

# Change-point marks by source (JSM 2026 talk): proposed method, comparator, both.
.bbq_cp_pal <- c(proposed = "#C62828", comparator = "#1565C0", both = "#6A1B9A")

# Time axis: a Date vector gets yearly breaks with rotated labels; any other x
# type is left to ggplot2's defaults.
.bbq_time_axis <- function(time, date_breaks = "12 months", date_labels = "%Y") {
  if (!inherits(time, "Date")) return(list())
  list(ggplot2::scale_x_date(date_breaks = date_breaks, date_labels = date_labels),
       ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)))
}

# Legend drawn inside the panel, anchored at `position` (npc), as in the talk.
.bbq_inset_legend <- function(position = c(0.995, 0.995)) {
  th <- ggplot2::theme(
    legend.justification = c(1, 1),
    legend.background = ggplot2::element_rect(fill = grDevices::adjustcolor("white", 0.78), colour = NA),
    legend.key = ggplot2::element_blank(),
    legend.text = ggplot2::element_text(size = 9),
    legend.key.size = ggplot2::unit(0.42, "cm"),
    legend.margin = ggplot2::margin(2, 4, 2, 4))
  if (as.package_version(getNamespaceVersion("ggplot2")) >= "3.5.0")
    th + ggplot2::theme(legend.position = "inside", legend.position.inside = position)
  else th + ggplot2::theme(legend.position = position)
}

# Block length recorded in a detection object (0 when it cannot be inferred).
.bbq_block_length <- function(detection) {
  os <- detection$detected_blocks$obs_start
  if (length(os) >= 2L) as.integer(stats::median(diff(os))) else 0L
}

# Quantile levels: prefer an explicit argument, else recover from the fit.
.bbq_taus <- function(fit, taus) {
  if (!is.null(taus)) return(taus)
  for (nm in c("taus", "tau_q", "tau")) if (!is.null(fit[[nm]])) return(fit[[nm]])
  stop("Quantile levels not found in fit; please pass `taus`.", call. = FALSE)
}

# Point-estimate fitted quantiles, n x m, reconstructed as X %*% beta + H %*% gamma.
# Uses the coherent point estimate (MAP mode under MAP, posterior median under
# MCMC) via .bbq_point_eta / .bbq_coefs, so the plotted central line matches the
# estimator used for detection and localization.
.bbq_fitted_quantiles <- function(fit, taus) {
  if (is.null(fit$map$par))
    stop("BBQ plots require a fit with fit$map$par (MAP mode or MCMC posterior median).",
         call. = FALSE)
  t(.bbq_point_eta(fit, taus))   # m x n -> n x m
}

# Significant blocks + their onset/localized observations from a
# detectChangepoints_gamma() result (block significant if any quantile cell is
# Holm-significant, matching the worked demo).
# Significant blocks + their onsets and localized change-points, resolved from the
# family/statistic/adjustment actually recorded in `detection`. Reading a hardcoded
# member (as this once did) silently returns nothing whenever the caller ran a basis
# other than "quantile", which is why a QSS-only detection drew no lines at all.
.bbq_sig_blocks <- function(detection, basis = NULL) {
  empty <- list(blocks = integer(0), onset = numeric(0), located = numeric(0))
  db <- detection$detected_blocks
  if (is.null(db)) return(empty)

  fam <- if (!is.null(detection$basis)) detection$basis else "quantile"
  # `basis` lets a caller pin the family so a plot marks the blocks flagged by
  # the SAME test its panel displays. Without it the old precedence applies and
  # quantile wins whenever it was run, which silently made an L-moment or QSS
  # panel carry quantile-basis change points.
  if (!is.null(basis) && length(basis) == 1L &&
      !is.null(detection$tests) && !is.null(detection$tests[[basis]])) {
    fam <- basis
  } else {
    fam <- if ("quantile" %in% fam) "quantile" else "qss"   # quantile wins if both were run
  }
  stat <- if (!is.null(detection$statistic)) detection$statistic else "ui"
  stat_name <- if (("hotelling_t2" %in% stat) && !("ui" %in% stat)) "hotelling_t2" else "ui"
  adjust <- if (!is.null(detection$adjust)) detection$adjust else "raw"

  sb <- detection$tests[[fam]][[stat_name]][[adjust]]
  if (is.null(sb)) {                                      # older detection objects
    nm <- if (fam == "qss") {
      if (stat_name == "hotelling_t2") paste0("significant_qss_t2_", adjust)
      else                             paste0("significant_qss_", adjust)
    } else {
      if (stat_name == "hotelling_t2") paste0("significant_wald_", adjust)
      else                             paste0("significant_", adjust)
    }
    sb <- detection[[nm]]
    if (is.null(sb) && !is.null(db[[nm]])) sb <- which(db[[nm]])
  }
  if (is.null(sb) || length(sb) == 0) return(empty)
  list(blocks = sb, onset = db$obs_start[sb], located = db$signal_obs[sb])
}


#' Plot the data process with fitted quantile bands over time
#'
#' Three-quantile fits show the median and IQR band, without outer-tail lines.
#'
#' Graph type 1: the observations with the five fitted quantile bands and, from a
#' \code{detectChangepoints_gamma()} result, the localized change-points as
#' circles on the series. Block-onset rules are off by default. When
#' \code{comparator} gives the change-points of another method, the circles are
#' colored by source: the proposed method, the comparator, or both when the two
#' fall within \code{match_tol} observations of each other, with a legend inside
#' the panel.
#'
#' @param fit A MAP fit from \code{getModel()}.
#' @param time Optional x-axis vector (default \code{seq_len(n)}). A \code{Date}
#'   vector gets yearly breaks with rotated labels (see \code{date_breaks}).
#' @param center,scale Map the fit-scale quantiles/data back to the display scale
#'   as \code{value * scale + center} (e.g. the standardization used before fitting).
#' @param detection Optional \code{detectChangepoints_gamma()} result. Its recorded
#'   \code{basis}, \code{statistic} and \code{adjust} select which blocks are
#'   flagged, and its \code{detected_blocks} supplies the block onset and the
#'   localized change-point \code{signal_obs}, the latter obtained under whichever
#'   \code{signal_position} was passed to \code{detectChangepoints_gamma()}.
#' @param show_onset,show_located Logical display toggles: draw the block-onset
#'   rules (default \code{FALSE}) and the localized change-points (default
#'   \code{TRUE}) recorded in \code{detection}. They only hide layers; the
#'   decisions themselves are made (and recorded) by
#'   \code{detectChangepoints_gamma()}.
#' @param basis Which detection family supplies the change-point marks. Default
#'   \code{NULL}: the family recorded in \code{detection}, quantile first.
#' @param comparator Optional integer vector of change-point positions
#'   (observation indices) from a comparator method, e.g.
#'   \code{changepoint::cpts()} of a binary-segmentation fit.
#' @param comparator_label,proposed_label Legend labels of the comparator and of
#'   the proposed method.
#' @param match_tol A proposed and a comparator change-point within this many
#'   observations of each other are drawn once, at the proposed location, as
#'   "Both". Default \code{NULL}: the block length recorded in \code{detection},
#'   since BBQ resolves a change only to its block.
#' @param legend_position Anchor of the inset legend (npc coordinates of its
#'   top-right corner); used only when \code{comparator} is supplied.
#' @param date_breaks,date_labels Breaks and label format of the x axis when
#'   \code{time} is a \code{Date} vector.
#' @param title Optional plot title.
#' @param xlab Label for the x axis (default \code{"time"}).
#' @param ylab Label for the y axis (default \code{"value"}). Set it to name the
#'   series actually plotted, e.g. \code{"residual"} when the model was fitted to
#'   regression residuals, or \code{"difference"} when it was fitted to their
#'   first-order differences. Use \code{NULL} to drop the label.
#' @return A ggplot object.
#' @export
plotQuantileProcess <- function(fit, time = NULL, center = 0, scale = 1,
                                detection = NULL, title = NULL,
                                xlab = "time", ylab = "value",
                                show_onset = FALSE, show_located = TRUE, basis = NULL,
                                comparator = NULL, comparator_label = "Comparator",
                                proposed_label = "Proposed (BBQ)", match_tol = NULL,
                                legend_position = c(0.995, 0.995),
                                date_breaks = "12 months", date_labels = "%Y") {
  .bbq_need_ggplot2()
  pal <- .bbq_pal
  taus <- if (!is.null(detection) && !is.null(detection$taus)) detection$taus
          else .bbq_taus(fit, NULL)
  m <- length(taus)
  q <- .bbq_fitted_quantiles(fit, taus) * scale + center
  n <- nrow(q)
  yv <- fit$y * scale + center
  if (is.null(time)) time <- seq_len(n)
  if (m == 3L) {
    if (!isTRUE(all.equal(as.numeric(taus), c(0.25, 0.5, 0.75))))
      stop("Three-quantile process plots require taus = c(0.25, 0.5, 0.75).", call. = FALSE)
    df <- data.frame(time = time, y = yv, q1 = q[, 1L], med = q[, 2L], q3 = q[, 3L])
    p <- ggplot2::ggplot(df, ggplot2::aes(x = time)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = q1, ymax = q3), fill = pal$crimson, alpha = 0.17) +
      ggplot2::geom_point(ggplot2::aes(y = y), color = pal$gray, alpha = 0.5, size = 1) +
      ggplot2::geom_line(ggplot2::aes(y = q1), color = pal$steel, linewidth = 0.5) +
      ggplot2::geom_line(ggplot2::aes(y = q3), color = pal$steel, linewidth = 0.5) +
      ggplot2::geom_line(ggplot2::aes(y = med), color = pal$ink, linewidth = 0.9)
  } else {
  df <- data.frame(time = time, y = yv,
                   lo = q[, 1], q1 = q[, 2], med = q[, ceiling(m / 2)], q3 = q[, m - 1], hi = q[, m])
  p <- ggplot2::ggplot(df, ggplot2::aes(x = time)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi), fill = pal$gray, alpha = 0.10) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = q1, ymax = q3), fill = pal$crimson, alpha = 0.17) +
    ggplot2::geom_point(ggplot2::aes(y = y), color = pal$gray, alpha = 0.5, size = 1) +
    ggplot2::geom_line(ggplot2::aes(y = lo), color = pal$brick, linetype = "dashed", linewidth = 0.4) +
    ggplot2::geom_line(ggplot2::aes(y = hi), color = pal$brick, linetype = "dashed", linewidth = 0.4) +
    ggplot2::geom_line(ggplot2::aes(y = q1), color = pal$steel, linewidth = 0.5) +
    ggplot2::geom_line(ggplot2::aes(y = q3), color = pal$steel, linewidth = 0.5) +
    ggplot2::geom_line(ggplot2::aes(y = med), color = pal$ink, linewidth = 0.9)
  }

  # localized change-points of the proposed method (from `detection`)
  lp <- integer(0)
  if (!is.null(detection)) {
    loc <- .bbq_sig_blocks(detection, basis = basis)
    if (isTRUE(show_onset) && length(loc$onset) > 0) {
      p <- p + ggplot2::geom_vline(xintercept = time[loc$onset], color = pal$crimson,
                                   linewidth = 0.6, alpha = 0.85)
    }
    lp <- loc$located[!is.na(loc$located) & loc$located >= 1 & loc$located <= n]
  }
  cp <- if (is.null(comparator)) integer(0) else as.integer(comparator)
  cp <- cp[!is.na(cp) & cp >= 1 & cp <= n]

  if (isTRUE(show_located)) {
    if (is.null(comparator)) {
      if (length(lp) > 0) {
        p <- p + ggplot2::geom_point(
          data = data.frame(x = time[lp], y = yv[lp]),
          ggplot2::aes(x = x, y = y), shape = 21, fill = pal$crimson, color = "black",
          size = 2.6, stroke = 0.8)
      }
    } else {
      # Comparator overlay (talk_figures_lmom.R): a proposed and a comparator point
      # within `match_tol` observations are one detection, drawn once at the
      # proposed location. All three legend keys are shown even when a source is
      # empty, so figures of different series read the same.
      tol <- if (!is.null(match_tol)) as.numeric(match_tol)
             else if (!is.null(detection)) .bbq_block_length(detection) else 0
      p_hit <- vapply(lp, function(b) any(abs(cp - b) <= tol), logical(1))
      c_hit <- vapply(cp, function(b) any(abs(lp - b) <= tol), logical(1))
      lv <- c(proposed_label, comparator_label, "Both")
      mk <- function(v, lab) data.frame(idx = as.integer(v), src = rep(lab, length(v)))
      d <- rbind(mk(lp[!p_hit], lv[1]), mk(cp[!c_hit], lv[2]), mk(lp[p_hit], lv[3]))
      miss <- setdiff(lv, unique(d$src))
      if (length(miss)) d <- rbind(d, data.frame(idx = NA_integer_, src = miss))
      d$src <- factor(d$src, levels = lv)
      d$x <- time[d$idx]; d$yv <- yv[d$idx]
      cols <- stats::setNames(unname(.bbq_cp_pal[c("proposed", "comparator", "both")]), lv)
      p <- p +
        ggplot2::geom_point(data = d, ggplot2::aes(x = x, y = yv, fill = src),
                            shape = 21, colour = "black", size = 2.6, stroke = 0.7,
                            inherit.aes = FALSE, na.rm = TRUE) +
        ggplot2::scale_fill_manual(values = cols, name = NULL, drop = FALSE) +
        ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 3.2)))
    }
  }
  p <- p + ggplot2::labs(x = xlab, y = ylab, title = title) + .bbq_theme() +
    .bbq_time_axis(time, date_breaks, date_labels)
  if (!is.null(comparator) && isTRUE(show_located)) p <- p + .bbq_inset_legend(legend_position)
  p
}


#' Plot the quantile-shape-statistic (QSS) process over time
#'
#' With taus = c(0.25, 0.5, 0.75), only Location and Scale panels are shown.
#'
#' Graph type 2: posterior Location, Scale, Skewness and Kurtosis over time, each
#' as a median line with a credible-band ribbon. Posterior draws are sorted within
#' each draw (non-crossing) before the shape statistics are formed, so the ratios
#' stay well defined.
#'
#' @param fit A MAP fit from \code{getModel()}.
#' @param eta Optional posterior predictive-quantile array from \code{getEta()};
#'   computed internally if not supplied.
#' @param H,X Optional design matrices passed to \code{getEta()} (default from fit).
#' @param time Optional x-axis vector (default \code{seq_len(n)}). A \code{Date}
#'   vector gets yearly breaks with rotated labels.
#' @param center,scale Map the fit-scale quantiles to the display scale.
#' @param level Credible-band level (default 0.95).
#' @param detection Optional \code{detectChangepoints_gamma()} result. Its recorded
#'   \code{basis} and \code{statistic} select which block test flags the blocks, and
#'   its \code{detected_blocks} supplies both the block onset (dashed line) and the
#'   localized change-point \code{signal_obs} (solid line), the latter obtained under
#'   whichever \code{signal_position} was passed to \code{detectChangepoints_gamma()}.
#' @param show_onset,show_located Logical display toggles (both default
#'   \code{FALSE}: the profile is drawn as bands only, as in the JSM 2026 talk):
#'   draw the block-onset rules and the localized change-point rules recorded in
#'   \code{detection}. They only hide layers; the decisions themselves are
#'   made (and recorded) by \code{detectChangepoints_gamma()}.
#' @param basis Which detection family supplies the change-point marks (default
#'   \code{NULL}: the family recorded in \code{detection}, quantile first).
#' @param seed Optional seed for \code{getEta()}.
#' @param title Optional plot title.
#' @param xlab Label for the x axis (default \code{"time"}).
#' @param ylab Label for the shared y axis. Default \code{NULL} (no label), since
#'   the panels are already named by their facet strips and each has its own
#'   free scale.
#' @param date_breaks,date_labels Breaks and label format of the x axis when
#'   \code{time} is a \code{Date} vector.
#' @return A ggplot object (two or four stacked, free-y facets).
#' @export
plotQSSProcess <- function(fit, eta = NULL, H = NULL, X = NULL, time = NULL,
                           center = 0, scale = 1, level = 0.95, detection = NULL,
                           seed = NULL, title = NULL,
                           xlab = "time", ylab = NULL,
                           show_onset = FALSE, show_located = FALSE, basis = NULL,
                           date_breaks = "12 months", date_labels = "%Y") {
  .bbq_need_ggplot2()
  pal <- .bbq_pal
  taus <- if (!is.null(detection) && !is.null(detection$taus)) detection$taus
          else .bbq_taus(fit, NULL)
  if (is.null(eta)) eta <- getEta(fit, H = H, X = X, seed = seed)
  eta <- eta * scale + center
  eta <- aperm(apply(eta, c(1, 3), sort), c(2, 1, 3))  # enforce non-crossing per draw
  qss <- getQSS(eta, taus = taus)                      # [iters, 4, n]
  a <- (1 - level) / 2
  qmid <- apply(qss, c(2, 3), stats::median, na.rm = TRUE)
  qlo  <- apply(qss, c(2, 3), stats::quantile, probs = a, na.rm = TRUE)
  qhi  <- apply(qss, c(2, 3), stats::quantile, probs = 1 - a, na.rm = TRUE)
  n <- dim(qss)[3]
  if (is.null(time)) time <- seq_len(n)
  labs4 <- dimnames(qss)[[2]]
  df <- do.call(rbind, lapply(seq_along(labs4), function(k) data.frame(
    time = time, stat = factor(labs4[k], levels = labs4),
    mid = qmid[k, ], lo = qlo[k, ], hi = qhi[k, ])))
  cols <- c(Location = pal$ink, Scale = pal$crimson, Skewness = pal$rose, Kurtosis = pal$brick)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = time))
  if (!is.null(detection)) {
    loc <- .bbq_sig_blocks(detection, basis = basis)
    ob <- loc$onset[!is.na(loc$onset) & loc$onset >= 1 & loc$onset <= n]
    if (isTRUE(show_onset) && length(ob) > 0)
      p <- p + ggplot2::geom_vline(xintercept = time[ob], color = pal$crimson,
                                   linetype = "dashed", linewidth = 0.35, alpha = 0.55)
    # localized change-point within each OOC block (Eq. 27): solid, so it is
    # distinguishable from the dashed block onset.
    lp <- loc$located[!is.na(loc$located) & loc$located >= 1 & loc$located <= n]
    if (isTRUE(show_located) && length(lp) > 0)
      p <- p + ggplot2::geom_vline(xintercept = time[lp], color = pal$crimson,
                                   linewidth = 0.6, alpha = 0.9)
  }
  p +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi, fill = stat), alpha = 0.22) +
    ggplot2::geom_line(ggplot2::aes(y = mid, color = stat), linewidth = 0.8) +
    ggplot2::scale_fill_manual(values = cols, guide = "none") +
    ggplot2::scale_color_manual(values = cols, guide = "none") +
    ggplot2::facet_wrap(~stat, ncol = 1, scales = "free_y", strip.position = "left") +
    ggplot2::labs(x = xlab, y = ylab, title = title) + .bbq_theme() +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "grey55", fill = NA, linewidth = 0.5),
      panel.spacing.y = ggplot2::unit(0.6, "lines")) +
    .bbq_time_axis(time, date_breaks, date_labels)
}

#' L-moment shape profile over time
#'
#' Three-quantile fits show only L-location and L-scale.
#'
#' The L-moment counterpart of \code{\link{plotQSSProcess}}: posterior median and
#' credible band for each of the four approximate L-moments over time, computed by
#' \code{\link{getLmom}} from the fitted quantiles.
#'
#' Change-point marks (off by default) use \code{basis = "lmom"}, so when they are
#' shown the vertical rules come from the SAME L-moment UI test that flags the
#' L-moment heatmap panel: an L-moment profile carrying quantile-basis change
#' points would be internally inconsistent.
#'
#' @inheritParams plotQSSProcess
#' @param basis Which detection family supplies the change-point marks. Defaults to
#'   \code{"lmom"}; pass another family name to override, or \code{NULL} for the
#'   package's legacy precedence (quantile wins when it was run).
#'
#' @note The panels are the L-moments \eqn{\lambda_r}, not the scale-free ratios
#'   \eqn{\tau_r}; see \code{\link{getLmom}}.
#'
#' @noRd
plotLmomProcess <- function(fit, eta = NULL, H = NULL, X = NULL, time = NULL,
                            center = 0, scale = 1, level = 0.95, detection = NULL,
                            seed = NULL, title = NULL,
                            xlab = "time", ylab = NULL,
                            show_onset = FALSE, show_located = FALSE,
                            basis = "lmom",
                            date_breaks = "12 months", date_labels = "%Y") {
  .bbq_need_ggplot2()
  pal <- .bbq_pal
  taus <- if (!is.null(detection) && !is.null(detection$taus)) detection$taus
          else .bbq_taus(fit, NULL)
  if (is.null(eta)) eta <- getEta(fit, H = H, X = X, seed = seed)
  eta <- eta * scale + center
  eta <- aperm(apply(eta, c(1, 3), sort), c(2, 1, 3))   # non-crossing per draw
  lm4 <- getLmom(eta, taus = taus)                      # [iters, 4, n]
  a <- (1 - level) / 2
  mid <- apply(lm4, c(2, 3), stats::median, na.rm = TRUE)
  lo  <- apply(lm4, c(2, 3), stats::quantile, probs = a, na.rm = TRUE)
  hi  <- apply(lm4, c(2, 3), stats::quantile, probs = 1 - a, na.rm = TRUE)
  n <- dim(lm4)[3]
  if (is.null(time)) time <- seq_len(n)
  labs4 <- dimnames(lm4)[[2]]
  df <- do.call(rbind, lapply(seq_along(labs4), function(k) data.frame(
    time = time, stat = factor(labs4[k], levels = labs4),
    mid = mid[k, ], lo = lo[k, ], hi = hi[k, ])))
  cols <- stats::setNames(c(pal$ink, pal$crimson, pal$rose, pal$brick)[seq_along(labs4)], labs4)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = time))
  if (!is.null(detection)) {
    loc <- .bbq_sig_blocks(detection, basis = basis)
    ob <- loc$onset[!is.na(loc$onset) & loc$onset >= 1 & loc$onset <= n]
    if (isTRUE(show_onset) && length(ob) > 0)
      p <- p + ggplot2::geom_vline(xintercept = time[ob], color = pal$crimson,
                                   linetype = "dashed", linewidth = 0.35, alpha = 0.55)
    lp <- loc$located[!is.na(loc$located) & loc$located >= 1 & loc$located <= n]
    if (isTRUE(show_located) && length(lp) > 0)
      p <- p + ggplot2::geom_vline(xintercept = time[lp], color = pal$crimson,
                                   linewidth = 0.6, alpha = 0.9)
  }
  p +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi, fill = stat), alpha = 0.22) +
    ggplot2::geom_line(ggplot2::aes(y = mid, color = stat), linewidth = 0.8) +
    ggplot2::scale_fill_manual(values = cols, guide = "none") +
    ggplot2::scale_color_manual(values = cols, guide = "none") +
    ggplot2::facet_wrap(~stat, ncol = 1, scales = "free_y", strip.position = "left") +
    ggplot2::labs(x = xlab, y = ylab, title = title) + .bbq_theme() +
    ggplot2::theme(
      panel.border = ggplot2::element_rect(color = "grey55", fill = NA, linewidth = 0.5),
      panel.spacing.y = ggplot2::unit(0.6, "lines")) +
    .bbq_time_axis(time, date_breaks, date_labels)
}


#' Plot block-shift coefficient diagnosis (heatmap)
#'
#' Graph type 3: block-shift heatmap(s) with black borders on the significant
#' cells. The function auto-detects, from the \code{detection} object, which test
#' family/families were run and shows the matching panel(s): a \strong{quantile}
#' panel (the quantile shift coefficients as a quantile-by-block map) and/or a \strong{QSS}
#' panel (the four QSS shift coefficients LS, ScS, SkS, KS by block). If only one
#' family was computed, only that panel is drawn; if both, both are stacked (via
#' \pkg{patchwork} when available). Borders follow the statistic recorded in
#' \code{detection}: per-cell where \eqn{|z|} exceeds the calibrated cell-max
#' threshold, or whole significant columns under the Hotelling \eqn{T^2}.
#'
#' @param fit A MAP fit from \code{getModel()}.
#' @param detection Optional \code{detectChangepoints_gamma()} result. Its
#'   \code{basis}/\code{statistic} fields drive which panels appear and how cells
#'   are bordered. Older results without those fields fall back to a quantile panel
#'   (Holm-bordered), plus a QSS panel if \code{z_qss} is present.
#' @param basis Optional subset of the families recorded in \code{detection} to
#'   draw, e.g. \code{"lmom"} for the L-moment panel alone (the talk's figure).
#'   Default \code{NULL}: every family the detection object carries.
#' @param block_labels Optional labels for the block (x) axis (default block index).
#'   Must be unique (they become factor levels). When supplied, the labels are
#'   rotated.
#' @param label_every Show every \code{label_every}-th block label on the x axis.
#'   Default \code{NULL}: all labels up to 12 blocks, about eight labels beyond.
#' @param note_clipping Logical; if \code{TRUE}, append to the panel label how many cells
#'   exceed \code{z_limit} and the largest |z| (default \code{FALSE}: the panel label is the
#'   basis name only, e.g. "Quantile Shift Coefficient"; author's rule, 2026-09-10).
#' @param title Optional plot title.
#' @param mark_cells Logical; when \code{TRUE} (default) the cells responsible for a
#'   OOC block are given a second, darker border. A block that the block-level
#'   rule flags is bordered in \code{block_color}; inside it, each cell whose
#'   cell-level posterior probability
#'   \eqn{p_{q,j} = 1 - P(\chi^2_1 \le \tilde z^2_{q,j})} falls below \code{alpha}
#'   is bordered in \code{cell_color}. This is the within-block localization step of
#'   Section 3.2: it runs only inside blocks that have already signaled, so it
#'   localizes a detected shift rather than adding a new family of tests, and the
#'   false alarm probability remains that of the block-level rule.
#' @param block_color Border color for cells of a OOC block (default grey).
#' @param cell_color Border color for the localized cells within a OOC block
#'   (default black).
#' @param xlab Label for the x axis (default \code{"block"}).
#' @param ylab Label for the y axis. Default \code{NULL} (no label), since the rows
#'   are already named by the quantile levels or the QSS contrasts. When both panels
#'   are stacked, \code{xlab}/\code{ylab} apply to each panel; \code{title} is
#'   applied once to the combined figure.
#' @param z_limit Positive scalar fixing the diverging fill scale to
#'   \code{c(-z_limit, z_limit)} whenever the fill is on a z scale -- the whitened
#'   \eqn{\tilde z} cells, or the studentized \code{z} fallback (default 3). A FIXED
#'   scale is the point: with a data-driven limit the same color means a different
#'   number in every figure, so two fits cannot be compared by eye. Cells beyond the
#'   limit are clipped to the end color (not dropped); the legend therefore labels its
#'   ends \code{"<= -z_limit"} and \code{">= z_limit"}, and the subtitle records that
#'   clipping occurred and how far out the extreme cell was. Set \code{NULL} to
#'   restore the old data-driven symmetric limit. Ignored when the fill is the raw
#'   posterior-mean \eqn{\gamma}, where a fixed \eqn{\pm 3} would be meaningless --
#'   those panels always use a data-driven limit.
#' @param pos_color,neg_color End colors of the diverging fill: positive
#'   \eqn{\tilde z} and negative \eqn{\tilde z} respectively, white at zero.
#'   Defaults are a muted, lightness-matched Morandi-style red and blue. There is
#'   no canonical hex for those, so pass your own values to match a house palette.
#' @return A ggplot object when one family is shown; a \pkg{patchwork} of two panels
#'   when both are shown (or a named list of ggplots if \pkg{patchwork} is absent).
#' @export
plotGammaHeatmap <- function(fit, detection = NULL, block_labels = NULL,
                             title = NULL, mark_cells = TRUE,
                             block_color = NULL, cell_color = "black",
                             pos_color = NULL, neg_color = NULL,
                             z_limit = 3,
                             xlab = "block", ylab = NULL,
                             basis = NULL, label_every = NULL, note_clipping = FALSE) {
  .bbq_need_ggplot2()
  pal <- .bbq_pal
  taus <- if (!is.null(detection) && !is.null(detection$taus)) detection$taus
          else .bbq_taus(fit, NULL)
  m <- length(taus)
  r <- if (!is.null(fit$H)) ncol(fit$H) else 0L
  if (r == 0) stop("No block-shift design (fit$H has no columns).", call. = FALSE)
  blk <- if (!is.null(block_labels)) block_labels else seq_len(r)
  if (length(blk) != r || anyDuplicated(blk))
    stop("block_labels must be ", r, " unique labels.", call. = FALSE)
  every <- if (!is.null(label_every)) max(1L, as.integer(label_every))
           else if (r <= 12L) 1L else ceiling(r / 8)
  x_breaks <- as.character(blk)[seq(1L, r, by = every)]

  ## ---- everything below renders the decisions RECORDED in `detection`: which
  ## families (basis), which block statistic, and which across-block rule
  ## (adjust) were chosen when detectChangepoints_gamma() was run. The plot has
  ## no decision arguments of its own. ----
  ALL_FAMS <- c("quantile", "qss", "lmom")
  fams <- if (is.null(detection)) "quantile"
          else if (!is.null(detection$basis)) detection$basis
          else c("quantile", if (!is.null(detection$z_qss)) "qss")
  fams <- intersect(ALL_FAMS, fams)
  if (!is.null(basis)) {
    bad <- setdiff(basis, fams)
    if (length(bad))
      stop("basis '", paste(bad, collapse = "', '"), "' was not run in this detection object.",
           call. = FALSE)
    fams <- intersect(fams, basis)
  }
  # keep only families the detection object actually carries cells for
  have <- function(f) switch(f,
    quantile = TRUE,
    qss      = !is.null(detection$z_qss),
    lmom     = !is.null(detection$z_lmom),
    FALSE)
  if (!is.null(detection)) fams <- fams[vapply(fams, have, logical(1))]
  if (length(fams) == 0) fams <- "quantile"
  stat <- if (!is.null(detection) && !is.null(detection$statistic)) detection$statistic else "ui"
  use_t2 <- ("hotelling_t2" %in% stat) && !("ui" %in% stat)   # UI wins if both were run
  stat_name <- if (use_t2) "hotelling_t2" else "ui"
  adjust <- if (!is.null(detection) && !is.null(detection$adjust)) detection$adjust else "raw"

  ## ---- OOC blocks under the recorded rule; older detection objects fall
  ## back to the flat calibrated aliases (calib only). ----
  get_sig <- function(fam) {
    if (is.null(detection)) return(integer(0))
    fam_res <- detection$tests[[fam]]
    if (!is.null(fam_res) && !is.null(fam_res[[stat_name]]) &&
        !is.null(fam_res[[stat_name]][[adjust]])) {
      return(fam_res[[stat_name]][[adjust]])
    }
    if (adjust != "calib") {
      warning("this detection object does not carry the '", adjust,
              "' adjustment member; using the calibrated flags.", call. = FALSE)
    }
    if (fam == "qss") {
      if (use_t2) detection$significant_qss_t2_calib else detection$significant_qss_calib
    } else {
      if (use_t2) detection$significant_wald_calib else detection$significant_calib
    }
  }

  ## ---- within-block localization (manuscript Sec 3.2): a OOC block is
  ## bordered in block_color; inside it, a cell is bordered in cell_color when
  ## its whitened statistic exceeds the SAME charting constant that flagged the
  ## block -- Eqs. (21)/(23) for the raw rule, (22)/(24) for the calibrated
  ## rule, and Eq. (25)'s adjusted constants under Bonferroni/Holm/BH. The flags
  ## are computed by detectChangepoints_gamma() ($cells); a cell exceedance
  ## implies the block exceedance, so cells localize the block decisions without
  ## adding a hypothesis family. ----
  blk_col   <- if (!is.null(block_color)) block_color else pal$gray
  pos_color <- if (!is.null(pos_color))   pos_color   else pal$morandi_red
  neg_color <- if (!is.null(neg_color))   neg_color   else pal$morandi_blue
  get_sig_cells <- function(fam, sig_cols, nr) {
    out <- matrix(FALSE, nr, r)
    if (!isTRUE(mark_cells) || is.null(detection) || length(sig_cols) == 0) return(out)
    cm <- detection$tests[[fam]][[stat_name]]$cells[[adjust]]
    if (is.null(cm) || !is.matrix(cm) || nrow(cm) != nr || ncol(cm) != r) {
      warning("cell-level flags are unavailable for the '", fam,
              "' family (older detection object?); blocks are bordered but ",
              "cells are not localized.", call. = FALSE)
      return(out)
    }
    cm
  }

  ## ---- single-panel builder; sig_cols = significant blocks (whole-column border),
  ## sig_cells = logical ncell x r matrix of localized cells within those blocks ----
  heat <- function(vals, rowlab, sig_cols, sig_cells, fill_lab, subtitle, diverging,
                   fixed_lim = NULL) {
    d <- expand.grid(ri = seq_len(nrow(vals)), bj = seq_len(r))
    d$val   <- vals[cbind(d$ri, d$bj)]
    d$sig   <- d$bj %in% sig_cols
    d$cell  <- sig_cells[cbind(d$ri, d$bj)]
    d$row   <- factor(rowlab[d$ri], levels = rowlab)
    d$block <- factor(blk[d$bj], levels = blk)

    # A FIXED fill scale is what makes two figures comparable; a data-driven one
    # silently redefines what a color means. Values beyond the fixed limit are
    # clipped to the end color rather than dropped (ggplot renders out-of-limits
    # as grey NA, which would read as "missing" instead of "extreme"), and the
    # clipping is disclosed in the subtitle.
    obs_max <- suppressWarnings(max(abs(vals), na.rm = TRUE))
    clipped <- 0L
    if (diverging && !is.null(fixed_lim) && is.finite(fixed_lim) && fixed_lim > 0) {
      lim <- fixed_lim
      clipped <- sum(abs(d$val) > lim, na.rm = TRUE)
      d$fill_val <- pmin(pmax(d$val, -lim), lim)
    } else {
      lim <- obs_max; if (!is.finite(lim) || lim == 0) lim <- 1
      d$fill_val <- d$val
    }
    if (clipped > 0L && isTRUE(note_clipping)) {
      subtitle <- paste0(subtitle, sprintf("  (fill fixed at +/-%g; %d cell%s clipped, max |z| = %.2f)",
                                           lim, clipped, if (clipped == 1L) "" else "s", obs_max))
    }

    g <- ggplot2::ggplot(d, ggplot2::aes(x = block, y = row)) +
      ggplot2::geom_tile(ggplot2::aes(fill = fill_val)) +
      ggplot2::geom_tile(data = d[d$sig, , drop = FALSE], fill = NA,
                         color = blk_col, linewidth = 0.6) +
      ggplot2::geom_tile(data = d[d$cell, , drop = FALSE], fill = NA,
                         color = cell_color, linewidth = 1.0) +
      ggplot2::scale_x_discrete(breaks = x_breaks)
    if (diverging) {
      # Legend for a FIXED scale: the end colors stand for "lim or more" / "-lim or less",
      # because everything beyond the limit is clipped to them (author, 2026-09-10).
      fixed <- !is.null(fixed_lim) && is.finite(fixed_lim) && fixed_lim > 0
      brks <- pretty(c(-lim, lim)); brks <- brks[brks > -lim & brks < lim]
      brks <- c(-lim, brks, lim)
      labs <- format(brks, trim = TRUE)
      if (fixed) { labs[1] <- paste0("\u2264 ", format(-lim)); labs[length(labs)] <- paste0("\u2265 ", format(lim)) }
      g <- g + ggplot2::scale_fill_gradient2(low = neg_color, mid = "white", high = pos_color,
                                             midpoint = 0, limits = c(-lim, lim),
                                             breaks = brks, labels = labs)
    } else {
      g <- g + ggplot2::scale_fill_gradient(low = "white", high = pal$crimson)
    }
    g <- g + ggplot2::labs(x = xlab, y = ylab, fill = fill_lab, subtitle = subtitle) + .bbq_theme()
    if (!is.null(block_labels))
      g <- g + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    g
  }

  panels <- list()

  ## ---- quantile panel ----
  if ("quantile" %in% fams) {
    if (!is.null(detection) && !is.null(detection$z_raw)) {
      if (!is.null(detection$z_white)) {
        vals <- detection$z_white
        flab <- expression(tilde(z)); sub <- "Quantile Shift Coefficient"; note <- ""
      } else {
        vals <- detection$z_raw
        flab <- "z"; sub <- "Quantile Shift Coefficient"
        note <- "  (studentized; whitened cells unavailable)"
      }
      rl <- rownames(vals); if (is.null(rl)) rl <- format(taus)
      zscale <- TRUE                  # z-tilde or studentized z: fix the fill scale
    } else {
      vals <- .bbq_coefs(fit, m, r)$gamma; rl <- format(taus)
      flab <- expression(gamma); sub <- "Quantile Shift Coefficient"
      note <- "  (posterior mean)"
      zscale <- FALSE                 # raw coefficients: +/-3 would be meaningless
    }
    # Panel label: the basis name only (author, 2026-09-10) -- no rule suffix, no note.
    sc <- get_sig("quantile")
    panels$quantile <- heat(vals, rl, sc, get_sig_cells("quantile", sc, nrow(vals)),
                            flab, sub, diverging = TRUE,
                            fixed_lim = if (zscale) z_limit else NULL)
  }

  ## ---- QSS panel ----
  if ("qss" %in% fams) {
    if (!is.null(detection$z_white_qss)) {
      vals <- detection$z_white_qss
      flab <- expression(tilde(z)); sub <- "QSS Shift Coefficient"; note <- ""
    } else {
      vals <- detection$z_qss
      flab <- "z"; sub <- "QSS Shift Coefficient"
      note <- "  (studentized; whitened cells unavailable)"
    }
    rl <- rownames(vals); if (is.null(rl)) rl <- c("LS", "ScS", "SkS", "KS")[seq_len(nrow(vals))]
    sc <- get_sig("qss")
    panels$qss <- heat(vals, rl, sc, get_sig_cells("qss", sc, nrow(vals)),
                       flab, sub, diverging = TRUE, fixed_lim = z_limit)
  }

  ## ---- alternative shape bases: identical rendering, different rotation ----
  shape_panel <- function(fam, zt, zs, default_rows, label) {
    if (!(fam %in% fams)) return(invisible(NULL))
    if (!is.null(zt)) {
      vals <- zt; flab <- expression(tilde(z)); note <- ""
    } else {
      vals <- zs; flab <- "z"; note <- "  (studentized; whitened cells unavailable)"
    }
    if (is.null(vals)) return(invisible(NULL))
    rl <- rownames(vals); if (is.null(rl)) rl <- default_rows[seq_len(nrow(vals))]
    sc <- get_sig(fam)
    panels[[fam]] <<- heat(vals, rl, sc, get_sig_cells(fam, sc, nrow(vals)),
                           flab, label,
                           diverging = TRUE, fixed_lim = z_limit)
  }
  shape_panel("lmom", detection$z_white_lmom, detection$z_lmom,
              c("L-location", "L-scale", "L-skewness", "L-kurtosis"),
              "L-moment Shift Coefficient")

  ## ---- return one panel, or stack them in the order requested ----
  panels <- panels[intersect(fams, names(panels))]
  if (length(panels) == 1L) {
    p <- panels[[1]]
    if (!is.null(title)) p <- p + ggplot2::labs(title = title)
    return(p)
  }
  if (requireNamespace("patchwork", quietly = TRUE)) {
    combo <- patchwork::wrap_plots(panels, ncol = 1L)
    if (!is.null(title)) combo <- combo + patchwork::plot_annotation(title = title)
    return(combo)
  }
  message(length(panels), " basis panels are present; install 'patchwork' to stack ",
          "them into one figure. Returning a named list of ggplot objects instead.")
  panels
}


#' Summary figure: quantile process, shape profile and shift heatmap
#'
#' A three-panel layout. Left:
#' the data with the fitted quantile bands and the localized change-points
#' (\code{\link{plotQuantileProcess}}, with the comparator overlay when
#' \code{comparator} is given). Right: the shape profile of one basis
#' (\code{\link{plotLmomProcess}} or \code{\link{plotQSSProcess}}, bands only)
#' above that basis's shift heatmap (\code{\link{plotGammaHeatmap}}). All three
#' panels read the same \code{detection} object, so the circles, the bordered
#' blocks and the heatmap subtitle follow one decision rule.
#'
#' @inheritParams plotQuantileProcess
#' @param detection A \code{detectChangepoints_gamma()} result whose \code{basis}
#'   includes the requested shape basis.
#' @param basis The shape basis shown on the right and used for the change-point
#'   marks: \code{"lmom"} (default) or \code{"qss"}.
#' @param eta,H,X Passed to the profile panel; see \code{\link{plotQSSProcess}}.
#' @param level Credible-band level of the profile (default 0.95).
#' @param block_labels Labels of the heatmap columns. Default \code{NULL}: the
#'   block start dates when \code{time} is a \code{Date} vector, else the block
#'   start indices.
#' @param label_every Passed to \code{\link{plotGammaHeatmap}}.
#' @param widths,heights Relative widths of the left and right columns, and
#'   relative heights of the profile and heatmap panels.
#' @param seed Optional seed for \code{getEta()}.
#' @return A \pkg{patchwork} object. The talk used 12.6 by 5.5 inches at 200 dpi.
#' @noRd
plotBBQSummary <- function(fit, detection, time = NULL, basis = c("lmom", "qss"),
                           eta = NULL, H = NULL, X = NULL, level = 0.95,
                           comparator = NULL, comparator_label = "Comparator",
                           proposed_label = "Proposed (BBQ)", match_tol = NULL,
                           ylab = "value", block_labels = NULL, label_every = NULL,
                           title = NULL, widths = c(1.05, 1), heights = c(1, 0.8),
                           date_breaks = "12 months", date_labels = "%Y", seed = NULL) {
  .bbq_need_ggplot2()
  if (!requireNamespace("patchwork", quietly = TRUE))
    stop("plotBBQSummary() needs the 'patchwork' package.", call. = FALSE)
  basis <- match.arg(basis)
  if (is.null(detection$tests[[basis]]))
    stop("`detection` does not carry the '", basis, "' basis; rerun ",
         "detectChangepoints_gamma() with it in `basis`.", call. = FALSE)
  n <- length(fit$y)
  if (is.null(time)) time <- seq_len(n)
  if (is.null(eta)) eta <- getEta(fit, H = H, X = X, seed = seed)
  if (is.null(block_labels)) {
    os <- detection$detected_blocks$obs_start
    block_labels <- if (inherits(time, "Date")) format(time[os]) else as.character(os)
  }
  no_x <- ggplot2::theme(axis.title.x = ggplot2::element_blank())
  p_q <- plotQuantileProcess(fit, time = time, detection = detection, basis = basis,
                             show_onset = FALSE, ylab = ylab,
                             comparator = comparator, comparator_label = comparator_label,
                             proposed_label = proposed_label, match_tol = match_tol,
                             date_breaks = date_breaks, date_labels = date_labels) + no_x
  prof <- if (basis == "lmom") plotLmomProcess else plotQSSProcess
  p_s <- prof(fit, eta = eta, H = H, X = X, time = time, level = level,
              detection = detection, basis = basis,
              show_onset = FALSE, show_located = FALSE,
              date_breaks = date_breaks, date_labels = date_labels) + no_x +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 7),
                   strip.text = ggplot2::element_text(size = 8))
  p_h <- plotGammaHeatmap(fit, detection = detection, block_labels = block_labels,
                          basis = basis, label_every = label_every,
                          note_clipping = FALSE) + no_x +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7),
                   axis.text.y = ggplot2::element_text(size = 7),
                   legend.key.width = ggplot2::unit(0.26, "cm"),
                   legend.key.height = ggplot2::unit(0.50, "cm"),
                   legend.text = ggplot2::element_text(size = 6),
                   legend.title = ggplot2::element_text(size = 7),
                   plot.subtitle = ggplot2::element_text(size = 9))
  right <- patchwork::wrap_elements(p_s) / patchwork::wrap_elements(p_h) +
    patchwork::plot_layout(heights = heights)
  fig <- (patchwork::wrap_elements(p_q) | right) + patchwork::plot_layout(widths = widths)
  if (!is.null(title)) fig <- fig + patchwork::plot_annotation(title = title)
  fig
}
