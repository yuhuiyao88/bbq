#' bbq: Bayesian Bi-Aspect Quantile Change-Point Detection
#'
#' Retrospective (Phase I) change-point detection in individual autocorrelated
#' observations by Bayesian multiple quantile regression at the quantile levels
#' (0.25, 0.50, 0.75), with location and scale inference through the median and
#' the interquartile range.
#'
#' The manuscript workflow is \code{\link{bbq}}, or its steps
#' \code{\link{bbqDesign}}, \code{\link{bbqTune}}, \code{\link{bbqFit}} and
#' \code{\link{bbqDetect}}. Figures use \code{\link{plotQuantileProcess}},
#' \code{\link{plotQSSProcess}} and \code{\link{plotGammaHeatmap}}.
#'
#' @importFrom stats dbeta dgamma dnorm plogis
#' @importFrom utils globalVariables
#' @keywords internal
"_PACKAGE"

# Column names used inside ggplot2 aesthetics.
globalVariables(c("block", "fill_val", "stat", "lo", "hi", "mid", "q1", "q3",
                         "med", "x", "y", "src"))
