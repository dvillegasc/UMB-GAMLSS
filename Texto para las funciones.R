#' Unit Maxwell-Boltzmann distribution
#'
#' @author David Villegas Ceballos, \email{david.villegas1@udea.edu.co}
#'
#' @description
#' These functions define the density, distribution function, quantile
#' function and random generation for the Unit Maxwell-Boltzmann distribution
#' with parameter \eqn{\mu}.
#'
#' @param x,q vector of (non-negative integer) quantiles.
#' @param p vector of probabilities.
#' @param mu vector of the mu parameter.
#' @param n number of random values to return.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are \eqn{P[X <= x]}, otherwise, \eqn{P[X > x]}.
#'
#' @references
#' Bicher, C., Bakouch, H. S., Biçor, H. D., Alomair, G., Hussain, T., y Almohisen, A. (2024). Unidad de Distribución Maxwell-Boltzmann y su Aplicación a Concentraciones de Datos Contaminantes. Axiomas, 13(4), 226.
#'
#' @seealso \link{UMB}.
#'
#' @details
#' The Unit Maxwell-Boltzmann distribution with parameter \eqn{\mu}
#' has a support in \eqn{(0, 1)} and density given by
#'
#' \eqn{f(x| \mu) = \frac{2 \mu}{(\mu+(2-\mu)x)^2} }
#'
#' for \eqn{0 < x < 1} and \eqn{\mu > 0}.
#'
#' @example examples/examples_dUHLG.R
#'
#' @export