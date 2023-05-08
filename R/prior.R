#' Gibbs fixed values object
#'
#' Specify hyperparameters \eqn{a_\sigma} and \eqn{b_\sigma} for Inverse Gamma
#' prior on \eqn{\sigma^2}.
#'
#' @param a_sigma A positive number
#' @param b_sigma A positive number
#'
#' @return An object of class \code{my_prior}.
#'
#' @export
get_prior = function(a_sigma, b_sigma)
{
	stopifnot(all(a_sigma > 0))
	stopifnot(all(b_sigma > 0))
	ret = list(a_sigma = a_sigma, b_sigma = b_sigma)
	class(ret) = "my_prior"
	return(ret)
}
