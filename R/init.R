#' Gibbs initial values object
#'
#' @param n Number of observations.
#' @param muX A vector of length \code{n}.
#' @param sigma2 A scalar
#'
#' @return An object of class \code{my_init}.
#'
#' @export
get_init = function(n, muX = NULL, sigma2 = NULL)
{
	if (is.null(muX)) {
		muX = numeric(n)
	}
	stopifnot(length(muX) == n)

	if (is.null(sigma2)) {
		sigma2 = 1
	}
	stopifnot(length(sigma2) == 1)

	ret = list(
		muX = muX,
		sigma2 = sigma2
	)
	class(ret) = "my_init"
	return(ret)
}
