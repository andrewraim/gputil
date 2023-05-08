#' Gibbs fixed values object
#'
#' Specify conditionals which should be kept fixed (i.e., not redrawn) during
#' sampling.
#'
#' @param muX Logical
#' @param sigma2 Logical
#'
#' @return An object of class \code{my_fixed}.
#'
#' @export
get_fixed = function(muX = FALSE, sigma2 = FALSE)
{
	ret = list(muX = muX, sigma2 = sigma2)
	class(ret) = "my_fixed"
	return(ret)
}
