#' Gibbs control object
#'
#' @param R Desired number of draws from the sampler.
#' @param burn Number of draws to burn before saving.
#' @param thin The period for saving samples after burn phase has ended.
#' @param report_period The period for printing progress.
#' @param save_latent Logical; if true, save draws of the latent variable
#' \eqn{mu(\bm{X})}.
#'
#' @return An object of class \code{my_control}.
#'
#' @export
get_control = function(R, burn = 0, thin = 1, report_period = R+1,
	save_latent = FALSE)
{
	ret = list(R = R, burn = burn, thin = thin, report_period = report_period,
		save_latent = save_latent)
	class(ret) = "my_control"
	return(ret)
}
