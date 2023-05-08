#' Gibbs sampler for Gaussian Process regression model
#'
#' @param y Vector of observations
#' @param K The matrix \eqn{k(X, X)}; see details.
#' @param init Initial values object from \code{get_init}
#' @param prior Hyperparameters object from \code{get_prior}
#' @param control Control object from \code{get_control}
#' @param fixed Object declaring fixed quantities from \code{get_fixed}
#' @param object An object of class \code{gputil_fit} returned by the Gibbs
#' sampler.
#' @param x An object of class \code{gputil_fit} returned by the Gibbs
#' sampler.
#' @param K22 The matrix \eqn{k(X_0, X_0)}; see details.
#' @param K21 The matrix \eqn{k(X_0, X)}; see details.
#'
#' @return
#' The \code{gibbs} function returns an object with sampler results. The
#' \code{predict} function returns an \eqn{R \times n_0} matrix of draws from
#' the posterior predictive distribution.
#'
#' @details
#' This is a Gibbs sampler for the model
#' \deqn{
#' \bm{y} = \mu(\bm{X}) + \bm{\epsilon}, \quad \bm{\epsilon} \sim
#' \text{N}(\bm{0}, \sigma^2 \bm{I}), \\
#' \mu(\bm{X}) \sim \text{N}(\bm{0}, k(\bm{X},\bm{X})) \\
#' \sigma^2 \sim \text{IG}(a_\sigma, b_\sigma)
#' }
#'
#' We have used the notation \eqn{\bm{X}} for the \eqn{n \times d} matrix of
#' \eqn{n} inputs whose corresponding observations are
#' \eqn{\bm{y} = (y_1, \ldots, y_n)}. We will also write \eqn{\bm{X}_0} as an
#' \eqn{n_0 \times d} matrix of inputs
#' \eqn{\bm{x}_{01}, \ldots, \bm{x}_{0 n_0}} on which to compute predictions.
#' The user may choose any covariance kernel \eqn{k(\bm{x}, \bm{x}')}, but it
#' is their responsibility to evaluate it on the inputs. Any parameters of the
#' kernel are considered fixed during MCMC.
#'
#' The \code{gibbs} functons requires the matrix \eqn{k(\bm{X}, \bm{X})} whose
#' \eqn{(i,j)}th element is \eqn{k(\bm{x}_i, \bm{x}_j)}. This is the argument
#' \code{K}. \code{K} may be provided as a matrix or the result of a spectral
#' decomposition; if \code{K} is a \code{list}, we assume eigenvalues are in
#' \code{K$values} and \code{K$vectors} is a matrix whose columns are the
#' corresponding eigenvectors. Otherwise, if \code{K} is not a \code{list}, it
#' is interpreted as a matrix.
#'
#' The \code{predict} function requires the \eqn{n_0 \times n_0} matrix
#' \eqn{k(\bm{X}_0, \bm{X}_0) = k(\bm{x}_{0i}, \bm{x}_{0j})} as argument
#' \code{K22}, and the matrix \eqn{k(\bm{X}_0, \bm{X}) = k(\bm{x}_{0i}, \bm{x}_j)}
#' as argument \code{K21}.
#'
#' @name gibbs
NULL

#' @name gibbs
#' @export
gibbs = function(y, K, init = get_init(), prior = get_prior(),
	control = get_control(), fixed = get_fixed())
{
	n = length(y)

	# Precompute quantities based on K
	if (is.list(K)) {
		U = K$vectors
		lambda = K$values
	} else {
		eig = eigen(K)
		U = eig$vectors
		lambda = eig$values
	}
	z = crossprod(U, y)

	stopifnot(class(control) == "my_control")
	R = control$R
	burn = control$burn
	thin = control$thin
	report_period = control$report_period
	save_latent = control$save_latent

	rep_keep = 0
	R_keep = ceiling((R - burn) / thin)

	# Set up histories
	sigma2_hist = numeric(R_keep)
	muX_hist = NULL
	if (save_latent) {
		muX_hist = matrix(NA, R_keep, n)
	}

	# Set up initial values.
	stopifnot(class(init) == "my_init")
	muX = init$muX
	sigma2 = init$sigma2

	# Set up prior
	stopifnot(class(prior) == "my_prior")
	a_sigma = prior$a_sigma
	b_sigma = prior$b_sigma

	# Set up fixed
	stopifnot(class(fixed) == "my_fixed")

	# Set up timers
	elapsed = list(muX = 0, sigma2 = 0)

	for (rep in 1:R) {
		if (rep %% report_period == 0) {
			raim::logger("Starting rep %d\n", rep)
		}

		# Draw [muX | rest]
		if (!fixed$muX) {
			st = Sys.time()

			# First draw from U^T muX, whose components are independent
			omega = 1/sigma2 + 1/lambda
			mean = 1/sigma2 * z / omega
			muX_tx = rnorm(n, mean, 1 / sqrt(omega))

			# Transform to muX
			muX = U %*% muX_tx
			elapsed$muX = elapsed$muX + as.numeric(Sys.time() - st, units = "secs")
		}

		# Draw [sigma2 | rest]
		if (!fixed$sigma2) {
			st = Sys.time()
			diff = y - muX
			norm2 = sum(diff^2)
			a_star = a_sigma + n / 2
			b_star = b_sigma + norm2 / 2
			sigma2 = raim::r_invgamma(1, a_star, b_star)
			elapsed$sigma2 = elapsed$sigma2 + as.numeric(Sys.time() - st, units = "secs")
		}

		if (rep > burn && rep %% thin == 0) {
			rep_keep = rep_keep + 1
			sigma2_hist[rep_keep] = sigma2
			if (save_latent) {
				muX_hist = matrix(NA, R_keep, n)
			}

		}
	}

	ret = list(muX_hist = muX_hist, sigma2_hist = sigma2_hist, R_keep = R_keep,
		elapsed = elapsed, R = R, burn = burn, thin = thin, fixed = fixed,
		K_eig = list(vectors = U, values = lambda), z = z, control = control)
	class(ret) = "gputil_fit"
	return(ret)
}

#' @name gibbs
#' @export
summary.gputil_fit = function(object, ...)
{
	pr = c(0.025, 0.5, 0.975)

	df_sigma2 = data.frame(
		mean(object$sigma2_hist),
		sd(object$sigma2_hist),
		t(quantile(object$sigma2_hist, probs = pr))
	)
	rownames(df_sigma2) = "sigma2"

	df = rbind(df_sigma2)
	colnames(df) = c("Mean", "SD", "2.5%", "50%", "97.5%")
	return(df)
}

#' @name gibbs
#' @export
print.gputil_fit = function(x, ...)
{
	cat("Summary of fit\n")
	print(summary(x))

	cat("----\n")
	raim::printf("Total iterations R: %d   burn: %d   thin: %d   Saved draws: %d\n",
		x$R, x$burn, x$thin, x$R_keep)

	kv = sprintf("%s: %g", names(x$elapsed), x$elapsed)
	raim::printf("Elapsed time [sec]  %s\n", paste(kv, collapse = "  "))

	is_fixed = unlist(x$fixed)
	if (any(is_fixed)) {
		idx = which(is_fixed)
		fixed_str = paste(names(x$fixed)[idx], collapse = ", ")
		raim::printf("Fixed: %s\n", fixed_str)
	}

}

#' @name gibbs
#' @export
predict.gputil_fit = function(object, K22, K21)
{
	z = object$z
	n = length(z)
	n0 = nrow(K22)
	stopifnot(n0 == nrow(K21))
	stopifnot(n == ncol(K21))
	stopifnot(n0 == ncol(K22))

	report_period = object$control$report_period

	U = object$K_eig$vectors
	lambda = object$K_eig$values
	K21U = K21 %*% U

	R = length(object$sigma2_hist)
	sigma2_draws = object$sigma2_hist
	mu0_mcmc = matrix(NA, R, n0)

	for (r in 1:R) {
		sigma2 = sigma2_draws[r]
		mu0 = K21U %*% (z / (sigma2 + lambda))
		Sigma0 = K22 - K21U %*% (1/(sigma2 + lambda) * t(K21U))
		mu0_mcmc[r,] = raim::r_mvnorm(1, mu0, Sigma0)
	}

	return(mu0_mcmc)
}
