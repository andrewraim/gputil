library(gputil)
library(ggplot2)

set.seed(1235)

mu_true = sin

sigma2_true = 0.05^2

n = 100
x = seq(-4*pi, 4*pi, length.out = n)
X = matrix(x)
muX_true = mu_true(X)
y = rnorm(n, muX_true, sqrt(sigma2_true))

dat_plot = data.frame(x = as.numeric(X), y = y)
gg = ggplot() +
	geom_point(data = dat_plot, aes(x, y)) +
	stat_function(fun = mu_true, n = 201, lty = 1) +
	ylab("mu(x)") +
	theme_minimal()
ggsave("data.pdf", gg)

# ----- Fit with sampler -----
init = get_init(n, muX = muX_true)
fixed = get_fixed(muX = FALSE)
prior = get_prior(a_sigma = 5, b_sigma = 10)
control = get_control(R = 5000, report_period = 1000, save_latent = FALSE)

K = exp(-plgp::distance(X,X))
K_eig = eigen(K)

gibbs_out = gibbs(y, K_eig, init, prior, control, fixed)
print(gibbs_out)

# Check trace plot and histogram of sigma2 draws
plot(gibbs_out$sigma2_hist, type = "l")
hist(gibbs_out$sigma2_hist)

# ----- Compute Predictions -----
R = gibbs_out$R_keep
n0 = 200
x0 = seq(-4*pi, 4*pi, length.out = n0)
X0 = matrix(x0)

if (FALSE) {
	# Draw from the full joint posterior predictive distribution
	K22 = exp(-plgp::distance(X0, X0))
	K21 = exp(-plgp::distance(X0, X))
	mu0_mcmc = predict(gibbs_out, K22, K21)
} else {
	# Draw from the marginal posterior predictive distribution for each input
	mu0_mcmc = matrix(NA, R, n0)
	for (i in 1:n0) {
		if (i %% 10 == 0) {
			raim::logger("Predicting input %d of %d\n", i, n0)
		}
		K22 = exp(-plgp::distance(X0[i,], X0[i,]))
		K21 = exp(-plgp::distance(X0[i,], X))
		mu0_mcmc[,i] = predict(gibbs_out, K22, K21)
	}
}

# Summarize with point estimates and interval bounds
mu0_hat = apply(mu0_mcmc, 2, mean)
mu0_lo = apply(mu0_mcmc, 2, quantile, prob = 0.05)
mu0_hi = apply(mu0_mcmc, 2, quantile, prob = 0.95)

# ----- Plot Predictions -----
dat_plot = data.frame(x = as.numeric(X), y = y)
dat0_plot = data.frame(x0 = x0, mu0_hat = mu0_hat, mu0_lo = mu0_lo, mu0_hi = mu0_hi)

gg = ggplot() +
	geom_point(data = dat_plot, aes(x, y)) +
	stat_function(fun = mu_true, n = 201, lty = 1) +
	geom_line(data = dat0_plot, aes(x = x0, y = mu0_hat), col = "blue", linewidth = 0.8) +
	geom_ribbon(data = dat0_plot,
		aes(x = x0, y = mu0_hat, ymin = mu0_lo, ymax = mu0_hi),
		fill = "blue", alpha = 0.25) +
	ylab("mu(x)") +
	theme_minimal()
ggsave("predict.pdf", gg)
