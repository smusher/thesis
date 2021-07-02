library(tidyverse)
library(ggdist)
library(distributional)
library(patchwork)
set.seed(2021)

#drc_priors
ec50 <- dist_normal(-4, 1)
hill <- dist_normal(1, 0.3)
emax <- dist_uniform(0, 1)
sigma <- dist_cauchy(0, 1)

(autoplot(ec50) + autoplot(hill)) / (autoplot(emax) + autoplot(sigma))

#pcf_priors
L <- dist_normal(log(1), log(5))
D <- dist_normal(log(0.1), log(5))
Ka <- dist_normal(log(1e4), log(5))
sigma <- dist_cauchy(0, 1)

(autoplot(L) + autoplot(D)) / (autoplot(Ka) + autoplot(sigma))

unfair_toss <-
	tibble(
		n_trials = c(10, 50, 200),
		n_success_fair = c(5, 25, 100)
		) %>%
	rowwise() %>%
	mutate(
		n_success = rbinom(1, n_trials, 0.8)
		) %>%
	expand(
		nesting(n_trials, n_success, n_success_fair), p_heads = seq(from = 0, to = 1, length.out = 101)
		) %>%
	group_by(n_trials) %>%
	mutate(
		prior = dbinom(
        	x = n_success_fair,
        	size = n_trials,
        	prob = p_heads
        	),
        likelihood = dbinom(
        	x = n_success,
        	size = n_trials,
        	prob = p_heads
        	),
        posterior = prior * likelihood
        ) %>%
  	mutate(
  		prior = prior / sum(prior),
        likelihood = likelihood / sum(likelihood),
        posterior = posterior / sum(posterior)
        ) %>%
  	pivot_longer(c(prior, likelihood, posterior))

ggplot(unfair_toss, aes(x = p_heads)) +
geom_line(aes(y = value, colour = name)) +
scale_x_continuous("proportion heads", breaks = c(0, .5, 1)) +
scale_y_continuous("plausibility", breaks = NULL) +
theme(panel.grid = element_blank()) +
facet_grid(cols = vars(n_trials), scales = "free")

tibble(x = dist_normal(0, 1)) %>%
  ggplot(aes(dist = x, y = "a")) +
  stat_dist_slab(aes(
    fill = stat(cut_cdf_qi(cdf, .width = c(.5, .8, .95, 1)))
  )) +
  scale_fill_brewer(direction = -1, na.translate = FALSE)

mixed_data <- tibble(
    n_trials = 50,
    n_batch = seq(1:9),
    p_heads = rnorm(9, 0.5, 0.1)
    ) %>%
  expand(nesting(n_trials, n_batch, p_heads), n_coin = seq(1:100)) %>%
  rowwise() %>%
  mutate(
    n_success = rbinom(1, n_trials, p_heads)
    )

ggplot(mixed_data) +
geom_dotplot(aes(x = n_success, y = stat(ncount), fill = n_batch), method="histodot", binwidth = 1, dotsize=2) +
facet_wrap(vars(n_batch)) +
scale_fill_distiller(palette = "YlGnBu") +
scale_x_continuous(limits = c(0, 50))

#coinflip_priors
alpha_j<- dist_normal(0, 1)

autoplot(alpha_j)


sigma_j <- dist_cauchy(0, 1)

autoplot(alpha_j) + autoplot(sigma_j)