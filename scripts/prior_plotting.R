library(tidyverse)
library(ggdist)
library(distributional)
library(patchwork)

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

#fair coin
fair_flips <- rbinom(10, 1, 0.5)

#unfair coin
unfair_flips <- rbinom(10, 1, 0.8)

fair_toss <-
	tibble(
		n_trials = c(10, 50, 100),
		) %>%
	rowwise() %>%
	mutate(
		n_success = rbinom(1, n_trials, 0.5)
		) %>%
	expand(
		nesting(n_trials, n_success), p_heads = seq(from = 0, to = 1, length.out = 51)
		) %>%
	group_by(n_trials) %>%
	mutate(
		prior = dunif(p_heads, 0, 1),
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

ggplot(fair_toss, aes(x = p_heads)) +
geom_line(aes(y = value, colour = name)) +
scale_x_continuous("proportion heads", breaks = c(0, .5, 1)) +
scale_y_continuous("plausibility", breaks = NULL) +
theme(panel.grid = element_blank()) +
facet_grid(rows = vars(name), cols = vars(n_trials))

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

unfair_toss %>%
filter(name == "posterior", n_trials == 50) 