library(brms)
library(tidyverse)
library(ggdist)
library(distributional)
library(patchwork)
library(tidybayes)
set.seed(2021)

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

model_1 <- brm(data = mixed_data, family = binomial,
      n_success | trials(n_trials) ~ 1 + (1 | n_batch),
      prior = c(prior(normal(0, 1), class = Intercept),
                prior(cauchy(0, 1), class = sd)),
      iter = 4000, warmup = 1000, chains = 4, cores = 4,
      seed = 12)

mixed_data %>%
expand(nesting(n_trials, n_batch, n_coin)) %>%
add_predicted_draws(model_1) -> plot_1

mixed_data %>%
expand(nesting(n_trials, n_batch, n_coin)) %>%
add_predicted_draws(model_1, re_formula = NA) -> plot_2

gather_draws(model_1, b_Intercept, r_n_batch[batch,v]) %>%
mutate(p_heads = inv_logit_scaled(.value)) -> plot_1

plot_1a <-
	plot_1 %>%
	ungroup() %>%
	filter(.variable == "b_Intercept") %>%
	select(p_heads)

plot_1b <-
	plot_1 %>%
	ungroup() %>%
	filter(.variable != "b_Intercept")

ggplot() +
stat_slab(data = plot_1a, aes(x = p_heads), fill = "grey50", colour = "grey50", alpha = 0.5) +
stat_slab(data = plot_1b, aes(x = p_heads, fill = batch, colour = batch), alpha = 0.5) +
facet_wrap(vars(batch)) +
scale_fill_distiller(palette = "YlGnBu", aesthetics = c("colour", "fill")) +
scale_x_continuous(limits = c(0, 1))