library(tidyverse)
library(broom)
library(brms)
library(tidybayes)
library(ggdist)
library(distributional)
library(ggbeeswarm)
brms_iter <- 45000
brms_warmup <- 5000
brms_chains <- 4
brms_thin <- 10
brms_seed <- 2021
interesting_constructs <- c("WT-GFP+SUR", "W311*-GFP+SUR")
interesting_nucleotides <- c("ATP", "TNP-ATP")
hill_equation <- function(log_concentration, ec50, hill, floor){
    floor + ((1 - floor) / (1 + 10^((ec50 - log_concentration) * -hill)))
}

concresp <-
    read_csv("data/combined_drc_data.csv") %>%
    mutate(log_concentration = log10(concentration)) %>%
    filter(log_concentration > -Inf) %>%
    group_by(unique_experiment_id) %>%
    mutate(observations = n()) %>%
    filter(observations >= 3) %>%
    ungroup() %>%
    drop_na()

readRDS("data/hill_fits/ua_free_fit_pop.rds") -> unadjusted_fits_free_floor_pop
readRDS("data/hill_fits/a_free_fit_pop.rds") -> adjusted_fits_free_floor_pop
readRDS("data/hill_fits/a_fixed_fit_pop.rds") -> adjusted_fits_fixed_floor_pop

readRDS("data/hill_fits/ua_free_fit.rds") -> unadjusted_fits_free_floor
readRDS("data/hill_fits/a_free_fit.rds") -> adjusted_fits_free_floor
readRDS("data/hill_fits/a_fixed_fit.rds") -> adjusted_fits_fixed_floor

nd_1 <- tibble(log_concentration = seq(-8, -2, length.out = 51))

unadjusted_fits_free_floor_pop %>%
mutate(augmented = map(fit, augment, newdata = nd_1)) %>%
unnest(augmented) -> ua_free_plot_pop

unadjusted_fits_free_floor %>%
mutate(augmented = map(fit, augment, newdata = nd_1)) %>%
unnest(augmented) -> ua_free_plot

adjusted_fits_fixed_floor_pop %>%
mutate(augmented = map(fit, augment, newdata = nd_1)) %>%
unnest(augmented) -> a_fixed_plot_pop

adjusted_fits_fixed_floor %>%
mutate(augmented = map(fit, augment, newdata = nd_1)) %>%
unnest(augmented) -> a_fixed_plot

ggplot() +
geom_line(
	data = ua_free_plot_pop %>% filter(measure == "current", construct %in% interesting_constructs, nucleotide %in% interesting_nucleotides),
	aes(x = log_concentration, y = .fitted, colour = construct)
	) +
geom_quasirandom(
	data = concresp %>% filter(measure == "current", construct %in% interesting_constructs, nucleotide %in% interesting_nucleotides),
	aes(x = log_concentration, y = response, fill = construct), shape = 21, width = 0.1
	) +
facet_grid(cols = vars(nucleotide))

ggplot() +
geom_line(
	data = ua_free_plot %>% filter(measure == "current", construct %in% interesting_constructs, nucleotide %in% interesting_nucleotides),
	aes(x = log_concentration, y = .fitted, colour = construct)
	) +
geom_point(
	data = concresp %>% filter(measure == "current", construct %in% interesting_constructs, nucleotide %in% interesting_nucleotides),
	aes(x = log_concentration, y = response, fill = construct), shape = 21
	) +
facet_wrap(vars(nucleotide, unique_experiment_id))

ggplot() +
geom_line(
	data = ua_free_plot_pop %>% filter(measure == "fluorescence", construct %in% interesting_constructs, nucleotide %in% interesting_nucleotides),
	aes(x = log_concentration, y = .fitted, colour = construct)
	) +
geom_quasirandom(
	data = concresp %>% filter(measure == "fluorescence", construct %in% interesting_constructs, nucleotide %in% interesting_nucleotides),
	aes(x = log_concentration, y = response, fill = construct), shape = 21, width = 0.1
	) +
facet_grid(cols = vars(method))

ggplot() +
geom_line(
	data = ua_free_plot %>% filter(measure == "fluorescence", construct %in% interesting_constructs, nucleotide %in% interesting_nucleotides),
	aes(x = log_concentration, y = .fitted, colour = construct)
	) +
geom_point(
	data = concresp %>% filter(measure == "fluorescence", construct %in% interesting_constructs, nucleotide %in% interesting_nucleotides),
	aes(x = log_concentration, y = response, fill = construct), shape = 21
	) +
facet_wrap(vars(method, unique_experiment_id))

ggplot() +
geom_line(
	data = a_fixed_plot_pop %>% filter(measure == "fluorescence", construct %in% interesting_constructs, nucleotide %in% interesting_nucleotides),
	aes(x = log_concentration, y = .fitted, colour = construct)
	) +
geom_quasirandom(
	data = concresp %>% filter(measure == "fluorescence", construct %in% interesting_constructs, nucleotide %in% interesting_nucleotides),
	aes(x = log_concentration, y = log2(response+1), fill = construct), shape = 21, width = 0.1
	) +
facet_grid(cols = vars(method))

ggplot() +
geom_line(
	data = ua_free_plot %>% filter(method == "unroofed", construct %in% interesting_constructs, nucleotide %in% interesting_nucleotides),
	aes(x = log_concentration, y = .fitted), colour = "black"
	) +
geom_line(
	data = a_fixed_plot %>% filter(method == "unroofed", construct %in% interesting_constructs, nucleotide %in% interesting_nucleotides),
	aes(x = log_concentration, y = .fitted), colour = "grey50"
	) +
geom_point(
	data = concresp %>% filter(method == "unroofed", construct %in% interesting_constructs, nucleotide %in% interesting_nucleotides),
	aes(x = log_concentration, y = response, fill = construct), shape = 21
	) +
geom_point(
	data = concresp %>% filter(method == "unroofed", construct %in% interesting_constructs, nucleotide %in% interesting_nucleotides),
	aes(x = log_concentration, y = log2(response+1), colour = construct), shape = 21
	) +
facet_wrap(vars(unique_experiment_id))

unadjusted_fits_free_floor %>%
mutate(tidied = map(fit, tidy)) %>%
unnest(tidied) %>%
select(construct, method, measure, nucleotide, unique_experiment_id, term, estimate) %>%
pivot_wider(names_from = term, values_from = estimate) -> ua_free_tidy

adjusted_fits_free_floor %>%
mutate(tidied = map(fit, tidy)) %>%
unnest(tidied) %>%
select(construct, method, measure, nucleotide, unique_experiment_id, term, estimate) %>%
pivot_wider(names_from = term, values_from = estimate) -> a_free_tidy

adjusted_fits_fixed_floor %>%
mutate(tidied = map(fit, tidy)) %>%
unnest(tidied) %>%
select(construct, method, measure, nucleotide, unique_experiment_id, term, estimate) %>%
pivot_wider(names_from = term, values_from = estimate) -> a_fixed_tidy

ranef_form_1 <- bf(
    ec50 ~ 0 + Intercept + (1|construct:method:measure:nucleotide),
    sigma ~ (1|construct:method:measure:nucleotide)
    )

priors <- c(
    prior(normal(-5, 3), class = b),
    prior(cauchy(0, 3), class = sd)
    )

brm(
    formula = ranef_form_1,
    family = gaussian(),
    data = ua_free_tidy,
    prior = priors,
    cores = getOption("mc.cores", 4),
    sample_prior = "yes",
    save_all_pars = TRUE,
    chains = brms_chains,
    iter = brms_iter,
    warmup= brms_warmup,
    thin  = brms_thin,
    seed = brms_seed,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    file = "/home/sam/thesis/data/hill_fits/drc_model_unadjusted_free.rds"
    ) -> model_1

brm(
    formula = ranef_form_1,
    family = gaussian(),
    data = a_free_tidy,
    prior = priors,
    cores = getOption("mc.cores", 4),
    sample_prior = "yes",
    save_all_pars = TRUE,
    chains = brms_chains,
    iter = brms_iter,
    warmup= brms_warmup,
    thin  = brms_thin,
    seed = brms_seed,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    file = "/home/sam/thesis/data/hill_fits/drc_model_adjusted_free.rds"
    ) -> model_2

brm(
    formula = ranef_form_1,
    family = gaussian(),
    data = a_fixed_tidy,
    prior = priors,
    cores = getOption("mc.cores", 4),
    sample_prior = "yes",
    save_all_pars = TRUE,
    chains = brms_chains,
    iter = brms_iter,
    warmup= brms_warmup,
    thin  = brms_thin,
    seed = brms_seed,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    file = "/home/sam/thesis/data/hill_fits/drc_model_adjusted_fixed.rds"
    ) -> model_3

ua_free_tidy %>%
expand(nesting(construct, method, measure, nucleotide)) %>%
add_predicted_draws(model_1) %>%
group_by(construct, method, measure, nucleotide) -> plot_1

a_fixed_tidy %>%
expand(nesting(construct, method, measure, nucleotide)) %>%
add_predicted_draws(model_3) %>%
group_by(construct, method, measure, nucleotide) -> plot_2

ggplot() +
stat_slab(
    data = plot_1 %>% filter(construct %in% interesting_constructs, nucleotide %in% interesting_nucleotides),
    aes(y = construct, x = .prediction, fill = interaction(measure, method), fill_ramp = stat(cut_cdf_qi(cdf, .width = c(.5, .8, .95), labels = scales::percent_format())))
    ) +
geom_point(
    data = ua_free_tidy %>% filter(construct %in% interesting_constructs, nucleotide %in% interesting_nucleotides),
    position = position_dodge(width=0.2), aes(y = construct, x = ec50, fill = interaction(measure, method)), shape=21, size = 2
    ) +
scale_fill_brewer(palette = "Pastel1", aesthetics = c("fill", "colour")) +
scale_fill_ramp_discrete(range = c(1, 0.2), na.translate = FALSE) +
facet_grid(cols = vars(nucleotide)) +
labs(fill_ramp = "Interval")

ggplot() +
stat_slab(
    data = plot_2 %>% filter(construct %in% interesting_constructs, nucleotide %in% interesting_nucleotides),
    aes(y = construct, x = .prediction, fill = interaction(measure, method), fill_ramp = stat(cut_cdf_qi(cdf, .width = c(.5, .8, .95), labels = scales::percent_format())))
    ) +
geom_point(
    data = a_fixed_tidy %>% filter(construct %in% interesting_constructs, nucleotide %in% interesting_nucleotides),
    position = position_dodge(width=0.2), aes(y = construct, x = ec50, fill = interaction(measure, method)), shape=21, size = 2
    ) +
scale_fill_brewer(palette = "Pastel1", aesthetics = c("fill", "colour")) +
scale_fill_ramp_discrete(range = c(1, 0.2), na.translate = FALSE) +
facet_grid(cols = vars(nucleotide)) +
labs(fill_ramp = "Interval")
