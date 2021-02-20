library(tidyverse)
library(broom)
library(nls.multstart)
library(brms)
library(tidybayes)
library(ggdist)
library(distributional)
brms_iter <- 45000
brms_warmup <- 5000
brms_chains <- 4
brms_thin <- 10
brms_seed <- 2021

concresp <-
    read_csv("data/combined_drc_data.csv") %>%
    mutate(log_concentration = log10(concentration)) %>%
    filter(log_concentration > -Inf) %>%
    group_by(unique_experiment_id) %>%
    mutate(observations = n()) %>%
    filter(observations >= 3) %>%
    ungroup() %>%
    drop_na()

currents <-
    concresp %>%
    filter(measure == "current")

fluorescence <-
    concresp %>%
    filter(measure == "fluorescence")

hill_equation <- function(log_concentration, ec50, hill, floor){
    floor + ((1 - floor) / (1 + 10^((ec50 - log_concentration) * -hill)))
}

unadjusted_fits_free_floor <-
    concresp %>%
    group_by(construct, method, measure, nucleotide, unique_experiment_id) %>%
    nest() %>%
    mutate(fit = map(data, ~ nls_multstart(
            response ~ hill_equation(log_concentration, ec50, hill, floor),
            data = data.frame(.),
            iter = c(10,10,10),
            lower = c(-Inf, 0, 0),
            start_lower = c(ec50 = -7, hill = 0, floor= 0),
            start_upper = c(ec50 = -1, hill = 2, floor = 1),
            na.action = na.omit
              )
    )
    )

unadjusted_fits_free_floor %>%
filter(map_lgl(fit, ~ .x$convInfo$isConv == TRUE)) %>%
mutate(tidied = map(fit, tidy)) %>%
unnest(tidied) %>%
select(construct, method, measure, nucleotide, unique_experiment_id, term, estimate) %>%
pivot_wider(names_from = term, values_from = estimate) -> ua_free_tidy

adjusted_fits_free_floor <-
    fluorescence %>%
    group_by(construct, method, unique_experiment_id) %>%
    nest() %>%
    mutate(fit = map(data, ~ nls_multstart(
            log2(response + 1) ~ hill_equation(log_concentration, ec50, hill, floor),
            data = data.frame(.),
            iter = c(10, 10, 10),
            lower = c(-Inf, 0, 0),
            start_lower = c(ec50 = -7, hill = 0, floor= 0),
            start_upper = c(ec50 = -1, hill = 2, floor = 1),
            na.action = na.omit
              )
    )
    )

adjusted_fits_free_floor %>%
filter(map_lgl(fit, ~ .x$convInfo$isConv == TRUE)) %>%
mutate(tidied = map(fit, tidy)) %>%
unnest(tidied) %>%
select(construct, method, unique_experiment_id, term, estimate) %>%
pivot_wider(names_from = term, values_from = estimate) %>%
mutate(measure = "fluorescence", nucleotide = "TNP-ATP") -> a_free_tidy

adjusted_fits_fixed_floor <-
    fluorescence %>%
    filter(!(unique_experiment_id %in% c(176,328,331))) %>%
    group_by(construct, method, unique_experiment_id) %>%
    nest() %>%
    mutate(fit = map(data, ~ nls_multstart(
            log2(response + 1) ~ 0.1 + ((1 - 0.1) / (1 + 10^((ec50 - log_concentration) * -hill))),
            data = data.frame(.),
            iter = c(15, 15),
            lower = c(-Inf, 0),
            start_lower = c(ec50 = -7, hill = 0),
            start_upper = c(ec50 = -2, hill = 2),
            na.action = na.omit
              )
    )
    )

adjusted_fits_fixed_floor %>%
filter(map_lgl(fit, ~ .x$convInfo$isConv == TRUE)) %>%
mutate(tidied = map(fit, tidy)) %>%
unnest(tidied) %>%
select(construct, method, unique_experiment_id, term, estimate) %>%
pivot_wider(names_from = term, values_from = estimate) %>%
mutate(measure = "fluorescence", nucleotide = "TNP-ATP")  %>%
filter(ec50 < 0) -> a_fixed_tidy

saveRDS(unadjusted_fits_free_floor, "data/hill_fits/ua_free_fit.rds")
saveRDS(adjusted_fits_free_floor, "data/hill_fits/a_free_fit.rds")
saveRDS(unadjusted_fits_free_floor, "data/hill_fits/a_fixed_fit.rds")

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

interesting_constructs <- c("WT-GFP+SUR", "W311*-GFP+SUR")
interesting_nucleotides <- c("ATP", "TNP-ATP")

ua_free_tidy %>%
expand(nesting(construct, method, measure, nucleotide)) %>%
add_predicted_draws(model_1) %>%
group_by(construct, method, measure, nucleotide) %>%
filter(.prediction > quantile(.prediction, 0.025) & .prediction < quantile(.prediction, 0.975)) -> plot_1

ggplot() +
stat_slab(
    data = plot_1 %>% filter(construct %in% interesting_constructs, nucleotide %in% interesting_nucleotides),
    aes(y = construct, x = .prediction, colour = interaction(method, measure)), fill = NA, size = 1
    ) +
geom_point(
    data = ua_free_tidy %>% filter(construct %in% interesting_constructs, nucleotide %in% interesting_nucleotides),
    position = position_dodge(width=0.2), aes(y = construct, x = ec50, fill = interaction(method, measure)), shape=21, size = 2
    ) +
scale_colour_brewer(palette = "Pastel1", aesthetics = c("fill", "colour")) +
facet_grid(cols = vars(nucleotide))

