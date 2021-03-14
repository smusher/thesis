library(tidyverse)
library(broom)
library(nls.multstart)
library(brms)
library(tidybayes)
library(ggbeeswarm)
library(ggdist)
library(distributional)

concresp <-
    read_csv("data/combined_drc_data.csv") %>%
    mutate(log_concentration = log10(concentration)) %>%
    filter(log_concentration > -Inf) %>%
    group_by(unique_experiment_id) %>%
    mutate(observations = n()) %>%
    filter(observations >= 3) %>%
    ungroup() %>%
    mutate(
        response = case_when(
            response < -0.1 ~ -0.1,
            response > 1.1 ~ 1.1,
            TRUE ~ response
            )
        ) %>%
    drop_na()

currents <-
    concresp %>%
    filter(measure == "current") %>%
    unite("construct_nucleotide", c(construct, nucleotide))

fluorescence <-
    concresp %>%
    filter(measure == "fluorescence") %>%
    unite("construct_method", c(construct, method))

# hill_equation <- function(log_concentration, ec50, hill, floor){
#     floor + ((1 - floor) / (1 + 10^((ec50 - log_concentration) * -hill)))
# }

# #population fits
# unadjusted_fits_free_floor_pop <-
#     concresp %>%
#     group_by(construct, method, measure, nucleotide) %>%
#     nest() %>%
#     mutate(fit = map(data, ~ nls_multstart(
#             response ~ hill_equation(log_concentration, ec50, hill, floor),
#             data = data.frame(.),
#             iter = c(10,10,10),
#             lower = c(-Inf, 0, 0),
#             start_lower = c(ec50 = -7, hill = 0, floor= 0),
#             start_upper = c(ec50 = -1, hill = 2, floor = 1),
#             na.action = na.omit
#               )
#     )
#     )

# adjusted_fits_free_floor_pop <-
#     fluorescence %>%
#     group_by(construct, method, measure, nucleotide) %>%
#     nest() %>%
#     mutate(fit = map(data, ~ nls_multstart(
#             log2(response + 1) ~ hill_equation(log_concentration, ec50, hill, floor),
#             data = data.frame(.),
#             iter = c(10, 10, 10),
#             lower = c(-Inf, 0, 0),
#             start_lower = c(ec50 = -7, hill = 0, floor= 0),
#             start_upper = c(ec50 = -1, hill = 2, floor = 1),
#             na.action = na.omit
#               )
#     )
#     )

# adjusted_fits_fixed_floor_pop <-
#     fluorescence %>%
#     group_by(construct, method, measure, nucleotide) %>%
#     nest() %>%
#     mutate(fit = map(data, ~ nls_multstart(
#             log2(response + 1) ~ 0.1 + ((1 - 0.1) / (1 + 10^((ec50 - log_concentration) * -hill))),
#             data = data.frame(.),
#             iter = c(15, 15),
#             lower = c(-Inf, 0),
#             start_lower = c(ec50 = -7, hill = 0),
#             start_upper = c(ec50 = -2, hill = 2),
#             na.action = na.omit
#               )
#     )
#     )

# saveRDS(unadjusted_fits_free_floor_pop, "data/hill_fits/ua_free_fit_pop.rds")
# saveRDS(adjusted_fits_free_floor_pop, "data/hill_fits/a_free_fit_pop.rds")
# saveRDS(adjusted_fits_fixed_floor_pop, "data/hill_fits/a_fixed_fit_pop.rds")

# unadjusted_fits_free_floor_pop %>%
# mutate(tidied = map(fit, tidy)) %>%
# unnest(tidied) %>%
# select(construct, method, measure, nucleotide, term, estimate) %>%
# pivot_wider(names_from = term, values_from = estimate) %>%
# rename(pop_ec50 = ec50) %>%
# right_join(concresp) -> ua_free_tidy_pop

# adjusted_fits_free_floor_pop %>%
# mutate(tidied = map(fit, tidy)) %>%
# unnest(tidied) %>%
# select(construct, method, measure, nucleotide, term, estimate) %>%
# pivot_wider(names_from = term, values_from = estimate) %>%
# rename(pop_ec50 = ec50) %>%
# right_join(fluorescence) -> a_free_tidy_pop

# adjusted_fits_fixed_floor_pop %>%
# mutate(tidied = map(fit, tidy)) %>%
# unnest(tidied) %>%
# select(construct, method, measure, nucleotide, term, estimate) %>%
# pivot_wider(names_from = term, values_from = estimate) %>%
# rename(pop_ec50 = ec50) %>%
# mutate(floor = 0.1) %>%
# right_join(fluorescence) -> a_fixed_tidy_pop

# #individual fits
# unadjusted_fits_free_floor <-
#     ua_free_tidy_pop %>%
#     group_by(construct, method, measure, nucleotide, unique_experiment_id) %>%
#     nest() %>%
#     mutate(fit = map(data, ~ nls_multstart(
#             response ~ hill_equation(log_concentration, ec50, hill, floor),
#             data = data.frame(.),
#             iter = c(100),
#             lower = c(-Inf),
#             start_lower = c(ec50 = -7),
#             start_upper = c(ec50 = -1),
#             na.action = na.omit
#               )
#     )
#     )

# adjusted_fits_free_floor <-
#     a_free_tidy_pop %>%
#     group_by(construct, method, measure, nucleotide, unique_experiment_id) %>%
#     nest() %>%
#     mutate(fit = map(data, ~ nls_multstart(
#             log2(response + 1) ~ hill_equation(log_concentration, ec50, hill, floor),
#             data = data.frame(.),
#             iter = c(100),
#             lower = c(-Inf),
#             start_lower = c(ec50 = -7),
#             start_upper = c(ec50 = -1),
#             na.action = na.omit
#               )
#     )
#     )

# adjusted_fits_fixed_floor <-
#     a_fixed_tidy_pop %>%
#     group_by(construct, method, measure, nucleotide, unique_experiment_id) %>%
#     nest() %>%
#     mutate(fit = map(data, ~ nls_multstart(
#             log2(response + 1) ~ hill_equation(log_concentration, ec50, hill, floor),
#             data = data.frame(.),
#             iter = c(100),
#             lower = c(-Inf),
#             start_lower = c(ec50 = -7),
#             start_upper = c(ec50 = -1),
#             na.action = na.omit
#               )
#     )
#     )

# saveRDS(unadjusted_fits_free_floor, "data/hill_fits/ua_free_fit.rds")
# saveRDS(adjusted_fits_free_floor, "data/hill_fits/a_free_fit.rds")
# saveRDS(adjusted_fits_fixed_floor, "data/hill_fits/a_fixed_fit.rds")

ggplot() +
geom_point(data = concresp, aes(x = log_concentration, y = response)) +
stat_dist_slab(aes(y = 0, dist = dist_normal(-4, 1)))

ggplot() +
geom_point(data = concresp, aes(x = log_concentration, y = response, colour = unique_experiment_id)) +
scale_color_viridis_c() +
facet_wrap(vars(construct,method,measure,nucleotide), labeller = label_wrap_gen(multi_line = FALSE))

fixed_brms_formula <-
    bf(
        log2(response + 1) ~ 0.1 + ((1 - 0.1) / (1 + 10^((ec50 - log_concentration) * -hill))),
        nl = TRUE,
        ec50 ~ 0 + construct_method + (construct_method||unique_experiment_id),
        hill ~ 0 + construct_method,
        sigma ~ (1||construct_method),
        family = gaussian()
        )

fixed_brms_priors <-
    c(
        #ec50 based on concentration range
        prior(normal(-4, 1), nlpar = ec50, class = b),
        #hill coefficient between 0 and 2
        prior(normal(1, 0.3), nlpar = hill, class = b),
        #standard cauchy prior for sigmas
        prior(cauchy(0, 1), nlpar = ec50, class = sd),
        prior(cauchy(0, 1), dpar = sigma, class = Intercept),
        prior(cauchy(0, 1), dpar = sigma, class = sd)
        )

free_brms_formula <-
    bf(
        response ~ floor + ((1 - floor) / (1 + 10^((ec50 - log_concentration) * -hill))),
        nl = TRUE,
        ec50 ~ 0 + construct_nucleotide + (construct_nucleotide||unique_experiment_id),
        floor + hill ~ 0 + construct_nucleotide,
        sigma ~ (1||construct_nucleotide),
        family = gaussian()
        )

free_brms_priors <-
    c(
        #floor prior uniform
        prior(uniform(0, 1), nlpar = floor, class = b, lb=0, ub=1),
        #ec50 based on concentration range
        prior(normal(-4, 1), nlpar = ec50, class = b),
        #hill coefficient between 0 and 2
        prior(normal(1, 0.3), nlpar = hill, class = b),
        #standard cauchy prior for sigmas
        prior(cauchy(0, 1), nlpar = ec50, class = sd),
        prior(cauchy(0, 1), dpar = sigma, class = Intercept),
        prior(cauchy(0, 1), dpar = sigma, class = sd)
        )

brms_iter <- 2000
brms_warmup <- 1000
brms_chains <- 4
brms_thin <- 1
brms_seed <- 2021

brm(
    formula = fixed_brms_formula,
    prior = fixed_brms_priors,
    data = fluorescence,
    chains = brms_chains,
    iter = brms_iter,
    warmup = brms_warmup,
    thin  = brms_thin,
    seed = brms_seed,
    control = list(adapt_delta = 0.95, max_treedepth = 10),
    cores = getOption("mc.cores", 4),
    file = "data/hill_fits/fixed_fluorescence_110321_short",
    sample_prior = "yes",
    save_all_pars = TRUE
    ) -> test_run

brm(
    formula = free_brms_formula,
    prior = free_brms_priors,
    data = currents,
    chains = brms_chains,
    iter = brms_iter,
    warmup = brms_warmup,
    thin  = brms_thin,
    seed = brms_seed,
    control = list(adapt_delta = 0.95, max_treedepth = 10),
    cores = getOption("mc.cores", 4),
    file = "data/hill_fits/free_currents_110321_short",
    sample_prior = "yes",
    save_all_pars = TRUE
    ) -> test_run_2

fluorescence %>%
ungroup() %>%
expand(construct_method, log_concentration = seq(-8, -2, length.out = 51)) %>%
add_fitted_draws(test_run, re_formula = NA) %>%
median_qi(.value, .width = .95) -> inferred_underlying

ggplot() +
geom_ribbon(data = inferred_underlying, aes(x = 10^log_concentration, ymin = .lower, ymax = .upper, colour = construct_method, fill = construct_method), alpha = 0.5) +
geom_quasirandom(data = fluorescence, aes(x = 10^log_concentration, y = log2(response+1), fill = construct_method), size = 3, shape = 21, width=0.1) +
scale_color_viridis_d(aesthetics = c("colour", "fill")) +
scale_x_log10() +
scale_y_continuous(breaks = c(0, 0.5, 1), minor_breaks = c(0.25, 0.75)) +
facet_wrap(vars(construct_method), labeller = label_wrap_gen(multi_line = FALSE)) +
coord_cartesian(ylim = c(-0.1, 1.2)) +
theme(legend.position = "none")

currents %>%
ungroup() %>%
expand(construct_nucleotide, log_concentration = seq(-8, -2, length.out = 51)) %>%
add_fitted_draws(test_run_2, re_formula = NA) %>%
median_qi(.value, .width = .95) -> inferred_underlying_2

ggplot() +
geom_ribbon(data = inferred_underlying_2, aes(x = 10^log_concentration, ymin = .lower, ymax = .upper, colour = construct_nucleotide, fill = construct_nucleotide), alpha = 0.5) +
geom_quasirandom(data = currents, aes(x = 10^log_concentration, y = response, fill = construct_nucleotide), size = 3, shape = 21, width=0.1) +
scale_color_viridis_d(aesthetics = c("colour", "fill")) +
scale_x_log10() +
scale_y_continuous(breaks = c(0, 0.5, 1), minor_breaks = c(0.25, 0.75)) +
facet_wrap(vars(construct_nucleotide), labeller = label_wrap_gen(multi_line = FALSE)) +
coord_cartesian(ylim = c(-0.1, 1.2)) +
theme(legend.position = "none")

fluorescence %>%
ungroup() %>%
expand(nesting(construct_method, unique_experiment_id), log_concentration = seq(-8, -2, length.out = 51)) %>%
add_fitted_draws(test_run) %>%
median_qi(.value, .width = .95) -> per_experiment

ggplot() +
geom_ribbon(data = per_experiment, aes(x = 10^log_concentration, ymin = .lower, ymax = .upper, colour = construct_method, fill = construct_method), alpha = 0.5) +
geom_point(data = fluorescence, aes(x = 10^log_concentration, y = response, fill = construct_method), size = 1.5, shape = 21) +
scale_color_viridis_d(aesthetics = c("colour", "fill")) +
scale_x_log10() +
scale_y_continuous(breaks = c(0, 0.5, 1), minor_breaks = c(0.25, 0.75)) +
facet_wrap(vars(construct_method, unique_experiment_id)) +
coord_cartesian(ylim = c(-0.1, 1.2)) +
theme(legend.position = "none")

currents %>%
ungroup() %>%
expand(nesting(construct_nucleotide, unique_experiment_id), log_concentration = seq(-8, -2, length.out = 51)) %>%
add_fitted_draws(test_run_2) %>%
median_qi(.value, .width = .95) -> per_experiment_2

ggplot() +
geom_ribbon(data = per_experiment_2, aes(x = 10^log_concentration, ymin = .lower, ymax = .upper, colour = construct_nucleotide, fill = construct_nucleotide), alpha = 0.5) +
geom_point(data = currents, aes(x = 10^log_concentration, y = response, fill = construct_nucleotide), size = 1.5, shape = 21) +
scale_color_viridis_d(aesthetics = c("colour", "fill")) +
scale_x_log10() +
scale_y_continuous(breaks = c(0, 0.5, 1), minor_breaks = c(0.25, 0.75)) +
facet_wrap(vars(construct_nucleotide, unique_experiment_id)) +
coord_cartesian(ylim = c(-0.1, 1.2)) +
theme(legend.position = "none")

test_run %>%
gather_draws(`b_.*`, regex = TRUE) %>%
separate(.variable, into = c(".variable", "construct_method"), sep = "_construct_method") %>%
drop_na() %>%
mutate(
    construct_method = str_replace_all(
        construct_method,
        c(
            "MUM" = "*-",
            "MUP" = "*+",
            "MU" = "*,",
            "GFPP" = "GFP+",
            "MGFP" = "-GFP",
            "SURM" = "SUR-",
            "SURP" = "SUR+"
            )
        )
    ) -> f_draws_1

test_run %>%
gather_draws(r_unique_experiment_id__ec50[unique_experiment_id,]) %>%
ungroup() %>%
select(-.variable) %>%
rename(ind_deviation = .value) %>%
left_join(fluorescence %>% select(construct_method, unique_experiment_id) %>% unique()) %>%
left_join(f_draws_1 %>% filter(.variable == "b_ec50") %>% select(-.variable)) %>%
mutate(ind_value = .value + ind_deviation) %>%
group_by(construct_method, unique_experiment_id) %>%
summarise(ind_value = median(ind_value)) -> f_draws_2

ggplot() +
stat_slab(
    data = f_draws_1 %>% filter(.variable == "b_ec50"),
    aes(y = construct_method, x = .value, fill = stat(cut_cdf_qi(cdf, .width = c(.5, .8, .95), labels = scales::percent_format())))
    ) +
geom_point(
    data = f_draws_2,
    aes(y = construct_method, x = ind_value),
    fill = "white",
    shape = 21,
    size = 2
    ) +
scale_fill_brewer(palette = "Blues", direction = -1) +
theme(legend.position = "none")

test_run_2 %>%
gather_draws(`b_.*`, regex = TRUE) %>%
separate(.variable, into = c(".variable", "construct_nucleotide"), sep = "_construct_nucleotide") %>%
drop_na() %>%
mutate(
    construct_nucleotide = str_replace_all(
        construct_nucleotide,
        c(
            "MUM" = "*-",
            "MUP" = "*+",
            "MU" = "*,",
            "GFPP" = "GFP+",
            "MGFP" = "-GFP",
            "SURM" = "SUR-",
            "SURP" = "SUR+",
            "TNPMATP" = "TNP-ATP",
            "MGMATP" = "MG-ATP",
            "TMD0M" = "TMD0-"
            )
        )
    ) -> c_draws_1

test_run_2 %>%
gather_draws(r_unique_experiment_id__ec50[unique_experiment_id,]) %>%
ungroup() %>%
select(-.variable) %>%
rename(ind_deviation = .value) %>%
left_join(currents %>% select(construct_nucleotide, unique_experiment_id) %>% unique()) %>%
left_join(c_draws_1 %>% filter(.variable == "b_ec50") %>% select(-.variable)) %>%
mutate(ind_value = .value + ind_deviation) %>%
group_by(construct_nucleotide, unique_experiment_id) %>%
summarise(ind_value = median(ind_value)) -> c_draws_2

ggplot() +
stat_slab(
    data = c_draws_1 %>% filter(.variable == "b_ec50"),
    aes(y = construct_nucleotide, x = .value, fill = stat(cut_cdf_qi(cdf, .width = c(.5, .8, .95), labels = scales::percent_format())))
    ) +
geom_point(
    data = c_draws_2,
    aes(y = construct_nucleotide, x = ind_value),
    fill = "white",
    shape = 21,
    size = 2
    ) +
scale_fill_brewer(palette = "Blues", direction = -1) +
theme(legend.position = "none")

f_draws_1 %>%
filter(construct_method == "W311*-GFP+SUR_unroofed") %>%
ungroup() %>%
select(-construct_method) %>%
rename(unroofed_contrast = .value) -> unroofed_wt

f_draws_1 %>%
filter(construct_method == "W311*-GFP+SUR_pcf") %>%
ungroup() %>%
select(-construct_method) %>%
rename(pcf_contrast = .value) -> pcf_wt

f_draws_1 %>%
left_join(unroofed_wt) %>%
filter(!is.na(construct_method), construct_method != "W311*-GFP+SUR_unroofed") %>%
mutate(contrast = .value - unroofed_contrast) -> unroofed_contrasts

f_draws_1 %>%
left_join(pcf_wt) %>%
filter(!is.na(construct_method), construct_method != "W311*-GFP+SUR_pcf") %>%
mutate(contrast = .value - pcf_contrast) -> pcf_contrasts

ggplot() +
stat_slab(
    data = unroofed_contrasts %>% filter(.variable == "b_ec50"),
    aes(y = construct_method, x = contrast, fill = stat(cut_cdf_qi(cdf, .width = c(.5, .8, .95), labels = scales::percent_format())))
    ) +
scale_fill_brewer(palette = "Blues", direction = -1) +
theme(legend.position = "none")

ggplot() +
stat_slab(
    data = pcf_contrasts %>% filter(.variable == "b_ec50"),
    aes(y = construct_method, x = contrast, fill = stat(cut_cdf_qi(cdf, .width = c(.5, .8, .95), labels = scales::percent_format())))
    ) +
scale_fill_brewer(palette = "Blues", direction = -1) +
theme(legend.position = "none")

c_draws_1 %>%
filter(construct_nucleotide == "WT-GFP+SUR_ATP") %>%
ungroup() %>%
select(-construct_nucleotide) %>%
rename(contrast = .value) -> atp_wt

c_draws_1 %>%
filter(construct_nucleotide == "WT-GFP+SUR_TNP-ATP") %>%
ungroup() %>%
select(-construct_nucleotide) %>%
rename(contrast = .value) -> tnpatp_wt

c_draws_1 %>%
filter(construct_nucleotide == "W311*-GFP+SUR_ATP") %>%
ungroup() %>%
select(-construct_nucleotide) %>%
rename(contrast = .value) -> atp_w311

c_draws_1 %>%
filter(construct_nucleotide == "W311*-GFP+SUR_TNP-ATP") %>%
ungroup() %>%
select(-construct_nucleotide) %>%
rename(contrast = .value) -> tnpatp_w311

c_draws_1 %>%
left_join(atp_wt) %>%
filter(!is.na(construct_nucleotide), construct_nucleotide != "WT-GFP+SUR_ATP") %>%
mutate(contrast = .value - contrast) -> atp_wt_contrasts

c_draws_1 %>%
left_join(pcf_wt) %>%
filter(!is.na(construct_nucleotide), construct_nucleotide != "WT-GFP+SUR_TNP-ATP") %>%
mutate(contrast = .value - contrast) -> tnpatp_wt_contrasts

c_draws_1 %>%
left_join(atp_w311) %>%
filter(!is.na(construct_nucleotide), construct_nucleotide != "W311*-GFP+SUR_ATP") %>%
mutate(contrast = .value - contrast) -> atp_w311_contrasts

c_draws_1 %>%
left_join(tnpatp_w311) %>%
filter(!is.na(construct_nucleotide), construct_nucleotide != "W311*-GFP+SUR_TNP-ATP") %>%
mutate(contrast = .value - contrast) -> tnpatp_w311_contrasts

ggplot() +
stat_slab(
    data = atp_wt_contrasts %>% filter(.variable == "b_ec50"),
    aes(y = construct_nucleotide, x = contrast, fill = stat(cut_cdf_qi(cdf, .width = c(.5, .8, .95), labels = scales::percent_format())))
    ) +
scale_fill_brewer(palette = "Blues", direction = -1) +
theme(legend.position = "none")

ggplot() +
stat_slab(
    data = atp_tnpatp_contrasts %>% filter(.variable == "b_ec50"),
    aes(y = construct_method, x = contrast, fill = stat(cut_cdf_qi(cdf, .width = c(.5, .8, .95), labels = scales::percent_format())))
    ) +
scale_fill_brewer(palette = "Blues", direction = -1) +
theme(legend.position = "none")

c_draws_1 %>%
group_by(.variable, construct_nucleotide) %>%
median_qi(.value, .width = .95) %>%
ungroup() %>%
filter(.variable == "b_ec50") %>%
select(construct_nucleotide, .value, .lower, .upper) %>%
mutate(across(starts_with("."), ~(10^.x)*1e6)) %>%
print(width=500, n=100)

f_draws_1 %>%
group_by(.variable, construct_method) %>%
median_qi(.value, .width = .95) %>%
ungroup() %>%
filter(.variable == "b_ec50") %>%
select(construct_method, .value, .lower, .upper) %>%
mutate(across(starts_with("."), ~(10^.x)*1e6)) %>%
print(width=500, n=100)
