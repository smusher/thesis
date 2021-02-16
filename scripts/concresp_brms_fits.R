library(tidyverse)
library(broom)
library(brms)
library(tidybayes)
library(bayesplot)
library(modelr)
library(patchwork)
library(ggbeeswarm)
brms_iter <- 5000
brms_warmup <- 2000
brms_chains <- 4
brms_thin <- 1
brms_seed <- 2019

currents_1 <-
    read_csv("data/electrophys_data.csv") %>%
    select(n, method, measure, construct, nucleotide, concentration, response) %>%
    mutate(experimenter = "Sam")

currents_2 <-
    read_csv("data/natascias_data.csv") %>%
    select(n, method, measure, construct, nucleotide, concentration, response) %>%
    mutate(experimenter = "Natascia")

currents <-
    bind_rows(currents_1, currents_2) %>%
    group_by(experimenter, n, method, measure, construct, nucleotide) %>%
    mutate(unique_experiment_id = cur_group_id())

fluorescence <-
    read_csv("data/unroofed_concresp_data.csv") %>%
    filter(dye < 480) %>%
    select(unique_experiment_id, method, measure, construct, concentration, response) %>%
    mutate(nucleotide = "TNP-ATP", experimenter = "Sam")

pcf_data <-
    read_csv("data/pcf_data.csv") %>%
    filter(dye < 480 | is.na(dye)) %>%
    select(unique_experiment_id, method, measure, construct, concentration, response) %>%
    mutate(nucleotide = "TNP-ATP", experimenter = "Sam")

concresp <-
    bind_rows(currents, fluorescence, pcf_data) %>%
    mutate(log_concentration = log10(concentration)) %>%
    filter(log_concentration > -Inf) %>%
    group_by(unique_experiment_id, experimenter, method, measure, construct, nucleotide) %>%
    mutate(n = n()) %>%
    filter(n > 3)

currents <-
    concresp %>%
    filter(measure == "current") %>%
    ungroup() %>%
    mutate(
        construct = fct_relevel(factor(construct), "WT-GFP+SUR"),
        nucleotide = fct_relevel(factor(nucleotide), "ATP")
        )

fluorescence <-
    concresp %>%
    filter(measure == "fluorescence") %>%
    ungroup() %>%
    mutate(
        construct = fct_relevel(factor(construct), "W311*-GFP+SUR")
        )

hill_priors <-
    c(
        prior(uniform(-8, 0), nlpar = "ec50", class = b, lb = -8, ub = 0),
        prior(uniform(0, 2), nlpar = "hill", class = b, lb = 0, ub = 2),
        prior(uniform(0, 1), nlpar = "floor", class = b, lb = 0, ub = 1),
        prior(cauchy(0, 5), class = sd, nlpar = "ec50"),
        prior(cauchy(0, 5), class = sd, nlpar = "hill"),
        prior(cauchy(0, 5), class = sd, nlpar = "floor"),
        prior(cauchy(0, 5), class = sigma)
    )

hill_priors_2 <-
    c(
        prior(uniform(-8, 0), nlpar = "ec50", class = b, lb = -8, ub = 0),
        prior(uniform(0, 2), nlpar = "hill", class = b, lb = 0, ub = 2),
        prior(cauchy(0, 5), class = sd, nlpar = "ec50"),
        prior(cauchy(0, 5), class = sd, nlpar = "hill"),
        prior(cauchy(0, 5), class = sigma)
    )

brm_hill_eq_1a <-
    bf(#fully specified model - no obvious correlations apart from maybe ec50 and floor
        response ~ floor + ((1 - floor) / (1 + 10^((ec50 - log_concentration) * -hill))),
        ec50 + hill + floor ~ 0 + construct:nucleotide + (1|construct:nucleotide),
        nl = TRUE
        )

brm_hill_eq_1b <-
    bf(#without correlations
        response ~ floor + ((1 - floor) / (1 + 10^((ec50 - log_concentration) * -hill))),
        ec50 + hill + floor ~ 0 + construct:nucleotide + (1||construct:nucleotide),
        nl = TRUE
        )

brm_hill_eq_1c <-
    bf(#group effects only for ec50
        response ~ floor + ((1 - floor) / (1 + 10^((ec50 - log_concentration) * -hill))),
        ec50 ~ 0 + construct:nucleotide + (1|n),
        floor + hill ~ 0 + construct:nucleotide,
        nl = TRUE
        )

brm_hill_eq_1d <-
    bf(#hill fixed between constructs
        response ~ floor + ((1 - floor) / (1 + 10^((ec50 - log_concentration) * -hill))),
        ec50 ~ 0 + construct:nucleotide + (1|n),
        floor ~ 0 + construct:nucleotide,
        hill ~ 1,
        nl = TRUE
        )

brm_hill_eq_2a <-
    bf(#fully specified model - no obvious correlations apart from maybe ec50 and floor
        response ~ (1 / (1 + 10^((ec50 - log_concentration) * -hill))),
        ec50 + hill ~ 0 + construct:method + (1|ID|n),
        nl = TRUE
        )

brm_hill_eq_2b <-
    bf(#without correlations
        response ~ (1 / (1 + 10^((ec50 - log_concentration) * -hill))),
        ec50 + hill ~ 0 + construct:method + (1|n),
        nl = TRUE
        )

brm_hill_eq_2c <-
    bf(#group effects only for ec50
        response ~ (1 / (1 + 10^((ec50 - log_concentration) * -hill))),
        ec50 ~ 0 + construct:method + (1|n),
        hill ~ 0 + construct:method,
        nl = TRUE
        )

brm_hill_eq_2d <-
    bf(#hill fixed between constructs
        response ~ (1 / (1 + 10^((ec50 - log_concentration) * -hill))),
        ec50 ~ 0 + construct:method + (1|n),
        hill ~ 1,
        nl = TRUE
        )

hill_brm_fit_1a <-
    brm(
        data = currents,
        family = gaussian(),
        prior = hill_priors,
        formula = brm_hill_eq_1a,
        control = list(adapt_delta = 0.99, max_treedepth=15),
        iter = brms_iter,
        warmup = brms_warmup,
        chains = brms_chains,
        thin = brms_thin,
        seed = brms_seed,
        file = "data/hill_fits/1a_currents.rds",
        sample_prior = "yes",
        save_all_pars = TRUE,
        cores = 4
        )

hill_brm_fit_1b <-
    brm(
        data = currents,
        family = gaussian(),
        prior = hill_priors,
        formula = brm_hill_eq_1b,
        control = list(adapt_delta = 0.99, max_treedepth=15),
        iter = brms_iter,
        warmup = brms_warmup,
        chains = brms_chains,
        thin = brms_thin,
        seed = brms_seed,
        file = "data/hill_fits/1b_currents.rds",
        sample_prior = "yes",
        save_all_pars = TRUE,
        cores = 4
        )

hill_brm_fit_1c <-
    brm(
        data = currents,
        family = gaussian(),
        prior = hill_priors,
        formula = brm_hill_eq_1c,
        iter = brms_iter,
        warmup = brms_warmup,
        chains = brms_chains,
        thin = brms_thin,
        seed = brms_seed,
        file = "data/hill_fits/1c_currents.rds",
        sample_prior = "yes",
        save_all_pars = TRUE,
        cores = 4
        )

hill_brm_fit_1d <-
    brm(
        data = currents,
        family = gaussian(),
        prior = hill_priors,
        formula = brm_hill_eq_1d,
        control = list(adapt_delta = 0.99, max_treedepth=15),
        iter = brms_iter,
        warmup = brms_warmup,
        chains = brms_chains,
        thin = brms_thin,
        seed = brms_seed,
        file = "data/hill_fits/1d_currents.rds",
        sample_prior = "yes",
        save_all_pars = TRUE,
        cores = 4
        )

hill_brm_fit_2a <-
    brm(
        data = fluorescence,
        family = gaussian(),
        prior = hill_priors_2,
        formula = brm_hill_eq_2a,
        control = list(adapt_delta = 0.99, max_treedepth=15),
        iter = brms_iter,
        warmup = brms_warmup,
        chains = brms_chains,
        thin = brms_thin,
        seed = brms_seed,
        file = "data/hill_fits/2a_fluorescence.rds",
        sample_prior = "yes",
        save_all_pars = TRUE,
        cores = 4
        )

hill_brm_fit_2b <-
    brm(
        data = fluorescence,
        family = gaussian(),
        prior = hill_priors_2,
        formula = brm_hill_eq_2b,
        control = list(adapt_delta = 0.99, max_treedepth=15),
        iter = brms_iter,
        warmup = brms_warmup,
        chains = brms_chains,
        thin = brms_thin,
        seed = brms_seed,
        file = "data/hill_fits/2b_fluorescence.rds",
        sample_prior = "yes",
        save_all_pars = TRUE,
        cores = 4
        )

hill_brm_fit_2c <-
    brm(
        data = fluorescence,
        family = gaussian(),
        prior = hill_priors_2,
        formula = brm_hill_eq_2c,
        control = list(adapt_delta = 0.99, max_treedepth=15),
        iter = brms_iter,
        warmup = brms_warmup,
        chains = brms_chains,
        thin = brms_thin,
        seed = brms_seed,
        file = "data/hill_fits/2c_fluorescence.rds",
        sample_prior = "yes",
        save_all_pars = TRUE,
        cores = 4
        )

hill_brm_fit_2d <-
    brm(
        data = fluorescence,
        family = gaussian(),
        prior = hill_priors_2,
        formula = brm_hill_eq_2d,
        control = list(adapt_delta = 0.99, max_treedepth=15),
        iter = brms_iter,
        warmup = brms_warmup,
        chains = brms_chains,
        thin = brms_thin,
        seed = brms_seed,
        file = "data/hill_fits/2d_fluorescence.rds",
        sample_prior = "yes",
        save_all_pars = TRUE,
        cores = 4
        )

posterior <- as.array(hill_brm_fit_1c)

dimnames(posterior)[[3]] %>%
stringr::str_replace_all(
    c(
        "construct" = "",
        "MUM" = "*-",
        "MU" = "*,",
        "M" = "-",
        "PP" = "P+",
        "nucleotide" = ""
        )
    ) -> dimnames(posterior)[[3]]

add_criterion(
    hill_brm_fit_1a,
    criterion = "loo",
    save_psis = TRUE,
    moment_match=TRUE,
    file = "data/hill_fits/1a_currents.rds"
    ) -> loo1
add_criterion(
    hill_brm_fit_1b,
    criterion = "loo",
    save_psis = TRUE,
    moment_match=TRUE,
    file = "data/hill_fits/1b_currents.rds") -> loo2
add_criterion(
    hill_brm_fit_1c,
    criterion = "loo",
    save_psis = TRUE,
    moment_match=TRUE,
    file = "data/hill_fits/1c_currents.rds") -> loo3
add_criterion(
    hill_brm_fit_1d,
    criterion = "loo",
    save_psis = TRUE,
    moment_match=TRUE,
    file = "data/hill_fits/1d_currents.rds") -> loo4
loo(hill_brm_fit_2a, save_psis = TRUE) -> loo5
loo(hill_brm_fit_2b, save_psis = TRUE) -> loo6
loo(hill_brm_fit_2c, save_psis = TRUE) -> loo7
loo(hill_brm_fit_2d, save_psis = TRUE) -> loo8

loo_compare(loo1, loo2, loo3, loo4)
loo_compare(loo5, loo6, loo7, loo8)

interesting_constructs <- 
    c(
    "WT-GFP+SUR", "W311*-GFP+SUR",
    "E179A-GFP+SUR", "W311*,E179A-GFP+SUR",
    "E179K-GFP+SUR", "W311*,E179K-GFP+SUR",
    "K39A-GFP+SUR", "W311*,K39A-GFP+SUR",
    "K39R-GFP+SUR", "W311*,K39R-GFP+SUR",
    "K39E-GFP+SUR", "W311*,K39E-GFP+SUR"
    )

colour_scheme <-
    c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3", "#a6d854", "#ffd92f") %>%
    rep(each = 2)
names(colour_scheme) <- interesting_constructs

currents_2 <-
    currents %>%
    filter(construct %in% interesting_constructs, nucleotide != "MG-ATP") %>%
    mutate(Anap = case_when(stringr::str_detect(construct, "W311*") == TRUE ~ TRUE, TRUE ~ FALSE))

currents_2 %>%
data_grid(
    log_concentration = seq_range(-7:-2, n = 101),
    construct,
    nucleotide
    ) %>%
add_fitted_draws(hill_brm_fit_1b, re_formula = NA) %>%
group_by(construct, nucleotide, log_concentration) %>%
point_interval(.value, .point = median, .interval = qi, .width = .95) %>%
mutate(Anap = case_when(stringr::str_detect(construct, "W311*") == TRUE ~ TRUE, TRUE ~ FALSE)) -> fitted_draws

w311 <-
    fitted_draws %>%
    filter(construct == "W311*-GFP+SUR") %>%
    ungroup() %>%
    select(-construct)

wt <-
    fitted_draws %>%
    filter(construct == "WT-GFP+SUR") %>%
    ungroup() %>%
    select(-construct)

ggplot() +
geom_ribbon(data = wt, aes(x = log_concentration, ymin = .lower, ymax = .upper), fill = "#636363", alpha = 0.5) +
geom_ribbon(data = fitted_draws %>% filter(Anap == FALSE), aes(x = log_concentration, ymin = .lower, ymax = .upper), fill = "#3182bd", alpha = 0.5) +
geom_line(data = fitted_draws %>% filter(Anap == FALSE), aes(x = log_concentration, y = .value)) +
geom_quasirandom(data = currents_2 %>% filter(Anap == FALSE), aes(x = log_concentration, y = response), shape = 21, fill = "white", size = 2.5, width = 0.2) +
facet_grid(rows = vars(construct), cols = vars(nucleotide)) +
cowplot::theme_cowplot(font_family = "IBM Plex Sans Condensed", font_size = 9) +
theme(strip.background  = element_rect(fill = "white"), legend.position = "none") +
scale_y_continuous(breaks = c(0.0, 0.5, 1.0)) +
labs(
    x = "[nucleotide] (M)",
    y = expression(I/I[max])
    ) -> test_1

ggplot() +
geom_ribbon(data = w311, aes(x = log_concentration, ymin = .lower, ymax = .upper), fill = "#636363", alpha = 0.5) +
geom_ribbon(data = fitted_draws %>% filter(Anap == TRUE), aes(x = log_concentration, ymin = .lower, ymax = .upper), fill = "#3182bd", alpha = 0.5) +
geom_line(data = fitted_draws %>% filter(Anap == TRUE), aes(x = log_concentration, y = .value)) +
geom_quasirandom(data = currents_2 %>% filter(Anap == TRUE), aes(x = log_concentration, y = response), shape = 21, fill = "white", size = 2.5, width = 0.2) +
facet_grid(rows = vars(construct), cols = vars(nucleotide)) +
cowplot::theme_cowplot(font_family = "IBM Plex Sans Condensed", font_size = 9) +
theme(strip.background  = element_rect(fill = "white"), legend.position = "none") +
scale_y_continuous(breaks = c(0.0, 0.5, 1.0)) +
labs(
    x = "[nucleotide] (M)",
    y = expression(I/I[max])
    ) -> test_2

fluorescence_2 <-
    fluorescence %>%
    filter(construct %in% interesting_constructs, method == "pcf")

fluorescence_2 %>%
data_grid(
    log_concentration = seq_range(-7:-2, n = 101),
    construct,
    method
    ) %>%
add_fitted_draws(hill_brm_fit_2b, re_formula = NA) %>%
group_by(construct, log_concentration) %>%
point_interval(.value, .point = median, .interval = qi, .width = .95) -> fitted_draws_2

w311 <-
    fitted_draws_2 %>%
    filter(construct == "W311*-GFP+SUR") %>%
    ungroup() %>%
    select(-construct)

ggplot() +
geom_ribbon(data = w311, aes(x = log_concentration, ymin = .lower, ymax = .upper), fill = "#636363", alpha = 0.5) +
geom_ribbon(data = fitted_draws_2, aes(x = log_concentration, ymin = .lower, ymax = .upper), fill = "#3182bd", alpha = 0.5) +
geom_line(data = fitted_draws_2, aes(x = log_concentration, y = .value)) +
geom_quasirandom(data = fluorescence_2, aes(x = log_concentration, y = response), shape = 21, fill = "white", size = 2.5, width = 0.2) +
facet_grid(cols = vars(construct)) +
cowplot::theme_cowplot(font_family = "IBM Plex Sans Condensed", font_size = 9) +
theme(strip.background  = element_rect(fill = "white"), legend.position = "none") +
scale_y_continuous(breaks = c(0.0, 0.5, 1.0)) +
labs(
    x = "[TNP-ATP] (M)",
    y = expression(F/F[max])
    ) -> test_3

(test_1 + test_2)/ test_3

frame <-
currents_2 %>%
ungroup() %>%
select(n, construct, nucleotide, Anap) %>%
unique()

hill_brm_fit_1b %>%
recover_types(currents) %>%
gather_draws(regex=TRUE, `b_.*`) %>%
separate(.variable, into = c(".variable", "construct"), sep = "_construct") %>%
separate(construct, into = c("construct", "nucleotide"), sep = ":nucleotide") %>%
mutate(
    construct = stringr::str_replace_all(construct,
    c(
        "MUM" = "*-",
        "MU" = "*,",
        "M" = "-",
        "PP" = "P+"
        )),
    nucleotide = stringr::str_replace_all(nucleotide, "M", "-"),
    .variable = stringr::str_replace_all(.variable, "b_", "")
    ) %>%
rename(pop_value = .value) %>%
right_join(frame) %>%
mutate(construct = ordered(construct, levels = rev(interesting_constructs))) -> pop_fits

wt <-
    pop_fits %>%
    filter(construct %in% c("W311*-GFP+SUR", "WT-GFP+SUR")) %>%
    group_by(Anap, construct, nucleotide, .variable) %>%
    point_interval(pop_value, .point = median, .interval = qi, .width = .95)

hill_brm_fit_1b %>%
recover_types(currents) %>%
gather_draws(c(r_n__ec50, r_n__hill, r_n__floor)[n,]) %>%
mutate(.variable = stringr::str_replace_all(.variable, "r_n__", "")) %>%
rename(group_deviation = .value) %>%
right_join(pop_fits) %>%
mutate(group_value = pop_value + group_deviation) %>%
group_by(n, Anap, construct, nucleotide, .variable) %>%
point_interval(group_value, .point = median, .interval = qi, .width = .95) %>%
mutate(construct = ordered(construct, levels = rev(interesting_constructs))) -> group_fits

fancy_scientific <- function(l) {
    l <- paste("10^", l, sep = "")
    parse(text=l)
}

ggplot() +
stat_slab(
    data = pop_fits %>% filter(.variable == "ec50", Anap == FALSE),
    aes(
        x = construct,
        y = pop_value,
        fill = construct,
        ),
    colour = "black", size = 0.5
    ) +
geom_quasirandom(
    data = group_fits %>% filter(.variable == "ec50", Anap == FALSE),
    aes(
        x = construct,
        y = group_value,
        fill = construct),
    shape = 21, width = 0.1, size = 2.5, stroke = 0.5
    ) +
geom_hline(
    data = wt %>% filter(.variable == "ec50", Anap == FALSE),
    aes(yintercept = .upper),
    linetype = 2
    ) +
geom_hline(
    data = wt %>% filter(.variable == "ec50", Anap == FALSE),
    aes(yintercept = .lower),
    linetype = 2
    ) +
facet_grid(rows = vars(nucleotide)) +
coord_flip() +
cowplot::theme_cowplot(font_family = "IBM Plex Sans Condensed", font_size = 9) +
theme(
    strip.background  = element_rect(fill = "white"),
    legend.position = "none",
    axis.title.y = element_blank()
    ) +
scale_colour_manual(values = colour_scheme, aesthetics = c('colour', 'fill')) +
scale_y_continuous(limits = c(-7, -2.5), labels = fancy_scientific) +
labs(y = "Fitted group-level and population-level IC50 parameter values (M)") -> test_4

ggplot() +
stat_slab(
    data = pop_fits %>% filter(.variable == "ec50", Anap == TRUE),
    aes(
        x = construct,
        y = pop_value,
        fill = construct,
        ),
    colour = "black", size = 0.5
    ) +
geom_quasirandom(
    data = group_fits %>% filter(.variable == "ec50", Anap == TRUE),
    aes(
        x = construct,
        y = group_value,
        fill = construct),
    shape = 21, width = 0.1, size = 2.5, stroke = 0.5
    ) +
geom_hline(
    data = wt %>% filter(.variable == "ec50", Anap == TRUE),
    aes(yintercept = .upper),
    linetype = 2
    ) +
geom_hline(
    data = wt %>% filter(.variable == "ec50", Anap == TRUE),
    aes(yintercept = .lower),
    linetype = 2
    ) +
facet_grid(rows = vars(nucleotide)) +
coord_flip() +
cowplot::theme_cowplot(font_family = "IBM Plex Sans Condensed", font_size = 9) +
theme(
    strip.background  = element_rect(fill = "white"),
    legend.position = "none",
    axis.title.y = element_blank()
    ) +
scale_colour_manual(values = colour_scheme, aesthetics = c('colour', 'fill')) +
scale_y_continuous(limits = c(-7, -2.5), labels = fancy_scientific) +
labs(y = "Fitted group-level and population-level IC50 parameter values (M)")  -> test_5

frame <-
fluorescence_2 %>%
ungroup() %>%
select(n, construct, method) %>%
unique()

hill_brm_fit_2b %>%
recover_types(fluorescence) %>%
gather_draws(regex=TRUE, `b_.*`) %>%
separate(.variable, into = c(".variable", "construct"), sep = "_construct") %>%
separate(construct, into = c("construct", "method"), sep = ":method") %>%
mutate(
    construct = stringr::str_replace_all(construct,
    c(
        "MUM" = "*-",
        "MU" = "*,",
        "M" = "-",
        "PP" = "P+"
        )),
    .variable = stringr::str_replace_all(.variable, "b_", "")
    ) %>%
rename(pop_value = .value) %>%
right_join(frame) %>%
mutate(construct = ordered(construct, levels = rev(interesting_constructs))) -> pop_fits

wt <-
    pop_fits %>%
    filter(construct %in% c("W311*-GFP+SUR", "WT-GFP+SUR")) %>%
    group_by(construct, method, .variable) %>%
    point_interval(pop_value, .point = median, .interval = qi, .width = .95)

hill_brm_fit_2b %>%
recover_types(fluorescence) %>%
gather_draws(r_n__ec50[n,]) %>%
mutate(.variable = stringr::str_replace_all(.variable, "r_n__", "")) %>%
rename(group_deviation = .value) %>%
right_join(pop_fits) %>%
mutate(group_value = pop_value + group_deviation) %>%
group_by(n, construct, method, .variable) %>%
point_interval(group_value, .point = median, .interval = qi, .width = .95) %>%
mutate(construct = ordered(construct, levels = rev(interesting_constructs))) -> group_fits

ggplot() +
stat_slab(
    data = pop_fits %>% filter(.variable == "ec50", method == "pcf"),
    aes(
        x = construct,
        y = pop_value,
        fill = construct,
        ),
    colour = "black", size = 0.5
    ) +
geom_quasirandom(
    data = group_fits %>% filter(.variable == "ec50", method == "pcf"),
    aes(
        x = construct,
        y = group_value,
        fill = construct),
    shape = 21, width = 0.1, size = 2.5, stroke = 0.5
    ) +
geom_hline(
    data = wt %>% filter(.variable == "ec50", method == "pcf"),
    aes(yintercept = .upper),
    linetype = 2
    ) +
geom_hline(
    data = wt %>% filter(.variable == "ec50", method == "pcf"),
    aes(yintercept = .lower),
    linetype = 2
    ) +
coord_flip() +
cowplot::theme_cowplot(font_family = "IBM Plex Sans Condensed", font_size = 9) +
theme(
    strip.background  = element_rect(fill = "white"),
    legend.position = "none",
    axis.title.y = element_blank()
    ) +
scale_colour_manual(values = colour_scheme, aesthetics = c('colour', 'fill')) +
scale_y_continuous(limits = c(-7, -2.5), labels = fancy_scientific) +
labs(y = "Fitted group-level and population-level EC50 parameter values (M)") -> test_6

test_4 + test_5 + (plot_spacer() / test_6)

frame <-
    expand.grid(
        construct = interesting_constructs,
        nucleotide = c("TNP-ATP", "ATP"),
        log_concentration = -8
        )

current_post_preds <-
    fitted(
        hill_brm_fit_1b,
        newdata = frame,
        nlpar = "ec50",
        summary = FALSE,
        re_formula = NA
    ) %>%
    as_tibble() %>%
    set_names(unite(frame, names, construct:nucleotide) %>% pull(names)) 

wt_atp_contrasts <-
    current_post_preds %>%
    select(!contains("TNP-ATP") & !contains("W311*")) %>%
    mutate(across(!`WT-GFP+SUR_ATP`, ~ .x - `WT-GFP+SUR_ATP`)) %>%
    select(-`WT-GFP+SUR_ATP`) %>%
    pivot_longer(everything(), names_to = "construct") %>%
    mutate(nucleotide = "ATP", construct = str_replace_all(construct, "_ATP", ""))

wt_tnp_atp_contrasts <-
    current_post_preds %>%
    select(contains("TNP-ATP") & !contains("W311*")) %>%
    mutate(across(!`WT-GFP+SUR_TNP-ATP`, ~ .x - `WT-GFP+SUR_TNP-ATP`)) %>%
    select(-`WT-GFP+SUR_TNP-ATP`) %>%
    pivot_longer(everything(), names_to = "construct") %>%
    mutate(nucleotide = "TNP-ATP", construct = str_replace_all(construct, "_TNP-ATP", ""))

wt_contrasts <-
    bind_rows(wt_atp_contrasts, wt_tnp_atp_contrasts)

anap_atp_contrasts <-
    current_post_preds %>%
    select(!contains("TNP-ATP") & contains("W311*")) %>%
    mutate(across(!`W311*-GFP+SUR_ATP`, ~ .x - `W311*-GFP+SUR_ATP`)) %>%
    select(-`W311*-GFP+SUR_ATP`) %>%
    pivot_longer(everything(), names_to = "construct") %>%
    mutate(nucleotide = "ATP", construct = str_replace_all(construct, "_ATP", ""))

anap_tnp_atp_contrasts <-
    current_post_preds %>%
    select(contains("TNP-ATP") & contains("W311*")) %>%
    mutate(across(!`W311*-GFP+SUR_TNP-ATP`, ~ .x - `W311*-GFP+SUR_TNP-ATP`)) %>%
    select(-`W311*-GFP+SUR_TNP-ATP`) %>%
    pivot_longer(everything(), names_to = "construct") %>%
    mutate(nucleotide = "TNP-ATP", construct = str_replace_all(construct, "_TNP-ATP", ""))

anap_contrasts <-
    bind_rows(anap_atp_contrasts, anap_tnp_atp_contrasts)

frame <-
    tibble(
        construct = str_subset(interesting_constructs, "W311*"),
        method = "pcf",
        log_concentration = -8
        )

fluorescence_post_preds <-
    fitted(
        hill_brm_fit_2b,
        newdata = frame,
        nlpar = "ec50",
        summary = FALSE,
        re_formula = NA,
        allow_new_levels = TRUE
    ) %>%
    as_tibble() %>%
    set_names(unite(frame, names, construct:method) %>% pull(names)) 

fluorescence_contrasts <-
    fluorescence_post_preds %>%
    mutate(across(!`W311*-GFP+SUR_pcf`, ~ .x - `W311*-GFP+SUR_pcf`)) %>%
    select(-`W311*-GFP+SUR_pcf`) %>%
    pivot_longer(everything(), names_to = "construct") %>%
    mutate(construct = str_replace_all(construct, "_pcf", ""))

ggplot() +
geom_vline(xintercept = 0) +
geom_histogram(data = wt_contrasts, aes(x = value, fill = construct, after_stat(density)), bins = 40) +
facet_grid(rows = vars(construct), cols = vars(nucleotide)) +
cowplot::theme_cowplot(font_family = "IBM Plex Sans Condensed", font_size = 9) +
theme(
    strip.background  = element_rect(fill = "white"),
    legend.position = "none"
    ) +
scale_colour_manual(values = colour_scheme, aesthetics = c('colour', 'fill')) +
labs(x = "Order of magnitude shift in IC50") -> wt_contrast_plot

ggplot() +
geom_vline(xintercept = 0) +
geom_histogram(data = anap_contrasts, aes(x = value, fill = construct, after_stat(density)), bins = 40) +
facet_grid(rows = vars(construct), cols = vars(nucleotide)) +
cowplot::theme_cowplot(font_family = "IBM Plex Sans Condensed", font_size = 9) +
theme(
    strip.background  = element_rect(fill = "white"),
    legend.position = "none"
    ) +
scale_colour_manual(values = colour_scheme, aesthetics = c('colour', 'fill')) +
labs(x = "Order of magnitude shift in IC50") -> anap_contrast_plot

ggplot() +
geom_vline(xintercept = 0) +
geom_histogram(data = fluorescence_contrasts, aes(x = value, fill = construct, after_stat(density)), bins = 40) +
facet_grid(rows = vars(construct)) +
cowplot::theme_cowplot(font_family = "IBM Plex Sans Condensed", font_size = 9) +
theme(
    strip.background  = element_rect(fill = "white"),
    legend.position = "none"
    ) +
scale_colour_manual(values = colour_scheme, aesthetics = c('colour', 'fill')) +
labs(x = "Order of magnitude shift in EC50") -> fluorescence_contrast_plot

wt_contrast_plot / anap_contrast_plot -> patch_1

patch_1 - (plot_spacer() / fluorescence_contrast_plot)
