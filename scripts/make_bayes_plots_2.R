library(tidyverse)
library(broom)
library(brms)
library(tidybayes)
library(cowplot)
library(ggbeeswarm)
library(ggforce)
library(knitr)
library(kableExtra)
extrafont::loadfonts()

filenames_new_priors <-
    list.files(path = "data/mwc_fits_new_priors", full.names = T)

models_new_priors <-
    filenames_new_priors %>%
    plyr::llply(readRDS)

names_new_priors <-
    filenames_new_priors %>%
    str_extract_all("(?<=mwc_fits_new_priors\\/)[^.]*") %>%
    unlist()

filenames_old_priors <-
    list.files(path = "data/mwc_fits", full.names = T)

models_old_priors <-
    filenames_old_priors %>%
    plyr::llply(readRDS)

names_old_priors <-
    filenames_old_priors %>%
    str_extract_all("(?<=mwc_fits\\/)[^.]*") %>%
    unlist()

names <-
    names_new_priors %>%
    append(names_old_priors)

models <-
    models_new_priors %>%
    append(models_old_priors)

model_fits <-
    tibble(name = names, fit = models) %>%
    separate(name, c("construct", "model", "hierarchy", "priors"), sep = "_") %>%
    filter(hierarchy == "population", model == "mwc") %>%
    mutate(priors = case_when(is.na(priors) ~ "original", TRUE ~ priors)) %>%
    select(-hierarchy)

posterior_spread <-
    model_fits %>%
    group_by(construct, model, priors) %>%
    mutate(fit = map(fit, spread_draws, `b_.*`, regex = T)) %>%
    unnest(fit) %>%
    rename(log_L = b_logL_Intercept, D = b_D_Intercept, log_Ka = b_logKa_Intercept) %>%
    mutate(probability = "posterior")

prior_spread <-
    model_fits %>%
    group_by(construct, model, priors) %>%
    mutate(fit = map(fit, spread_draws, `prior_b_.*`, regex = T)) %>%
    unnest(fit) %>%
    rename(log_L = prior_b_logL, D = prior_b_D, log_Ka = prior_b_logKa) %>%
    mutate(probability = "prior")

spreaded_draws <-
    bind_rows(posterior_spread, prior_spread) %>%
    ungroup() %>%
    select(construct, model, probability, priors, D, log_Ka, log_L)

ggplot(spreaded_draws %>% filter(probability == "posterior", construct == "W311*,C166S-GFP+SUR")) +
geom_autodensity(colour = "black", position = "identity", alpha = 0.5, aes(fill = priors)) +
stat_density_2d(geom = "contour", aes(x = .panel_x, y = .panel_y, colour = priors)) +
facet_matrix(vars(D:log_L), layer.lower = 2, layer.diag = 1, layer.upper = 3) +
theme_cowplot(font_family = "IBM Plex Sans Condensed", font_size = 9) +
theme(strip.background  = element_rect(fill = "white"), legend.position = "top")

gathered_draws <-
    spreaded_draws %>%
    pivot_longer(cols = D:log_L, names_to = "parameter", values_to =  "estimate")

ggplot(gathered_draws %>% filter(construct == "W311*-GFP+SUR" | probability == "prior")) +
geom_density(aes(estimate, colour = priors, linetype = probability), fill = NA) +
facet_wrap(vars(parameter), scales = "free") +
theme_cowplot(font_family = "IBM Plex Sans Condensed", font_size = 9) +
theme(strip.background  = element_rect(fill = "white"), legend.position = "top")

gathered_draws_2 <-
    gathered_draws %>%
    filter(construct == "W311*-GFP+SUR") %>%
    select(-priors)

gathered_draws_3 <-
    gathered_draws %>%
    filter(construct != "W311*-GFP+SUR")

ggplot() +
geom_density(data = gathered_draws %>% filter(probability == "posterior"), aes(estimate, fill = construct), colour = "white", alpha = 0.5) +
facet_wrap(vars(parameter, priors), scales = "free", ncol = 3) +
theme_cowplot(font_family = "IBM Plex Sans Condensed", font_size = 9) +
theme(strip.background  = element_rect(fill = "white"), legend.position = "top")

ggplot() +
geom_density(data = gathered_draws_2 %>% filter(probability == "posterior"), aes(estimate), colour = "white", fill = "grey", alpha = 0.5) +
geom_density(data = gathered_draws_3 %>% filter(probability == "posterior"), aes(estimate, fill = construct), colour = "white", alpha = 0.5) +
facet_grid(cols = vars(parameter), rows = vars(priors), scales = "free") +
theme_cowplot(font_family = "IBM Plex Sans Condensed", font_size = 9) +
theme(strip.background  = element_rect(fill = "white"), legend.position = "top")

pcf_data <-
    read_csv("data/pcf_data.csv") %>%
    filter(dye < 480 | is.na(dye), concentration > 0) %>%
    select(unique_experiment_id, construct, measure, method, concentration, response)

pcf_data_summary <-
    pcf_data %>%
    group_by(construct, measure, concentration) %>%
    summarise(
        se = sd(response)/sqrt(length(response)),
        response = mean(response)
    )

x <- c(10^seq(-7, -2, length.out = 51))
z <- c(0, 1)

frame <-
    expand.grid(
        concentration = x,
        binding_mask = z,
        unique_experiment_id = 1
    ) %>%
    as_tibble()

fits <-
    model_fits %>%
    mutate(fit = map(fit, fitted_draws, newdata = frame, re_formula = NA, allow_new_levels = TRUE)) %>%
    unnest(fit) %>%
    group_by(construct, model, priors, binding_mask, concentration) %>%
    point_interval(.value, .point = median, .interval = qi) %>%
    mutate(
        measure = case_when(binding_mask == 1 ~ "fluorescence", TRUE ~ "current"),
        response = case_when(measure == "fluorescence" ~ 1 - .value, TRUE ~ .value),
        .lower = case_when(measure == "fluorescence" ~ 1 - .lower, TRUE ~ .lower),
        .upper = case_when(measure == "fluorescence" ~ 1 - .upper, TRUE ~ .upper)
    )

fancy_scientific <- function(l) {
    l <- format(l, scientific = TRUE)
    l <- gsub("^(.*)e", "10^", l)
    parse(text=l)
}

split_construct <- function(l) {
    gsub("[+]", "\n +", l)
}

fits <-
    fits %>%
    ungroup %>%
    mutate(
    	model = case_when(model == "coop" ~ "Negative Cooperativity", model == "mwc" ~ "Full MWC"),
        priors = case_when(is.na(priors) ~ "original", TRUE ~ priors)
    	)

ggplot() +
geom_ribbon(data = fits, aes(x = concentration, ymin = .lower, ymax = .upper, fill = measure), alpha = 0.25) +
geom_quasirandom(data = pcf_data, aes(x = concentration, y = response, fill = measure), shape = 21, alpha = 0.25, size = 1.5, width = 0.2) +
geom_line(data = fits, aes(x = concentration, y = response, colour = measure), size = 1) +
geom_errorbar(data = pcf_data_summary, aes(x = concentration, ymin = response - se, ymax = response + se), width = 0.2) +
geom_point(data = pcf_data_summary, aes(x = concentration, y = response, fill = measure), shape = 21, size = 3) +
facet_grid(rows = vars(construct), cols = vars(priors), labeller = labeller(construct = split_construct)) +
theme_cowplot(font_family = "IBM Plex Sans Condensed", font_size = 12) +
theme(strip.background  = element_rect(fill = "white"), legend.position = "top") +
scale_colour_manual(values = c("#1E88E5", "#FB8C00"), aesthetics = c("colour", "fill")) +
scale_x_log10(labels = fancy_scientific) +
scale_y_continuous(breaks = c(0.0, 0.5, 1.0)) +
labs(
    x = "[TNP-ATP] (M)",
    y = expression(F/F[max] ~ or ~ I/I[max])
    ) -> fit_plot_models

fits_2 <-
    fits %>%
    filter(priors == "original", construct != "W311*-GFP+SUR") %>%
    select(-priors)

fits_3 <-
    fits %>%
    filter(priors != "original", construct != "W311*-GFP+SUR")

ggplot() +
geom_ribbon(data = fits_3 %>% filter(model == "Full MWC"), aes(x = concentration, ymin = .lower, ymax = .upper, fill = measure), alpha = 0.25) +
geom_quasirandom(data = pcf_data %>% filter(construct != "W311*-GFP+SUR"), aes(x = concentration, y = response, fill = measure), shape = 21, alpha = 0.25, size = 1.5, width = 0.2) +
geom_line(data = fits_3 %>% filter(model == "Full MWC"), aes(x = concentration, y = response, colour = measure), size = 1) +
geom_line(data = fits_2 %>% filter(model == "Full MWC"), aes(x = concentration, y = response, colour = measure), size = 1, linetype = 2) +
geom_errorbar(data = pcf_data_summary %>% filter(construct != "W311*-GFP+SUR"), aes(x = concentration, ymin = response - se, ymax = response + se), width = 0.2) +
geom_point(data = pcf_data_summary %>% filter(construct != "W311*-GFP+SUR"), aes(x = concentration, y = response, fill = measure), shape = 21, size = 3) +
facet_grid(cols = vars(construct), rows = vars(priors), labeller = labeller(construct = split_construct)) +
theme_cowplot(font_family = "IBM Plex Sans Condensed", font_size = 12) +
theme(strip.background  = element_rect(fill = "white"), legend.position = "top") +
scale_colour_manual(values = c("#1E88E5", "#FB8C00"), aesthetics = c("colour", "fill")) +
scale_x_log10(labels = fancy_scientific) +
scale_y_continuous(breaks = c(0.0, 0.5, 1.0)) +
labs(
    x = "[TNP-ATP] (M)",
    y = expression(F/F[max] ~ or ~ I/I[max])
    ) -> fit_plot_mwc

ggplot() +
geom_ribbon(data = fits %>% filter(model == "Negative Cooperativity"), aes(x = concentration, ymin = .lower, ymax = .upper, fill = measure), alpha = 0.25) +
geom_quasirandom(data = pcf_data, aes(x = concentration, y = response, fill = measure), shape = 21, alpha = 0.25, size = 1.5, width = 0.2) +
geom_line(data = fits %>% filter(model == "Negative Cooperativity"), aes(x = concentration, y = response, colour = measure), size = 1) +
geom_errorbar(data = pcf_data_summary, aes(x = concentration, ymin = response - se, ymax = response + se), width = 0.2) +
geom_point(data = pcf_data_summary, aes(x = concentration, y = response, fill = measure), shape = 21, size = 3) +
facet_grid(cols = vars(construct), rows = vars(priors), labeller = labeller(construct = split_construct)) +
theme_cowplot(font_family = "IBM Plex Sans Condensed", font_size = 12) +
theme(strip.background  = element_rect(fill = "white"), legend.position = "top") +
scale_colour_manual(values = c("#1E88E5", "#FB8C00"), aesthetics = c("colour", "fill")) +
scale_x_log10(labels = fancy_scientific) +
scale_y_continuous(breaks = c(0.0, 0.5, 1.0)) +
labs(
    x = "[TNP-ATP] (M)",
    y = expression(F/F[max] ~ or ~ I/I[max])
    ) -> fit_plot_coop


gathered_draws %>%
filter(parameter == "log_L", priors == "broader", model == "mwc", probability == "posterior") %>%
group_by(construct) %>%
mutate(rank = percent_rank(estimate)) %>%
filter(rank < 0.51) %>%
filter(estimate == max(estimate))