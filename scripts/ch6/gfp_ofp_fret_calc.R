library(tidyverse)
library(brms)
library(tidybayes)
library(ggdist)
library(distributional)

import <- function(directory, pattern) {

    oldwd <- getwd()

    setwd(directory)

    filenames <-
        list.files(pattern = pattern)

    tempdf <-
        filenames %>%
        plyr::llply(
            read_csv, col_types = cols(
                construct = col_character(),
                date = col_date(),
                cumulative_exposure = col_integer(),
                excitation_wavelength = col_double(),
                wavelength = col_double(),
                raw_spectra = col_double(),
                background_corrected_spectra = col_double()
            )
        ) %>%
        bind_rows(.id = "n") %>%
        type_convert()

    setwd(oldwd)

    return(tempdf)
}

spectra <-
    import("/home/sam/current_analysis/kir_sur_fret/whole_cell", "spectra") %>%
    filter(wavelength > 450, wavelength < 700) %>%
    select(n, construct, wavelength, background_corrected_spectra, excitation_wavelength)

ggplot(spectra, aes(x = wavelength, y = background_corrected_spectra, colour = factor(excitation_wavelength), group = interaction(excitation_wavelength, n))) +
geom_line() +
facet_wrap(vars(construct), scales = "free")

gfp_peak <-
    spectra %>%
    filter(wavelength > 516, wavelength < 518) %>%
    pivot_wider(names_from = excitation_wavelength, values_from = background_corrected_spectra) %>%
    group_by(n, construct) %>%
    summarise(gfp_intensity_max = mean(`490`))

mo_peak <-
    spectra %>%
    filter(wavelength > 573, wavelength < 575) %>%
    pivot_wider(names_from = excitation_wavelength, values_from = background_corrected_spectra) %>%
    group_by(n, construct) %>%
    summarise(mo_intensity_max = mean(`565`))

max_intensities <-
    left_join(gfp_peak, mo_peak) %>%
    mutate(max_ratio = gfp_intensity_max / mo_intensity_max)

spectra_summarised <-
    spectra %>%
    pivot_wider(names_from = excitation_wavelength, values_from = background_corrected_spectra) %>%
    left_join(gfp_peak) %>%
    group_by(n, construct) %>%
    mutate(
        norm_385 = `385` / gfp_intensity_max,
        norm_490 = `490` / gfp_intensity_max,
        norm_565 = `565` / gfp_intensity_max
        ) %>%
    select(construct, n, wavelength, norm_385, norm_490, norm_565) %>%
    pivot_longer(cols = c(norm_385, norm_490, norm_565), names_to = "excitation_wavelength", values_to = "intensity")

ggplot(spectra_summarised %>% filter(excitation_wavelength != "norm_385")) +
geom_line(aes(x = wavelength, y = intensity, colour = factor(excitation_wavelength), group = interaction(excitation_wavelength, n))) +
facet_wrap(vars(construct), scales = "free")

gfp_idealised <-
    spectra_summarised %>%
    filter(construct == "WT-GFP+SUR", excitation_wavelength == "norm_490") %>%
    ungroup() %>%
    select(wavelength, intensity) %>%
    group_by(wavelength) %>%
    summarise(
        gfp_intensity = mean(intensity)
        )

spectra_subtracted <-
    left_join(spectra_summarised, gfp_idealised) %>%
    filter(excitation_wavelength == "norm_490") %>%
    group_by(construct, n) %>%
    mutate(
        intensity = intensity - gfp_intensity
        ) %>%
    bind_rows(spectra_summarised %>% filter(excitation_wavelength == "norm_565")) %>%
    filter(wavelength > 525)

ggplot(spectra_subtracted) +
geom_line(aes(x = wavelength, y = intensity, colour = factor(excitation_wavelength), group = interaction(excitation_wavelength, n))) +
facet_wrap(vars(construct), scales = "free")

mo_ratio <-
    spectra_subtracted %>%
    ungroup() %>%
    select(construct, n, wavelength, excitation_wavelength, intensity) %>%
    filter(wavelength > 580, wavelength < 610) %>%
    pivot_wider(names_from = excitation_wavelength, values_from = intensity) %>%
    group_by(construct, n, wavelength) %>%
    summarise(
        ratio = norm_490 / norm_565
        ) %>%
    group_by(construct, n) %>%
    summarise(
        stdev = sd(ratio),
        ratio = mean(ratio)
        ) %>%
    filter(construct != "WT-GFP+SUR", n!=3)

control <-
    mo_ratio %>%
    filter(construct == "SUR_mOrange") %>%
    summarise(
        control_mean = mean(ratio)
        )

mo_ratio_centered <-
    mo_ratio %>%
    filter(construct != "SUR_mOrange") %>%
    left_join(max_intensities) %>%
    mutate(
        ratio_centered = (ratio - control %>% pull(control_mean)) / max_ratio,
        log_ratio_centered = log(ratio_centered)
        )

centered_prior_calc <-
    mo_ratio_centered %>%
    ungroup() %>%
    summarise(mu =  mean(log_ratio_centered), sigma = sd(log_ratio_centered))

ggplot() +
geom_point(data = mo_ratio_centered, aes(x = construct, y = log_ratio_centered, colour = factor(n))) +
stat_dist_slab(data = centered_prior_calc, aes(dist = dist_normal(mu, sigma))) +
coord_flip()

priors <- c(
    prior(normal(0.34, 0.25), class = b),
    prior(cauchy(0, 1), dpar = sigma, class = Intercept),
    prior(cauchy(0, 1), dpar = sigma, class = sd)
    )

form_1 <-
    bf(
        log_ratio_centered ~ 0 + construct,
        sigma ~ (1|control)
        )

brm(
    formula = form_1,
    family = gaussian(),
    data = mo_ratio_centered,
    prior = priors,
    cores = getOption("mc.cores", 4),
    sample_prior = "yes",
    save_all_pars = TRUE,
    chains = 4,
    iter = 40000,
    warmup= 20000,
    thin  = 10,
    control = list(adapt_delta = 0.99, max_treedepth = 15),
    file = "/home/sam/thesis/data/other_fits/gfp_ofp_fret.rds"
    ) -> model_1

mo_ratio_centered %>%
expand(construct, control) %>%
add_predicted_draws(model_1) %>%
ungroup() -> plot_1

ggplot() +
stat_slab(data = plot_1, aes(x = construct, y = .prediction, colour = construct), fill = NA, size = 1) +
geom_point(data = mo_ratio_centered, position = position_dodge(width=0.2), aes(x = construct, y = log_ratio_centered, fill = construct), shape=21, size = 2) +
coord_flip() +
scale_fill_brewer(palette = "Pastel1", aesthetics = c("colour", "fill"))

contrasts <-
    plot_1 %>%
    filter(construct == "SUR_mOrange") %>%
    select(.draw, .prediction) %>%
    rename(contrast_ctrl = ".prediction") %>%
    left_join(plot_1 %>% filter(construct != "SUR_mOrange")) %>%
    mutate(contrast = .prediction)

ggplot() +
stat_slab(data = contrasts, aes(x = construct, y = exp(contrast), fill = stat(cut_cdf_qi(cdf, .width = c(.5, .8, .95), labels = scales::percent_format())))) +
geom_point(data = mo_ratio_centered, position = position_dodge(width=0.2), aes(x = construct, y = exp(log_ratio_centered)), shape=21, size = 2) +
coord_flip() +
scale_fill_brewer(palette = "Blues", direction = -1, na.translate = FALSE) +
labs(fill = "Interval", x = "Construct", y = "Fold increase in FRET") +
scale_y_log10()