library(tidyverse)
library(broom)
library(nls.multstart)
library(patchwork)
library(ggbeeswarm)

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
                concentration = col_double(),
                dye = col_double(),
                response = col_double(),
                raw_intensity = col_double(),
                corrected_intensity = col_double()
            )
        ) %>%
        bind_rows(.id = "n") %>%
        type_convert()

    setwd(oldwd)

    return(tempdf)
}

data_1 <-
    import("/home/sam/current_analysis/paper/pcf/W311-GFP+SUR/reanalyse", "intensities") %>%
    mutate(method = "pcf")

data_2 <-
    import("/home/sam/current_analysis/paper/tnp-atp_doseresponse_unroofed/W311-GFP+SUR", "intensities") %>%
    mutate(method = "unroofed")

w311_fluorescence <-
    bind_rows(data_1, data_2) %>%
    filter(dye < 480) %>%
    group_by(method, n) %>%
    mutate(
        normalised_intensity = raw_intensity / max(raw_intensity),
        cumulative_exposure = cumulative_exposure / 1000
        ) %>%
    select(method, n, cumulative_exposure, concentration, raw_intensity, normalised_intensity)

bleaching_equation <- function(time, a, k){
    a * exp(time * -k) + (1 - a)
}

pcf_bleaching_fits <-
	w311_fluorescence %>%
	filter(concentration == 0) %>%
	group_by(method, n) %>%
	nest() %>%
	mutate(fit = map(data, ~ nls_multstart(
            normalised_intensity ~ bleaching_equation(cumulative_exposure, a, k),
            data = data.frame(.),
            iter = 500,
            lower = c(a = 0, k = 0),
            start_lower = c(a = 0, k = 0.01),
            start_upper = c(a = 2, k = 1),
            na.action = na.omit
              )
    )
    )

nd_1 <-
	tibble(
		cumulative_exposure = seq(min(w311_fluorescence$cumulative_exposure), max(w311_fluorescence$cumulative_exposure), length.out = 51)
	)

pcf_bleaching_plots <-
	pcf_bleaching_fits %>%
	mutate(fit = map(fit, augment, newdata = nd_1)) %>%
	unnest(fit)

fancy_scientific <- function(l) {
    l <- format(l, scientific = TRUE)
    l <- gsub("^(.*)e", "10^", l)
    l <- gsub("0.001", "10^-3", l)
    parse(text=l)
}

ggplot() +
geom_point(data = w311_fluorescence, aes(x = cumulative_exposure, y = normalised_intensity, fill = factor(concentration)), shape = 21, size = 3) +
geom_line(data = pcf_bleaching_plots , aes(x = cumulative_exposure, y = .fitted), linetype = 2, size = 1) +
scale_fill_brewer(
	palette = "PuOr",
	direction = -1,
	labels = fancy_scientific,
    aesthetics = c("fill", "colour")
	) +
facet_wrap(vars(method, n)) +
labs(
	x = "Cumulative exposure (seconds)",
	y = "Normalised fluorescence intensity",
	fill = "[TNPATP] (M)"
	) +
theme(legend.position = "bottom")

pcf_bleaching_tidy <-
    pcf_bleaching_fits %>%
    mutate(fit = map(fit, tidy)) %>%
    unnest(fit)

ggplot() +
geom_quasirandom(data = pcf_bleaching_tidy, aes(x = term, y = estimate, fill = term), shape = 21, size = 3, width = 0.2) +
scale_fill_brewer(palette = "Set2") +
theme(legend.position = "bottom") +
facet_wrap(vars(method)) +
scale_y_log10()

pcf_data_min <-
    w311_fluorescence %>%
    group_by(method, n) %>%
    filter(concentration == 0) %>%
    filter(cumulative_exposure == max(cumulative_exposure)) %>%
    group_by(method) %>%
    mutate(average = mean(normalised_intensity))

ggplot() +
geom_quasirandom(data = pcf_data_min, aes(x = concentration, y = normalised_intensity), shape = 21, size = 3, width = 0.1, fill = "white") +
scale_fill_brewer(palette = "Set2") +
theme(legend.position = "bottom") +
coord_cartesian(ylim = c(0, 1)) +
facet_wrap(vars(method)) +
geom_hline(data = pcf_data_min, aes(yintercept = average))

mean(pcf_data_min$normalised_intensity)