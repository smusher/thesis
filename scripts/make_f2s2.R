library(tidyverse)
library(broom)
library(ggbeeswarm)

currents <-
	read_csv("data/tea_block.csv") %>%
    mutate(log_concentration = log10(tea_concentration))

currents_summary <-
    currents %>%
    group_by(construct, tea_concentration, log_concentration) %>%
    summarise(
        se = sd(response)/sqrt(length(response)),
        response = mean(response)
    )

x <- seq(-4, 1, length.out = 51)
frame <- tibble(log_concentration = x)

hill_equation <- function(log_concentration, ec50, hill, floor){
    floor + ((1 - floor) / (1 + 10^((ec50 - log_concentration) * -hill)))
}

fits <-
    currents %>%
    group_by(construct) %>%
    nest() %>%
    mutate(fit = map(data, ~ nls(
            response ~ hill_equation(log_concentration, ec50, hill, 0),
            data = data.frame(.),
            algorithm = "port",
            start = c(ec50 = -1, hill = 1),
            lower = c(ec50 = -6, hill = 0)
              )
    )
    )

plot_fits <-
	fits %>%
	mutate(fit = map(fit, augment, newdata = frame)) %>%
    unnest(fit) %>%
    mutate(tea_concentration = 10^log_concentration)

tidy_fits <-
	fits %>%
	mutate(fit = map(fit, tidy)) %>%
    unnest(fit)

fancy_scientific <- function(l) {
    l <- format(l, scientific = TRUE)
    l <- gsub("^(.*)e", "10^", l)
    parse(text=l)
}

ggplot() +
geom_quasirandom(data = currents, aes(x = tea_concentration/3, y = response, fill = construct), shape = 21, alpha = 0.5, size = 1.5, width = 0.2) +
geom_line(data = plot_fits, aes(x = tea_concentration/3, y = .fitted, colour = construct), size = 1) +
geom_errorbar(data = currents_summary, aes(x = tea_concentration/3, ymin = response - se, ymax = response + se), width = 0.2) +
geom_point(data = currents_summary, aes(x = tea_concentration/3, y = response, fill = construct), shape = 21, size = 2.5) +
hrbrthemes::theme_ipsum_ps() +
scale_colour_brewer(palette = "Set2", aesthetics = c("colour", "fill")) +
scale_x_log10(labels = fancy_scientific) +
scale_y_continuous(breaks = c(0.0, 0.5, 1.0)) +
labs(
    x = "[Triethylamine] (M)",
    y = expression(I/I[max])
    )

test <- tibble(log_concentration = 1e-3)