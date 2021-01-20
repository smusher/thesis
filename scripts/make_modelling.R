library(tidyverse)
library(patchwork)

binding_func <- function(x, L, D, Ka) {
    (Ka * x * (1 + Ka * x)^3 + L * D * Ka * x * (1 + D * Ka * x)^3) /
    ((1 + Ka * x)^4 + L * (1 + D * Ka * x)^4)
}

gating_func <- function(x, L, D, Ka) {
    (L * (1 + D * Ka * x)^4) /
    ((1 + Ka * x)^4 + L * (1 + D * Ka * x)^4)
}

concentration <- c(10 ^ seq(-7, -2, by = 0.1))
L = c(0.02, 0.2, 2)
D = c(0.8, 0.4, 0.04)
Ka = c(1.5e4)

model_frame <-
    expand.grid(
        x = concentration,
        Ka = Ka,
        L = L,
        D = D
        ) %>%
    as_tibble() %>%
    mutate(
        binding = 1 - binding_func(x, L, D, Ka),
        gating = gating_func(x, L, D, Ka) / (L / (1+L))
        )

ggplot() +
scale_x_log10() +
scale_colour_brewer(palette = "RdYlBu") +
scale_linetype_manual(
    name = "modelling",
    values = c(
        "bound fraction" = 1,
        "F / Fmax" = 1,
        "open fraction" = 2,
        "I / Imax" = 2,
        "limits" = 3
        )
    ) +
guides(linetype = guide_legend(nrow = 3), colour = guide_legend(nrow = 2)) +
geom_line(data = model_frame, aes(x, binding, linetype = "F / Fmax", colour = factor(L))) +
geom_line(data = model_frame, aes(x, gating, linetype = "I / Imax", colour = factor(L))) +
facet_grid(. ~ D, labeller = label_both) +
hrbrthemes::theme_ipsum_ps() -> Dplot

ggplot() +
scale_x_log10() +
scale_colour_brewer(palette = "RdYlBu") +
scale_linetype_manual(
    name = "modelling",
    values = c(
        "bound fraction" = 1,
        "F / Fmax" = 1,
        "open fraction" = 2,
        "I / Imax" = 2,
        "limits" = 3
        )
    ) +
guides(linetype = guide_legend(nrow = 3), colour = guide_legend(nrow = 2)) +
geom_line(data = model_frame, aes(x, binding, linetype = "F / Fmax", colour = factor(D))) +
geom_line(data = model_frame, aes(x, gating, linetype = "I / Imax", colour = factor(D))) +
facet_grid(. ~ L, labeller = label_both) +
hrbrthemes::theme_ipsum_ps() -> Lplot

Dplot / Lplot