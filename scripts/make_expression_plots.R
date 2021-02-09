library(tidyverse)
library(cowplot)
library(ggbeeswarm)
library(knitr)
library(kableExtra)
extrafont::loadfonts()

surface_expression <-
    bind_rows(
        read_csv("data/anap_construct_surface_expression.csv"),
        read_csv("data/tmd0_surface_expression_one.csv"),
        read_csv("data/tmd0_surface_expression_two.csv")
        )

ggplot() +
geom_point(
    data = surface_expression,
    shape = 21,
    stroke = 0.5,
    size = 2,
    colour = "black",
    aes(x = construct, y = counts, fill = anap_present),
    position = "dodge"
) +
scale_y_log10(labels = fancy_scientific) +
coord_flip() +
scale_fill_brewer(palette = "Set2", guide = guide_legend(title = "ANAP", title.position = "top")) +
theme_cowplot(font_family = "IBM Plex Sans Condensed Light", font_size = 9) +
labs(y = "Arbitrary Light Units (ALU)") +
theme(axis.title.y = element_blank(), legend.position = "right", strip.background  = element_rect(fill = "white")) -> kir_comparison_plot

ggplot() +
geom_quasirandom(
    data = sur_construct_comparison,
    shape = 21,
    stroke = 0.5,
    size = 2,
    colour = "black",
    aes(x = kir_construct, y = counts, fill = sur_construct)
) +
scale_y_log10() +
coord_flip() +
scale_fill_brewer(palette = "Set2", na.value = "white", guide = guide_legend(title = "SUR Construct", title.position = "top")) +
labs(y = "Arbitrary Light Units (ALU)") +
theme_cowplot(font_family = "IBM Plex Sans Condensed Light", font_size = 9) +
theme(axis.title.y = element_blank(), legend.position = "right", strip.background  = element_rect(fill = "white")) -> sur_comparison_plot

ggplot() +
geom_quasirandom(
    data = tolbutamide_currents,
    shape = 21,
    stroke = 0.5,
    size = 2,
    colour = "black",
    aes(x = construct, y = tolbutamide_response, fill = construct)
) +
scale_y_log10() +
coord_flip() +
scale_fill_brewer(palette = "Set2", guide = guide_legend(title = "Construct", title.position = "top")) +
labs(y = "Current fraction remaining with 100 uM Tolbutamide") +
theme_cowplot(font_family = "IBM Plex Sans Condensed Light", font_size = 9) +
theme(axis.title.y = element_blank(), legend.position = "none", strip.background  = element_rect(fill = "white")) -> tolbutamide_plot
