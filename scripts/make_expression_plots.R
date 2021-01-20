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

kir_construct_comparison <-
    surface_expression %>%
    filter(
        kir_construct %in% c("WT-GFP", "WT", "WT-HA-GFP", "WT-HA", "W311*-GFP", "W311*", "W311*-HA", "W311*-HA-GFP"),
        sur_construct %in% c("SUR", NA),
        chaperone %in% c(NA, "Tolbutamide")
        ) %>%
    unite(construct, kir_construct, sur_construct, sep = "+")

sur_construct_comparison <-
    surface_expression %>%
    filter(
        kir_construct %in% c("WT-GFP", "WT-HA-GFP", "W311*-GFP", "W311*-HA-GFP"),
        sur_construct %in% c("SUR", "TMD0_195", "TMD0_232", NA),
        chaperone %in% c(NA, "Tolbutamide"),
        !(kir_construct %in% c("W311*-GFP", "W311*-HA-GFP") & anap_present == FALSE)
        )

tolbutamide_currents <-
    read_csv("data/electrophys_data.csv") %>%
    select(unique_experiment_id, construct, tolbutamide_response) %>%
    group_by(construct, unique_experiment_id) %>%
    summarise(tolbutamide_response = mean(tolbutamide_response)) %>%
    drop_na()

fancy_scientific <- function(l) {
    l <- format(l, scientific = TRUE)
    l <- gsub("^(.*)e", "10^", l)
    parse(text=l)
}

split_construct <- function(l) {
    gsub("[+]", "\n +", l)
}

ggplot() +
geom_point(
    data = kir_construct_comparison,
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
