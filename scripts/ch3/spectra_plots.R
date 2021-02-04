library(tidyverse)
library(patchwork)
source("/home/sam/thesis/scripts/plot_theme.R")

spectra <-
	read_csv("/home/sam/current_analysis/paper/pcf/W311-GFP+SUR_ATP/W311-GFP+SUR_ATP_2017-11-30_spectra.csv") %>%
	filter(dye < 480)
images <-
	read_csv("/home/sam/current_analysis/paper/pcf/W311-GFP+SUR_ATP/W311-GFP+SUR_ATP_2017-11-30_images.csv") %>%
	mutate(
		intensity = case_when(bleaching_corrected_intensity < 0 ~ 0, TRUE ~ bleaching_corrected_intensity),
		wavelength = trunc(wavelength)
		) %>%
	group_by(cumulative_exposure, y, wavelength) %>%
	summarise(intensity = mean(intensity))

ggplot() +
geom_line(data=spectra, aes(x = wavelength, y = bleaching_corrected_spectra, colour = factor(concentration), group = cumulative_exposure)) +
theme_thesis() +
coord_cartesian(xlim = c(420, 650), ylim = c(-20, 200)) +
scale_colour_brewer(palette = "Set2") -> a

ggplot() +
geom_raster(data=images, aes(x = wavelength, y = y, fill = intensity)) +
theme_thesis() +
scale_fill_viridis_c(option = "inferno") +
facet_grid(rows = vars(cumulative_exposure)) +
coord_cartesian(xlim = c(420, 650)) -> b

ggsave("/home/sam/thesis/figures/ch3/atp_quenching.svg", a/b)