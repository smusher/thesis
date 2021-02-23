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

spectra <-
	read_csv("/home/sam/current_analysis/paper/tnp-atp_doseresponse_unroofed/W311-GFP+SUR/W311-GFP+SUR_2019-06-17_5_spectra.csv") %>%
	filter(dye < 480)
images <-
	read_csv("/home/sam/current_analysis/paper/tnp-atp_doseresponse_unroofed/W311-GFP+SUR/W311-GFP+SUR_2019-06-17_5_images.csv") %>%
	mutate(
		intensity = case_when(bleaching_corrected_intensity < 0 ~ 0, TRUE ~ bleaching_corrected_intensity),
		wavelength = trunc(wavelength)
		) %>%
	group_by(cumulative_exposure, y, wavelength) %>%
	summarise(intensity = mean(intensity))
intensities <-
	read_csv("/home/sam/current_analysis/paper/tnp-atp_doseresponse_unroofed/W311-GFP+SUR/W311-GFP+SUR_2019-06-17_5_intensities.csv") %>%
	filter(dye < 480)

ggplot() +
geom_line(data=spectra, aes(x = wavelength, y = bleaching_corrected_spectra, colour = factor(concentration), group = cumulative_exposure)) +
theme_thesis() +
coord_cartesian(xlim = c(420, 650), ylim = c(-20, 150)) +
scale_colour_brewer(palette = "PuOr", direction = -1) -> a

ggplot() +
geom_raster(data=images, aes(x = wavelength, y = y, fill = intensity)) +
theme_thesis() +
scale_fill_viridis_c(option = "inferno") +
facet_grid(rows = vars(cumulative_exposure)) +
coord_cartesian(xlim = c(420, 650)) -> b

ggplot() +
geom_point(data = intensities, aes(x = concentration, y = 1 - response, fill = factor(concentration)), shape = 21, size = 3) +
scale_fill_brewer(palette = "PuOr", direction = -1) +
scale_x_log10() +
theme_thesis() -> c

(b + a) / c