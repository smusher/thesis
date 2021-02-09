library(tidyverse)
source("/home/sam/thesis/scripts/plot_theme.R")

#read in spectra for donor and acceptor
anap <- read_csv("scripts/ch3/anap.csv")
tnpatp <- read_csv("scripts/ch3/tnp-atp.csv")
J <-
	left_join(anap, tnpatp) %>%
	mutate(overlap = donor_emission_normalised * molar_extinction_coefficient * wavelength^4) %>%
	summarise(overlap_sum = sum(overlap)) %>%
	pull(overlap_sum)

qy <- 0.22 #Anap quantum yield (Zagotta 2016)
sc <- 0.00008785 #Scale constant
ri <- 1.33 #Refractive index of the medium (water in this case)
k <- 2/3 #Kappa-squared - 2/3 assumes freely rotating donor and acceptor.

R0 <- (sc * k * qy * J * (ri^-4))^(1/6)

frame <-
	tibble(distance = seq(0, 80, length.out = 51)) %>%
	mutate(efficiency = 1/(1+((distance / R0)^6)))

ggplot() +
	geom_line(data = frame, aes(x = distance, y = efficiency), size = 1) +
	theme_thesis() +
	geom_vline(xintercept = R0, linetype = 2) +
	xlab("Distance (Å)") +
	ylab("FRET efficiency") -> fret_plot

ggsave("/home/sam/thesis/figures/ch3/fret_efficiency.svg",
	plot = fret_plot,
	device = "svg",
	width = 10, height = 6, units = "cm")
