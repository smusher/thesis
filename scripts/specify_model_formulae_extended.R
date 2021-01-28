K1 <- 10^-4
DK1 <- K1*(1/0.04)
A1 <- F1/K1
C1 <- K1/DK1
L <- 1

K2 <- 10^-4
DK2 <- K1*(1/5)
A2 <- F2/K2
C2 <- K2/DK2

F1 <- 10^seq(-7, -2, length.out = 16)
F2 <- 10^seq(-7, -2, length.out = 16)

expand.grid(F1, F2) %>%
tibble() %>%
rename(ATP = Var1, PIP2 = Var2) -> frame

calc_p <- function(L, F1, K1, DK1, F2, K2, DK2){
	A1 <- F1/K1
	C1 <- K1/DK1
	A2 <- F2/K2
	C2 <- K2/DK2
	(L*(1 + C1*A1)*(1 + C2*A2)) / ((1 + A1) * (1 + A2))
}

calc_open <- function(p){
	p^4/ (1+p^4)
}

calc_bound <- function(F, K, DK, p){
	A <- F/K
	C <- K/DK
	((A/(1+A)) + ((p^4) * ((C*A)/(1 + C*A)))) / (1 + (p^4))
}

frame %>%
mutate(
	p = calc_p(L, ATP, K1, DK1, PIP2, K2, DK2),
	open_frac = calc_open(p),
	atp_bound_frac = calc_bound(ATP, K1, DK1, p),
	pip_bound_frac = calc_bound(PIP2, K2, DK2, p)
	) %>%
pivot_longer(open_frac:pip_bound_frac, names_to = "measure", values_to = "value") -> toplot

ggplot(toplot) +
geom_line(aes(x = ATP, y = value, linetype = measure)) +
scale_linetype_manual(values = c("open_frac" = 1, "atp_bound_frac" = 2, "pip_bound_frac" = 3)) +
facet_wrap(vars(PIP2), labeller = labeller(PIP2 = label_both)) +
scale_x_log10()
