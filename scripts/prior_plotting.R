library(tidyverse)
library(ggdist)
library(distributional)
library(patchwork)

#drc_priors
ec50 <- dist_normal(-4, 1)
hill <- dist_normal(1, 0.3)
emax <- dist_uniform(0, 1)
sigma <- dist_cauchy(0, 1)

(autoplot(ec50) + autoplot(hill)) / (autoplot(emax) + autoplot(sigma))

#pcf_priors
L <- dist_normal(log(1), log(5))
D <- dist_normal(log(0.1), log(5))
Ka <- dist_normal(log(1e4), log(5))
sigma <- dist_cauchy(0, 1)

(autoplot(L) + autoplot(D)) / (autoplot(Ka) + autoplot(sigma))
