library(tidyverse)

subunit <- c(1,2,3,4)
site <- c("A", "B")
state <- c("open", "closed")
bound <-c(TRUE, FALSE)

expand.grid("subunit"=subunit, "site"=site, "state"=state, "bound"=bound) %>%
unite(sites, subunit:site) %>%
pivot_wider(names_from = sites, values_from = bound) %>%
unnest() %>%
expand(state,`1_A`,`2_A`,`3_A`,`4_A`,`1_B`,`2_B`,`3_B`,`4_B`) %>%
mutate(state_id = row_number()) %>%
pivot_longer(`1_A`:`4_B`, names_to = c("subunit", "site"), values_to = "bound", names_sep = "_") %>%
pivot_wider(names_from = site, values_from = bound) %>%
rowwise() %>%
mutate(
	AB = case_when(A==TRUE & B==TRUE ~ TRUE, TRUE~FALSE)
	) %>%
group_by(state, state_id) %>%
mutate(
	AB_bound = sum(AB),
	A_bound = sum(A) - AB_bound,
	B_bound = sum(B) - AB_bound
	) %>%
ungroup() %>%
select(state, A_bound, B_bound, AB_bound) %>%
unique()