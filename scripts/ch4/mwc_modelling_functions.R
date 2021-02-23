s4u <- "1"

s3u1a <- "4*Ka*Fa*1"
s2u2a <- "3/2*Ka*Fa*4*Ka*Fa*1"
s1u3a <- "2/3*Ka*Fa*3/2*Ka*Fa*4*Ka*Fa*1"
s4a <- "1/4*Ka*Fa*2/3*Ka*Fa*3/2*Ka*Fa*4*Ka*Fa*1"

s3u1b <- "4*Kb*Fb*1"
s2u2b <- "3/2*Kb*Fb*4*Kb*Fb*1"
s1u3b <- "2/3*Kb*Fb*3/2*Kb*Fb*4*Kb*Fb*1"
s4b <- "1/4*Kb*Fb*2/3*Kb*Fb*3/2*Kb*Fb*4*Kb*Fb*1"

s2u1a1b <- "3*Kb*Fb*4*Ka*Fa*1"
s1u1a2b <- "1*Kb*Fb*3*Kb*Fb*4*Ka*Fa*1"
s1a3b <- "1/3*Kb*Fb*1*Kb*Fb*3*Kb*Fb*4*Ka*Fa*1"

s1u2a1b <- "1*Ka*Fa*3*Kb*Fb*4*Ka*Fa*1"
s2a2b <- "1/2*Kb*Fb*1*Ka*Fa*3*Kb*Fb*4*Ka*Fa*1"

s3a1b <- "1/3*Ka*Fa*1*Ka*Fa*3*Kb*Fb*4*Ka*Fa*1"

#closed states, monomers bound to multiple ligands
s3u1ab <- "C*Kb*Fb*4*Ka*Fa*1"

s2u1a1ab <- "3*Ka*Fa*C*Kb*Fb*4*Ka*Fa*1"
s1u2a1ab <- "Ka*Fa*3*Ka*Fa*C*Kb*Fb*4*Ka*Fa*1"
s3a1ab <- "1/3*Ka*Fa*Ka*Fa*3*Ka*Fa*C*Kb*Fb*4*Ka*Fa*1"

s2u1b1ab <- "3*Kb*Fb*C*Kb*Fb*4*Ka*Fa*1"
s1u2b1ab <- "Kb*Fb*3*Kb*Fb*C*Kb*Fb*4*Ka*Fa*1"
s3b1ab <- "1/3*Kb*Fb*Kb*Fb*3*Kb*Fb*C*Kb*Fb*4*Ka*Fa*1"

s1u1a1b1ab <- "2*Ka*Fa*3*Kb*Fb*C*Kb*Fb*4*Ka*Fa*1"
s2a1b1ab <- "1/2*Ka*Fa*2*Ka*Fa*3*Kb*Fb*C*Kb*Fb*4*Ka*Fa*1"
s1a2b1ab <- "1/2*Kb*Fb*2*Ka*Fa*3*Kb*Fb*C*Kb*Fb*4*Ka*Fa*1"

s2u2ab <- "1/2*C*Kb*Fb*3*Ka*Fa*C*Kb*Fb*4*Ka*Fa*1"
s1u1a2ab <- "2*Ka*Fa*1/2*C*Kb*Fb*3*Ka*Fa*C*Kb*Fb*4*Ka*Fa*1"
s2a2ab <- "1/2*Ka*Fa*2*Ka*Fa*1/2*C*Kb*Fb*3*Ka*Fa*C*Kb*Fb*4*Ka*Fa*1"

s1a1b2ab <- "Kb*Fb*2*Ka*Fa*1/2*C*Kb*Fb*3*Ka*Fa*C*Kb*Fb*4*Ka*Fa*1"
s1u1b2ab <- "2*Kb*Fb*1/2*C*Kb*Fb*3*Ka*Fa*C*Kb*Fb*4*Ka*Fa*1"
s2b2ab <- "1/2*Kb*Fb*2*Kb*Fb*1/2*C*Kb*Fb*3*Ka*Fa*C*Kb*Fb*4*Ka*Fa*1"

s1u3ab <- "1/3*C*Kb*Fb*2*Ka*Fa*1/2*C*Kb*Fb*3*Ka*Fa*C*Kb*Fb*4*Ka*Fa*1"
s1a3ab <- "Ka*Fa*1/3*C*Kb*Fb*2*Ka*Fa*1/2*C*Kb*Fb*3*Ka*Fa*C*Kb*Fb*4*Ka*Fa*1"
s1b3ab <- "Kb*Fb*1/3*C*Kb*Fb*2*Ka*Fa*1/2*C*Kb*Fb*3*Ka*Fa*C*Kb*Fb*4*Ka*Fa*1"
s4ab <- "1/4*C*Kb*Fb*Ka*Fa*1/3*C*Kb*Fb*2*Ka*Fa*1/2*C*Kb*Fb*3*Ka*Fa*C*Kb*Fb*4*Ka*Fa*1"

#open states, each monomer bound to one ligand only
s4u_o <- "L*1"

s3u1a_o <- "Da*L*4*Ka*Fa*1"
s2u2a_o <- "Da^2*L*3/2*Ka*Fa*4*Ka*Fa*1"
s1u3a_o <- "Da^3*L*2/3*Ka*Fa*3/2*Ka*Fa*4*Ka*Fa*1"
s4a_o <- "Da^4*L*1/4*Ka*Fa*2/3*Ka*Fa*3/2*Ka*Fa*4*Ka*Fa*1"

s3u1b_o <- "Db*L*4*Kb*Fb*1"
s2u2b_o <- "Db^2*L*3/2*Kb*Fb*4*Kb*Fb*1"
s1u3b_o <- "Db^3*L*2/3*Kb*Fb*3/2*Kb*Fb*4*Kb*Fb*1"
s4b_o <- "Db^4*L*1/4*Kb*Fb*2/3*Kb*Fb*3/2*Kb*Fb*4*Kb*Fb*1"

s2u1a1b_o <- "Da*Db*L*3*Kb*Fb*4*Ka*Fa*1"
s1u1a2b_o <- "Da*Db^2*L*1*Kb*Fb*3*Kb*Fb*4*Ka*Fa*1"
s1a3b_o <- "Da*Db^3*L*1/3*Kb*Fb*1*Kb*Fb*3*Kb*Fb*4*Ka*Fa*1"

s1u2a1b_o <- "Da^2*Db*L*1*Ka*Fa*3*Kb*Fb*4*Ka*Fa*1"
s2a2b_o <- "Da^2*Db^2*L*1/2*Kb*Fb*1*Ka*Fa*3*Kb*Fb*4*Ka*Fa*1"

s3a1b_o <- "Da^3*Db*L*1/3*Ka*Fa*1*Ka*Fa*3*Kb*Fb*4*Ka*Fa*1"

#open states, monomers bound to multiple ligands
s3u1ab_o <- "Da*Db*L*C*Kb*Fb*4*Ka*Fa*1"

s2u1a1ab_o <- "Da^2*Db*L*3*Ka*Fa*C*Kb*Fb*4*Ka*Fa*1"
s1u2a1ab_o <- "Da^3*Db*L*Ka*Fa*3*Ka*Fa*C*Kb*Fb*4*Ka*Fa*1"
s3a1ab_o <- "Da^4*Db*L*1/3*Ka*Fa*Ka*Fa*3*Ka*Fa*C*Kb*Fb*4*Ka*Fa*1"

s2u1b1ab_o <- "Da*Db^2*L*3*Kb*Fb*C*Kb*Fb*4*Ka*Fa*1"
s1u2b1ab_o <- "Da*Db^3*L*Kb*Fb*3*Kb*Fb*C*Kb*Fb*4*Ka*Fa*1"
s3b1ab_o <- "Da*Db^3*L*1/3*Kb*Fb*Kb*Fb*3*Kb*Fb*C*Kb*Fb*4*Ka*Fa*1"

s1u1a1b1ab_o <- "Da^2*Db^2*L*2*Ka*Fa*3*Kb*Fb*C*Kb*Fb*4*Ka*Fa*1"
s2a1b1ab_o <- "Da^3*Db^2*L*1/2*Ka*Fa*2*Ka*Fa*3*Kb*Fb*C*Kb*Fb*4*Ka*Fa*1"
s1a2b1ab_o <- "Da^2*Db^3*L*1/2*Kb*Fb*2*Ka*Fa*3*Kb*Fb*C*Kb*Fb*4*Ka*Fa*1"

s2u2ab_o <- "Da^2*Db^2*L*1/2*C*Kb*Fb*3*Ka*Fa*C*Kb*Fb*4*Ka*Fa*1"
s1u1a2ab_o <- "Da^3*Db^2*L*2*Ka*Fa*1/2*C*Kb*Fb*3*Ka*Fa*C*Kb*Fb*4*Ka*Fa*1"
s2a2ab_o <- "Da^4*Db^2*L*1/2*Ka*Fa*2*Ka*Fa*1/2*C*Kb*Fb*3*Ka*Fa*C*Kb*Fb*4*Ka*Fa*1"

s1a1b2ab_o <- "Da^3*Db^3*L*Kb*Fb*2*Ka*Fa*1/2*C*Kb*Fb*3*Ka*Fa*C*Kb*Fb*4*Ka*Fa*1"
s1u1b2ab_o <- "Da^2*Db^3*L*2*Kb*Fb*1/2*C*Kb*Fb*3*Ka*Fa*C*Kb*Fb*4*Ka*Fa*1"
s2b2ab_o <- "Da^2*Db^4*L*1/2*Kb*Fb*2*Kb*Fb*1/2*C*Kb*Fb*3*Ka*Fa*C*Kb*Fb*4*Ka*Fa*1"

s1u3ab_o <- "Da^3*Db^3*L*1/3*C*Kb*Fb*2*Ka*Fa*1/2*C*Kb*Fb*3*Ka*Fa*C*Kb*Fb*4*Ka*Fa*1"
s1a3ab_o <- "Da^4*Db^3*L*Ka*Fa*1/3*C*Kb*Fb*2*Ka*Fa*1/2*C*Kb*Fb*3*Ka*Fa*C*Kb*Fb*4*Ka*Fa*1"
s1b3ab_o <- "Da^3*Db^4*L*Kb*Fb*1/3*C*Kb*Fb*2*Ka*Fa*1/2*C*Kb*Fb*3*Ka*Fa*C*Kb*Fb*4*Ka*Fa*1"
s4ab_o <- "Da^4*Db^4*L*1/4*C*Kb*Fb*Ka*Fa*1/3*C*Kb*Fb*2*Ka*Fa*1/2*C*Kb*Fb*3*Ka*Fa*C*Kb*Fb*4*Ka*Fa*1"

closed_states <-
	c(
	s4u, s3u1a, s2u2a, s1u3a, s4a, s3u1b, s2u2b, s1u3b, s4b,
	s2u1a1b, s1u2a1b, s1u1a2b, s3a1b, s1a3b, s2a2b,
	s3u1ab, s2u1a1ab, s1u2a1ab, s3a1ab, s2u1b1ab, s1u2b1ab,
	s3b1ab, s1u1a1b1ab, s2a1b1ab, s1a2b1ab, s2u2ab, s1u1a2ab,
	s2a2ab, s1a1b2ab, s1u1b2ab, s2b2ab, s1u3ab, s1a3ab, s1b3ab, s4ab
		)

open_states <-
	c(
	s4u_o, s3u1a_o, s2u2a_o, s1u3a_o, s4a_o, s3u1b_o, s2u2b_o, s1u3b_o, s4b_o,
	s2u1a1b_o, s1u2a1b_o, s1u1a2b_o, s3a1b_o, s1a3b_o, s2a2b_o,
	s3u1ab_o, s2u1a1ab_o, s1u2a1ab_o, s3a1ab_o, s2u1b1ab_o, s1u2b1ab_o,
	s3b1ab_o, s1u1a1b1ab_o, s2a1b1ab_o, s1a2b1ab_o, s2u2ab_o, s1u1a2ab_o,
	s2a2ab_o, s1a1b2ab_o, s1u1b2ab_o, s2b2ab_o, s1u3ab_o, s1a3ab_o, s1b3ab_o, s4ab_o
		)

all_states <-
	c(closed_states, open_states)

a_bound_states <-
	c(
	s3u1a, paste0("2*", s2u2a), paste0("3*", s1u3a), paste0("4*", s4a), s2u1a1b, paste0("2*", s1u2a1b), s1u1a2b, paste0("3*", s3a1b),
	s1a3b, paste0("2*", s2a2b), s3u1ab, paste0("2*", s2u1a1ab), paste0("3*", s1u2a1ab), paste0("4*", s3a1ab), s2u1b1ab, s1u2b1ab,
	s3b1ab, paste0("2*", s1u1a1b1ab), paste0("3*", s2a1b1ab), paste0("2*", s1a2b1ab), paste0("2*", s2u2ab), paste0("3*", s1u1a2ab),
	paste0("4*", s2a2ab), paste0("3*", s1a1b2ab), paste0("2*", s1u1b2ab), paste0("2*", s2b2ab), paste0("3*", s1u3ab), paste0("4*", s1a3ab), paste0("3*", s1b3ab), paste0("4*", s4ab),
	s3u1a_o, paste0("2*", s2u2a_o), paste0("3*", s1u3a_o), paste0("4*", s4a_o), s2u1a1b_o, paste0("2*", s1u2a1b_o), s1u1a2b_o, paste0("3*", s3a1b_o),
	s1a3b_o, paste0("2*", s2a2b_o), s3u1ab_o, paste0("2*", s2u1a1ab_o), paste0("3*", s1u2a1ab_o), paste0("4*", s3a1ab_o), s2u1b1ab_o, s1u2b1ab_o,
	s3b1ab_o, paste0("2*", s1u1a1b1ab_o), paste0("3*", s2a1b1ab_o), paste0("2*", s1a2b1ab_o), paste0("2*", s2u2ab_o), paste0("3*", s1u1a2ab_o),
	paste0("4*", s2a2ab_o), paste0("3*", s1a1b2ab_o), paste0("2*", s1u1b2ab_o), paste0("2*", s2b2ab_o), paste0("3*", s1u3ab_o), paste0("4*", s1a3ab_o), paste0("3*", s1b3ab_o), paste0("4*", s4ab_o)
		)

binding_eq <-
    paste(
        "((",
        paste(a_bound_states, collapse = ")+("),
        "))/(4*((",
        paste(all_states, collapse = ")+("),
        ")))",
        sep = ""
    )

gating_eq <-
    paste(
        "((",
        paste(open_states, collapse = ")+("),
        "))/((",
        paste(all_states, collapse = ")+("),
        "))",
        sep = ""
    )

norm_gating_eq <-
	paste("(", gating_eq, ")/(", str_replace_all(gating_eq, "Fa", "0"), ")", sep = "")

mwc_formula_string <-
    paste(
        "response ~ (binding_mask*(",
        binding_eq,
        "))+((1-binding_mask)*(",
        gating_eq,
        "))",
        sep = ""
    )
