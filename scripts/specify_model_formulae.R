library(tidyverse)

all_states_1 <-
    list(
        "(1)",
        "(4 * Ka * concentration)",
        "(6 * Ka^2 * concentration^2)",
        "(4 * Ka^3 * concentration^3)",
        "(Ka^4 * concentration^4)",
        "(L)",
        "(L * D * 4 * Ka * concentration)",
        "(L * D^2 * 6 * Ka^2 * concentration^2)",
        "(L * D^3 * 4 * Ka^3 * concentration^3)",
        "(L * D^4 * Ka^4 * concentration^4)"
    )

open_states_1 <-
    list(
        "(L)",
        "(L * D * 4 * Ka * concentration)",
        "(L * D^2 * 6 * Ka^2 * concentration^2)",
        "(L * D^3 * 4 * Ka^3 * concentration^3)",
        "(L * D^4 * Ka^4 * concentration^4)"
    )

bound_states_1 <-
    list(
        "(1 * 4 * Ka * concentration)",
        "(2 * 6 * Ka^2 * concentration^2)",
        "(3 * 4 * Ka^3 * concentration^3)",
        "(4 * Ka^4 * concentration^4)",
        "(1 * L * D * 4 * Ka * concentration)",
        "(2 * L * D^2 * 6 * Ka^2 * concentration^2)",
        "(3 * L * D^3 * 4 * Ka^3 * concentration^3)",
        "(4 * L * D^4 * Ka^4 * concentration^4)"
    )

binding_eq_1 <-
    paste(
        "(",
        paste(bound_states_1, collapse = " + "),
        ") / (4 * (",
        paste(all_states_1, collapse = " + "),
        "))",
        sep = ""
    )

gating_eq_1 <-
    paste(
        "((",
        paste(open_states_1, collapse = " + "),
        ") / (",
        paste(all_states_1, collapse = " + "),
        ") / (L / (L + 1)))",
        sep = ""
    )

mwc_formula_full <-
    paste(
        "response ~ (binding_mask * (",
        binding_eq_1,
        ")) + ((1 - binding_mask) * (",
        gating_eq_1,
        "))",
        sep = ""
    )

all_states_2 <-
    list(
        "(1)",
        "(4 * Ka * concentration)",
        "(6 * Ka^2 * concentration^2)",
        "(4 * Ka^3 * concentration^3)",
        "(Ka^4 * concentration^4)",
        "(L)",
        "(L * D * 4 * Ka * concentration)",
        "(L * D * 6 * Ka^2 * concentration^2)",
        "(L * D * 4 * Ka^3 * concentration^3)",
        "(L * D * Ka^4 * concentration^4)"
    )

open_states_2 <-
    list(
        "(L)",
        "(L * D * 4 * Ka * concentration)",
        "(L * D * 6 * Ka^2 * concentration^2)",
        "(L * D * 4 * Ka^3 * concentration^3)",
        "(L * D * Ka^4 * concentration^4)"
    )

bound_states_2 <-
    list(
        "(1 * 4 * Ka * concentration)",
        "(2 * 6 * Ka^2 * concentration^2)",
        "(3 * 4 * Ka^3 * concentration^3)",
        "(4 * Ka^4 * concentration^4)",
        "(1 * L * D * 4 * Ka * concentration)",
        "(2 * L * D * 6 * Ka^2 * concentration^2)",
        "(3 * L * D * 4 * Ka^3 * concentration^3)",
        "(4 * L * D * Ka^4 * concentration^4)"
    )

binding_eq_2 <-
    paste(
        "(",
        paste(bound_states_2, collapse = " + "),
        ") / (4 * (",
        paste(all_states_2, collapse = " + "),
        "))",
        sep = ""
    )

gating_eq_2 <-
    paste(
        "((",
        paste(open_states_2, collapse = " + "),
        ") / (",
        paste(all_states_2, collapse = " + "),
        ") / (L / (L + 1)))",
        sep = ""
    )

mwc_formula_single_site <-
    paste(
        "response ~ (binding_mask * (",
        binding_eq_2,
        ")) + ((1 - binding_mask) * (",
        gating_eq_2,
        "))",
        sep = ""
    )

all_states_3 <-
    list(
        "(1)",
        "(4 * Ka * concentration)",
        "(6 * Ka^2 * concentration^2 * c)",
        "(4 * Ka^3 * concentration^3 * c^2)",
        "(Ka^4 * concentration^4 * c^3)",
        "(L)",
        "(L * D * 4 * Ka * concentration)",
        "(L * D^2 * 6 * Ka^2 * concentration^2 * c)",
        "(L * D^3 * 4 * Ka^3 * concentration^3 * c^2)",
        "(L * D^4 * Ka^4 * concentration^4 * c^3)"
    )

open_states_3 <-
    list(
        "(L)",
        "(L * D * 4 * Ka * concentration)",
        "(L * D^2 * 6 * Ka^2 * concentration^2 * c)",
        "(L * D^3 * 4 * Ka^3 * concentration^3 * c^2)",
        "(L * D^4 * Ka^4 * concentration^4 * c^3)"
    )

bound_states_3 <-
    list(
        "(1 * 4 * Ka * concentration)",
        "(2 * 6 * Ka^2 * concentration^2 * c)",
        "(3 * 4 * Ka^3 * concentration^3 * c^2)",
        "(4 * Ka^4 * concentration^4 * c^3)",
        "(1 * L * D * 4 * Ka * concentration)",
        "(2 * L * D^2 * 6 * Ka^2 * concentration^2 * c)",
        "(3 * L * D^3 * 4 * Ka^3 * concentration^3 * c^2)",
        "(4 * L * D^4 * Ka^4 * concentration^4 * c^3)"
    )

binding_eq_3 <-
    paste(
        "(",
        paste(bound_states_3, collapse = " + "),
        ") / (4 * (",
        paste(all_states_3, collapse = " + "),
        "))",
        sep = ""
    )

gating_eq_3 <-
    paste(
        "((",
        paste(open_states_3, collapse = " + "),
        ") / (",
        paste(all_states_3, collapse = " + "),
        ") / (L / (L + 1)))",
        sep = ""
    )

mwc_formula_cooperative <-
    paste(
        "response ~ (binding_mask * (",
        binding_eq_3,
        ")) + ((1 - binding_mask) * (",
        gating_eq_3,
        "))",
        sep = ""
    )
