inhibition_bound <-
"((Ka * Fa * ((1 + Ka * Fa) ^ 3)) + (L * D * Ka * Fa * ((1 + D * Ka * Fa) ^ 3)) / (((1 + Ka * Fa) ^ 4) + (L * ((1 + D * Ka * Fa) ^ 4)))"

inhibition_open <-
"((L * ((1 + D * Ka * Fa) ^ 4)) / (((1 + Ka * Fa) ^ 4) + (L * ((1 + D * Ka * Fa) ^ 4)))) / (L / (L + 1))"

activation_bound <-
"((Kb * Fb * ((1 + Kb * Fb) ^ 3)) + (L * E * Kb * Fb * ((1 + E * Kb * Fb) ^ 3)) / (((1 + Kb * Fb) ^ 4) + (L * ((1 + E * Kb * Fb) ^ 4)))"

activation_open <-
"(((L * ((1 + E * Kb * Fb) ^ 4)) / (((1 + Kb * Fb) ^ 4) + (L * ((1 + E * Kb * Fb) ^ 4)))) - (L / (L + 1))) / (((L * E^4) / (1 + L * E^4)) - (L / (L + 1)))"
