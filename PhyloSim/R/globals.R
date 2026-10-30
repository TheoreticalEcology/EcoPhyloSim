# These variables are referenced in the Leipzig branch of runSimulation()
# (type = "Leipzig"), which is known non-functional legacy code awaiting
# Florian's decision on whether to fix, keep, or remove. globalVariables()
# here only silences the R CMD check NOTE — it does not make this code work.
utils::globalVariables(c("callLeipzig", "x", "y", "nSpec", 
                         "specRate", "runs", "densityStrength", "outVec"))