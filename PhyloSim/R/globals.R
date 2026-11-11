# Variables used in the type = "Leipzig" branch of runSimulation(), which is
# currently not functional (see #112). Declared here only to avoid the
# R CMD check NOTE about undefined global variables.
utils::globalVariables(c("callLeipzig", "x", "y", "nSpec", 
                         "specRate", "runs", "densityStrength", "outVec"))

# Loop variable created by foreach(), which R CMD check cannot see.
# Keep it when the Leipzig variables above are removed (#112).
utils::globalVariables("i")
