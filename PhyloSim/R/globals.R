# Variables used in the type = "Leipzig" branch of runSimulation(), which is
# currently not functional (see #112). Declared here only to avoid the
# R CMD check NOTE about undefined global variables.
utils::globalVariables(c("callLeipzig", "x", "y", "nSpec", 
                         "specRate", "runs", "densityStrength", "outVec"))

# Loop variable of foreach() in runPhylosimAnalysis.R and runSimulationBatch.R.
utils::globalVariables("i")
