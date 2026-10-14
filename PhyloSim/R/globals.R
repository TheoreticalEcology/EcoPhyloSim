# Suppress R CMD check notes for variables used in .C() calls
utils::globalVariables(c("callLeipzig", "x", "y", "nSpec", 
                         "specRate", "runs", "densityStrength", "outVec"))