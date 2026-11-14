## Neutral model with global dispersal
simu <- NeutralMod(x = 20, y = 20, runs = 100, specRate = 2, seed = 1)
image(simu$Output[[1]]$specMat)
specRich(simu)

## Local dispersal with exponential distances (mean 2 cells), saving
## generations 50 and 100, and conspecific negative density dependence
simu <- NeutralMod(x = 20, y = 20, runs = c(50, 100), specRate = 2, seed = 1,
                   dispersal = function(n) stats::rexp(n, rate = 1/2),
                   negativeDD = 0.5, negativeDDCut = 2)
rac(simu, which.result = "all")
