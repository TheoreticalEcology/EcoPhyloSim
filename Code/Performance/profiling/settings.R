# Named settings (createCompletePar() arguments, added to x = 50, y = 50, specRate = 2), used by profile-run.R.
MR <- c("mortality", "recruitment")
settings <- list(
  neutralGlobal    = list(dispersal = "global"),
  neutralLocal1    = list(dispersal = 1),
  neutralLocal3    = list(dispersal = 3),
  env              = list(dispersal = 1, environment = 1),
  dd1              = list(dispersal = 1, negativeDD = 1, negativeDDCut = 1),
  dd2              = list(dispersal = 1, negativeDD = 1, negativeDDCut = 2),
  dd3              = list(dispersal = 1, negativeDD = 1, negativeDDCut = 3),
  dd2gauss         = list(dispersal = 1, negativeDD = 1, negativeDDCut = 2, negativeDDKernel = "gaussian"),
  dd2species       = list(dispersal = 1, negativeDD = 1, negativeDDCut = 2, negativeDDKernel = "species"),
  envDD2           = list(dispersal = 1, environment = 1, negativeDD = 1, negativeDDCut = 2),
  posDD            = list(dispersal = 1, negativeDD = 1, negativeDDCut = 2, positiveDD = 0.1, positiveDDCut = 3),
  dd2recruit       = list(dispersal = 1, negativeDD = 1, negativeDDCut = 2, fitnessActsOn = MR),
  dd2recruitOnly   = list(dispersal = 1, negativeDD = 1, negativeDDCut = 2, fitnessActsOn = "recruitment"),
  dd2recruitGlobal = list(dispersal = "global", negativeDD = 1, negativeDDCut = 2, fitnessActsOn = MR),
  envRecruit       = list(dispersal = 1, environment = 1, fitnessActsOn = "recruitment"),
  envRecruit2      = list(dispersal = 2, environment = 1, fitnessActsOn = "recruitment"),
  envRecruitGlobal = list(dispersal = "global", environment = 1, fitnessActsOn = "recruitment"),
  envRecruitNarrow = list(dispersal = 1, environment = 1, nicheWidth = 0.01, fitnessActsOn = "recruitment"),
  speciesRecruit   = list(dispersal = 1, negativeDD = 1, negativeDDCut = 2, negativeDDKernel = "species", fitnessActsOn = "recruitment"),
  repro            = list(dispersal = 1, environment = 1, fitnessActsOn = "reproduction"),
  reproDD          = list(dispersal = 1, negativeDD = 1, environment = 1, fitnessActsOn = c("reproduction", "mortality")),
  ddGlobal         = list(dispersal = "global", negativeDD = 1, negativeDDCut = 2),
  envDDGlobal      = list(dispersal = "global", negativeDD = 1, negativeDDCut = 2, environment = 1),
  fission1         = list(dispersal = 1, fission = 1, x = 30, y = 30),
  fissionDD        = list(dispersal = 1, fission = 1, negativeDD = 1, x = 30, y = 30),
  protracted       = list(dispersal = 1, protracted = 5),
  redQueen         = list(dispersal = 1, redQueen = 0.5, redQueenStrength = 2)
)
# complete parameter list for a setting
settingPar <- function(name, ...) do.call(createCompletePar, modifyList(list(x = 50, y = 50, specRate = 2),
                                                                         modifyList(settings[[name]], list(...))))
