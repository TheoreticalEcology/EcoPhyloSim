
simu <- runSimulation(createCompletePar(x=50, y=50, runs=1000, specRate=2))
\donttest{
summStat <- calculateSummaryStatistics(simu)
}
