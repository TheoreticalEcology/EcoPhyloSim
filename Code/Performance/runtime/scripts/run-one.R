# Runs one simulation and appends its CPU time to a results file (semicolon-separated).
# usage: Rscript run-one.R <tag> <library folder> <parameters> <seed> <results file> <cpu label>
# <parameters>: createCompletePar() arguments as R code, e.g. "x=100,y=100,dispersal=1,runs=500";
# the abbreviations G, SP, GA, MR stand for "global", "species", "gaussian", c("mortality", "recruitment").
# Writes a header line if the file does not exist yet. Columns: tag; parameters; seed; cpu label; CPU seconds of runSimulation(); elapsed seconds; deaths; time.
# Libraries of versions before the renaming (#164) are loaded as PhyloSim; old argument names
# (density, densityCut) must then be used in <parameters>.
a <- commandArgs(TRUE)
lib <- a[2]
pkg <- if (dir.exists(file.path(lib, "PhyloSim"))) "PhyloSim" else "EcoPhyloSim"
suppressMessages(library(pkg, lib.loc = lib, character.only = TRUE))
G <- "global"; SP <- "species"; GA <- "gaussian"; MR <- c("mortality", "recruitment")
p <- eval(parse(text = paste0("list(", a[3], ")")))
p <- modifyList(list(x = 50, y = 50, specRate = 2), p)
p$seed <- as.integer(a[4])
par <- do.call(createCompletePar, p)
Sys.sleep(1)  # time for the caller to set affinity and priority
t <- system.time(invisible(capture.output(sim <- runSimulation(par))))
if (!file.exists(a[5])) cat("tag;parameters;seed;cpu;cpu_s;elapsed_s;deaths;time\n", file = a[5])
cat(sprintf("%s;%s;%s;%s;%.3f;%.3f;%.0f;%s\n", a[1], a[3], a[4], a[6], t[["user.self"]], t[["elapsed"]],
            p$x * p$y * max(p$runs), format(Sys.time(), "%Y-%m-%d %H:%M:%S")), file = a[5], append = TRUE)
