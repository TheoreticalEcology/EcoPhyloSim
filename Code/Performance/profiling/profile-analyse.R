# usage: Rscript profile-analyse.R <dir with unzipped .sleepy> <EcoPhyloSim.dll>
# Self and inclusive time shares per function, only for samples inside callModel.
# EcoPhyloSim addresses are resolved with addr2line (DWARF, incl. inlined functions), because Very Sleepy
# assigns addresses to the nearest preceding symbol.
bin <- Sys.getenv("RTOOLS_BIN", "C:/rtools45/x86_64-w64-mingw32.static.posix/bin")  # nm and addr2line
args <- commandArgs(TRUE); d <- args[1]; dll <- args[2]
sym <- read.table(file.path(d, "Symbols.txt"), quote = "\"", col.names = c("addr", "mod", "fun", "file", "line"),
                  colClasses = "character", comment.char = "")
modName <- sub("[.]dll$", "", basename(dll)); eco <- sym$mod == modName
hex <- function(x) as.numeric(paste0("0x", sub("^0x", "", x)))  # exact below 2^53
toHex <- function(x) vapply(x, function(v) { s <- ""; for (k in 1:16) { s <- paste0(c(0:9, letters[1:6])[v %% 16 + 1], s); v <- v %/% 16 }; paste0("0x", s) }, "")
# load base of the DLL: the most-sampled function's lowest address minus its file address, rounded to 64 KiB
nmTab <- read.table(text = system2(file.path(bin, "nm.exe"), dll, stdout = TRUE), fill = TRUE,
                    col.names = c("a", "t", "s"), colClasses = "character")
top <- names(sort(table(sym$fun[eco]), decreasing = TRUE))[1]
a0 <- hex(nmTab$a[nmTab$s == paste0("_", top)][1])
r0 <- min(hex(sym$addr[eco & sym$fun == top]))
delta <- floor((r0 - a0) / 65536) * 65536
fileAddr <- toHex(hex(sym$addr[eco]) - delta)
a2l <- system2(file.path(bin, "addr2line.exe"), c("-f", "-i", "-C", "-a", "-e", shQuote(dll)), input = fileAddr, stdout = TRUE)
# parse: "0xADDR" then pairs (function, file:line), innermost first
chains <- list(); cur <- NULL
i <- 1
while (i <= length(a2l)) {
  if (grepl("^0x[0-9a-f]+$", a2l[i])) { cur <- a2l[i]; chains[[cur]] <- character(0); i <- i + 1; next }
  f <- sub("[(].*", "", a2l[i]); loc <- sub(".*/", "", a2l[i + 1])
  chains[[cur]] <- c(chains[[cur]], if (length(chains[[cur]]) == 0) paste0(f, " @", loc) else f)
  i <- i + 2
}
names(chains) <- sym$addr[eco][match(names(chains), fileAddr)]
frames <- setNames(as.list(paste0("[", sym$mod, "]")), sym$addr)
frames[sym$addr[eco]] <- chains[sym$addr[eco]]
cs <- strsplit(readLines(file.path(d, "Callstacks.txt")), " ")
tm <- vapply(cs, function(x) as.numeric(x[1]), 0)
fr <- lapply(cs, function(x) unlist(frames[x[-1]], use.names = FALSE))
keep <- vapply(fr, function(f) any(f == "callModel"), TRUE)
tm <- tm[keep]; fr <- fr[keep]; tot <- sum(tm)
leafLine <- vapply(fr, `[`, "", 1)
leafFun <- sub(" @.*", "", leafLine)
cleaned <- lapply(fr, function(f) unique(sub(" @.*", "", f)))
incl <- sort(tapply(rep(tm, lengths(cleaned)), unlist(cleaned), sum), decreasing = TRUE)
self <- tapply(tm, leafFun, sum)
res <- data.frame(fun = names(incl), inclusive = round(100 * incl / tot, 1),
                  self = round(100 * ifelse(is.na(self[names(incl)]), 0, self[names(incl)]) / tot, 1))
cat(sprintf("samples in callModel: %.2f s (%d stacks)\n", tot, length(tm)))
res$fun <- substr(gsub("std::|, std::allocator<[^>]*> ?|Rcpp::", "", res$fun), 1, 70)
options(width = 200)
cat("
Model functions (inclusive / self, % of callModel time):
")
print(res[grepl("^(Landscape|LocalEnvironment|GlobalEnvironment|Individual|Species|Phylogeny|RandomGen|PhylSimModel)::", res$fun) & res$inclusive >= 0.5, ], row.names = FALSE)
cat("
Self time by function (>= 1%):
")
print(head(res[order(-res$self), ][res$self[order(-res$self)] >= 1, c("fun", "self")], 25), row.names = FALSE)
cat("\nTop source lines (self time):\n")
ll <- sort(tapply(tm, leafLine, sum), decreasing = TRUE)
print(head(data.frame(line = substr(names(ll), 1, 90), self = round(100 * ll / tot, 1)), 15), row.names = FALSE)
