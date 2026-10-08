# Simulation code for validating use of one skyline interval.
# Based on https://github.com/jugne/sRanges-material/blob/main/validation/estimate_all_params_ext_dna/simulate.R
# Edits made by Kate Truman using ChatGPT and Claude

rm(list = ls())
options(digits = 16)

library(ape)
library(phytools)
library(ggplot2)
library(stringr)
library(beastio)
library(FossilSim)
library(R.utils)


## helper function
"%notin%" <- Negate("%in%")
setwd("/home/ket581/skyline/SUMMER/skyline-SRFBD-validation/posteriors_for_tree_summaries")

## change paths if needed
templates_dir <- "templates/"
beast_dir <- "beast/bin/beast"
run_count <- 0
while (dir.exists(paste("simulated_trees/", run_count, sep = ""))) {
  run_count <- run_count + 1
}
out_dir <- paste("simulated_trees/", run_count, "/", sep = "")
dir.create(paste("simulated_trees/", run_count, sep = ""))
file.create(paste(out_dir, "out", sep = ""))


set.seed(5647829)


### Parameter transformation functions
f_lambda <- function(diversification, turnover) {
  return(diversification / (1 - turnover))
}

f_mu <- function(turnover, lambda) {
  return(turnover * lambda)
}

f_psi <- function(mu, sampling_prop) {
  return(mu * sampling_prop / (1 - sampling_prop))
}

## Set initial params to null, origin will be fixed (but still estimated)
div_rate <- NULL
turnover <- NULL
sampling_prop <- NULL
sampl_extant_prob <- NULL
lambda <- NULL
mu <- NULL
psi <- NULL
redraws <- 0
origin <- 4

# Boundaries of skyline intervals. This should always include zero, and not the origin.
times <- c(0)
ints <- length(times)

## function to redraw parameters from their posteriors (uniform in all cases here)
redraw_params <- function() {
  div_rate <<- runif(ints, 0.7, 0.9)
  turnover <<- runif(ints, 0.2, 0.8)
  sampling_prop <<-  runif(ints, 0.2, 0.8)
  sampl_extant_prob <<- runif(1, 0.7, 1.) # 1

  lambda <<- f_lambda(div_rate, turnover)
  mu <<- f_mu(turnover, lambda)
  psi <<- f_psi(mu, sampling_prop)
}

# number of trees to simulate
ntrees <- 200


## The following min and max values are set to ensure a small enough tree.
## If they are violated, the simulation will be rejected and new parameter values redrawn.
## This is a bit problematic, since we do not account in any way for rejected simulations and parameter combinations when doing inference.
## In inference we only condition on survival (so not on min or max number of samples).

# min for extant samples and max for total nodes
min_ext_samples <- 5
max_nodes <- 1000

# min and max value for fossils
min_fossils <- 0
max_fossils <- 10000


# This structure is getting rather large as the number of intervals increase - it would be worth finding a dynamic way to define this if there are significantly more intervals.
true_rates <- data.frame(
  div_rate_1 = numeric(),
  turnover_1 = numeric(),
  sampling_prop_1 = numeric(),
  rho = numeric(), birth_1 = numeric(),
  death_1 = numeric(),
  sampling_1 = numeric(),
  times = numeric(), origin = numeric(), mrca = numeric(), tree = character(),
  n_samples = numeric(), n_extant = numeric(), n_ranges = numeric(), draws = numeric()
)

# simulating trees and fossils with parameters
trees <- list()
beast_trees <- list()
fossils <- list()
taxonomy <- list()
samp_trees <- list()
while (length(trees) < ntrees) {
  redraw_params()
  redraws <- redraws + 1
  # Simulate tree using skyline parameters
  tree_tmp <- TreeSim::sim.rateshift.age(age = origin, numbsim = 1, lambda = lambda, mu = mu, times = times, mrca = FALSE, complete = TRUE)
  
  if (length(tree_tmp[[1]]) == 1 || tree_tmp[[1]]$Nnode > max_nodes) {
    next # reject and redraw parameters
  }
  mrca <- max(ape::node.depth.edgelength(tree_tmp[[1]]))
  if (length(which((mrca - ape::node.depth.edgelength(tree_tmp[[1]])) < 1e-7)) < min_ext_samples) {
    next # reject and redraw parameters
  }
  t <- tree_tmp[1][[1]]
    print(ape::write.tree(t))
  newick <- ape::write.tree(t)
  origin <- tree.max(as.phylo(t))
  # Call to FossilSim needs to include the origin in the interval input
  horizons <- c(times, origin)
  taxonomy_tmp <- sim.taxonomy(tree_tmp[[1]], beta = 0, lambda.a = 0)
  # Simulate fossils using skyline rates.
  fossils_tmp <- FossilSim::sim.fossils.intervals(rates = psi, taxonomy = taxonomy_tmp, interval.ages = horizons)
  plot(fossils_tmp, tree=tree_tmp[[1]],taxonomy = taxonomy_tmp, show.taxonomy = TRUE, show.tip.label=TRUE)
  beast_tree_tmp = beast.fbd.format(tree_tmp[[1]], fossils_tmp, rho=sampl_extant_prob, digits=16)

  tree_after_rho <- ape::read.tree(text = beast_tree_tmp)
  mrca <- max(ape::node.depth.edgelength(tree_after_rho))
  n_ext <- length(which((mrca - ape::node.depth.edgelength(tree_after_rho)) < 1e-7))
  n_fossils <- length(tree_after_rho$tip.label) - n_ext
  if (n_ext < min_ext_samples || n_fossils < min_fossils || n_fossils > max_fossils) {
    next # reject and redraw parameters
  }
  write(beast_tree_tmp, file = paste(out_dir, "out", sep = ""), append = TRUE)
  trees <- c(trees, tree_tmp)
  beast_trees[[length(trees)]] <- beast_tree_tmp
  fossils[[length(trees)]] <- fossils_tmp
  taxonomy <- c(taxonomy, taxonomy_tmp)
  true_rates_tmp <- data.frame(div_rate[1], turnover[1],
    sampling_prop[1],
    rho = sampl_extant_prob, lambda[1], mu[1], psi[1], times = paste(times, collapse = " "), origin, mrca = mrca, tree = beast_tree_tmp,
    n_samples = n_fossils + n_ext, n_extant = n_ext, n_ranges = 0, draws = redraws
  )
  true_rates <- rbind(true_rates, true_rates_tmp)
  redraws <- 0
  i <- length(trees)
  print(i)
}

# Now create simulation xmls for each tree. We don't currently include DNA or morphological characters.

for (i in 1:ntrees) {
  beast_tree = beast_trees[[i]]
  true_rates$tree[i] = beast_tree
  tmp_tree = ape::read.tree(text=beast_tree)
  sample_times <- ape::node.depth.edgelength(tmp_tree)
sample_times <- max(sample_times) - sample_times
sample_times_round <- round(sample_times, 15)
sample_times_round[sample_times_round < 1e-10] <- 0
mrca <- max(ape::node.depth.edgelength(tmp_tree))


# 1. rename: youngest -> _first, oldest -> _last
prefix  <- sub("_[^_]*$", "", tmp_tree$tip.label)    # drops the trailing counter
new_lab <- tmp_tree$tip.label
for (tip in unique(prefix)) {
  idx <- which(prefix == tip)
  idx <- idx[order(sample_times_round[idx], idx)]    # youngest first
  new_lab[idx[1]] <- paste0(tip, "_first")
  if (length(idx) > 1) new_lab[idx[length(idx)]] <- paste0(tip, "_last")
  if (length(idx) > 2) warning(tip, " has ", length(idx), " tips; intermediates not renamed")
}
tmp_tree$tip.label <- new_lab

# 2. taxon list, extant list, ranges from the renamed tree
taxon <- list(); taxon_extant <- list(); strat_ranges <- list(); strat_ranges_refs <- list()
j <- 0
for (tip in unique(prefix)) {
  first <- paste0(tip, "_first"); last <- paste0(tip, "_last")
  has_last <- last %in% tmp_tree$tip.label
  taxon <- append(taxon, first)
  if (has_last) taxon <- append(taxon, last)
  end_lab <- if (has_last) last else first
  strat_ranges <- append(strat_ranges, paste0(
    '<stratigraphicRange id="r', j, '" spec="StratigraphicRange" firstOccurrence="@',
    first, '" lastOccurrence="@', end_lab, '"/>'))
  strat_ranges_refs <- append(strat_ranges_refs, paste0('<stratigraphicRange idref="r', j, '"/>'))
  for (lab in c(first, if (has_last) last))
    if (sample_times_round[match(lab, tmp_tree$tip.label)] == 0) taxon_extant <- append(taxon_extant, lab)
  j <- j + 1
}

# 3. strings, with a hard stop on missing labels
taxon_extant_str <- paste0(
  vapply(taxon_extant, function(tx) paste0("<sequence spec='Sequence' taxon='", tx, "' value='?'/>"), ""),
  collapse = "\n\t\t\t")

taxon_str <- c(); taxon_set_str <- c(); taxa_age_str <- c()
for (tx in taxon) {
  ti <- match(tx, tmp_tree$tip.label)
  if (is.na(ti)) stop("taxon ", tx, " not found among tree tips")
  taxon_str     <- c(taxon_str, paste0("<sequence spec='Sequence' taxon='", tx, "' value='?'/>"))
  taxon_set_str <- c(taxon_set_str, paste0("<taxon spec='Taxon' id='", tx, "'/>"))
  taxa_age_str  <- c(taxa_age_str, paste0(tx, "=", sample_times_round[ti]))
}
taxon_str     <- paste0(taxon_str, collapse = "\n\t\t\t")
taxon_set_str <- paste0(taxon_set_str, collapse = "\n\t\t\t\t\t\t")
taxa_age_str  <- paste0(taxa_age_str, collapse = ", ")

  ####### now create inference xmls, with the previosuly simulated data
  sim <- readLines(paste0(templates_dir, "ssRanges_inference_template_one_int.xml"))

  sim <- gsub(
    pattern = "<inputTaxa/>",
    replace = taxon_set_str, x = sim
  )

  sim <- gsub(
    pattern = "<intervalTimes/>",
    replace = paste(times, collapse = " "), x = sim
  )

  sim <- gsub(
    pattern = "<inputTaxaAge/>",
    replace = taxa_age_str, x = sim
  )

  sim <- gsub(
    pattern = "<inputStratRanges/>",
    replace = paste0(strat_ranges, collapse = "\n\t\t\t\t"), x = sim
  )
  sim <- gsub(
    pattern = "<inputStratRangesRef/>",
    replace = paste0(strat_ranges_refs, collapse = "\n\t\t\t"), x = sim
  )

  rnd_origin <- runif(1, mrca * 2, 500) # just so origin is not smaller than mrca
  rnd_div_rate <- runif(ints, 0.7, 0.9)
  rnd_turnover <- runif(ints, 0.2, 0.8)
  rnd_sampling_prop <- runif(ints, 0.2, 0.8)
  rnd_sampl_extant_prob <- runif(1, 0.7, 1.)

  sim <- gsub(
    pattern = "<initOrigin/>",
    replace = paste0("<parameter id='origin' lower='0.0'
                                name='stateNode'>", rnd_origin, "</parameter>"),
    x = sim
  )
  sim <- gsub(
    pattern = "<initDiversificationRate/>",
    replace = paste0("<parameter id='netDiversification' lower='0.0'
                                name='stateNode'>", paste(rnd_div_rate, collapse = " "), "</parameter>"),
    x = sim
  )
  sim <- gsub(
    pattern = "<initTurnover/>",
    replace = paste0("<parameter id='turnOver' lower='0.' upper = '1.'
                                name='stateNode'>", paste(rnd_turnover, collapse = " "), "</parameter>"),
    x = sim
  )
  sim <- gsub(
    pattern = "<initSamplingProportion/>",
    replace = paste0("<parameter id='samplingProportion' lower='0.0'
                                name='stateNode'>", paste(rnd_sampling_prop, collapse = " "), "</parameter>"),
    x = sim
  )

  sim <- gsub(
    pattern = "<initSamplingAtPresentProb/>",
    replace = paste0("<parameter id='rho' lower='0.0'
                                name='stateNode'>", rnd_sampl_extant_prob, "</parameter>"),
    x = sim
  )
  cat("mrca", mrca, "\n")
  index <- i - 1
  inf_dir <- paste0("inf/", index)
  inf_full <- paste0("inf/", (i - 1))
  dir.create(inf_dir, recursive = T)

sim <- gsub(
  pattern = "<logname/>",
  replace = "<logger logEvery=\"5000\" fileName=\"sRanges.$(seed).log\">",
  x = sim
)

sim <- gsub(
  pattern = "<logtreename/>",
  replace = "<logger logEvery=\"5000\" fileName=\"sRanges.$(seed).trees\">",
  x = sim
)

  writeLines(sim, con = paste(inf_dir, "/sRanges_inference.xml", sep = ""))
}
save.image(file = "simulation.RData")

write.csv(true_rates, "true_rates.csv")
file.create(paste(out_dir, "true_rates.csv", sep = ""))
write.csv(true_rates, file = paste(out_dir, "true_rates.csv", sep = ""))
