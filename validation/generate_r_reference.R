# Generate reference outputs from R-QUEFTS (Rquefts) for validating PyQUEFTS.
#
# Usage:  Rscript generate_r_reference.R [n_scenarios] [output_dir]
#
# Writes:
#   r_reference_model.csv   randomized model runs (inputs + 12 outputs)
#   r_reference_supply.csv  randomized nutSupply1 / nutSupply2 runs
#   r_reference_examples.csv documented examples (reagro tutorial, man pages)
#   r_session_info.txt

suppressPackageStartupMessages(library(Rquefts))

args <- commandArgs(trailingOnly = TRUE)
n_scen <- if (length(args) >= 1) as.integer(args[1]) else 5000
outdir <- if (length(args) >= 2) args[2] else "."
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

set.seed(20251112)
crops <- c("", "Barley", "Cassava", "Chickpea", "Chillies", "Cotton", "Cowpea", "Field bean",
           "Groundnut", "Jute", "Kenaf", "Lentil", "Maize", "Millet (bulrush)", "Mungbean", "Onion",
           "Pigeonpea", "Potato", "Rapeseed", "Rice", "Sesame", "Sorghum", "Soyabean", "Sugarbeet",
           "Sugarcane", "Sunflower", "Sweetpotato", "Tobacco", "Wheat")

runif_mix <- function(n, lo, hi, p_zero = 0.1, p_extreme = 0.05) {
  x <- runif(n, lo, hi)
  z <- runif(n) < p_zero
  x[z] <- 0
  e <- runif(n) < p_extreme
  x[e] <- runif(sum(e), hi, hi * 5)
  x
}

# ---------------------------------------------------------------- model runs
inp <- data.frame(
  crop = sample(crops, n_scen, replace = TRUE),
  N_base_supply = runif_mix(n_scen, 0, 200),
  P_base_supply = runif_mix(n_scen, 0, 40),
  K_base_supply = runif_mix(n_scen, 0, 300),
  N_recovery = runif(n_scen, 0.2, 0.8),
  P_recovery = runif(n_scen, 0.05, 0.4),
  K_recovery = runif(n_scen, 0.2, 0.8),
  N = runif_mix(n_scen, 0, 250, p_zero = 0.3),
  P = runif_mix(n_scen, 0, 80, p_zero = 0.3),
  K = runif_mix(n_scen, 0, 200, p_zero = 0.3),
  store_att = runif(n_scen, 500, 15000),
  SeasonLength = sample(c(60, 90, 110, 120, 150, 180, 240, 400), n_scen, replace = TRUE),
  stringsAsFactors = FALSE
)
# organ branches: some runs with tiny storage organ (organs 1 / 2)
k <- sample(n_scen, round(0.08 * n_scen))
inp$store_att[k] <- runif(length(k), 0, 200)
inp$leaf_att <- inp$store_att * runif(n_scen, 0.1, 1.2)
inp$stem_att <- inp$store_att * runif(n_scen, 0.1, 1.2)
k2 <- sample(k, round(length(k) / 2))
inp$stem_att[k2] <- runif(length(k2), 0, 200)   # organs == 1 branch
inp$leaf_att[k2] <- runif(length(k2), 0, 300)

outnames <- names(run(quefts()))
out <- matrix(NA_real_, nrow = n_scen, ncol = length(outnames))
colnames(out) <- outnames
soil0 <- quefts_soil()
for (i in seq_len(n_scen)) {
  crop <- quefts_crop(inp$crop[i])
  soil <- soil0
  soil$N_base_supply <- inp$N_base_supply[i]
  soil$P_base_supply <- inp$P_base_supply[i]
  soil$K_base_supply <- inp$K_base_supply[i]
  soil$N_recovery <- inp$N_recovery[i]
  soil$P_recovery <- inp$P_recovery[i]
  soil$K_recovery <- inp$K_recovery[i]
  fert <- list(N = inp$N[i], P = inp$P[i], K = inp$K[i])
  biom <- list(leaf_att = inp$leaf_att[i], stem_att = inp$stem_att[i],
               store_att = inp$store_att[i], SeasonLength = inp$SeasonLength[i])
  q <- quefts(soil, crop, fert, biom)
  out[i, ] <- run(q)
}
write.csv(cbind(inp, out), file.path(outdir, "r_reference_model.csv"), row.names = FALSE)

# ---------------------------------------------------------------- supply runs
ns <- 2000
sup <- data.frame(
  temp = runif(ns, 10, 35),
  pH = runif(ns, 3.5, 9.5),
  SOC = runif(ns, 0, 80),
  Kex = runif(ns, 0, 40),
  Polsen = runif(ns, 0, 40),
  Ptotal = runif(ns, 0, 2000)
)
s1 <- nutSupply1(sup$pH, sup$SOC, sup$Kex, sup$Polsen)
colnames(s1) <- paste0("s1_", colnames(s1))
s1t <- nutSupply1(sup$pH, sup$SOC, sup$Kex, sup$Polsen, sup$Ptotal)
colnames(s1t) <- paste0("s1tot_", colnames(s1t))
s2 <- nutSupply2(sup$temp, sup$pH, sup$SOC, sup$Kex, sup$Polsen, sup$Ptotal)
colnames(s2) <- paste0("s2_", colnames(s2))
write.csv(cbind(sup, s1, s1t, s2), file.path(outdir, "r_reference_supply.csv"), row.names = FALSE)

# ---------------------------------------------------------------- documented examples
ex <- list()
ex$default <- run(quefts())
q <- quefts(quefts_soil(), quefts_crop("Barley"), list(N = 100, P = 0, K = 0),
            list(leaf_att = 2200, stem_att = 2700, store_att = 4800, SeasonLength = 110))
ex$barley_N100 <- run(q)
q$SeasonLength <- 162; q$leaf_att <- 2651; q$stem_att <- 5053; q$store_att <- 8208
q$N <- 100; q$P <- 50; q$K <- 50
ex$barley_long <- run(q)
rice <- quefts_crop("Rice"); q <- quefts(crop = rice)
q$leaf_att <- 2651; q$stem_att <- 5053; q$store_att <- 8208
fert(q) <- list(N = 0, P = 0, K = 0)
for (N in c(1, 51, 101, 191)) { q["N"] <- N; ex[[paste0("rice_N", N)]] <- run(q) }
exdf <- do.call(rbind, ex)
write.csv(data.frame(case = rownames(exdf), exdf), file.path(outdir, "r_reference_examples.csv"), row.names = FALSE)

writeLines(c(R.version.string, paste("Rquefts", as.character(packageVersion("Rquefts"))),
             paste("generated", format(Sys.time()))), file.path(outdir, "r_session_info.txt"))
cat("done:", n_scen, "model scenarios,", ns, "supply scenarios\n")
