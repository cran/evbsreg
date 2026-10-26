################################################################
## script_05_simulation_scenario3.R
##
## Reproduces the corresponding simulation tables of the paper (appendix):
## Monte Carlo simulation - Scenario 3 (heterogeneous dispersion: 10% of the
## observations are generated with alpha reduced by 0.5, i.e. 0.7 instead of 1.2;
## with alpha = 0.5, as before, the contaminated cases would have alpha = 0).
################################################################

library(evbsreg)

m <- 5000               # change to 500 for a quick check
gama_vals <- c(-0.20, 0, 0.20)   # as in the paper's tables
n_vals    <- c(60, 120, 180)

cat(sprintf("Scenario 3 (alpha contamination): m=%d, n=%s, gama=%s\n",
            m, paste(n_vals, collapse=","), paste(gama_vals, collapse=",")))

results <- list()
for (g in gama_vals) {
  for (n in n_vals) {
    key <- sprintf("gama=%g_n=%d", g, n)
    cat(sprintf("\n--- %s ---\n", key))
    t0 <- Sys.time()
    res <- evbsreg.fit.mc(m = m, n = n,
                          beta0 = 1, beta1 = 2, alpha = 1.2, gama = g,   # 10% of cases get alpha - 0.5 = 0.7
                          scenario = "robustness", semente = 2023)
    print(res)
    cat(sprintf("Elapsed: %.1f s\n",
                as.numeric(difftime(Sys.time(), t0, units = "secs"))))
    results[[key]] <- res
  }
}

save(results, file = file.path(tempdir(), "results_scenario3.RData"))
cat("\nResults saved to results_scenario3.RData\n")
