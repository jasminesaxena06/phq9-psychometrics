# =============================================================================
# Psychometric evaluation of the PHQ-9 and a simulated computer adaptive test
# Data: NHANES 2017-2018 (US adults, publicly available)
#
# Steps
#   1. Load and clean the data
#   2. Descriptives
#   3. Reliability (alpha, item-total correlations)
#   4. Confirmatory factor analysis (is the PHQ-9 one factor?)
#   5. Measurement invariance by sex
#   6. Item response theory (graded response model)
#   7. Simulated computer adaptive test (CAT)
#
# Run the whole file from the project folder. Results are written to output/.
# =============================================================================

# ---- 0. Setup ---------------------------------------------------------------

required <- c("haven", "psych", "lavaan", "mirt", "mirtCAT", "latticeExtra")
to_install <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(to_install) > 0) install.packages(to_install)

library(psych)
library(lavaan)
library(mirt)
library(mirtCAT)

set.seed(2026)
dir.create("data", showWarnings = FALSE)
dir.create("output", showWarnings = FALSE)

# Close anything left open by an earlier run that stopped part-way
while (sink.number() > 0) sink()
graphics.off()

# Everything printed from here on is also saved to output/results.txt
sink("output/results.txt", split = TRUE)

# ---- 1. Load and clean the data ---------------------------------------------

base_url <- "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2017/DataFiles/"

get_nhanes <- function(file) {
  path <- file.path("data", file)
  if (!file.exists(path)) {
    message("Downloading ", file, " ...")
    ok <- tryCatch({
      download.file(paste0(base_url, file), path, mode = "wb", quiet = TRUE)
      TRUE
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (!ok) {
      unlink(path)
      stop("Could not download ", file, ". Download it manually from\n  ",
           base_url, file, "\nand save it in the data folder, then run again.")
    }
  }
  haven::read_xpt(path)
}

dpq  <- get_nhanes("DPQ_J.xpt")   # PHQ-9 depression screener
demo <- get_nhanes("DEMO_J.xpt")  # demographics (sex, age)

nhanes_items <- sprintf("DPQ%03d", seq(10, 90, by = 10))
item_names <- c("anhedonia", "low_mood", "sleep", "fatigue", "appetite",
                "worthless", "concentration", "psychomotor", "suicidal")

d <- merge(as.data.frame(dpq)[, c("SEQN", nhanes_items)],
           as.data.frame(demo)[, c("SEQN", "RIAGENDR", "RIDAGEYR")],
           by = "SEQN")

# NHANES stores zeros as tiny numbers (e.g. 5.4e-79), so round.
# 7 = refused and 9 = don't know become missing.
d[nhanes_items] <- lapply(d[nhanes_items], function(x) {
  x <- round(as.numeric(x))
  x[x %in% c(7, 9)] <- NA
  x
})
names(d)[match(nhanes_items, names(d))] <- item_names

n_screened <- nrow(d)
d <- d[complete.cases(d[item_names]), ]
d$sex <- factor(d$RIAGENDR, levels = 1:2, labels = c("Male", "Female"))
X <- d[item_names]

# ---- 2. Descriptives --------------------------------------------------------

cat("\n=== 2. DESCRIPTIVES ===\n")
total <- rowSums(X)
cat(sprintf("Respondents with PHQ-9 data: %d; complete on all 9 items: %d\n",
            n_screened, nrow(X)))
cat(sprintf("Mean age %.1f (SD %.1f); %.1f%% female\n",
            mean(d$RIDAGEYR), sd(d$RIDAGEYR), 100 * mean(d$sex == "Female")))
cat(sprintf("PHQ-9 total: mean %.2f, median %d, SD %.2f\n",
            mean(total), as.integer(median(total)), sd(total)))
cat(sprintf("Scoring 10 or more (usual cut-off for moderate depression): %.1f%%\n",
            100 * mean(total >= 10)))

item_dist <- t(sapply(X, function(x) prop.table(table(factor(x, levels = 0:3)))))
colnames(item_dist) <- c("Not at all", "Several days", "Over half the days", "Nearly every day")
cat("\nProportion choosing each response, by item:\n")
print(round(item_dist, 3))
write.csv(round(item_dist, 4), "output/item_distributions.csv")

# ---- 3. Reliability ---------------------------------------------------------

cat("\n=== 3. RELIABILITY ===\n")
rel <- psych::alpha(X)
poly <- psych::polychoric(X)$rho
rel_ordinal <- psych::alpha(poly, n.obs = nrow(X))
cat(sprintf("Cronbach's alpha: %.3f\n", rel$total$raw_alpha))
cat(sprintf("Ordinal alpha (from polychoric correlations): %.3f\n",
            rel_ordinal$total$raw_alpha))
cat("\nCorrected item-total correlations:\n")
print(round(setNames(rel$item.stats$r.drop, item_names), 3))

# ---- 4. Confirmatory factor analysis ----------------------------------------

cat("\n=== 4. CONFIRMATORY FACTOR ANALYSIS (one factor, WLSMV) ===\n")
model <- paste("depression =~", paste(item_names, collapse = " + "))
fit <- cfa(model, data = X, ordered = item_names, estimator = "WLSMV")

fm <- fitMeasures(fit, c("chisq.scaled", "df", "pvalue.scaled", "cfi.scaled",
                         "tli.scaled", "rmsea.scaled", "rmsea.ci.lower.scaled",
                         "rmsea.ci.upper.scaled", "srmr"))
print(round(fm, 3))
cat("Rough guide: CFI/TLI > .95, RMSEA < .06, SRMR < .08 suggest good fit.\n",
    "(The chi-square is almost always significant in samples this large.)\n")

std <- standardizedSolution(fit)
loadings <- std[std$op == "=~", "est.std"]
names(loadings) <- item_names
cat("\nStandardised factor loadings:\n")
print(round(loadings, 3))

omega <- sum(loadings)^2 / (sum(loadings)^2 + sum(1 - loadings^2))
cat(sprintf("\nOrdinal omega (from CFA loadings): %.3f\n", omega))

# ---- 5. Measurement invariance by sex ---------------------------------------

cat("\n=== 5. MEASUREMENT INVARIANCE BY SEX ===\n")
cat("Does the PHQ-9 measure depression the same way in men and women?\n")
cat("Configural model: same one-factor structure, all parameters free in each group.\n")
cat("Full invariance model: loadings and thresholds equal across groups.\n")
cat("(For ordinal items, loadings and thresholds are tested together; see Wu & Estabrook, 2016.)\n\n")
mi <- tryCatch({
  fit_group <- function(equal = character(0)) {
    cfa(model, data = d, group = "sex", ordered = item_names,
        estimator = "WLSMV", parameterization = "theta", group.equal = equal)
  }
  configural <- fit_group()
  full       <- fit_group(c("loadings", "thresholds"))

  mi_fits <- sapply(list(configural = configural, full_invariance = full),
                    fitMeasures, fit.measures = c("cfi.scaled", "rmsea.scaled", "srmr"))
  print(round(t(mi_fits), 3))
  cat(sprintf("Change in CFI: %.3f\n",
              mi_fits["cfi.scaled", "full_invariance"] - mi_fits["cfi.scaled", "configural"]))
  cat("A drop in CFI of more than .01 is a common sign that invariance does not hold.\n",
      "(The chi-square test below is very sensitive in samples this large.)\n\n")
  print(lavTestLRT(configural, full))
  TRUE
}, error = function(e) {
  cat("Invariance models could not be fitted:", conditionMessage(e), "\n")
  FALSE
})

# ---- 6. Item response theory ------------------------------------------------

cat("\n=== 6. ITEM RESPONSE THEORY (graded response model) ===\n")
irt <- mirt(X, model = 1, itemtype = "graded", verbose = FALSE)

pars <- coef(irt, IRTpars = TRUE, simplify = TRUE)$items
cat("Item parameters. a = discrimination (higher = more informative);\n",
    "b1-b3 = severity at which a higher response becomes more likely than a lower one.\n")
print(round(pars, 3))
write.csv(round(pars, 4), "output/irt_item_parameters.csv")

cat("\nItem fit (S-X2; small p-values flag possible misfit):\n")
print(itemfit(irt))

cat("\nOverall fit (C2 statistic):\n")
print(M2(irt, type = "C2"))

# Where on the depression scale does the PHQ-9 measure precisely?
theta_grid <- matrix(seq(-3, 4, by = 0.01))
info <- testinfo(irt, theta_grid)
precise <- theta_grid[info >= 1 / 0.4^2]   # SE <= 0.4, roughly reliability >= .84
cat(sprintf("\nTest information peaks at theta = %.2f.\n", theta_grid[which.max(info)]))
if (length(precise) > 0) {
  cat(sprintf("Measurement is precise (SE <= 0.4) between theta = %.2f and %.2f.\n",
              min(precise), max(precise)))
  cat(sprintf("%.1f%% of respondents score below that range, where the PHQ-9 measures imprecisely.\n",
              100 * mean(fscores(irt, method = "MAP")[, 1] < min(precise))))
}

png("output/test_information.png", width = 1600, height = 1100, res = 200)
print(plot(irt, type = "infoSE", main = "PHQ-9 test information and standard error"))
invisible(dev.off())

png("output/item_information.png", width = 1600, height = 1100, res = 200)
print(plot(irt, type = "infotrace", facet_items = FALSE,
           main = "Information provided by each item",
           auto.key = list(space = "right", lines = TRUE, points = FALSE),
           par.settings = list(superpose.line = list(col = hcl.colors(9, "Dark 3"),
                                                     lty = rep(1:3, 3), lwd = 2))))
invisible(dev.off())

# ---- 7. Simulated computer adaptive test ------------------------------------

cat("\n=== 7. SIMULATED COMPUTER ADAPTIVE TEST ===\n")
cat("Real response patterns are fed through a CAT that picks the most informative\n",
    "next item and stops once the score is precise enough (standard error target).\n\n")

n_sim <- min(500, nrow(X))
sim_idx <- sample(nrow(X), n_sim)
patterns <- as.matrix(X[sim_idx, ])
theta_full <- fscores(irt, method = "MAP", response.pattern = patterns)[, "F1"]
severity <- cut(rowSums(patterns), breaks = c(-1, 4, 9, 14, 27),
                labels = c("Minimal (0-4)", "Mild (5-9)", "Moderate (10-14)", "Severe (15+)"))

run_cat <- function(se_target) {
  res <- mirtCAT(mo = irt, local_pattern = patterns, criteria = "MI",
                 start_item = "MI", method = "MAP",
                 design = list(min_SEM = se_target))
  data.frame(se_target = se_target,
             n_items   = sapply(res, function(r) length(r$items_answered)),
             theta_cat = sapply(res, function(r) r$thetas[1]),
             theta_full = theta_full,
             severity  = severity)
}

cat_results <- do.call(rbind, lapply(c(0.3, 0.4, 0.5), run_cat))

summary_tab <- do.call(rbind, lapply(split(cat_results, cat_results$se_target), function(s) {
  data.frame(se_target        = s$se_target[1],
             mean_items       = round(mean(s$n_items), 2),
             pct_fewer_than_9 = round(100 * mean(s$n_items < 9), 1),
             r_with_full_test = round(cor(s$theta_cat, s$theta_full), 3))
}))
cat("Items needed and agreement with the full 9-item score:\n")
print(summary_tab, row.names = FALSE)

cat("\nMean number of items used, by severity (SE target 0.4):\n")
s04 <- cat_results[cat_results$se_target == 0.4, ]
print(round(tapply(s04$n_items, s04$severity, mean), 2))

write.csv(cat_results, "output/cat_simulation.csv", row.names = FALSE)

png("output/cat_items_by_severity.png", width = 1600, height = 1100, res = 200)
plot(jitter(s04$theta_full), jitter(s04$n_items, amount = 0.15),
     xlab = "Depression severity (theta, full PHQ-9)",
     ylab = "Items needed by the CAT",
     main = "How many items the CAT needed (SE target 0.4)",
     pch = 16, col = adjustcolor("steelblue", 0.5), ylim = c(1, 9))
invisible(dev.off())

# ---- Done -------------------------------------------------------------------

cat("\n=== SESSION INFO ===\n")
print(sessionInfo())
sink()
message("Finished. Results are in the output folder.")
