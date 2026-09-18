# Week 4 Q4 - Variant prioritization (filled thresholds)
# Data: ../data/variants_q4.tsv (synthetic teaching table).

# install.packages(c("readr", "dplyr"))  # run once if needed

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

variants_path <- file.path("data", "variants_q4.tsv")
variants <- read_tsv(variants_path, comment = "#", show_col_types = FALSE)

cat("=== Raw input (12 variants) ===\n")
print(variants, n = 12)

# ============================================================
# Filtering thresholds (filled)
# ============================================================
# Technical quality
require_pass <- TRUE    # keep only FILTER == "PASS" (drops LowQual, FAIL)
min_dp       <- 20       # depth floor (drops DP < 20)
min_gq       <- 20       # genotype quality floor (drops GQ < 20)

# Population frequency
max_af       <- 0.001    # rare-variant ceiling (drops common alleles > 0.1%)

# Functional consequence (impactful)
impact_consequences <- c(
  "stop_gained",
  "frameshift_variant",
  "splice_acceptor_variant",
  "splice_donor_variant",
  "missense_variant"
)

# Hard excludes (always drop)
filtered <- variants %>%
  filter(
    if (require_pass) FILTER == "PASS" else TRUE,
    if (!is.na(min_dp)) DP >= min_dp else TRUE,
    if (!is.na(min_gq)) GQ >= min_gq else TRUE,
    if (!is.na(max_af)) AF <= max_af else TRUE,
    CONSEQUENCE %in% impact_consequences
  )

cat("\n=== After hard excludes ===\n")
cat(sprintf("Kept %d of %d variants\n", nrow(filtered), nrow(variants)))
print(filtered)

# Soft ranking: ClinVar Pathogenic > Likely_pathogenic > Conflicting > Uncertain > Likely_benign > Benign
# + tie-break: lower AF first, higher DP first
filtered <- filtered %>%
  mutate(
    clinvar_rank = case_when(
      CLINVAR_SIG == "Pathogenic"                                ~ 1L,
      grepl("Likely_pathogenic", CLINVAR_SIG, ignore.case=TRUE)  ~ 2L,
      grepl("Conflicting",        CLINVAR_SIG, ignore.case=TRUE)  ~ 3L,
      grepl("Uncertain",          CLINVAR_SIG, ignore.case=TRUE)  ~ 4L,
      grepl("Likely_benign",      CLINVAR_SIG, ignore.case=TRUE)  ~ 5L,
      grepl("Benign",             CLINVAR_SIG, ignore.case=TRUE)  ~ 6L,
      TRUE                                                         ~ 7L
    )
  ) %>%
  arrange(clinvar_rank, AF, desc(DP)) %>%
  select(-clinvar_rank)

cat("\n=== Top shortlist (ranked) ===\n")
print(filtered, n = nrow(filtered))

# Save shortlist
dir.create("results", showWarnings = FALSE)
write_tsv(filtered, "results/q4_shortlist.tsv")
cat("\nWrote results/q4_shortlist.tsv\n")

# ============================================================
# Sanity report: which hard excludes dropped which variants
# ============================================================
cat("\n=== Hard-exclude report (why each dropped variant was dropped) ===\n")
dropped <- variants %>% anti_join(filtered, by = c("CHROM","POS","REF","ALT"))
for (i in seq_len(nrow(dropped))) {
  v <- dropped[i, ]
  reasons <- c()
  if (v$FILTER != "PASS") reasons <- c(reasons, sprintf("FILTER=%s", v$FILTER))
  if (!is.na(min_dp) && v$DP < min_dp) reasons <- c(reasons, sprintf("DP=%d<%d", v$DP, min_dp))
  if (!is.na(min_gq) && v$GQ < min_gq) reasons <- c(reasons, sprintf("GQ=%d<%d", v$GQ, min_gq))
  if (!is.na(max_af) && v$AF > max_af) reasons <- c(reasons, sprintf("AF=%.5f>%.5f", v$AF, max_af))
  if (!(v$CONSEQUENCE %in% impact_consequences)) reasons <- c(reasons, sprintf("CONSEQUENCE=%s (not impactful)", v$CONSEQUENCE))
  cat(sprintf("  DROP %s:%d %s>%s  [%s]  reasons: %s\n",
              v$CHROM, v$POS, v$REF, v$ALT,
              ifelse(nchar(v$GENE) > 0, v$GENE, "intergenic"),
              paste(reasons, collapse = "; ")))
}

# ============================================================
# sessionInfo for the log
# ============================================================
sink("logs/q4_sessionInfo.txt")
print(sessionInfo())
sink()
cat("\nWrote logs/q4_sessionInfo.txt\n")