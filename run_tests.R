#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(testthat)
  # Don't load spatialGE here - individual tests use devtools::load_all()
  # to load the development version
})

# Run all tests and capture results
results <- test_dir(
  "tests/testthat",
  reporter = SummaryReporter$new()
)

# Print detailed results
cat("\n\n=== TEST SUMMARY ===\n\n")
cat(sprintf("Total tests: %s\n", results[[1]]$total))
cat(sprintf("Passed: %s\n", results[[1]]$passed))
cat(sprintf("Failed: %s\n", results[[1]]$failed))
cat(sprintf("Skipped: %s\n", results[[1]]$skipped))

if (results[[1]]$failed > 0) {
  cat("\n=== FAILED TESTS ===\n")
  for (f in results[[1]]$failed_tests) {
    cat(sprintf("- %s: %s\n", f$file, f$test))
  }
}