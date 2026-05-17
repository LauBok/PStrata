# Build internal package data (R/sysdata.rda) from inst/ config files
# Run this script from the package root directory:
#   source("data-raw/build_sysdata.R")

.family_info_lines <- readLines("inst/family_info.txt")
.link_info_lines <- readLines("inst/link_info.txt")
.function_implement_lines <- readLines("inst/function_implement.txt")

save(
  .family_info_lines,
  .link_info_lines,
  .function_implement_lines,
  file = "R/sysdata.rda"
)

cat("Saved R/sysdata.rda\n")
