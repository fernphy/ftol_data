# R script for releasing new versions of FTOL data
# requires [github cli](https://cli.github.com/) to be installed
# automatically updates CFF file

library(gert)

# Specify version and release notes
new_ver <- "v1.5.1"
notes <- paste(
  "Built with DNA sequences in [GenBank](https://ftp.ncbi.nlm.nih.gov/genbank/) release 258 (cutoff date 2023-10-15)"
)

# Format CFF
cff <- glue::glue('
cff-version: 1.1.0
authors:
- name: "FTOL working group"
title: "Fern Tree of Life (FTOL) data"
type: data
version: {new_ver}
date-released: {Sys.Date()}')

# Write and commit CFF
readr::write_lines(cff, "CITATION.cff")
git_add("CITATION.cff")
git_commit("Update CITATION.cff")

if (nrow(git_status()) > 0) stop("Must have clean git repo before releasing")

# Push release
git_push()
system(glue::glue(
  'gh release create {new_ver} --title "{new_ver}" --notes "{notes}" --latest'
))
