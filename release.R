# R script for releasing new versions of FTOL data
# requires [github cli](https://cli.github.com/) to be installed
# automatically updates CFF file

# new_ver/notes are now derived automatically (docs/updating.md's semver
# rule: 3rd digit if the GenBank release is unchanged from the last release,
# 2nd digit if it's new) instead of hand-typed. First-digit (breaking
# format) bumps are never inferred -- if one is warranted, set new_ver by
# hand below instead of relying on this block.

library(gert)
library(targets)

ftol_repo <- "../ftol" # sibling checkout of the main pipeline repo

gb_release <- targets::tar_read(
  gb_release, store = file.path(ftol_repo, "_targets")
)
date_cutoff <- as.character(
  targets::tar_read(date_cutoff, store = file.path(ftol_repo, "_targets"))
)

last_tag <- system("git tag --sort=-v:refname", intern = TRUE)[1]
last_notes <- system(
  glue::glue('gh release view {last_tag} --json body -q .body'),
  intern = TRUE
)
last_gb_release <- stringr::str_match(last_notes, "release (\\d+)")[, 2]

ver_parts <- as.integer(strsplit(sub("^v", "", last_tag), "\\.")[[1]])
if (identical(as.character(gb_release), last_gb_release)) {
  ver_parts[3] <- ver_parts[3] + 1L # tree-only change, same GenBank release
} else {
  ver_parts[2] <- ver_parts[2] + 1L # new GenBank release
  ver_parts[3] <- 0L
}
new_ver <- paste0("v", paste(ver_parts, collapse = "."))

notes <- glue::glue(
  "Built with DNA sequences in [GenBank](https://ftp.ncbi.nlm.nih.gov/genbank/) ",
  "release {gb_release} (cutoff date {date_cutoff})"
)

message(glue::glue(
  "Last release: {last_tag} (GenBank {last_gb_release}). ",
  "Proposed release: {new_ver}.\n{notes}"
))

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

# Push the CFF-bump commit (mechanical; safe to automate)
git_push()

# Creating the GitHub release is public and hard to reverse -- confirm
# new_ver/notes above, then run this line yourself (or have Claude run it
# after you say go):
message(glue::glue(
  'gh release create {new_ver} --title "{new_ver}" --notes "{notes}" --latest'
))
