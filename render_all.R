# render_all.R
# Run from the biomarker_reports/ directory:
#   Rscript render_all.R          # renders full Quarto website
#   Rscript render_all.R single   # renders each TA report standalone (embed-resources)

args <- commandArgs(trailingOnly = TRUE)
mode <- if (length(args) > 0) args[1] else "website"

ta_files <- c(
  "autoimmune_fcrn.qmd",
  "oncology_checkpoint.qmd",
  "cardiovascular_pcsk9.qmd",
  "neurology_alzheimers.qmd"
)

if (mode == "website") {
  message("==> Rendering full Quarto website project ...")
  system("quarto render")

} else {
  # Standalone mode: render each file with embed-resources = true
  for (f in ta_files) {
    message(sprintf("\n==> Rendering %s ...", f))
    cmd <- sprintf(
      'quarto render %s --to html -M embed-resources:true -M self-contained:true',
      f
    )
    ret <- system(cmd)
    if (ret != 0) warning(sprintf("Non-zero exit for %s", f))
  }
  message("\nDone. Standalone HTML files written alongside each .qmd")
}
