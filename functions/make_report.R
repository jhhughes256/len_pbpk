# Clean up the environment
rm(list = ls())

# Load the libraries
library(knitr)
library(rmarkdown)

# Set the root dir because my rmds live in rmds/ subfolder
opts_knit$set(root.dir = '../.')

# Main function
report <- function(file, n.file = "", report.dir = "reports") {
  # Set the name of the report file
  base.name <- sub(pattern = ".Rmd", replacement = "", x = basename(file))

  # Make nfiles with always 2 digits
  n.file <- ifelse(as.integer(n.file) < 10, paste0("0", n.file), n.file)
  file.name <- paste0(n.file, "-", base.name, ".html")

  # Render
  render(
    input = file,
    output_format = html_document(
      toc = TRUE,
      toc_depth = 1,
      code_folding = "hide"
    ),
    output_file = file.name,
    output_dir = report.dir,
    envir = new.env()
  )
}
