library(quarto)

input_dir <- "paper/results"
output_dir <- input_dir


files <- list.files(
  input_dir,
  pattern = "^Supp.+\\.csv$",
  full.names = TRUE
)

for (i in seq_along(files)) {
  file = files[i]
  name <- tools::file_path_sans_ext(basename(file))

  message("Rendering: ", name)
  tmp <- readLines("paper/table_template.qmd")
  tmp <- gsub("{title}", gsub(".csv", "", basename(file), fixed = T), tmp, fixed = TRUE)
  writeLines(tmp, "paper/temp.qmd")
  quarto::quarto_render(
    input = "paper/temp.qmd",
    output_file = paste0(name, ".pdf"),
    execute_params = list(
      csv = normalizePath(file),
      caption = gsub("_", " ", name)
    ),
    quiet = F
  )
  file.copy(paste0(name, ".pdf"), to = file.path(dirname(input_dir) , paste0(name, ".pdf")))
  file.remove(paste0(name, ".pdf"))
}
