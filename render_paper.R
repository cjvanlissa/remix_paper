library(quarto)
library(rmarkdown)

rmarkdown::render("manuscript.Rmd")
tmp <- readLines("manuscript.Rmd")
tmp[grepl("floatsintext", tmp)] <- "floatsintext      : yes"
tmp <- tmp[-do.call(`:`, as.list(grep("HIGHLIGHTS", tmp, fixed = TRUE)))]
tmp <- tmp[-(grep("**Original Publication**", tmp, fixed = TRUE):grep("**Funding statement**", tmp, fixed = TRUE))]
writeLines(tmp, "preprint.Rmd")
rmarkdown::render("preprint.Rmd")

files <- list.files(
  "results",
  pattern = "^Supp.+\\.csv$",
  full.names = TRUE
)

for (i in seq_along(files)) {
  file = files[i]
  name <- tools::file_path_sans_ext(basename(file))

  message("Rendering: ", name)
  tmp <- readLines("table_template.qmd")
  tmp <- gsub("{title}", gsub(".csv", "", basename(file), fixed = T), tmp, fixed = TRUE)
  writeLines(tmp, "temp.qmd")
  quarto::quarto_render(
    input = "temp.qmd",
    output_file = paste0(name, ".pdf"),
    execute_params = list(
      csv = normalizePath(file),
      caption = gsub("_", " ", name)
    ),
    quiet = F
  )
  #file.copy(paste0(name, ".pdf"), to = file.path(dirname(input_dir) , paste0(name, ".pdf")))
  #file.remove(paste0(name, ".pdf"))
}
