#' @description open/create a qmd or Rmd file and switch to the document directory
myopen <- function(file) {
  if (!tools::file_ext(file) %in% c("Rmd", "qmd")) {
    stop(cat("\nWrong file type. Only open/create .Rmd or .qmd files\n"))
  }
  
  current_path = dirname(".")
  doc_path = dirname(file)
  
  if (!file.exists(file)) {
    if (!dir.exists(doc_path)) dir.create(doc_path)
    file.create(file)
  }
  
  rstudioapi::navigateToFile(file)
  
  if (current_path != doc_path) {
    setwd(doc_path)
    cat(paste0("Set working dir to: \n", doc_path, "\n"))
  }
}

