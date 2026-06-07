Sys.setenv(RSTUDIO_PANDOC = "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools")

rmarkdown::render(
  input      = "d:/OneDrive/\U0001F4CA_R_Statistics/Repurp/vignettes/repurp-intro.Rmd",
  output_dir = "d:/OneDrive/\U0001F4CA_R_Statistics/Repurp/vignettes"
)
