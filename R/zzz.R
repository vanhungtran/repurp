.onLoad <- function(libname, pkgname) {
  if (.Platform$OS.type == "windows") {
    grDevices::windowsFonts(sans = grDevices::windowsFont("Arial"))
  }
}
