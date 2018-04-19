# welcome message on package attaching
.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to TDAstats")
}
