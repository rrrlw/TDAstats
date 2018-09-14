# welcome message on package attaching
.onAttach <- function(libname, pkgname) {
}

# unload C++ DLL for proper cleanup
.onUnload <- function (libpath) {
  library.dynam.unload("TDAstats", libpath)
}
