context("TDAstats should load and unload cleanly")
library("TDAstats")

test_that("Loading and unload happen cleanly", {
  # get dlls after loading
  loadNamespace("TDAstats")
  dll.val <- getLoadedDLLs()
  
  # get dlls after unloading
  unloadNamespace("TDAstats")
  dll.val.1 <- getLoadedDLLs()
  
  # check if everything is fine
  expect_true("TDAstats" %in% names(dll.val))
  expect_false("TDAstats" %in% names(dll.val.1))
  expect_equal(length(dll.val), length(dll.val.1) + 1)
})