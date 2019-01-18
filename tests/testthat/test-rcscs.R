context("test-rcscs")

setup({
  features <- read.table("data/small_GNPS_buckettable.tsv", sep = "\t", comment.char = "", header = T, row.names = 1)
  cssraw <- read.table("data/small_GNPS_edges.tsv", sep = "\t", header = T)
  css <- prepare_css("data/small_GNPS_edges.tsv")
})


test_that("prepare_css sum", {
  expect_equal(sum(prepare_css("data/small_GNPS_edges.tsv")), 23.82165, tolerance = 10^5)
})

test_that("GNPS directory parsing", {
  expect_equal(length(read_GNPS_dir("tests/testthat/data/dir_example/")), 3)
})

test_that("Compute weighted CSCS distance", {
  library(foreach)
  features <- read.table("data/small_GNPS_buckettable.tsv", sep = "\t", comment.char = "", header = T, row.names = 1)
  css <- prepare_css("data/small_GNPS_edges.tsv")
  expect_equal(sum(as.matrix(cscs(features, css))), 25.14, tolerance = 10^-2)
})


