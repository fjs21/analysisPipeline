test_that("conversion from Human to Mouse GeneID works", {
  # a vector of human GeneIDs
  HsEG = c('5156','4155')
  MmEG = c('18595','17196')
  names(MmEG) = HsEG
  expect_equal(convertHumanToMouseGeneID(HsEG), MmEG)
})
