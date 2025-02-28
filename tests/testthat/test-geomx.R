
# geomx function testing
#
testthat::test_that("import_geomx_samplesheet", {
   # use test file
   samplefile <- system.file("data", "SampleSheet.csv", package="platjam")

   # test list of vectors
   exp_lens <- setNames(c(4, 4, 1, 1, 65, 3, 65),
      c("Header", "Reads", "Sequencing_Settings",
         "BCLConvert_Settings",
         "BCLConvert_Data",
         "DragenGeoMxNGS_Settings",
         "DragenGeoMxNGS_Data"))
   exp_class1 <- setNames(rep("character", 7),
      names(exp_lens))
   samplelist <- import_geomx_samplesheet(samplefile, return_type="list")
   testthat::expect_equal(
      lengths(samplelist),
      exp_lens)
   testthat::expect_equal(
      sapply(samplelist, class),
      exp_class1)

   exp_class2 <- setNames(rep("data.frame", 7),
      names(exp_lens))
   exp_ncols <- setNames(c(2, 2, 2, 2, 4, 2, 1),
      names(exp_lens))
   # test list of data.frame
   sampledfs <- import_geomx_samplesheet(samplefile, return_type="dflist")
   testthat::expect_equal(
      sapply(sampledfs, class),
      exp_class2)
   testthat::expect_equal(
      sapply(sampledfs, ncol),
      exp_ncols)

   # test indices
   sampleindices <- import_geomx_samplesheet(samplefile, return_type="indices")
   expected_ind5_1 <- c("GCACTCGA+GTTACTCT",
      "ACGAGATG+CAAGCCAC",
      "TTCCGTTA+ATCCAGAC",
      "TCGTCGAC+CGGATTGG",
      "TCTCTAAC+ATAGCGTC",
      "GAATTGGT+CGCGTTGA")
   testthat::expect_equal(
      length(sampleindices),
      64)
   testthat::expect_equal(
      sampleindices[c(1:5, 64)],
      expected_ind5_1)

   # test indices, no revcomp
   sampleindices_norc <- import_geomx_samplesheet(samplefile,
      return_type="indices",
      do_revcomp=FALSE)
   expected_ind5_1_norc <- c("GCACTCGA+AGAGTAAC",
      "ACGAGATG+GTGGCTTG",
      "TTCCGTTA+GTCTGGAT",
      "TCGTCGAC+CCAATCCG",
      "TCTCTAAC+GACGCTAT",
      "GAATTGGT+TCAACGCG")
   testthat::expect_equal(
      length(sampleindices_norc),
      64)
   testthat::expect_equal(
      sampleindices_norc[c(1:5, 64)],
      expected_ind5_1_norc)
})

# revcomp
testthat::test_that("revcomp", {
   # test sequences
   test_dna <- c("AGAGTAAC", "GTGGCTTG", "GTCTGGAT",
      "CCAATCCG", "GACGCTAT", "TCAACGCG")
   exp_revcomp <- c("GTTACTCT", "CAAGCCAC", "ATCCAGAC",
      "CGGATTGG", "ATAGCGTC", "CGCGTTGA")

   testthat::expect_equal(
      revcomp(test_dna),
      exp_revcomp)
   testthat::expect_equal(
      revcomp(exp_revcomp),
      test_dna)
   testthat::expect_equal(
      revcomp(revcomp(test_dna)),
      test_dna)
})

testthat::test_that("import_geomx_GNP", {
   # use test file
   samplefile <- system.file("data", "GNP_config.ini", package="platjam")
   samplelist <- import_geomx_samplesheet(samplefile, return_type="dflist")
   # test list of vectors
   exp_lens <- setNames(c(2, 2, 4, 2),
      c("Sequencing",
         "Processing_v2",
         "AOI_List",
         "Targets"))
   exp_nrows <- setNames(c(2, 12, 64, 74),
      names(exp_lens))
   testthat::expect_equal(
      lengths(samplelist),
      exp_lens)
   testthat::expect_equal(
      sapply(samplelist, nrow),
      exp_nrows)
})
