
context("curate_to_df_by_pattern")

# define test data
df <- data.frame(
   pattern=c("NOV14_p2w5_VEH",
      "NS644_p2w5VEH",
      "NOV14_p4w4_VEH",
      "NOV14_UL3_VEH",
      "NS644_UL3VEH",
      "NS50644_UL3VEH"),
   batch=c("NOV14",
      "NS644",
      "NOV14",
      "NOV14",
      "NS644",
      "NS50644"),
   group=c("p2w5_Veh",
      "p2w5_Veh",
      "p4w4_Veh",
      "UL3_Veh",
      "UL3_Veh",
      "UL3_Veh")
);
x <- c(
   "NOV14_p2w5_VEH_25_v2_CoordSort_deduplicated_SingleFrag_38to100.bam",
   "NOV14_p4w4_VEHrep1_25_v2_CoordSort_deduplicated_SingleFrag_38to100.bam",
   "NOV14_UL3_VEH_25_v2_CoordSort_deduplicated_SingleFrag_38to100.bam",
   "NS644_UL3VEH_25_v3_CoordSort_deduplicated_SingleFrag_38to100.bam",
   "NOV14_p2w5_VEH_50_v2_CoordSort_dedup_singleFragment.bam",
   "NOV14_UL3_VEH_50_v2_CoordSort_dedup_singleFragment.bam",
   "NS50644_UL3VEH_25_v3_CoordSort_deduplicated_SingleFrag.bam",
   "NS644_p2w5VEH_12p5_v3_CoordSort_deduplicated_SingleFrag_38to100.bam")

k <- c(1, 1, 2, 3, 4, 4, 5, 6)
dflabel <- paste0(df$batch[k], "_",
   df$group[k], "_rep", c(1,2,1,1,1,2,1,1));
df_out <- data.frame(
   Pattern=df$pattern[k],
   Batch=df$batch[k],
   Group=df$group[k],
   Label=dflabel,
   row.names=dflabel,
   Filename=x[c(1, 5, 8, 2, 3, 6, 4, 7)])
k2 <- c(1, 4, 5, 7, 2, 6, 8, 3);


testthat::test_that("curate_order_df", {
   testthat::expect_equal(
      curate_to_df_by_pattern(x, df),
      df_out)
})

testthat::test_that("curate_order_x", {
   testthat::expect_equal(
      curate_to_df_by_pattern(x, df, order_priority="x"),
      df_out[k2,])
})
