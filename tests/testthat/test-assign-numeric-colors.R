
# context("assign_numeric_colors")


# utility function only used when creating tests below
# - it prints a string to assign to a variable via copy and paste
# - it returns the color vector to view with jamba::showColors()
expected_colors <- function(x, ramp="RdBu_r", values=TRUE, ...) {
   color_values <- assign_numeric_colors(x, color=ramp, ...)(x);
   names(color_values) <- x;
   expected_string <- paste0("c(",
      jamba::cPaste(
         paste0('"', color_values, '"'),
         sep=", "),
      ")\n");
   cat(expected_string);
   if (TRUE %in% values) {
      return(color_values)
   }
   invisible(expected_string)
}
## Example usage:
# expected_colors(-3:3, "RdBu_r", color_max=2, lens=2)
# jamba::showColors(expected_colors(-3:4, "RdBu_r", color_max=3, lens=2))

testthat::test_that("assign_numeric_colors", {
   # create expected color vector, then run test
   expected_x3 <- c("#053061FF", "#2D73AEFF", "#8ABEDAFF",
      "#F7F7F7FF", "#ED9E81FF", "#BA3237FF", "#67001FFF")
   x <- -3:3;
   testthat::expect_equal(
      assign_numeric_colors(x, color="RdBu_r")(x),
      expected_x3)

   expected_x34 <- c("#376191FF", "#6098C4FF", "#AFC6DEFF", "#F7F7F7FF",
      "#ECB1A6FF", "#D56B5BFF", "#9D3A3CFF", "#67001FFF")
   x <- -3:4;
   testthat::expect_equal(
      assign_numeric_colors(x, color="RdBu_r")(x),
      expected_x34)

   expected_x34_max2 <- c("#053061FF", "#053061FF", "#6098C4FF",
      "#F7F7F7FF", "#D56B5BFF", "#67001FFF", "#67001FFF", "#67001FFF")
   x <- -3:4;
   testthat::expect_equal(
      assign_numeric_colors(x, color="RdBu_r", color_max=2)(x),
      expected_x34_max2)

   expected_x34_max3_lens2 <- c("#053061FF", "#2569A6FF", "#77B3D4FF",
      "#F7F7F7FF", "#E98E70FF", "#B02631FF", "#67001FFF", "#67001FFF")
   x <- -3:4;
   testthat::expect_equal(
      assign_numeric_colors(x, color="RdBu_r", color_max=3, lens=2)(x),
      expected_x34_max3_lens2)
})
