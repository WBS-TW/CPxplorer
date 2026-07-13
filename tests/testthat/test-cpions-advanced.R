test_that("CPions advanced settings fixture matches expected output", {
    actual <- advanced_fixture()

    expect_gt(nrow(actual), 0)
    expect_true(all(c(
        "Molecule_Formula", "Molecule_Halo_perc", "Compound_Class", "TP",
        "Charge", "Adduct_Annotation", "Adduct_Isotopologue",
        "Adduct_Formula", "Isotope_Formula", "m/z", "Rel_ab"
    ) %in% names(actual)))

    expect_snapshot_value(actual, style = "json2")
})
