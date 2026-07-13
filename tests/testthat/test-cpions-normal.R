test_that("CPions normal settings fixture matches expected output", {
    actual <- normal_fixture()

    expect_gt(nrow(actual), 0)
    expect_true(all(c(
        "Molecule_Formula", "Compound_Class", "Halo_perc", "Charge",
        "Adduct", "Adduct_Formula", "Isotopologue", "Isotope_Formula",
        "m/z", "Rel_ab"
    ) %in% names(actual)))

    expect_snapshot_value(actual, style = "json2")
})
