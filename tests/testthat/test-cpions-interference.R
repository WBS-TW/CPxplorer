test_that("CPions interfering ions fixture matches expected output", {
    actual <- interference_fixture()

    expect_true(all(c(
        "difflag", "difflead", "reslag", "reslead", "interference"
    ) %in% names(actual)))
    expect_type(actual$interference, "character")

    expect_snapshot_value(
        dplyr::select(
            actual,
            Molecule_Formula,
            Adduct,
            `m/z`,
            Rel_ab,
            difflag,
            difflead,
            reslag,
            reslead,
            interference
        ),
        style = "json2"
    )
})
