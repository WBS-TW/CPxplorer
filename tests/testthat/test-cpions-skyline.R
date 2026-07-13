test_that("CPions Skyline fixture matches expected output", {
    actual <- skyline_fixture()

    expect_identical(
        names(actual),
        c(
            "Molecule List Name",
            "Molecule Name",
            "Precursor Charge",
            "Label Type",
            "Precursor m/z",
            "Explicit Retention Time",
            "Explicit Retention Time Window",
            "Ion Selection Rank",
            "Original Rel_ab Rank",
            "Selected Reason",
            "Interfering Candidate Count",
            "Closest Interference m/z",
            "Closest Interference Molecule",
            "Resolution Needed",
            "Interference at MS Res?",
            "Note"
        )
    )
    expect_true(all(actual$`Label Type` %in% c("Quan", "Qual")))
    expect_equal(sum(actual$`Label Type` == "Quan"), 2L)
    expect_equal(sum(actual$`Label Type` == "Qual"), 4L)
    expect_true(all(actual$`Interference at MS Res?` == "NO" | is.na(actual$`Interference at MS Res?`)))
    expect_false(any(grepl("[INTERFERENCE]", actual$Note, fixed = TRUE)))
})

test_that("Skyline Note warns when final exported transition interferes", {
    source_data <- dplyr::bind_rows(
        CPxplorer:::getAdduct_normal("[PCA+Cl]-", 10:10, 3:3, 3L, 5L),
        CPxplorer:::getAdduct_normal("[PCO+Cl]-", 10:10, 3:3, 3L, 5L)
    )

    actual <- CPxplorer:::build_skyline_transition_list(
        source_data,
        mode = "normal",
        quant_ion = "Most intense",
        ms_resolution = 20000L,
        strategy = "balanced",
        preferred_qual_n = 2L
    )

    final_interference <- CPxplorer:::compute_transition_interference(actual, 20000L)
    warned_rows <- grepl("[INTERFERENCE]", actual$Note, fixed = TRUE)

    expect_identical(actual$`Interference at MS Res?`, final_interference$interference)
    expect_identical(warned_rows, actual$`Interference at MS Res?` == "YES")
    expect_false(any(warned_rows))
})

test_that("Filtered Skyline selection avoids interfered Quan ions", {
    actual <- skyline_dense_fixture()

    quan_rows <- dplyr::filter(actual, `Label Type` == "Quan")

    expect_gt(nrow(quan_rows), 0)
    expect_false(any(grepl("[INTERFERENCE]", quan_rows$Note, fixed = TRUE)))
    expect_true(all(table(actual$`Molecule Name`, actual$`Label Type`)[, "Quan"] == 1))
})

test_that("Filtered Skyline selection yields a final non-conflicting transition set", {
    actual <- skyline_dense_fixture()

    final_interference <- actual |>
        CPxplorer:::compute_transition_interference(ms_resolution = 20000L)

    non_interfered <- dplyr::filter(final_interference, interference == "NO")

    expect_equal(nrow(non_interfered), nrow(actual))
    expect_true(all(actual$`Label Type` %in% c("Quan", "Qual")))
    expect_true(all(actual$`Interference at MS Res?` == "NO"))
})

test_that("Filtered Skyline selection limits Qual ions per molecule", {
    actual <- skyline_dense_fixture()

    qual_counts <- actual |>
        dplyr::filter(`Label Type` == "Qual") |>
        dplyr::count(`Molecule Name`, name = "n_qual")

    expect_true(all(qual_counts$n_qual <= 2L))
})

test_that("Abundance Skyline strategy selects highest abundance Quan ions", {
    source_data <- tibble::tibble(
        Molecule_Formula = c("C10H17Cl3", "C10H17Cl3", "C11H19Cl3", "C11H19Cl3"),
        Adduct = "[PCA+Cl]-",
        Compound_Class = "PCA",
        Charge = -1L,
        `m/z` = c(100, 101, 200, 201),
        Rel_ab = c(40, 80, 60, 20)
    )

    actual <- CPxplorer:::build_skyline_transition_list(
        source_data,
        mode = "normal",
        quant_ion = "Interference-filtered",
        ms_resolution = 20000L,
        strategy = "abundance",
        preferred_qual_n = 0L
    )

    quan_rows <- dplyr::filter(actual, `Label Type` == "Quan")

    expect_equal(nrow(quan_rows), 2L)
    expect_true(all(quan_rows$`Original Rel_ab Rank` == 1L))
    expect_setequal(quan_rows$Note, c("{[PCA+Cl]-}{80}", "{[PCA+Cl]-}{60}"))
})

test_that("Interference Skyline strategy prefers non-interfering Quan ions", {
    source_data <- tibble::tibble(
        Molecule_Formula = c("C10H17Cl3", "C10H17Cl3", "C11H19Cl3", "C11H19Cl3"),
        Adduct = "[PCA+Cl]-",
        Compound_Class = "PCA",
        Charge = -1L,
        `m/z` = c(100, 101, 200, 201),
        Rel_ab = c(90, 80, 70, 60),
        interference = c("YES", "NO", "YES", "NO")
    )

    actual <- CPxplorer:::build_skyline_transition_list(
        source_data,
        mode = "normal",
        quant_ion = "Interference-filtered",
        ms_resolution = 20000L,
        strategy = "interference",
        preferred_qual_n = 0L
    )

    quan_rows <- dplyr::filter(actual, `Label Type` == "Quan")

    expect_equal(nrow(quan_rows), 2L)
    expect_true(all(quan_rows$`Original Rel_ab Rank` == 2L))
    expect_true(all(quan_rows$`Selected Reason` == "Least interference: higher-abundance ion interfered"))
})

test_that("Interference Skyline strategy falls back to least-interfering Quan ions", {
    source_data <- tibble::tibble(
        Molecule_Formula = c("C10H17Cl3", "C10H17Cl3", "C11H19Cl3"),
        Adduct = "[PCA+Cl]-",
        Compound_Class = "PCA",
        Charge = -1L,
        `m/z` = c(100.000, 300.000, 100.001),
        Rel_ab = c(90, 80, 70),
        interference = "YES"
    )

    actual <- CPxplorer:::build_skyline_transition_list(
        source_data,
        mode = "normal",
        quant_ion = "Interference-filtered",
        ms_resolution = 20000L,
        strategy = "interference",
        preferred_qual_n = 0L
    )

    target_quan <- actual |>
        dplyr::filter(`Molecule Name` == "C10H17Cl3", `Label Type` == "Quan")

    expect_equal(target_quan$`Precursor m/z`, 300)
    expect_equal(target_quan$`Original Rel_ab Rank`, 2L)
    expect_equal(target_quan$`Selected Reason`, "Least interference: all ions interfere")
})
