test_that("CPquant imports GC and LC demo Skyline reports", {
    for (mode in c("GC", "LC")) {
        df <- cpquant_demo_base_df(mode)

        expect_s3_class(df, "tbl_df")
        expect_gt(nrow(df), 0)
        expect_true(all(c(
            "Replicate_Name",
            "Sample_Type",
            "Molecule_List",
            "Molecule",
            "Area",
            "Isotope_Label_Type",
            "Transition_Note",
            "Rel_Ab",
            "Quantification_Group"
        ) %in% names(df)))
        expect_true(any(df$Sample_Type == "Standard"))
        expect_true(any(df$Sample_Type %in% c("Unknown", "Blank")))
        expect_true(any(df$Isotope_Label_Type == "Quan"))
        expect_true(all(is.finite(df$Area)))
        expect_false(any(is.na(df$Rel_Ab[!(df$Molecule_List %in% c("IS", "RS", "VS"))])))
    }
})

test_that("CPquant helper plots and isotope QC work with GC and LC demo data", {
    for (mode in c("GC", "LC")) {
        df <- cpquant_demo_base_df(mode, add_interference_note = TRUE)

        expect_s3_class(CPxplorer:::plot_skyline_output(df), "plotly")
        expect_s3_class(CPxplorer:::plot_quanqualratio(df), "plotly")
        expect_s3_class(CPxplorer:::plot_meas_vs_theor_ratio(df), "plotly")

        qc <- CPxplorer:::prepare_isotope_pattern_qc(df)
        expect_named(qc, c("ion_qc", "summary_qc"))
        expect_gt(nrow(qc$ion_qc), 0)
        expect_gt(nrow(qc$summary_qc), 0)
        expect_true(all(qc$summary_qc$qc_flag %in% c("Pass", "Review", "Fail", "No signal")))
        expect_false(any(is.infinite(qc$summary_qc$max_ion_ratio_error)))
    }
})

test_that("CPquant prepares calibration and sample data from GC and LC demos", {
    for (mode in c("GC", "LC")) {
        df <- cpquant_demo_base_df(mode)
        standards <- cpquant_prepare_standards(df)
        samples <- cpquant_prepare_samples(df)

        expect_gt(nrow(standards), 0)
        expect_gt(nrow(samples), 0)
        expect_true(all(c("RF", "intercept", "cal_rsquared", "Sum_RF_group") %in% names(standards)))
        expect_true(all(is.finite(standards$RF)))
        expect_true(all(is.finite(standards$Sum_RF_group)))
        expect_true(any(standards$RF > 0))
        expect_true(all(samples$Sample_Type %in% c("Unknown", "Blank")))
        expect_true(all(is.finite(samples$Relative_Area)))
    }
})

test_that("CPquant deconvolution helper accepts demo-shaped matrices", {
    for (mode in c("GC", "LC")) {
        df <- cpquant_demo_base_df(mode)
        standards <- cpquant_prepare_standards(df)
        samples <- cpquant_prepare_samples(df)

        standard_matrix <- standards |>
            dplyr::filter(RF > 0) |>
            dplyr::select(Molecule, Batch_Name, RF) |>
            dplyr::group_by(Molecule, Batch_Name) |>
            dplyr::summarise(RF = sum(RF, na.rm = TRUE), .groups = "drop") |>
            tidyr::pivot_wider(names_from = Batch_Name, values_from = RF, values_fill = 0) |>
            dplyr::arrange(Molecule)

        sample_data <- samples |>
            dplyr::filter(Sample_Type == "Unknown") |>
            dplyr::select(Replicate_Name, Molecule, Relative_Area) |>
            dplyr::group_by(Replicate_Name, Molecule) |>
            dplyr::summarise(Relative_Area = sum(Relative_Area, na.rm = TRUE), .groups = "drop") |>
            tidyr::nest(data = c(Molecule, Relative_Area)) |>
            dplyr::slice(1) |>
            dplyr::pull(data)

        sample_vector <- sample_data[[1]] |>
            dplyr::right_join(dplyr::select(standard_matrix, Molecule), by = "Molecule") |>
            dplyr::arrange(Molecule) |>
            dplyr::mutate(Relative_Area = tidyr::replace_na(Relative_Area, 0)) |>
            dplyr::select(Relative_Area)

        matrix_input <- standard_matrix |>
            tibble::column_to_rownames("Molecule") |>
            as.matrix()
        matrix_input <- sweep(matrix_input, 2, colSums(matrix_input), "/")
        matrix_input[!is.finite(matrix_input)] <- 0

        sum_rf <- matrix(colSums(matrix_input), nrow = 1)
        result <- CPxplorer:::perform_deconvolution(
            sample_vector,
            matrix_input,
            sum_rf,
            sample_name = paste(mode, "demo")
        )

        expect_named(result, c("sum_deconv_RF", "deconv_coef", "deconv_resolved", "deconv_rsquared"))
        expect_length(result$deconv_coef, ncol(matrix_input))
        expect_length(result$deconv_resolved, nrow(matrix_input))
        expect_true(all(is.finite(result$deconv_coef)))
        expect_true(all(is.finite(result$deconv_resolved)))
    }
})
