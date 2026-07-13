cpquant_demo_path <- function(mode) {
    file <- switch(
        mode,
        GC = file.path("GC", "SkylineReportGC_NCI_Orbi_demo.xlsx"),
        LC = file.path("LC", "SkylineReportLC_APCI_Orbi_demo.xlsx"),
        stop("Unknown CPquant demo mode")
    )

    system.file("cpxplorer-demo_2026-04-01", file, package = "CPxplorer", mustWork = TRUE)
}

cpquant_demo_base_df <- function(mode, add_interference_note = FALSE) {
    df <- readxl::read_excel(
        cpquant_demo_path(mode),
        guess_max = 5000,
        na = c("", "NA", "#N/A", "N/A")
    ) |>
        dplyr::rename(Molecule_List = tidyr::any_of(c("Molecule List", "MoleculeList"))) |>
        dplyr::rename(Transition_Note = tidyr::any_of(c("Transition Note", "TransitionNote")))

    if (add_interference_note) {
        first_cp_note <- which(!(df$Molecule_List %in% c("IS", "RS", "VS")))[[1]]
        df$Transition_Note[[first_cp_note]] <- paste0(df$Transition_Note[[first_cp_note]], "[INTERFERENCE]")
    }

    CPxplorer:::normalize_cpquant_skyline(df)
}

cpquant_prepare_standards <- function(df, min_rsquared = 0.9) {
    df |>
        dplyr::filter(
            Sample_Type == "Standard",
            !Molecule_List %in% c("IS", "RS"),
            Isotope_Label_Type == "Quan",
            Batch_Name != "NA"
        ) |>
        dplyr::mutate(
            C_range = stringr::str_extract_all(Quantification_Group, "\\d+"),
            C_min = as.numeric(purrr::map_chr(C_range, ~ .x[1])),
            C_max = as.numeric(purrr::map_chr(C_range, function(x) if (length(x) > 1) x[2] else x[1]))
        ) |>
        dplyr::select(-C_range) |>
        dplyr::group_by(Batch_Name, Sample_Type, Molecule, Molecule_List, C_number, Cl_number, PCA, Quantification_Group, C_min, C_max) |>
        tidyr::nest() |>
        dplyr::filter(C_number >= C_min & C_number <= C_max) |>
        dplyr::mutate(models = purrr::map(data, ~ stats::lm(Area ~ Analyte_Concentration, data = .x))) |>
        dplyr::mutate(coef = purrr::map(models, stats::coef)) |>
        dplyr::mutate(RF = purrr::map_dbl(models, ~ stats::coef(.x)["Analyte_Concentration"])) |>
        dplyr::mutate(intercept = purrr::map(coef, purrr::pluck("(Intercept)"))) |>
        dplyr::mutate(cal_rsquared = purrr::map(models, summary)) |>
        dplyr::mutate(cal_rsquared = purrr::map(cal_rsquared, purrr::pluck("r.squared"))) |>
        dplyr::select(-coef) |>
        tidyr::unnest(c(RF, intercept, cal_rsquared)) |>
        dplyr::mutate(RF = dplyr::if_else(RF < 0, 0, RF)) |>
        dplyr::mutate(cal_rsquared = ifelse(is.nan(cal_rsquared), 0, cal_rsquared)) |>
        dplyr::mutate(RF = dplyr::if_else(cal_rsquared < min_rsquared, 0, RF)) |>
        dplyr::ungroup() |>
        dplyr::group_by(Batch_Name) |>
        dplyr::mutate(Sum_RF_group = sum(RF, na.rm = TRUE)) |>
        dplyr::ungroup()
}

cpquant_prepare_samples <- function(df) {
    df |>
        dplyr::filter(
            Sample_Type %in% c("Unknown", "Blank"),
            !Molecule_List %in% c("IS", "RS", "VS"),
            Isotope_Label_Type == "Quan"
        ) |>
        dplyr::group_by(Replicate_Name) |>
        dplyr::mutate(Relative_Area = Area / sum(Area, na.rm = TRUE)) |>
        dplyr::ungroup() |>
        dplyr::select(-Mass_Error_PPM, -Isotope_Label_Type, -Chromatogram_Precursor_MZ, -Analyte_Concentration, -Batch_Name) |>
        dplyr::mutate(dplyr::across(Relative_Area, ~ replace(., is.nan(.), 0)))
}
