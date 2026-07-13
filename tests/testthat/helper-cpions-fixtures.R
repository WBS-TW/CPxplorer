normalize_tbl <- function(x) {
    order_cols <- intersect(c(
        "Molecule_Formula",
        "Compound_Class",
        "TP",
        "Adduct",
        "Adduct_Annotation",
        "Isotopologue",
        "m/z"
    ), names(x))

    if (length(order_cols) == 0) {
        return(x)
    }

    dplyr::arrange(x, dplyr::across(dplyr::all_of(order_cols)))
}

normal_fixture <- function() {
    CPxplorer:::getAdduct_normal(
        adduct_ions = "[PCA+Cl]-",
        C = 10:10,
        Cl = 3:3,
        Clmax = 3L,
        threshold = 5L
    ) |>
        normalize_tbl()
}

advanced_fixture <- function() {
    CPxplorer:::getAdduct_advanced(
        Class = "PCA",
        Adduct_Ion = "+Cl",
        TP = "None",
        Charge = "-",
        C = 10:10,
        Cl = 3:3,
        Clmax = 3L,
        Br = 0:0,
        Brmax = 0L,
        threshold = 5L
    ) |>
        normalize_tbl()
}

interference_fixture <- function() {
    dplyr::bind_rows(
        CPxplorer:::getAdduct_normal("[PCA+Cl]-", 10:10, 3:3, 3L, 5L),
        CPxplorer:::getAdduct_normal("[PCO+Cl]-", 10:10, 3:3, 3L, 5L)
    ) |>
        CPxplorer:::compute_interference(ms_resolution = 20000L) |>
        normalize_tbl()
}

skyline_fixture <- function() {
    interference_fixture() |>
        CPxplorer:::build_skyline_transition_list(
            mode = "normal",
            quant_ion = "Interference-filtered",
            ms_resolution = 20000L,
            strategy = "balanced",
            preferred_qual_n = 2L
        )
}

skyline_dense_fixture <- function() {
    dplyr::bind_rows(
        CPxplorer:::getAdduct_normal("[PCA+Cl]-", 10:12, 3:5, 5L, 5L),
        CPxplorer:::getAdduct_normal("[PCO+Cl]-", 10:12, 3:5, 5L, 5L)
    ) |>
        CPxplorer:::compute_interference(ms_resolution = 20000L) |>
        CPxplorer:::build_skyline_transition_list(
            mode = "normal",
            quant_ion = "Interference-filtered",
            ms_resolution = 20000L,
            strategy = "balanced",
            preferred_qual_n = 2L
        )
}
