test_that("CPions app tab flow renders expected outputs", {
    skip_if_not_installed("shiny")

    app <- shiny::testServer(CPxplorer:::CPions_server, {
        session$setInputs(
            Cmin = 10,
            Cmax = 10,
            Clmin = 3,
            Clmax = 3,
            Brmin = 1,
            Brmax = 1,
            Adducts = "[PCA+Cl]-",
            threshold = 5,
            ISRS_input = "",
            go1 = 1
        )

        session$flushReact()
        expect_true(!is.null(output$Table_norm))

        session$setInputs(
            Cmin_adv = 10,
            Cmax_adv = 10,
            Clmin_adv = 3,
            Clmax_adv = 3,
            Brmin_adv = 1,
            Brmax_adv = 1,
            Compclass_adv = "PCA",
            Adducts_adv = "+Cl",
            Charge_adv = "-",
            TP_adv = "None",
            threshold_adv = 5,
            ISRS_input_adv = "",
            go_adv = 1
        )

        session$flushReact()
        expect_true(!is.null(output$Table_adv))

        session$setInputs(MSresolution = 20000, interfere_mode = "normal", go2 = 1)
        session$flushReact()
        expect_true(!is.null(output$Plotly))
        expect_true(!is.null(output$Plotly2))
        expect_true(!is.null(output$Table2))

        session$setInputs(
            QuantIon = "Interference-filtered",
            skyline_strategy = "Balanced",
            skyline_qual_n = 2,
            skylineoutput = "mz",
            skyline_mode = "normal",
            go3 = 1
        )
        session$flushReact()
        expect_true(!is.null(output$Table3))
    })
})
