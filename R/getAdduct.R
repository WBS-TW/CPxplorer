
getAdduct_normal <- function(adduct_ions, C, Cl, Clmax, threshold) {
        ####################################################################         
        # Regex to extract strings
        #################################################################### 
        ion_modes <- str_extract(adduct_ions, "(?<=\\]).{1}") # Using lookbehind assertion to extract ion mode
        fragment_ions <- str_extract(adduct_ions, "(?<=.{4}).+?(?=\\])") # extract after the 3rd character and before ]
        group <- str_extract(adduct_ions, "[^\\[].{2}") # Using positive lookbehind for [)
        
        if (group == "PCA") {
                data <- crossing(C, Cl) %>% #set combinations of C and Cl
                        filter(C >= Cl) %>% # filter so Cl dont exceed C atoms
                        filter(Cl <= Clmax) %>% # limit chlorine atoms. 
                        mutate(H = 2*C+2-Cl) %>% # add H atoms
                        mutate(Molecule_Formula = paste0("C", C, "H", H, "Cl", Cl)) %>% #add chemical formula
                        select(Molecule_Formula, C, H, Cl) # move Formula to first column
        } else if (group == "PCO") {
                data <- crossing(C, Cl) %>% 
                        filter(C >= Cl) %>% 
                        filter(Cl <= Clmax) %>% 
                        mutate(H = 2*C-Cl) %>% 
                        mutate(Molecule_Formula = paste0("C", C, "H", H, "Cl", Cl)) %>% 
                        select(Molecule_Formula, C, H, Cl) 
        }  else {
                print("Input not correct, only PCA or PCO is allowed")
        }
        
        
        # adding ion modes to the data frame to be inserted to isopattern, only -1 or +1 allowed
        if (ion_modes == "-") {
                data <- data %>%
                        mutate(Charge = as.integer(-1))
        }else if (ion_modes == "+") {
                data <- data %>%
                        mutate(Charge = as.integer(1))
        }
        
        ####################################################################       
        ####### generate input data for envipat based on fragment_ions
        ####################################################################         
        
        data <- generateInput_Envipat_normal(data = data, group = group, adduct_ions = adduct_ions, fragment_ions = fragment_ions)
        
        # Remove formula without Cl after adduct formations
        data <- data %>%
                filter(Cl > 0)
        
        # Create empty list for all ion formulas
        CP_allions <- list()
        data_ls <- list()
        
        #################################################################### 
        # function to get isotopic patterns for all PCAs. 
        # data("isotopes") needs to be loaded in app.R
        #################################################################### 
        getisotopes <- function(x) {enviPat::isopattern(isotopes = isotopes, 
                                                        chemforms = x, 
                                                        threshold = threshold, 
                                                        emass = 0.00054857990924,
                                                        plotit = FALSE, 
                                                        charge = Charge)}
        
        if (fragment_ions == "+Br") { #this is for Br adduct
                for (j in seq_along(data$Adduct_Formula)) {
                        Adduct_Formula <- data$Adduct_Formula[j]
                        Molecule_Formula <- data$Molecule_Formula[j]
                        Charge <- data$Charge[j]
                        Halo_perc <- data$Halo_perc[j]
                        dat <- getisotopes(x = as.character(data$Adduct_Formula[j]))
                        dat <- as.data.frame(dat[[1]])
                        dat <- dat %>% 
                                mutate(abundance = round(abundance, 1)) %>%
                                mutate(`m/z` = round(`m/z`, 6)) %>%
                                mutate(Isotope_Formula = paste0("[12C]", `12C`, "[13C]", `13C`, "[1H]", `1H`, "[2H]", `2H`, "[35Cl]", `35Cl`, "[37Cl]", `37Cl`, "[79Br]", `79Br`, "[81Br]", `81Br`)) %>%
                                mutate(Molecule_Formula = Molecule_Formula) %>%
                                mutate(Halo_perc = Halo_perc) %>%
                                mutate(Adduct_Formula =  Adduct_Formula) %>%
                                mutate(Charge = Charge) %>%
                                mutate(Isotopologue = case_when(
                                        `13C` + (`37Cl`+`81Br`)*2 == 0 ~ "",
                                        `13C` + (`37Cl`+`81Br`)*2 == 1 ~ "+1",
                                        `13C` + (`37Cl`+`81Br`)*2 == 2 ~ "+2",
                                        `13C` + (`37Cl`+`81Br`)*2 == 3 ~ "+3",
                                        `13C` + (`37Cl`+`81Br`)*2 == 4 ~ "+4",
                                        `13C` + (`37Cl`+`81Br`)*2 == 5 ~ "+5",
                                        `13C` + (`37Cl`+`81Br`)*2 == 6 ~ "+6",
                                        `13C` + (`37Cl`+`81Br`)*2 == 7 ~ "+7",
                                        `13C` + (`37Cl`+`81Br`)*2 == 8 ~ "+8",
                                        `13C` + (`37Cl`+`81Br`)*2 == 9 ~ "+9",
                                        `13C` + (`37Cl`+`81Br`)*2 == 10 ~ "+10",
                                        `13C` + (`37Cl`+`81Br`)*2 == 11 ~ "+11",
                                        `13C` + (`37Cl`+`81Br`)*2 == 12 ~ "+12",
                                        `13C` + (`37Cl`+`81Br`)*2 == 13 ~ "+13",
                                        `13C` + (`37Cl`+`81Br`)*2 == 14 ~ "+14",
                                        `13C` + (`37Cl`+`81Br`)*2 == 15 ~ "+15",
                                        `13C` + (`37Cl`+`81Br`)*2 == 16 ~ "+16",
                                        `13C` + (`37Cl`+`81Br`)*2 == 17 ~ "+17",
                                        `13C` + (`37Cl`+`81Br`)*2 == 18 ~ "+18",
                                        `13C` + (`37Cl`+`81Br`)*2 == 19 ~ "+19",
                                        `13C` + (`37Cl`+`81Br`)*2 == 20 ~ "+20")) %>%
                                mutate(Adduct = paste0(adduct_ions, " ", Isotopologue)) %>%
                                rename(Rel_ab = abundance) %>%
                                select(Molecule_Formula, Halo_perc, Charge, Adduct, Adduct_Formula, Isotopologue, Isotope_Formula, `m/z`, Rel_ab, `12C`, `13C`, `1H`, `2H`, `35Cl`, `37Cl`, `79Br`, `81Br`)
                        data_ls[[j]] <- dat
                }
        }else { # for other adducts
                for (j in seq_along(data$Adduct_Formula)) {
                        Adduct_Formula <- data$Adduct_Formula[j]
                        Molecule_Formula <- data$Molecule_Formula[j]
                        Charge <- data$Charge[j]
                        Halo_perc <- data$Halo_perc[j]
                        dat <- getisotopes(x = as.character(data$Adduct_Formula[j]))
                        dat <- as.data.frame(dat[[1]])
                        dat <- dat %>% 
                                mutate(abundance = round(abundance, 1)) %>%
                                mutate(`m/z` = round(`m/z`, 6)) %>%
                                mutate(Isotope_Formula = paste0("[12C]", `12C`, "[13C]", `13C`, "[1H]", `1H`, "[2H]", `2H`, "[35Cl]", `35Cl`, "[37Cl]", `37Cl`)) %>%
                                mutate(Molecule_Formula = Molecule_Formula) %>%
                                mutate(Halo_perc = Halo_perc) %>%
                                mutate(Adduct_Formula =  Adduct_Formula) %>%
                                mutate(Charge = Charge) %>%
                                mutate(Isotopologue = case_when(
                                        `13C` + (`37Cl`)*2 == 0 ~ "",
                                        `13C` + (`37Cl`)*2 == 1 ~ "+1",
                                        `13C` + (`37Cl`)*2 == 2 ~ "+2",
                                        `13C` + (`37Cl`)*2 == 3 ~ "+3",
                                        `13C` + (`37Cl`)*2 == 4 ~ "+4",
                                        `13C` + (`37Cl`)*2 == 5 ~ "+5",
                                        `13C` + (`37Cl`)*2 == 6 ~ "+6",
                                        `13C` + (`37Cl`)*2 == 7 ~ "+7",
                                        `13C` + (`37Cl`)*2 == 8 ~ "+8",
                                        `13C` + (`37Cl`)*2 == 9 ~ "+9",
                                        `13C` + (`37Cl`)*2 == 10 ~ "+10",
                                        `13C` + (`37Cl`)*2 == 11 ~ "+11",
                                        `13C` + (`37Cl`)*2 == 12 ~ "+12",
                                        `13C` + (`37Cl`)*2 == 13 ~ "+13",
                                        `13C` + (`37Cl`)*2 == 14 ~ "+14",
                                        `13C` + (`37Cl`)*2 == 15 ~ "+15",
                                        `13C` + (`37Cl`)*2 == 16 ~ "+16",
                                        `13C` + (`37Cl`)*2 == 17 ~ "+17",
                                        `13C` + (`37Cl`)*2 == 18 ~ "+18",
                                        `13C` + (`37Cl`)*2 == 19 ~ "+19",
                                        `13C` + (`37Cl`)*2 == 20 ~ "+20")) %>%
                                mutate(Adduct = paste0(adduct_ions, " ", Isotopologue)) %>%
                                rename(Rel_ab = abundance) %>%
                                select(Molecule_Formula, Halo_perc, Charge, Adduct, Adduct_Formula, Isotopologue, Isotope_Formula, `m/z`, Rel_ab, `12C`, `13C`, `1H`, `2H`, `35Cl`, `37Cl`)
                        data_ls[[j]] <- dat
                }
        }
        
        # combine all elements in list list to get dataframe
        data_ls <- do.call(rbind, data_ls)
        
        
        # combine both all adduct ions
        CP_allions <- rbind(CP_allions, data_ls)
        return(CP_allions)
        
}



#----------------------------------------------------------------------------------------------------

getAdduct_BCA <- function(adduct_ions, C, Cl, Br, Clmax, Brmax, threshold) {
        ####################################################################         
        # Regex to extract strings
        #################################################################### 
        ion_modes <- str_extract(adduct_ions, "(?<=\\]).{1}") # Using lookbehind assertion to extract ion mode
        fragment_ions <- str_extract(adduct_ions, "(?<=.{4}).+?(?=\\])") # extract after the 3rd character and before ]
        group <- str_extract(adduct_ions, "[^\\[].{2}") # Using positive lookbehind for [)
        
        if (group == "BCA") {
                data <- crossing(C, Cl, Br) %>% #get combinations of C, Cl, Br
                        filter(C >= Cl) %>% # filter so Cl dont exceed C atoms
                        filter(Cl <= Clmax) %>% # limit chlorine atoms.
                        filter(Br <= Brmax) %>%
                        filter(Br + Cl <= C) %>%
                        mutate(H = 2*C+2-Cl-Br) %>% # add H atoms
                        mutate(Molecule_Formula = paste0("C", C, "H", H, "Cl", Cl, "Br", Br)) %>% #add chemical formula
                        select(Molecule_Formula, C, H, Cl, Br) # move Formula to first column
        }  else {
                print("Input not correct, only BCA is allowed")
        }
        
        # check chem_forms
        # if (any(check_chemform(isotopes = isotopes, chemforms = data$Formula)$warning == TRUE)) {print("Warning: incorrect formula")} else {"All correct"}
        
        # adding ion modes to the data frame to be inserted to isopattern, only -1 or +1 allowed
        if (ion_modes == "-") {
                data <- data %>%
                        mutate(Charge = as.integer(-1))
        }else if (ion_modes == "+") {
                data <- data %>%
                        mutate(Charge = as.integer(1))
        }
        
        ####################################################################       
        ####### generate input data for envipat based on fragment_ions
        #################################################################### 
        
        
        data <- generateInput_Envipat_BCA(data = data, group = group, adduct_ions = adduct_ions, fragment_ions = fragment_ions)      
        
        
        
        # Remove formula without Cl after adduct formations
        data <- data %>%
                filter(Cl > 0)
        
        # Create empty list for all ion formulas
        CP_allions <- list()
        data_ls <- list()
        
        
        # function to get isotopic patterns for all PCAs. Threshold based on the app, neutral form. data("isotopes") needs to be loaded first
        getisotopes <- function(x) {enviPat::isopattern(isotopes = isotopes, 
                                                        chemforms = x, 
                                                        threshold = threshold, 
                                                        emass = 0.00054857990924,
                                                        plotit = FALSE, 
                                                        charge = Charge)}
        
        
        for (j in seq_along(data$Adduct_Formula)) {
                Adduct_Formula <- data$Adduct_Formula[j]
                Molecule_Formula <- data$Molecule_Formula[j]
                Charge <- data$Charge[j]
                Halo_perc <- data$Halo_perc[j]
                dat <- getisotopes(x = as.character(data$Adduct_Formula[j]))
                dat <- as.data.frame(dat[[1]])
                dat <- dat %>% 
                        mutate(abundance = round(abundance, 1)) %>%
                        mutate(`m/z` = round(`m/z`, 6)) %>%
                        mutate(Isotope_Formula = paste0("[12C]", `12C`, "[13C]", `13C`, "[1H]", `1H`, "[2H]", `2H`, "[35Cl]", `35Cl`, "[37Cl]", `37Cl`, "[79Br]", `79Br`, "[81Br]", `81Br`)) %>%
                        mutate(Molecule_Formula = Molecule_Formula) %>%
                        mutate(Halo_perc = Halo_perc) %>%
                        mutate(Adduct_Formula = Adduct_Formula) %>%
                        mutate(Charge = Charge) %>%
                        mutate(Isotopologue = case_when(
                                `13C` + (`37Cl`+`81Br`)*2 == 0 ~ "",
                                `13C` + (`37Cl`+`81Br`)*2 == 1 ~ "+1",
                                `13C` + (`37Cl`+`81Br`)*2 == 2 ~ "+2",
                                `13C` + (`37Cl`+`81Br`)*2 == 3 ~ "+3",
                                `13C` + (`37Cl`+`81Br`)*2 == 4 ~ "+4",
                                `13C` + (`37Cl`+`81Br`)*2 == 5 ~ "+5",
                                `13C` + (`37Cl`+`81Br`)*2 == 6 ~ "+6",
                                `13C` + (`37Cl`+`81Br`)*2 == 7 ~ "+7",
                                `13C` + (`37Cl`+`81Br`)*2 == 8 ~ "+8",
                                `13C` + (`37Cl`+`81Br`)*2 == 9 ~ "+9",
                                `13C` + (`37Cl`+`81Br`)*2 == 10 ~ "+10",
                                `13C` + (`37Cl`+`81Br`)*2 == 11 ~ "+11",
                                `13C` + (`37Cl`+`81Br`)*2 == 12 ~ "+12",
                                `13C` + (`37Cl`+`81Br`)*2 == 13 ~ "+13",
                                `13C` + (`37Cl`+`81Br`)*2 == 14 ~ "+14",
                                `13C` + (`37Cl`+`81Br`)*2 == 15 ~ "+15",
                                `13C` + (`37Cl`+`81Br`)*2 == 16 ~ "+16",
                                `13C` + (`37Cl`+`81Br`)*2 == 17 ~ "+17",
                                `13C` + (`37Cl`+`81Br`)*2 == 18 ~ "+18",
                                `13C` + (`37Cl`+`81Br`)*2 == 19 ~ "+19",
                                `13C` + (`37Cl`+`81Br`)*2 == 20 ~ "+20")) %>%
                        mutate(Adduct = paste0(adduct_ions, " ", Isotopologue)) %>%
                        rename(Rel_ab = abundance) %>%
                        select(Molecule_Formula, Halo_perc, Charge, Adduct, Adduct_Formula, Isotopologue, Isotope_Formula, `m/z`, Rel_ab, `12C`, `13C`, `1H`, `2H`, `35Cl`, `37Cl`, `79Br`, `81Br`)
                data_ls[[j]] <- dat
        }
        
        
        # combine all elements in list list to get dataframe
        data_ls <- do.call(rbind, data_ls)
        
        
        # combine both all adduct ions
        CP_allions <- rbind(CP_allions, data_ls)
        return(CP_allions)
        
}

#----------------------------------------------------------------------------------------------------

getAdduct_advanced <- function(Compounds, Adduct_Ion, TP, Charge, C, Cl, Clmax, Br, Brmax, threshold) {
        ####################################################################         
        # Regex to extract strings
        #################################################################### 
        
        
        if (Compounds == "PCA") {
                data <- crossing(C, Cl) |> #set combinations of C and Cl
                        filter(C >= Cl) |> # filter so Cl dont exceed C atoms
                        filter(Cl <= Clmax) |> # limit chlorine atoms. 
                        mutate(H = case_when(# add H atoms
                                TP == "None" ~ 2*C+2-Cl,
                                TP == "+OH" ~ 2*C+2-Cl,
                                TP == "+2OH" ~ 2*C+2-Cl,
                                TP == "-Cl+OH" ~ 2*C+2-Cl+1,
                                TP == "-2Cl+2OH" ~ 2*C+2-Cl+2,
                                TP == "+SO4" ~ 2*C+2-Cl-1))  |> 
                        mutate(Cl = case_when(
                                TP == "-Cl+OH" ~ Cl-1,
                                TP == "-2Cl+2OH" ~ Cl-2,
                                .default = Cl)) |> 
                        mutate(Molecule_Formula = paste0("C", C, "H", H, "Cl", Cl)) |> 
                        mutate(Molecule_Formula = case_when( #DOUBLE CHECK THE FORMULA IS CORRECT!!!!!
                                TP == "None" ~ paste0("C", C, "H", H, "Cl", Cl),
                                TP == "+OH" ~ paste0("C", C, "H", H, "Cl", Cl, "O"),
                                TP == "+2OH" ~ paste0("C", C, "H", H, "Cl", Cl, "O2"),
                                TP == "-Cl+OH" ~ paste0("C", C, "H", H, "Cl", Cl, "O"),
                                TP == "-2Cl+2OH" ~ paste0("C", C, "H", H, "Cl", Cl, "O2"),
                                TP == "+SO4" ~ paste0("C", C, "H", H, "Cl", Cl, "SO4"))) 

        } else if (Compounds == "PCO") {
                data <- crossing(C, Cl) |> 
                        filter(C >= Cl) |> 
                        filter(Cl <= Clmax) |> 
                        mutate(H = case_when(# add H atoms. 
                                TP == "None" ~ 2*C-Cl, #PCO general formula 2*C-Cl
                                TP == "+OH" ~ 2*C-Cl,
                                TP == "+2OH" ~ 2*C-Cl,
                                TP == "-Cl+OH" ~ 2*C-Cl+1,
                                TP == "-2Cl+2OH" ~ 2*C-Cl+2,
                                TP == "+SO4" ~ 2*C-Cl-1))  |> 
                        mutate(Cl = case_when(
                                TP == "-Cl+OH" ~ Cl-1,
                                TP == "-2Cl+2OH" ~ Cl-2,
                                .default = Cl)) |> 
                        mutate(Molecule_Formula = paste0("C", C, "H", H, "Cl", Cl)) |> 
                        mutate(Molecule_Formula = case_when( #DOUBLE CHECK THE FORMULA IS CORRECT!!!!!
                                TP == "None" ~ paste0("C", C, "H", H, "Cl", Cl),
                                TP == "+OH" ~ paste0("C", C, "H", H, "Cl", Cl, "O"),
                                TP == "+2OH" ~ paste0("C", C, "H", H, "Cl", Cl, "O2"),
                                TP == "-Cl+OH" ~ paste0("C", C, "H", H, "Cl", Cl, "O"),
                                TP == "-2Cl+2OH" ~ paste0("C", C, "H", H, "Cl", Cl, "O2"),
                                TP == "+SO4" ~ paste0("C", C, "H", H, "Cl", Cl, "SO4"))) 
                
        } else if (Compounds == "BCA") {
                data <- crossing(C, Cl, Br) |>  #get combinations of C, Cl, Br
                        filter(C >= Cl) |>  # filter so Cl dont exceed C atoms
                        filter(Cl <= Clmax) |>  # limit chlorine atoms.
                        filter(Br <= Brmax) |> 
                        filter(Br + Cl <= C) |> 
                        mutate(H = case_when(# add H atoms. 
                                TP == "None" ~ 2*C+2-Cl-Br, #BCA general formula
                                TP == "+OH" ~ 2*C+2-Cl-Br,
                                TP == "+2OH" ~ 2*C+2-Cl-Br,
                                TP == "-Cl+OH" ~ 2*C+2-Cl-Br+1,
                                TP == "-2Cl+2OH" ~ 2*C+2-Cl-Br+2,
                                TP == "+SO4" ~ 2*C+2-Cl-Br-1))  |> 
                        
                        mutate(Molecule_Formula = case_when( #DOUBLE CHECK THE FORMULA IS CORRECT!!!!!
                                TP == "None" ~ paste0("C", C, "H", H, "Cl", Cl, "Br", Br),
                                TP == "+OH" ~ paste0("C", C, "H", H, "Cl", Cl, "Br", Br, "O"),
                                TP == "+2OH" ~ paste0("C", C, "H", H, "Cl", Cl, "Br", Br, "O2"),
                                TP == "-Cl+OH" ~ paste0("C", C, "H", H, "Cl", Cl, "Br", Br, "O"),
                                TP == "-2Cl+2OH" ~ paste0("C", C, "H", H, "Cl", Cl, "Br", Br, "O2"),
                                TP == "+SO4" ~ paste0("C", C, "H", H, "Cl", Cl, "Br", Br, "SO4"))) 
                
                
        }
        
        
        
        # adding ion modes to the data frame to be inserted to isopattern, only -1 or +1 allowed
        if (Charge == "-") {
                data <- data |>
                        mutate(Charge = as.integer(-1))
        }else if (Charge == "+") {
                data <- data |>
                        mutate(Charge = as.integer(1))
        }
        
        ####################################################################       
        ####### generate input data for envipat based on adduct ions
        ####################################################################         
        
        data <- generateInput_Envipat_advanced(data = data, Compounds = Compounds, Adduct_Ion = Adduct_Ion, TP = TP, Charge = Charge)
        
        # Remove formula without Cl after adduct formations
        data <- data |>
                filter(Cl > 0)
        
        # Create empty list for all ion formulas
        CP_allions <- list()
        data_ls <- list()
        
        #################################################################### 
        # function to get isotopic patterns for all PCA/PCO/BCA. 
        # data("isotopes") needs to be loaded in app.R
        #################################################################### 
        getisotopes <- function(x) {enviPat::isopattern(isotopes = isotopes, 
                                                        chemforms = x, 
                                                        threshold = threshold, 
                                                        emass = 0.00054857990924,
                                                        plotit = FALSE, 
                                                        charge = Charge)}
        
        for (j in seq_along(data$Molecule_Formula)) {
                Adduct_Formula <- data$Adduct_Formula[j]
                Molecule_Formula <- data$Molecule_Formula[j]
                Compound_Class <- data$Compound_Class[j]
                TP <- data$TP[j]
                Charge <- data$Charge[j]
                Molecule_Halo_perc <- data$Molecule_Halo_perc[j]
                Adduct_Annotation <- data$Adduct_Annotation[j]
                dat <- getisotopes(x = as.character(data$Adduct_Formula[j]))
                dat <- as.data.frame(dat[[1]])
                
                dat <- dat |> 
                        mutate(abundance = round(abundance, 1)) |>
                        mutate(`m/z` = round(`m/z`, 6)) 
                
                dat <- create_elements(dat) |> 
                        mutate(Isotope_Formula = create_formula_isotope(`12C`,`13C`, `1H`, `2H`, `35Cl`, `37Cl`, `79Br`, `81Br`, 
                                                                        `16O`, `17O`, `18O`, `32S`, `33S`, `34S`, `36S`)) |> 
                        mutate(Molecule_Formula = Molecule_Formula) |>
                        mutate(Molecule_Halo_perc = Molecule_Halo_perc) |>
                        mutate(Compound_Class = Compound_Class) |> 
                        mutate(TP = TP) |> 
                        mutate(Adduct_Annotation =  Adduct_Annotation) |>
                        mutate(Adduct_Formula =  Adduct_Formula) |>
                        mutate(Charge = Charge) |>
                        mutate(Isotopologue = case_when(
                                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 0 ~ "",
                                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 1 ~ "+1",
                                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 2 ~ "+2",
                                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 3 ~ "+3",
                                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 4 ~ "+4",
                                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 5 ~ "+5",
                                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 6 ~ "+6",
                                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 7 ~ "+7",
                                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 8 ~ "+8",
                                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 9 ~ "+9",
                                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 10 ~ "+10",
                                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 11 ~ "+11",
                                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 12 ~ "+12",
                                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 13 ~ "+13",
                                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 14 ~ "+14",
                                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 15 ~ "+15",
                                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 16 ~ "+16",
                                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 17 ~ "+17",
                                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 18 ~ "+18",
                                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 19 ~ "+19",
                                `13C` + (`37Cl`+`81Br` + `18O` + `34S`)*2 == 20 ~ "+20")) |> 
                        mutate(Adduct_Isotopologue = paste0(Adduct_Ion, " ", Isotopologue)) |>
                        rename(Rel_ab = abundance) |>
                        select(Molecule_Formula, Molecule_Halo_perc, Compound_Class, TP, Charge, Adduct_Annotation, Adduct_Isotopologue, Adduct_Formula, Isotopologue, Isotope_Formula, `m/z`, Rel_ab, `12C`, `13C`, `1H`, `2H`, `35Cl`, `37Cl`, everything())
                data_ls[[j]] <- dat
        }
        
        
        # combine all elements in list list to get dataframe
        data_ls <- do.call(rbind, data_ls)
        
        
        # combine both all adduct ions
        CP_allions <- rbind(CP_allions, data_ls)
        return(CP_allions)
        
}
