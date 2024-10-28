# ========== VARIABLES ==========

# harmonization method levels and labels
harm_method_levels <- c(
    "suvr",
    "cl",
    "combat_nocovar",
    "combat_linear",
    "combat_gam",
    "peace__nocovar",
    "peace__age_sex_apoe"
)
harm_method_labels <- c(
    "unharmonized",
    "Centiloid, scaled",
    "ComBat, no covariates",
    "ComBat + age, sex, APOE",
    "GAM-ComBat + age, sex, APOE",
    "PEACE, no covariates",
    "PEACE + age, sex, APOE"
)
harm_method_dict <- harm_method_labels
names(harm_method_dict) <- harm_method_levels

# ggseg plot layout
ggseg_layout <- "
ABCD
EFGH
"
ggseg_layout_all <- "
ABCD##
EFGHIJ
"

# ========== FUNCTIONS ==========

mutate_combat <- function(.data) {

    # convert datatypes of certain columns (sex, apoe, clinical group, amyloid-positive)

    .data <- .data %>%
        mutate(
            across(
                all_of(c("sex", "apoe", "clinical_group", "clinical_group_extended")),
                factor
            ),
            across(
                starts_with("age"),
                round_age
            )
        )

    return(.data)

}

harm_method_to_factor <- function(.data, harm_method_col = "harmonization_method") {

    # convert to factor col
    .data <- .data %>%
        mutate(
            across(all_of(harm_method_col), ~ factor(.x, harm_method_levels, harm_method_labels) %>% droplevels)
        )

    return(.data)

}

round_age <- function(age, digits = 2) {

    return(round(age, digits = digits))

}