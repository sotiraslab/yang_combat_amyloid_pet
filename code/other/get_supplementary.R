rm(list = ls())

library(tidyverse)

roi_df <- read_rds("data/csv/roi.RDS")
supp_df <- roi_df %>%
    filter(study == "oasis", feature_type == "suvr") %>%
    select(fs_label, roi_type) %>%
    mutate(
        roi_type = if_else(
            str_detect(fs_label, "cerebellum.cortex"),
            "Reference ROI", roi_type
        )
    ) %>%
    group_by(roi_type) %>%
    arrange(fs_label, .by_group = TRUE) %>%
    ungroup() %>%
    rename(
        `FreeSurfer ROI` = fs_label,
        `ROI subgroup` = roi_type
    )

supp_df %>%
    write_csv("tables/roi_list.csv")
