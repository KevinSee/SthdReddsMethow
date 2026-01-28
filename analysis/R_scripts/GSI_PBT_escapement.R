# Author: Kevin See
# Purpose: create table of escapement estimates in the Methow, by PBT / GSI assignment
# Created: 1/28/26
# Last Modified: 1/28/26
# Notes:

#-----------------------------------------------------------------
# load needed libraries
library(tidyverse)
# library(PITcleanr)
library(janitor)
# library(readxl)
library(writexl)
# library(magrittr)
library(msm)
library(here)
# library(DescTools)
library(sroem)



#-----------------------------------------------------------------
# set year
yr = 2025


# get escapement estimates by origin
esc_est <-
  query_dabom_results(dabom_dam_nm = "RockIsland",
                      query_year = yr,
                      result_type = "escape_summ")

# set use black box estimates at intermediate sites (with sites upstream)
bb_esc_est <-
  esc_est |>
  mutate(final_site = str_remove(location, "_bb$")) |>
  group_by(final_site) |>
  mutate(n_est = n_distinct(location),
         keep = case_when(n_est == 1 ~ T,
                          n_est == 2 &
                            str_detect(location, "_bb") ~ T,
                          .default = F)) |>
  ungroup() |>
  filter(keep) |>
  select(-c(n_est,
            keep,
            location)) |>
  relocate(final_site,
           .after = "spawn_year")

# get tag summary, and filter for tags in the Methow
tag_df <-
  query_dabom_results(dabom_dam_nm = "RockIsland",
                      query_year = yr,
                      result_type = "tag_summ") |>
  filter(group == "Methow") |>
  mutate(final_site = str_remove(final_node, "_U$"),
         across(final_site,
                ~ str_remove(., "_D$")))

# look at proportion of tags in each "patch" by PBT / GSI assignment
# start with hatchery tags by PBT assignment
met_escp_prop_df <-
  tag_df |>
  filter(origin == "H") |>
  mutate(across(popname,
                ~ str_split_i(., "_", 1)),
         across(popname,
                ~ case_when(origin == "W" ~ "---",
                            origin == "H" &
                              is.na(.) ~ "---",
                            .default = .))) |>
  count(final_site,
        origin,
        popname,
        name = "n_tags") |>
  add_column(data_source = "PBT",
             .before = "popname") |>
  # add wild tags with GSI assignments
  bind_rows(
    tag_df |>
      filter(origin == "W") |>
      count(final_site,
            origin,
            popname = gsi_assignment,
            name = "n_tags") |>
      add_column(data_source = "GSI",
                 .before = "popname")
  ) |>
  arrange(final_site,
          origin,
          popname) |>
  group_by(final_site,
           origin) |>
  mutate(prop_tags = n_tags / sum(n_tags),
         prop_tags_se = sqrt((prop_tags * (1 - prop_tags)) / sum(n_tags))) |>
  ungroup() |>
  # add escapement estimates
  left_join(bb_esc_est |>
              select(final_site,
                     origin,
                     tot_escp = median,
                     tot_escp_se = sd),
            by = join_by(final_site,
                         origin)) |>
  mutate(grp_escp = prop_tags * tot_escp) |>
  rowwise() |>
  mutate(grp_escp_se = msm::deltamethod(~ x1 * x2,
                                        mean = c(prop_tags,
                                                 tot_escp),
                                        cov = diag(c(prop_tags_se,
                                                     tot_escp_se)^2))) |>
  ungroup() |>
  group_by(final_site) |>
  mutate(prop_tot_escp = grp_escp / sum(grp_escp)) |>
  ungroup() |>
  mutate(across(origin,
                ~ case_match(.,
                             "H" ~ "HOR",
                             "W" ~ "NOR",
                             .default = .)))

list("Methow" = met_escp_prop_df) |>
  write_xlsx(path = here("analysis/data/derived_data",
                         "Methow_PBT_GSI_escapement.xlsx"))
