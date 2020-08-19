# summarizes distributional statistics related to allocaiton results.
# For feasible allocation, what are the averages attributes for each cell

# # Marginal Propensity to Consume Check Tables

# Files:
# source('fs_opti_support.R')
ls_output <- fs_opti_support()
srt_simu_path <- ls_output$srt_simu_path
srt_csv_path_root <- ls_output$srt_csv_path_root

# st_ds_file_name <- 'snwx_ds_small.csv'
st_ds_file_name <- 'snwx_ds_docdense.csv'
# Load Raw File Here Once (to save time)
ar_svr_csv <- c('marital', 'kids', 'income', 'mass', 'age', 'savings', 'pop', 'educ', 'headinc', 'spouseinc')

if (!exists('df_mass_states')) {
  tm_start_csv <- proc.time()
  print('start loading df_mass_states')

  df_mass_states <- as_tibble(read.csv(paste0(srt_simu_path, st_ds_file_name), header=FALSE)) %>%
    rename_all(~c(ar_svr_csv)) %>%
    mutate(marital_group = factor(marital), kids_group = factor(kids))

  tm_csv <- proc.time() - tm_start_csv
  print(paste0('df_mass_states loaded, tm_csv:', tm_csv))
} else {
  print('already loaded df_mass_states previously')
}

# summarize
# REconTools::ff_summ_percentiles(df_mass_states)

# Generate Marital, Kids, Income, Age, Savings, educ groups
fl_multiple <- 58056
ar_svr_groups <- c('marital_group', 'kids_group', 'ymin_group')
# ar_svr_groups_noage <- c('marital_group', 'kids_group', 'ymin_group')
ar_svr_values <- c('marital', 'kids', 'income', 'age', 'savings', 'pop', 'educ', 'headinc', 'spouseinc')

# Loop Over possibly different Max Ages
ar_it_min_age <- c(18, 22)
ar_it_max_age <- c(64, 69, 79, 100)
# ar_it_min_age <- c(18)
# ar_it_max_age <- c(64)
# it_min_age <- 18
# it_max_age <- 64
# bl_full <- 0
for (it_min_age in ar_it_min_age) {
  for (it_max_age in ar_it_max_age) {
    for (bl_full in c(0, 1)) {

      if (bl_full == 0) {
        ar_agecut = c(17, seq(19, it_max_age, by=5))
        # ar_ycut <- c(0, 0.3797, 0.4747, 0.5696, 0.725, 1.044, 1.441, 2.106, 7, 100)
        ar_ycut <- c(seq(0, 8, by=10000/58056), 100)
        st_suffix <- '_dist_grouped_type_10k'
      }
      if (bl_full == 1) {
        ar_agecut = c(17, seq(19, it_max_age, by=5))
        # ar_ycut <- c(0, 0.3797, 0.4747, 0.5696, 0.725, 1.044, 1.441, 2.106, 7, 100)
        ar_ycut <- c(seq(0, 8, by=20000/58056), 100)
        st_suffix <- '_dist_grouped_type_20k'
      }
      # if (bl_full == 1) {
      #   ar_ycut <- c(seq(0, 7, by=1), 100)
      # }

      # Save Path, save to csv
      srt_csv_path <- paste0(srt_csv_path_root, 'nsw_a', it_min_age, 't', it_max_age, st_suffix, '/')
      dir.create(file.path(srt_csv_path), showWarnings = FALSE, recursive = TRUE)

      # Cut by Age
      df_mass_states_working <- df_mass_states %>%
        filter(age <= it_max_age) %>%
        filter(age >= it_min_age)

      # Income Groups
      df_mass_states_working <- df_mass_states_working %>%
        mutate(ymin_group = (cut(income, ar_ycut))) %>%
        group_by(marital_group, kids_group, ymin_group) %>%
        summarize(marital = sum(marital*mass)/sum(mass),
                  kids = sum(kids*mass)/sum(mass),
                  income = sum(income*mass)/sum(mass),
                  age = sum(age*mass)/sum(mass),
                  savings = sum(savings*mass)/sum(mass),
                  pop = sum(pop*mass)/sum(mass),
                  educ = sum(educ*mass)/sum(mass),
                  headinc = sum(headinc*mass)/sum(mass),
                  spouseinc = sum(spouseinc*mass)/sum(mass),
                  mass = sum(mass)) %>%
        ungroup()

      # # Age Groups
      # df_mass_states_working <- df_mass_states_working %>%
      #   mutate(age_group = (cut(age, ar_agecut))) %>%
      #   group_by(marital_group, kids_group, age_group, ymin_group) %>%
      #   summarize(marital = sum(marital*mass)/sum(mass),
      #             kids = sum(kids*mass)/sum(mass),
      #             income = sum(income*mass)/sum(mass),
      #             age = sum(age*mass)/sum(mass),
      #             savings = sum(savings*mass)/sum(mass),
      #             pop = sum(pop*mass)/sum(mass),
      #             educ = sum(educ*mass)/sum(mass),
      #             headinc = sum(headinc*mass)/sum(mass),
      #             spouseinc = sum(spouseinc*mass)/sum(mass),
      #             mass = sum(mass)) %>%
      #   ungroup()

      # Analyze results:
      REconTools::ff_summ_percentiles(df_mass_states_working,FALSE)
      print(paste0('csv export folder: ', srt_csv_path))

      # Export
        df_alloc_combine_all_mean <- df_mass_states_working %>%
          ungroup() %>% group_by(!!!syms(ar_svr_groups)) %>%
          summarize(marital = sum(marital*mass)/sum(mass),
                    kids = sum(kids*mass)/sum(mass),
                    income = sum(income*mass)/sum(mass),
                    age = sum(age*mass)/sum(mass),
                    savings = sum(savings*mass)/sum(mass),
                    pop = sum(pop*mass)/sum(mass),
                    educ = sum(educ*mass)/sum(mass),
                    headinc = sum(headinc*mass)/sum(mass),
                    spouseinc = sum(spouseinc*mass)/sum(mass),
                    mass = sum(mass))
      write.csv(df_alloc_combine_all_mean,
                paste0(srt_csv_path, "df_dist_marry_kids_ymin.csv"),
                row.names = TRUE)


      # Export Grouped Average Statistics ------
      ls_svr_groups <- ar_svr_groups
      for (svr_group in ls_svr_groups) {

        # group mean
        df_alloc_combine_group_mean <- df_mass_states_working %>%
          ungroup() %>% group_by(!!sym(svr_group)) %>%
          summarize(marital = sum(marital*mass)/sum(mass),
                    kids = sum(kids*mass)/sum(mass),
                    income = sum(income*mass)/sum(mass),
                    age = sum(age*mass)/sum(mass),
                    savings = sum(savings*mass)/sum(mass),
                    pop = sum(pop*mass)/sum(mass),
                    educ = sum(educ*mass)/sum(mass),
                    headinc = sum(headinc*mass)/sum(mass),
                    spouseinc = sum(spouseinc*mass)/sum(mass),
                    mass = sum(mass))

        # Export
        write.csv(df_alloc_combine_group_mean,
                  paste0(srt_csv_path, "df_dist_",svr_group, ".csv"),
                  row.names = TRUE)

        # All but Group mean
        ls_svr_groups_oneless <- ls_svr_groups[ls_svr_groups != svr_group]
        df_alloc_combine_group_mean_oneless <- df_mass_states_working %>%
          ungroup() %>% group_by(!!!syms(ls_svr_groups_oneless)) %>%
          summarize(marital = sum(marital*mass)/sum(mass),
                    kids = sum(kids*mass)/sum(mass),
                    income = sum(income*mass)/sum(mass),
                    age = sum(age*mass)/sum(mass),
                    savings = sum(savings*mass)/sum(mass),
                    pop = sum(pop*mass)/sum(mass),
                    educ = sum(educ*mass)/sum(mass),
                    headinc = sum(headinc*mass)/sum(mass),
                    spouseinc = sum(spouseinc*mass)/sum(mass),
                    mass = sum(mass))

        # Export
        write.csv(df_alloc_combine_group_mean_oneless,
                  paste0(srt_csv_path, "df_dist_",svr_group,"_without.csv"),
                  row.names = TRUE)

      }

      # # Average Statistics without Age ------
      # ls_svr_groups <- ar_svr_groups_noage
      # for (svr_group in ls_svr_groups) {

      #   # All but Group mean
      #   ls_svr_groups_oneless <- ls_svr_groups[ls_svr_groups != svr_group]
      #   df_alloc_combine_group_mean_oneless <- df_mass_states_working %>%
      #     ungroup() %>% group_by(!!!syms(ls_svr_groups_oneless)) %>%
      #    summarize(marital = sum(marital*mass)/sum(mass),
      #              kids = sum(kids*mass)/sum(mass),
      #              income = sum(income*mass)/sum(mass),
      #              age = sum(age*mass)/sum(mass),
      #              savings = sum(savings*mass)/sum(mass),
      #              pop = sum(pop*mass)/sum(mass),
      #              educ = sum(educ*mass)/sum(mass),
      #              headinc = sum(headinc*mass)/sum(mass),
      #              spouseinc = sum(spouseinc*mass)/sum(mass),
      #              mass = sum(mass))

      #   # Export
      #   write.csv(df_alloc_combine_group_mean_oneless,
      #             paste0(srt_csv_path, "df_dist_",svr_group,"_without_noage.csv"),
      #             row.names = TRUE)
      # }
    }
  }
}
