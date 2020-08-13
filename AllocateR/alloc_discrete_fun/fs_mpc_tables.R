# # Marginal Propensity to Consume Check Tables

# Files:
# source('fs_opti_support.R')
ls_output <- fs_opti_support()
st_b0b1 <- ls_output$st_b0b1
st_file_type_withspouse_shock <- ls_output$st_file_type_withspouse_shock
st_file_type_withspouse_shock <- ls_output$st_file_type_withspouse_shock
snm_simu_csv_withspouse_shock <- ls_output$snm_simu_csv_withspouse_shock
srt_simu_path <- ls_output$srt_simu_path
srt_csv_path_root <- ls_output$srt_csv_path_root
ar_rho <- ls_output$ar_rho
### Variable Names and Paths
## Common Parameters Across

# Max Phase Out given 1200*2 + 500*4 = 4400
fl_max_phaseout = 238000
it_bin_dollar_before_phaseout = 500
# Dollar Per Check
fl_percheck_dollar = 100
# Meaning of Ymin Ymax simulated interval of 1
fl_multiple = 58056
# Number of Max Checks
it_max_checks_1st = 44
it_max_checks_2nd = 88
# Number of Tax Paying Households
fl_tax_hh = 128580000
# Number of Income Groups to Use: use 25 for 10,000 = 1
# Age Conditions

# Load Raw File Here Once (to save time)
ar_svr_csv <- c('age', 'marital', 'kids', 'checks', 'ymin', 'mass', 'survive', 'vtilde', 'ctilde')
df_plan_v_tilde_full <- as_tibble(read.csv(paste0(srt_simu_path, snm_simu_csv_withspouse_shock), header=FALSE)) %>%
  rename_all(~c(ar_svr_csv))

ar_it_max_age <- c(64, 69, 79, 100)
for (it_max_age in ar_it_max_age) {

  for (bl_full in c(0)) {

    if (bl_full == 1) {
      st_suffix <- '_mpc_all_y_ages'
      ar_agecut = seq(17, it_max_age, by=1)
      ar_ycut <- NULL
    }

    if (bl_full == 0) {
      ar_agecut = c(17, seq(19, it_max_age, by=5))
      ar_ycut <- c(0, 0.3797, 0.4747, 0.5696, 0.725, 1.044, 1.441, 2.106, 7)
      st_suffix <- '_mpc_grouped_y_ages'
    }

    srt_csv_path <- paste0(srt_csv_path_root, st_b0b1, '_a', it_max_age, st_suffix, '/')
    dir.create(file.path(srt_csv_path), showWarnings = FALSE, recursive = TRUE)

    # Parameters -------
    # it_max_age = 79
    it_min_age = 18

    # Image Save Suffix
    st_img_suf_age_ybin <- paste0(it_min_age, 't', it_max_age)
    print(paste0('st_img_suf_age_ybin:', st_img_suf_age_ybin, ' started'))

    # File Path
    # srt_simu_path <- 'C:/Users/fan/Documents/Dropbox (UH-ECON)/PrjNygaardSorensenWang/Output/'
    # Column Names
    # Variables That Identify Individual Types
    ar_svr_groups <- c('marital', 'kids', 'age_group', 'ymin_group')
    ar_svr_groups_noage <- c('marital', 'kids', 'ymin_group')
    ar_svr_groups_stats <- c('mass', 'survive')
    # Number of Checks and Planner Value
    svr_checks <- 'checks'
    svr_v_value <- 'vtilde'
    svr_c_value <- 'ctilde'
    svr_mass <- 'mass'

    # Call Processor Function ----
    ### Process
    # Load dataset and average by
    ls_prc_outputs_core <- ffp_snw_process_inputs_core(
      srt_simu_path = srt_simu_path,
      snm_simu_csv = snm_simu_csv_withspouse_shock,
      df_plan_v_tilde_full = df_plan_v_tilde_full,
      fl_max_phaseout = fl_max_phaseout,
      it_bin_dollar_before_phaseout = it_bin_dollar_before_phaseout,
      fl_percheck_dollar = fl_percheck_dollar,
      fl_multiple = fl_multiple,
      it_max_checks = it_max_checks_1st,
      fl_tax_hh = fl_tax_hh,
      it_max_age = it_max_age,
      it_min_age = it_min_age,
      ar_ycut = ar_ycut,
      ar_agecut = ar_agecut,
      ar_svr_csv = ar_svr_csv,
      ar_svr_groups = ar_svr_groups,
      ar_svr_groups_stats = ar_svr_groups_stats,
      svr_checks = svr_checks,
      svr_v_value = svr_v_value,
      svr_c_value = svr_c_value,
      svr_mass = svr_mass,
      ar_rho = ar_rho,
      bl_non_inc_adjust = TRUE,
      bl_print = FALSE,
      bl_print_verbose = FALSE)

    print(paste0('core data loaded, ', st_img_suf_age_ybin))

    ### Outputs
    # Mass Numbers
    mass_sum=ls_prc_outputs_core$mass_sum
    mass_sum_covid_vox_actual=ls_prc_outputs_core$mass_sum_covid_vox_actual
    # Dataframes
    df_input_il_noninc=ls_prc_outputs_core$df_input_il_noninc
    df_id=ls_prc_outputs_core$df_id
    # Merge input_il and df id and include only the first check
    df_main_mpc <- df_input_il_noninc %>% filter(D_il == 1) %>%
      left_join(df_id, by = 'id_i') %>%
      ungroup %>% mutate(mpc = c_alpha_il/(fl_percheck_dollar/fl_multiple))

    # Analyze results:
    REconTools::ff_summ_percentiles(df_main_mpc,FALSE)
    print(paste0('csv export folder: ', srt_csv_path))

    # Export Grouped Average Statistics ------
    ls_svr_groups <- ar_svr_groups
    for (svr_group in ls_svr_groups) {

      # group mean
      df_alloc_combine_group_mean <- df_main_mpc %>%
        ungroup() %>% group_by(!!sym(svr_group)) %>%
        summarize(mpc_1stcheck = sum(mpc*mass)/sum(mass),
                  c_1stcheck = sum(c_alpha_il*mass)/sum(mass),
                  mass_sum = sum(mass))

      # Export
      write.csv(df_alloc_combine_group_mean,
                paste0(srt_csv_path, "df_mpc_",svr_group, ".csv"),
                row.names = TRUE)

      # All but Group mean
      ls_svr_groups_oneless <- ls_svr_groups[ls_svr_groups != svr_group]
      df_alloc_combine_group_mean_oneless <- df_main_mpc %>%
        ungroup() %>% group_by(!!!syms(ls_svr_groups_oneless)) %>%
        summarize(mpc_1stcheck = sum(mpc*mass)/sum(mass),
                  c_1stcheck = sum(c_alpha_il*mass)/sum(mass),
                  mass_sum = sum(mass))

      # Export
      write.csv(df_alloc_combine_group_mean_oneless,
                paste0(srt_csv_path, "df_mpc_",svr_group,"_without.csv"),
                row.names = TRUE)

    }

    # Average Statistics without Age ------
    ls_svr_groups <- ar_svr_groups_noage
    for (svr_group in ls_svr_groups) {

      # All but Group mean
      ls_svr_groups_oneless <- ls_svr_groups[ls_svr_groups != svr_group]
      df_alloc_combine_group_mean_oneless <- df_main_mpc %>%
        ungroup() %>% group_by(!!!syms(ls_svr_groups_oneless)) %>%
        summarize(mpc_1stcheck = sum(mpc*mass)/sum(mass),
                  c_1stcheck = sum(c_alpha_il*mass)/sum(mass),
                  mass_sum = sum(mass))
      # Export
      write.csv(df_alloc_combine_group_mean_oneless,
                paste0(srt_csv_path, "df_mpc_",svr_group,"_without_noage.csv"),
                row.names = TRUE)
    }
  }
}
