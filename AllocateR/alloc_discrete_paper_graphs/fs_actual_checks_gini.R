# NA source?
# D:/Dropbox (UH-ECON)/PrjNygaardSorensenWang/Results202103/capi_rh1/20210225_bidenchk_tiny_b1_xi0_manna_168/csv/b1_a64_14ca14ck_18t64
# df_alloc_all_optimalg47_14ca14ck.csv

for (it_file in c(1, 2)) {

  if (it_file == 1) {
    srt_results_root <- 'D:/Dropbox (UH-ECON)/PrjNygaardSorensenWang/Results202103/capi_rh1/20210225_bidenchk_tiny_b1_xi0_manna_168/csv/b1_a64_14ca14ck_18t64'
    spn_csv_full_path <- file.path(srt_results_root, 'df_alloc_all_optimalg47_14ca14ck.csv')
    spn_csv_full_path <- file.path(srt_results_root, 'df_alloc_all_feasible_14ca14ck.csv')
    df_alloc_cur <- as_tibble(read.csv(spn_csv_full_path))
  }
  if (it_file == 2) {
    srt_results_root <- 'D:/Dropbox (UH-ECON)/PrjNygaardSorensenWang/Results202103/capi_rh1/20210225_bidenchk_a65zh266zs5_b1_xi0_manna_168/csv/b1_a64_14ca14ck_18t64'
    spn_csv_full_path <- file.path(srt_results_root, 'df_alloc_all_optimalg47_14ca14ck.csv')
    spn_csv_full_path <- file.path(srt_results_root, 'df_alloc_all_feasible_14ca14ck.csv')
    df_alloc_cur <- as_tibble(read.csv(spn_csv_full_path))
  }

  ar_mass_wgt <- df_alloc_cur$mass
  marital <- df_alloc_cur$marital
  kids <- df_alloc_cur$kids
  ar_household_size <- (kids + marital + 1)
  ar_mass_wgt_hh <- ar_mass_wgt*ar_household_size
  ar_mass_wgt_hh <- ar_mass_wgt_hh/sum(ar_mass_wgt_hh)

  ar_checks_actual <- df_alloc_cur$actual
  ar_mass_wgt <- ar_mass_wgt_hh/sum(ar_mass_wgt_hh)
  fl_gini <- ff_dist_gini_random_var(ar_checks_actual, ar_mass_wgt)
  print(paste0('fl_gini=', fl_gini))

  c_no_allocation <- df_alloc_cur$optimal_c_1st_c_no_allocation
  c_no_allocation_mean <- sum(ar_mass_wgt_hh*log(c_no_allocation))
  sd_log_c <- sqrt(sum(ar_mass_wgt_hh*(log(c_no_allocation)-c_no_allocation_mean)^2))
  print(paste0('sd_log_c=', sd_log_c))

}
