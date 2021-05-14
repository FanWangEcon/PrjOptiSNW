# fs_mpc_tables_increments_202103 is almost the same as fs_mpc_tables_increments_202102, except now
# rather than having a fixed 20k bin, finer bins can be specified
# can specify bin width, as well as final bin point.
# This helps with drawing a scatter-plot that illustrats the relationship between
# A and alpha for the allocation problem.

library(tidyverse)
library(REconTools)
# library(PrjOptiAlloc)

library(forcats)
library(foreach)
library(doParallel)

# it_no_cores <- detectCores(logical = TRUE)
# cl <- makeCluster(5)
# registerDoParallel(cl)

# Number of ways to cut income bins
ls_it_income_cuts <- c(1,2,3,4)
ls_it_income_cuts <- c(1,3)

# Types of allocation files to consider
ls_st_file_suffix_trumpchk <-
  c('snwx_trumpchk_moredense_a65zh266zs5_b1_xi0_manna_168_bt95',
    'snwx_trumpchk_moredense_a65zh266zs5_b1_xi0_manna_168_bt60',
    'snwx_trumpchk_moredense_a65zh266zs5_b1_xi0_manna_168_married',
    'snwx_trumpchk_moredense_a65zh266zs5_b1_xi0_manna_168_unmarried',
    'snwx_trumpchk_moredense_a65zh266zs5_b1_xi0_manna_168')
ls_st_file_suffix_trumpchk <- rev(ls_st_file_suffix_trumpchk)

ls_st_file_suffix_bidenchk <-
  c('snwx_bidenchk_moredense_a65zh266zs5_b1_xi0_manna_168_bt95',
    'snwx_bidenchk_moredense_a65zh266zs5_b1_xi0_manna_168_bt60',
    'snwx_bidenchk_moredense_a65zh266zs5_b1_xi0_manna_168_married',
    'snwx_bidenchk_moredense_a65zh266zs5_b1_xi0_manna_168_unmarried',
    'snwx_bidenchk_moredense_a65zh266zs5_b1_xi0_manna_168')
ls_st_file_suffix_bidenchk <- rev(ls_st_file_suffix_bidenchk)

ls_st_file_suffix_bidenchk_mixturealter <-
  c('snwx_bidenchk_moredense_a65zh266zs5_b1_xi0_manna_168_3o6',
    'snwx_bidenchk_moredense_a65zh266zs5_b1_xi0_manna_168_4o6',
    'snwx_bidenchk_moredense_a65zh266zs5_b1_xi0_manna_168_8o6',
    'snwx_bidenchk_moredense_a65zh266zs5_b1_xi0_manna_168_9o6')
ls_st_file_suffix_bidenchk_mixturealter <- rev(ls_st_file_suffix_bidenchk_mixturealter)

ls_st_file_suffix_bchkbnoui <-
  c('snwx_bchknoui_moredense_a65zh266zs5_b1_xi0_manna_168')

ls_st_file_suffix_bchklock <-
  c('snwx_bchklock_moredense_a65zh266zs5_b1_xi0_manna_168_bt95',
    'snwx_bchklock_moredense_a65zh266zs5_b1_xi0_manna_168_bt60',
    'snwx_bchklock_moredense_a65zh266zs5_b1_xi0_manna_168_married',
    'snwx_bchklock_moredense_a65zh266zs5_b1_xi0_manna_168_unmarried',
    'snwx_bchklock_moredense_a65zh266zs5_b1_xi0_manna_168')
ls_st_file_suffix_bchklock <- rev(ls_st_file_suffix_bchklock)

ls_st_file_suffix_bidenchk_betaedu <-
  c('snwx_bidenchk_docdense_e2hm2_b1_xi0_manna_88_bt99',
    'snwx_bidenchk_docdense_e2hm2_b1_xi0_manna_88_bt97',
    'snwx_bidenchk_docdense_e2hm2_b1_xi0_manna_88_bt95',
    'snwx_bidenchk_docdense_e2hm2_b1_xi0_manna_88_bt90',
    'snwx_bidenchk_docdense_e2hm2_b1_xi0_manna_88_bt80',
    'snwx_bidenchk_docdense_e2hm2_b1_xi0_manna_88_bt70',
    'snwx_bidenchk_docdense_e2hm2_b1_xi0_manna_88_bt60',
    'snwx_bidenchk_docdense_e2hm2_b1_xi0_manna_88_bt50',
    'snwx_bidenchk_docdense_e1lm2_b1_xi0_manna_88_bt99',
    'snwx_bidenchk_docdense_e1lm2_b1_xi0_manna_88_bt97',
    'snwx_bidenchk_docdense_e1lm2_b1_xi0_manna_88_bt95',
    'snwx_bidenchk_docdense_e1lm2_b1_xi0_manna_88_bt90',
    'snwx_bidenchk_docdense_e1lm2_b1_xi0_manna_88_bt80',
    'snwx_bidenchk_docdense_e1lm2_b1_xi0_manna_88_bt70',
    'snwx_bidenchk_docdense_e1lm2_b1_xi0_manna_88_bt60',
    'snwx_bidenchk_docdense_e1lm2_b1_xi0_manna_88_bt50')

ls_st_file_suffix_trumpchk_betaedu <-
  c('snwx_trumpchk_docdense_e2hm2_b1_xi0_manna_88_bt99',
    'snwx_trumpchk_docdense_e2hm2_b1_xi0_manna_88_bt97',
    'snwx_trumpchk_docdense_e2hm2_b1_xi0_manna_88_bt95',
    'snwx_trumpchk_docdense_e2hm2_b1_xi0_manna_88_bt90',
    'snwx_trumpchk_docdense_e2hm2_b1_xi0_manna_88_bt80',
    'snwx_trumpchk_docdense_e2hm2_b1_xi0_manna_88_bt70',
    'snwx_trumpchk_docdense_e2hm2_b1_xi0_manna_88_bt60',
    'snwx_trumpchk_docdense_e2hm2_b1_xi0_manna_88_bt50',
    'snwx_trumpchk_docdense_e1lm2_b1_xi0_manna_88_bt99',
    'snwx_trumpchk_docdense_e1lm2_b1_xi0_manna_88_bt97',
    'snwx_trumpchk_docdense_e1lm2_b1_xi0_manna_88_bt95',
    'snwx_trumpchk_docdense_e1lm2_b1_xi0_manna_88_bt90',
    'snwx_trumpchk_docdense_e1lm2_b1_xi0_manna_88_bt80',
    'snwx_trumpchk_docdense_e1lm2_b1_xi0_manna_88_bt70',
    'snwx_trumpchk_docdense_e1lm2_b1_xi0_manna_88_bt60',
    'snwx_trumpchk_docdense_e1lm2_b1_xi0_manna_88_bt50')


ls_st_file_suffix_test <- c('snwx_bidenchk_tiny_b1_xi0_manna_168')

bl_per_capita <- TRUE
fl_rho <- 1

# list to run
ls_st_file_suffix <- c(ls_st_file_suffix_bidenchk,
                       ls_st_file_suffix_bidenchk_mixturealter,
                       ls_st_file_suffix_trumpchk,
                       ls_st_file_suffix_bchklock,
                       ls_st_file_suffix_bchkbnoui)

ls_st_file_suffix <- c(ls_st_file_suffix_bchkbnoui)

for (it_income_cuts in ls_it_income_cuts) {
  if (it_income_cuts == 1) {
    # 20k interval, between 0 and 100k and 100 million
    # generate income cut-offs
    fl_bin_start <- 0
    # width equal to 20,000
    fl_bin_width <- 2e4
    # final point is 100 million
    fl_bin_final_end <- 1e8
    # final segment starting point is 100,000 dollars
    fl_bin_final_start <- 1e5
    # File name
    st_snm_suffix_save <- '_0t100k20kbin'

  } else if (it_income_cuts == 2) {
    # 10k interval, between 0 and 150k and 100 million
    # generate income cut-offs
    fl_bin_start <- 0
    # width equal to 20,000
    fl_bin_width <- 1e4
    # final point is 100 million
    fl_bin_final_end <- 1e8
    # final segment starting point is 100,000 dollars
    fl_bin_final_start <- 1.5e5
    # File name
    st_snm_suffix_save <- '_0t150k10kbin'

  } else if (it_income_cuts == 3) {
    # 10k interval, between 0 and 200k and 100 million
    # generate income cut-offs
    fl_bin_start <- 0
    # width equal to 20,000
    fl_bin_width <- 1e4
    # final point is 100 million
    fl_bin_final_end <- 1e8
    # final segment starting point is 100,000 dollars
    fl_bin_final_start <- 2e5
    # File name
    st_snm_suffix_save <- '_0t200k10kbin'

  } else if (it_income_cuts == 4) {
    # 5k interval, between 0 and 200k and 100 million
    # generate income cut-offs
    fl_bin_start <- 0
    # width equal to 20,000
    fl_bin_width <- 5e3
    # final point is 100 million
    fl_bin_final_end <- 1e8
    # final segment starting point is 100,000 dollars
    fl_bin_final_start <- 2e5
    # File name
    st_snm_suffix_save <- '_0t200k5kbin'

  }
  ar_income_bins <- c(seq(fl_bin_start, fl_bin_final_start, by=fl_bin_width),
                      fl_bin_final_end)
  # Solve iteratively
  for (st_which_solu in ls_st_file_suffix) {
  # foreach (st_which_solu=ls_st_file_suffix) %dopar% {
    rm(df_plan_v_tilde_full)

    # Files:
    # source('fs_opti_support.R')
    # st_which_solu <- 'b1_manna'
    # st_which_solu <- paste0('b1_xi0_manna_88_', st_file_suffix)
    ls_output <- fs_opti_support_202103(st_which_solu, bl_per_capita=bl_per_capita, fl_rho=fl_rho)
    st_b0b1 <- ls_output$st_b0b1
    st_file_type_withspouse_shock <- ls_output$st_file_type_withspouse_shock
    st_file_type_withspouse_shock <- ls_output$st_file_type_withspouse_shock
    snm_simu_csv_withspouse_shock <- ls_output$snm_simu_csv_withspouse_shock

    srt_simu_path <- ls_output$srt_simu_path
    srt_csv_path_root <- ls_output$srt_csv_path_root
    srt_imgcsv_mpcapc_root <- ls_output$srt_imgcsv_mpcapc_root
    ar_rho <- ls_output$ar_rho
    ### Variable Names and Paths
    ## Common Parameters Across

    # Max Phase Out given 1200*2 + 500*4 = 4400
    fl_max_phaseout = 200000
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
    if (!exists('df_plan_v_tilde_full')) {
      print('start loading df_plan_v_tilde_full')
      df_plan_v_tilde_full <- as_tibble(
        read.csv(file.path(srt_simu_path, snm_simu_csv_withspouse_shock), header=FALSE)) %>%
        rename_all(~c(ar_svr_csv))
    } else {
      print('already loaded df_plan_v_tilde_full previously')
    }

    for (it_age_type in c(1, 2)) {
      if (it_age_type == 1) {
        it_max_age = 64
        it_min_age = 18
      }
      if (it_age_type == 2) {
        it_max_age = 99
        it_min_age = 18
      }
      st_img_suf_age_ybin <- paste0(it_min_age, 't', it_max_age)

      # Various Parameters
      it_max_checks <- 168
      ar_agecut <- c(it_min_age-1, it_max_age)
      ar_svr_groups <- c('marital', 'kids', 'age_group', 'ymin_group')
      svr_v_value <- 'vtilde'
      svr_c_value <- 'ctilde'
      svr_group_id <- 'group_id'
      svr_checks <- 'checks'
      svr_mass <- 'mass'
      svr_c_zerochk <- 'c_avg_chk0_usd'

      # Select Subset
      df_plan_v_tilde <- df_plan_v_tilde_full %>%
        filter(vtilde != 0) %>%
        filter(checks <= it_max_checks) %>%
        filter(age <= it_max_age) %>%
        filter(age >= it_min_age)

      # Average over all ages, since table does not show age as group
      df_plan_v_tilde_yjm <- df_plan_v_tilde %>%
        mutate(age_group = cut(age, ar_agecut)) %>%
        group_by(marital, kids, checks, ymin, age_group) %>%
        summarize(vtilde = sum(vtilde*mass)/sum(mass),
                  ctilde = sum(ctilde*mass)/sum(mass),
                  mass = sum(mass),
                  survive = mean(survive)) %>%
        ungroup()

      # Cutting Income Grid, want 0 to 20k, every 20k, than over 100k groups
      # ar_ycut_usd <- c(0, 20000, 40000, 60000, 80000, 100000, 100000000)
      ar_ycut_usd <- ar_income_bins
      ar_ycut <- ar_ycut_usd/fl_multiple
      it_inc_groups = length(ar_ycut)
      fl_inc_gap = (ar_ycut[3]-ar_ycut[2])*fl_multiple
      fl_inc_min = min(df_plan_v_tilde_yjm %>% pull(ymin))*fl_multiple
      subtitle = paste0('1 unit along x-axis = $', round(fl_inc_gap),
                        ', x-axis min = $', round(fl_inc_min),
                        ', x-axis final group >= $', round(ar_ycut[length(ar_ycut)-1]*fl_multiple))
      print(subtitle)

      # Cut income Groups
      df_plan_v_tilde_ygrpjm <- df_plan_v_tilde_yjm %>%
        mutate(ymin_group = (cut(ymin, ar_ycut))) %>%
        group_by(marital, kids, checks, age_group, ymin_group) %>%
        summarize(vtilde = sum(vtilde*mass)/sum(mass),
                  ctilde = sum(ctilde*mass)/sum(mass),
                  mass = sum(mass),
                  survive = mean(survive)) %>%
        ungroup()

      # Grouping
      ls_svr_group_vars <- ar_svr_groups
      df_il <- df_plan_v_tilde_ygrpjm %>%
        arrange(!!!syms(ls_svr_group_vars)) %>%
        group_by(!!!syms(ls_svr_group_vars)) %>%
        mutate(!!sym(svr_group_id) := (row_number()==1)*1) %>%
        ungroup() %>%
        rowid_to_column(var = "id") %>%
        mutate(!!sym(svr_group_id) := cumsum(!!sym(svr_group_id)))

      # Cleaning
      df_il <- df_il %>%
        rename(id_i = !!sym(svr_group_id)) %>%
        mutate(id_il = row_number()) %>%
        arrange(id_i, svr_checks) %>% group_by(id_i) %>%
        mutate(D_max_i = max(!!sym(svr_checks))) %>%
        mutate(D_il = !!sym(svr_checks)) %>%
        mutate(beta_i = 1/n())

      # Consumption without checks: store in USD
      df_il_checkzero <- df_il %>%
        ungroup() %>%
        filter(checks == 0) %>%
        select(id_i, ctilde) %>%
        mutate(ctilde = ctilde*fl_multiple) %>%
        rename(!!sym(svr_c_zerochk) := ctilde)

      # Marginal Consumption Changes
      df_il_U <- df_il %>%
        mutate(c_alpha_il = lead(!!sym(svr_c_value)) - (!!sym(svr_c_value)),
               v_alpha_il = lead(!!sym(svr_v_value)) - (!!sym(svr_v_value))) %>%
        rename(c_A_il = !!sym(svr_c_value),
               v_A_il = !!sym(svr_v_value)) %>%
        ungroup()

      # Drop max check
      df_il_U <- df_il_U %>%
        filter(D_il != max(df_il$D_il)) %>%
        mutate(D_il = D_il + 1)

      # Smooth out MPC
      # Function copied from PrjOptiAlloc\R\ffp_snw_welfarechecks_input.R-L1099
      ffi_alpha_non_increasing_adj_flatten <- function(ar_alpha, fl_min_inc_bd = 1e-20) {
        # Flattening, flatten so that the marginal value of the next check must be
        # lower than the marginal value of the check prior
        ar_cur <- ar_alpha
        ar_cur_adj <- rep(NA, length(ar_cur))
        ar_cur_adj[1] <- ar_cur[1]
        for (it_ctr in 2:length(ar_cur)) {
          fl_cur_val <- ar_cur[it_ctr]
          fl_last_val <- ar_cur_adj[it_ctr-1]
          if (fl_cur_val > fl_last_val) {
            ar_cur_adj[it_ctr] = fl_last_val
          } else {
            ar_cur_adj[it_ctr] = fl_cur_val
          }
        }
        # return
        return(ar_cur_adj)
      }

      # Run Smoothing Function
      df_il_U <- df_il_U %>%
        arrange((D_il)) %>%
        group_by(id_i) %>%
        do(c_alpha_il_noninc = ffi_alpha_non_increasing_adj_flatten(.$c_alpha_il)) %>%
        unnest(c(c_alpha_il_noninc)) %>%
        group_by(id_i) %>%
        mutate(D_il = row_number()) %>%
        left_join(df_il_U, by=(c('id_i'='id_i', 'D_il'='D_il')))

      # MPC calculation
      df_il_U <- df_il_U %>%
        mutate(MPC_smooth = (c_alpha_il_noninc/(fl_percheck_dollar/fl_multiple))) %>%
        mutate(MPC_raw = (c_alpha_il/(fl_percheck_dollar/fl_multiple)))

      for (MPC_type in c(1,2,3,4)) {

        snm_save_csv = ""

        if (MPC_type == 1 | MPC_type == 3) {
          # Select subset
          df_il_U_select <- df_il_U %>% select(id_i, marital, kids, ymin_group, checks, mass, MPC_smooth) %>%
            mutate(MPC=MPC_smooth)
          if (MPC_type == 1) {
            snm_save_csv <- 'mpc_smooth'
          }
          if (MPC_type == 3) {
            snm_save_csv <- 'apc_smooth'
          }
        }
        if (MPC_type == 2 | MPC_type == 4) {
          # Select subset
          df_il_U_select <- df_il_U %>% select(id_i, marital, kids, ymin_group, checks, mass, MPC_raw) %>%
            mutate(MPC=MPC_raw)
          if (MPC_type == 2) {
            snm_save_csv <- 'mpc_raw'
          }
          if (MPC_type == 4) {
            snm_save_csv <- 'apc_raw'
          }
        }
        snm_save_csv = paste0(snm_save_csv,
                              '_bykidsmarital20k_allchecks_',
                              st_img_suf_age_ybin,
                              st_snm_suffix_save, '.csv')

        if (MPC_type == 1 | MPC_type == 2) {
          # Average MPCs by group
          df_MPC_results <- df_il_U_select %>%
            select(id_i, marital, kids, ymin_group, checks, mass, MPC) %>%
            pivot_wider(names_from = checks,
                        values_from = MPC)
        }
        if (MPC_type == 3 | MPC_type == 4) {
          # Average MPCs by group
          df_MPC_results <- df_il_U_select %>%
            select(marital, kids, ymin_group, checks, mass, MPC) %>%
            arrange(marital, kids, ymin_group, checks) %>%
            group_by(marital, kids, ymin_group, isna = is.na(MPC)) %>%
            mutate(APC = ifelse(isna, NA, cummean(MPC))) %>%
            ungroup() %>%
            select(id_i, marital, kids, ymin_group, checks, mass, APC) %>%
            pivot_wider(names_from = checks,
                        values_from = APC)
        }

        # Merge results with additional column to show consumption at zero
        df_MPC_results <- df_MPC_results %>% ungroup() %>%
          left_join(df_il_checkzero, by="id_i") %>%
          select(marital, kids, ymin_group, mass, !!sym(svr_c_zerochk), everything()) %>%
          select(-id_i)

        # CSV Save
        write.csv(df_MPC_results,
                  paste0(srt_imgcsv_mpcapc_root, snm_save_csv),
                  row.names = TRUE)

        # view(df_MPC_results)

        # summarize
        # REconTools::ff_summ_percentiles(df_il_U_select, bl_statsasrows = FALSE)
      }
    }
  }
}
# stopCluster(cl)
