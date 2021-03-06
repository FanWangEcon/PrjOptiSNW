# # Marginal Propensity to Consume Check Tables
# Using the full output of all checks, summarize MPCs for different types of households by child count, marital status, and income bins at 20k intervals. Summarize MPC for each check.

library(tidyverse)
library(REconTools)
# library(PrjOptiAlloc)
library(forcats)

# it_no_cores <- detectCores(logical = TRUE)
# cl <- makeCluster(6)
# registerDoParallel(cl)


# list of files to solve/loop generate MPC/APC over
it_which_list <- 4
if (it_which_list == 1) {
  ls_st_file_suffix <- c('gm300', 'gm105', 'bt95', 'bt99', 'rr2', 'rr1')
} else if (it_which_list == 2) {
  ls_st_file_suffix <- c('bt50', 'bt60', 'bt70', 'bt80', 'bt90')
} else if (it_which_list == 3) {
  ls_st_file_suffix <- c('b1_manna', 'b1_xi0_tax', 'b0_xi0p25_manna')
} else if (it_which_list == 4) {
  # common string components
  st_base_name <- 'snwx_v_planner_docdense'
  st_base_middle <- 'b1_xi0_manna_88'
  # numeric values to loop over
  # ar_st_beta_val <- c('bt50', 'bt60', 'bt70', 'bt80', 'bt90', 'bt95', 'bt97', 'bt97')
  ar_st_beta_val <- c('bt99')
  ar_st_edu_type <- c('e1lm2', 'e2hm2')
  # initialize string list
  ls_st_file_suffix <- vector(mode = "list", length = length(ar_st_beta_val)*length(ar_st_edu_type))

  # generate list
  it_ctr = 0
  for (st_beta_val in ar_st_beta_val) {
    for (st_edu_type in ar_st_edu_type) {
      it_ctr = it_ctr + 1
      # snm_file_name <- 'snwx_v_planner_docdense_e2hm2_b1_xi0_manna_88_bt90'
      snm_file_name <- paste(st_base_name, st_edu_type, st_base_middle, st_beta_val, sep ='_')
      ls_st_file_suffix[it_ctr] <- snm_file_name
    }
  }
}

# Solve iteratively
for (st_which_solu in ls_st_file_suffix) {
# foreach (st_file_suffix=ls_st_file_suffix) %dopar% {
  rm(df_plan_v_tilde_full)

  # Files:
  # source('fs_opti_support.R')
  # st_which_solu <- 'b1_manna'
  # st_which_solu <- paste0('b1_xi0_manna_88_', st_file_suffix)
  ls_output <- fs_opti_support(st_which_solu)
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
  if (!exists('df_plan_v_tilde_full')) {
    print('start loading df_plan_v_tilde_full')
    df_plan_v_tilde_full <- as_tibble(read.csv(paste0(srt_simu_path, snm_simu_csv_withspouse_shock), header=FALSE)) %>%
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
    it_max_checks <- 245
    ar_agecut <- c(it_min_age-1, it_max_age)
    ar_svr_groups <- c('marital', 'kids', 'age_group', 'ymin_group')
    svr_v_value <- 'vtilde'
    svr_c_value <- 'ctilde'
    svr_group_id <- 'group_id'
    svr_checks <- 'checks'
    svr_mass <- 'mass'

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
    ar_ycut_usd <- c(0, 20000, 40000, 60000, 80000, 100000, 100000000)
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
      snm_save_csv = paste0(snm_save_csv, '_bykidsmarital20k_allchecks_', st_img_suf_age_ybin, '.csv')


      if (MPC_type == 1 | MPC_type == 2) {
        # Average MPCs by group
        df_MPC_results <- df_il_U_select %>%
          select(marital, kids, ymin_group, checks, mass, MPC) %>%
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
          select(marital, kids, ymin_group, checks, mass, APC) %>%
          pivot_wider(names_from = checks,
                      values_from = APC)
      }

      # CSV Save
      write.csv(df_MPC_results,
                paste0(srt_csv_path_root, snm_save_csv),
                row.names = TRUE)

      # view(df_MPC_results)

      # summarize
      # REconTools::ff_summ_percentiles(df_il_U_select, bl_statsasrows = FALSE)
    }
  }
}
# stopCluster(cl)
