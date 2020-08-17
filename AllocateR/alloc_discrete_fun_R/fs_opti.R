# This is a R file, not RMD file, the loops through allocation problems to generate results.
library(tidyverse)
library(REconTools)
# library(PrjOptiAlloc)
library(forcats)

# library(foreach)
# library(doParallel)

# File Names, Paths Etc -------
ls_output <- fs_opti_support()
st_file_type_withspouse_shock <- ls_output$st_file_type_withspouse_shock
snm_simu_csv_withspouse_shock <- ls_output$snm_simu_csv_withspouse_shock
srt_simu_path <- ls_output$srt_simu_path
bl_save_img <- ls_output$bl_save_img
spt_img_save <- ls_output$spt_img_save
srt_csv_path <- ls_output$srt_csv_path
ar_rho <- ls_output$ar_rho
it_max_age <- ls_output$it_max_age

# Graph Title -----------
slb_add_title_round1 = 'Round One Policy, '
slb_add_title_round2 = 'Round Two Policy, '
slb_subtitle_stack_2nd = 'optimal2nd_total = 1st actual + 2nd round optimal (given 1st actual) '
slb_subtitle_joint_1st2nd = 'optimal 1st and 2nd rounds (given 1st actual) '

# Image Sizes ----------------
it_img_width=300
it_img_height=180
st_img_units='mm'
it_img_res=300
it_img_pointsize=3

# Locally set parameters -------------------
# Max Phase Out given 1200*2 + 500*4 = 4400
fl_max_phaseout = 238000
it_bin_dollar_before_phaseout = 2500
# Dollar Per Check
fl_percheck_dollar = 100
# Meaning of Ymin Ymax simulated interval of 1
fl_multiple = 58056
# Number of Tax Paying Households
fl_tax_hh = 128580000
# Number of Income Groups to Use: use 25 for 10,000 = 1
# Age Conditions
it_age_bins = 1

# Variable Names -----------------------
# Variables That Identify Individual Types
ar_svr_groups <- c('marital', 'kids', 'age_group', 'ymin_group')
ar_svr_groups_noage <- c('marital', 'kids', 'ymin_group')
ar_svr_groups_stats <- c('mass', 'survive')
# Number of Checks and Planner Value
svr_checks <- 'checks'
svr_v_value <- 'vtilde'
svr_c_value <- 'ctilde'
svr_mass <- 'mass'

# Parameters that vary ------------------
# Age Limits

# Load File ---------------------
ar_svr_csv <- c('age', 'marital', 'kids', 'checks', 'ymin', 'mass', 'survive', 'vtilde', 'ctilde')
if (!exists('df_plan_v_tilde_full')) {
  tm_start_csv <- proc.time()
  print('start loading df_plan_v_tilde_full')
  df_plan_v_tilde_full <- as_tibble(
    read.csv(paste0(srt_simu_path, snm_simu_csv_withspouse_shock), header=FALSE)) %>%
    rename_all(~c(ar_svr_csv))
  # df_plan_v_tilde_full <- vroom(paste0(srt_simu_path, snm_simu_csv_withspouse_shock), col_names=ar_svr_csv)
  tm_csv <- proc.time() - tm_start_csv
  print(paste0('df_plan_v_tilde_full loaded, tm_csv:', tm_csv))
} else {
  print('already loaded df_plan_v_tilde_full previously')
}

# Numberings
# Overall we are doing 14 allocation problems now, 7 for each round, with the second conditional on first round allocation
# Threshold allocation
# 4400 max Feasible no age
# 4400 max Optimal 4 Age group
# 4400 max Optimal 47 age groups
# 20k max Feasible no age
# 20k max Optimal 4 age group
# 20k max Optimal 47 age groups

# parallel
# it_no_cores <- detectCores(logical = TRUE)
# cl <- makeCluster(2)
# registerDoParallel(cl)

# Loop Over Soluiton Files
# it_solu_type <- 1
# it_age_type <- 1
# it_solu_round <- 1

# foreach (it_solu_type=c(1,2,3,5,6,4,7)) %dopar% {
# for (it_solu_type in c(4)) {
# for (it_age_type in c(1)) {
for (it_solu_type in c(1,2,3,5,6,4,7)) {
  for (it_solu_type in c(1,4)) {
    if (it_age_type == 1) {
      it_max_age = 64
      it_min_age = 18
    }
    if (it_age_type == 2) {
      it_max_age = 69
      it_min_age = 18
    }
    if (it_age_type == 3) {
      it_max_age = 74
      it_min_age = 18
    }
    if (it_age_type == 4) {
      it_max_age = 64
      it_min_age = 22
    }
    if (it_age_type == 5) {
      it_max_age = 69
      it_min_age = 22
    }
    if (it_age_type == 6) {
      it_max_age = 74
      it_min_age = 22
    }
    # Image Save Suffix
    st_img_suf_age_ybin <- paste0(it_min_age, 't', it_max_age)

    for (it_solu_round in c(1, 2)) {

      it_solu_group <- it_solu_round*100 + it_solu_type
      print(paste0('solution group: ', it_solu_group, '=it_solu_group started'))

      # Start timer
      tm_start_solu <- proc.time()

      # Base settings
      bl_graph <- TRUE
      bl_threshold <- FALSE
      it_age_bins <- 1
      ar_rho <- 1 - (10^(c(seq(-2,2, length.out=8))))
      ar_rho <- unique(c(1,ar_rho))
      # ar_rho <- c(1)
      ar_rho_g47 <- c(1)
      ar_rho_g47 <- unique(c(1,ar_rho))
      bl_optimal_with_age_groups <- FALSE

      # Number of Checks to Consider
      if ((it_solu_group >= 101 & it_solu_group <=104) | (it_solu_group >= 201 & it_solu_group <=204)) {
        it_max_checks <- 44
        spt_img_save_use <- paste0(spt_img_save, '_4400_', st_img_suf_age_ybin,'/')
        srt_csv_path_use <- paste0(srt_csv_path, '_4400_', st_img_suf_age_ybin,'/')
        st_solu_name_check_count <- '44 checks'
        st_image_name <- ''
      } else {
        it_max_checks <- 200
        spt_img_save_use <- paste0(spt_img_save, '_20k_', st_img_suf_age_ybin,'/')
        srt_csv_path_use <- paste0(srt_csv_path, '_20k_', st_img_suf_age_ybin,'/')
        st_solu_name_check_count <- '200 checks'
        st_image_name <- '_20k'
      }

      # Generate Path if does not exist
      dir.create(file.path(spt_img_save_use), showWarnings = FALSE, recursive = TRUE)
      dir.create(file.path(srt_csv_path_use), showWarnings = FALSE, recursive = TRUE)

      # Seven Types of Allocations
      if (it_solu_group == 101) {
        bl_threshold <- TRUE
        st_tfo_method <- 'thresold'
        st_suffix_csvimg <- paste0(st_tfo_method, st_image_name)
        st_img_suf <- paste0('_round1_', st_suffix_csvimg)
        st_solu_name <- paste('1st threshold', it_max_checks)
      }
      if (it_solu_group == 102 | it_solu_group == 105) {
        st_tfo_method <- 'feasible'
        st_suffix_csvimg <- paste0(st_tfo_method, st_image_name)
        st_img_suf <- paste0('_round1_', st_suffix_csvimg)
        st_solu_name <- paste('1st feasible', it_max_checks)
      }
      if (it_solu_group == 103 | it_solu_group == 106) {
        bl_optimal_with_age_groups <- TRUE
        it_age_bins <- 4
        st_tfo_method <- 'optimalg4'
        st_suffix_csvimg <- paste0(st_tfo_method, st_image_name)
        st_img_suf <- paste0('_round1_', st_suffix_csvimg)
        st_solu_name <- paste('1st optimal g4', it_max_checks)
      }
      if (it_solu_group == 104 | it_solu_group == 107) {
        bl_optimal_with_age_groups <- TRUE
        bl_graph <- FALSE
        it_age_bins <- it_max_age - it_min_age + 1
        st_tfo_method <- 'optimalg47'
        st_suffix_csvimg <- paste0(st_tfo_method, st_image_name)
        st_img_suf <- paste0('_round1_', st_suffix_csvimg)
        st_solu_name <- paste('1st optimal g4', it_max_checks)
        ar_rho <- ar_rho_g47
      }

      # Seven Types of Allocations 2nd round
      if (it_solu_group == 201) {
        bl_threshold <- TRUE
        st_tfo_method <- 'thresold'
        st_suffix_csvimg <- paste0(st_tfo_method, st_image_name)
        st_img_suf <- paste0('_round2_', st_suffix_csvimg)
        st_solu_name <- paste('2nd threshold', it_max_checks)
      }
      if (it_solu_group == 202 | it_solu_group == 205) {
        st_tfo_method <- 'feasible'
        st_suffix_csvimg <- paste0(st_tfo_method, st_image_name)
        st_img_suf <- paste0('_round2_', st_suffix_csvimg)
        st_solu_name <- paste('2nd feasible', it_max_checks)
      }
      if (it_solu_group == 203 | it_solu_group == 206) {
        bl_optimal_with_age_groups <- TRUE
        it_age_bins <- 4
        st_tfo_method <- 'optimalg4'
        st_suffix_csvimg <- paste0(st_tfo_method, st_image_name)
        st_img_suf <- paste0('_round2_', st_suffix_csvimg)
        st_solu_name <- paste('2nd optimal g4', it_max_checks)
      }
      if (it_solu_group == 204 | it_solu_group == 207) {
        bl_optimal_with_age_groups <- TRUE
        bl_graph <- FALSE
        it_age_bins <- it_max_age - it_min_age + 1
        st_tfo_method <- 'optimalg47'
        st_suffix_csvimg <- paste0(st_tfo_method, st_image_name)
        st_img_suf <- paste0('_round2_', st_suffix_csvimg)
        st_solu_name <- paste('2nd optimal g4', it_max_checks)
        ar_rho <- ar_rho_g47
      }

      # First Round Allocation Solve ------------------------
      if (it_solu_group >= 101 & it_solu_group <= 107) {

        print(paste0(it_solu_group, ': First Round Allocation Started'))

        ls_prc_outputs_zs5_1st <- ffp_snw_process_inputs(
          srt_simu_path = srt_simu_path,
          snm_simu_csv = snm_simu_csv_withspouse_shock,
          df_plan_v_tilde_full = df_plan_v_tilde_full,
          fl_max_phaseout = fl_max_phaseout,
          it_bin_dollar_before_phaseout = it_bin_dollar_before_phaseout,
          fl_percheck_dollar = fl_percheck_dollar,
          fl_multiple = fl_multiple,
          it_max_checks = it_max_checks,
          fl_tax_hh = fl_tax_hh,
          it_max_age = it_max_age,
          it_min_age = it_min_age,
          it_age_bins = it_age_bins,
          ar_svr_csv = ar_svr_csv,
          ar_svr_groups = ar_svr_groups,
          ar_svr_groups_stats = ar_svr_groups_stats,
          svr_checks = svr_checks,
          svr_v_value = svr_v_value,
          svr_c_value = svr_c_value,
          svr_mass = svr_mass,
          ar_rho = ar_rho,
          bl_threshold = bl_threshold,
          bl_non_inc_adjust = TRUE,
          bl_print = FALSE,
          bl_print_verbose = FALSE)

        print(paste0(it_solu_group, ': First Round Allocation Finished'))

        # First Round Solution Results
        tb_rho_rev_c=ls_prc_outputs_zs5_1st$tb_rho_rev_c
        tb_rho_rev_v=ls_prc_outputs_zs5_1st$tb_rho_rev_v

        df_input_il_noninc_covar_zs5_1st=ls_prc_outputs_zs5_1st$df_input_il_noninc_covar
        df_alloc_i_long_covar_c_zs5_1st=ls_prc_outputs_zs5_1st$df_alloc_i_long_covar_c
        df_alloc_i_long_covar_v_zs5_1st=ls_prc_outputs_zs5_1st$df_alloc_i_long_covar_v

        stg_subtitle=ls_prc_outputs_zs5_1st$stg_subtitle
        stg_caption=ls_prc_outputs_zs5_1st$stg_caption

        # Graphing
        if (bl_graph) {
          allocate_type_levels <- c(Actual = "actual", Optimal = "optimal")
          # Generic non feasible specific Figures
          st_img_suffix <- paste0(st_file_type_withspouse_shock, '_', st_img_suf_age_ybin, st_img_suf)
          ls_pl_generic <- ffp_snw_graph_feasible(
            ar_rho=ar_rho,
            df_input_il_noninc_covar=df_input_il_noninc_covar_zs5_1st,
            df_alloc_i_long_covar_c=
              df_alloc_i_long_covar_c_zs5_1st %>%
              mutate(allocate_type = fct_recode(allocate_type, !!!allocate_type_levels)),
            df_alloc_i_long_covar_v=
              df_alloc_i_long_covar_v_zs5_1st %>%
              mutate(allocate_type = fct_recode(allocate_type, !!!allocate_type_levels)),
            ls_st_gen_imgs = c('mass', 'mc', 'mv', 'mlogc'),
            ls_st_save_imgs = c('mass', 'mc', 'mv', 'mlogc'),
            bl_optimal_with_age_groups = bl_optimal_with_age_groups,
            st_file_type = st_file_type_withspouse_shock,
            slb_add_title = slb_add_title_round1,
            stg_subtitle=stg_subtitle, stg_caption=stg_caption,
            it_img_width=it_img_width, it_img_height=it_img_height,
            st_img_units='mm', it_img_res=it_img_res, it_img_pointsize=it_img_pointsize,
            st_img_suffix=st_img_suffix, bl_save_img=bl_save_img,
            spt_img_save=spt_img_save_use)

          # Generate Graphs
          if ((it_solu_group >= 101 & it_solu_group <=104) | (it_solu_group >= 201 & it_solu_group <=204)) {
            ls_st_checks_cv <- c('checks_cjv', 'checks_c', 'checks_v')
          } else {
            ls_st_checks_cv <- c('checks_cjv', 'checks_cjv_log', 'checks_c', 'checks_c_log', 'checks_v', 'checks_v_log')
          }
          st_img_suffix <- paste0(st_file_type_withspouse_shock, '_', st_img_suf_age_ybin, st_img_suf)
          ls_pl <- ffp_snw_graph_feasible(
            ar_rho=ar_rho,
            df_input_il_noninc_covar=df_input_il_noninc_covar_zs5_1st,
            df_alloc_i_long_covar_c=
              df_alloc_i_long_covar_c_zs5_1st %>%
              mutate(allocate_type = fct_recode(allocate_type, !!!allocate_type_levels)),
            df_alloc_i_long_covar_v=
              df_alloc_i_long_covar_v_zs5_1st %>%
              mutate(allocate_type = fct_recode(allocate_type, !!!allocate_type_levels)),
            ls_st_gen_imgs = ls_st_checks_cv,
            ls_st_save_imgs = ls_st_checks_cv,
            bl_optimal_with_age_groups = bl_optimal_with_age_groups,
            st_file_type = st_file_type_withspouse_shock,
            slb_add_title = slb_add_title_round1,
            stg_subtitle=stg_subtitle, stg_caption=stg_caption,
            it_img_width=it_img_width, it_img_height=it_img_height,
            st_img_units='mm', it_img_res=it_img_res, it_img_pointsize=it_img_pointsize,
            st_img_suffix=st_img_suffix, bl_save_img=bl_save_img,
            spt_img_save=spt_img_save_use)
        }
      }

      # Second Round Allocation Solve ------------------------
      if (it_solu_group >= 201 & it_solu_group <= 207) {

        print(paste0(it_solu_group, ': Second Round Allocation Started'))

        bl_given_firstcheck <- TRUE
        ls_prc_outputs_zs5_2nd <- ffp_snw_process_inputs(
          srt_simu_path = srt_simu_path,
          snm_simu_csv = snm_simu_csv_withspouse_shock,
          df_plan_v_tilde_full = df_plan_v_tilde_full,
          fl_max_phaseout = fl_max_phaseout,
          it_bin_dollar_before_phaseout = it_bin_dollar_before_phaseout,
          fl_percheck_dollar = fl_percheck_dollar,
          fl_multiple = fl_multiple,
          it_max_checks = it_max_checks,
          fl_tax_hh = fl_tax_hh,
          it_max_age = it_max_age,
          it_min_age = it_min_age,
          it_age_bins = it_age_bins,
          ar_svr_csv = ar_svr_csv,
          ar_svr_groups = ar_svr_groups,
          ar_svr_groups_stats = ar_svr_groups_stats,
          svr_checks = svr_checks,
          svr_v_value = svr_v_value,
          svr_c_value = svr_c_value,
          svr_mass = svr_mass,
          ar_rho = ar_rho,
          bl_given_firstcheck = bl_given_firstcheck,
          bl_threshold = bl_threshold,
          bl_non_inc_adjust = TRUE,
          bl_print = FALSE,
          bl_print_verbose = FALSE)

        tb_rho_rev_c=ls_prc_outputs_zs5_2nd$tb_rho_rev_c
        tb_rho_rev_v=ls_prc_outputs_zs5_2nd$tb_rho_rev_v

        df_input_il_noninc_covar_zs5_2nd=ls_prc_outputs_zs5_2nd$df_input_il_noninc_covar
        df_alloc_i_long_covar_c_zs5_2nd=ls_prc_outputs_zs5_2nd$df_alloc_i_long_covar_c
        df_alloc_i_long_covar_v_zs5_2nd=ls_prc_outputs_zs5_2nd$df_alloc_i_long_covar_v

        print(paste0(it_solu_group, ': Second Round Allocation Finished'))

        if (bl_graph) {
          # 2nd MV and MC ------
          stg_subtitle_stack=paste0(slb_subtitle_stack_2nd, stg_subtitle)
          st_img_suffix <- paste0(st_file_type_withspouse_shock, '_', st_img_suf_age_ybin, st_img_suf)
          ls_pl_r2 <- ffp_snw_graph_feasible(ar_rho=ar_rho,
                                             df_input_il_noninc_covar=df_input_il_noninc_covar_zs5_2nd,
                                             df_alloc_i_long_covar_c=df_alloc_i_long_covar_c_zs5_2nd,
                                             df_alloc_i_long_covar_v=df_alloc_i_long_covar_v_zs5_2nd,
                                             ls_st_gen_imgs = c('mv', 'mc', 'mlogc'),
                                             ls_st_save_imgs = c('mv', 'mc', 'mlogc'),
                                             bl_optimal_with_age_groups = bl_optimal_with_age_groups,
                                             st_file_type = st_file_type_withspouse_shock,
                                             slb_add_title = slb_add_title_round2,
                                             stg_subtitle=stg_subtitle_stack, stg_caption=stg_caption,
                                             it_img_width=it_img_width, it_img_height=it_img_height,
                                             st_img_units='mm', it_img_res=it_img_res, it_img_pointsize=it_img_pointsize,
                                             st_img_suffix=st_img_suffix, bl_save_img=bl_save_img,
                                             spt_img_save=spt_img_save_use)

          # 2nd Allocation Figures ------
          # Generate Graphs
          if ((it_solu_group >= 101 & it_solu_group <=104) | (it_solu_group >= 201 & it_solu_group <=204)) {
            ls_st_checks_cv <- c('checks_c', 'checks_v')
          } else {
            ls_st_checks_cv <- c('checks_c', 'checks_c_log', 'checks_v', 'checks_v_log')
          }
          allocate_type_levels <- c(Actual = "actual", Optimal = "optimal")
          st_img_suffix <- paste0(st_file_type_withspouse_shock, '_', st_img_suf_age_ybin, st_img_suf)
          ls_pl_2nd <- ffp_snw_graph_feasible(ar_rho=ar_rho,
                                              df_input_il_noninc_covar=df_input_il_noninc_covar_zs5_2nd,
                                              df_alloc_i_long_covar_c=
                                                df_alloc_i_long_covar_c_zs5_2nd %>%
                                                mutate(allocate_type = fct_recode(allocate_type, !!!allocate_type_levels)),
                                              df_alloc_i_long_covar_v=
                                                df_alloc_i_long_covar_v_zs5_2nd %>%
                                                mutate(allocate_type = fct_recode(allocate_type, !!!allocate_type_levels)),
                                              ls_st_gen_imgs = ls_st_checks_cv,
                                              ls_st_save_imgs = ls_st_checks_cv,
                                              bl_optimal_with_age_groups = bl_optimal_with_age_groups,
                                              st_file_type = st_file_type_withspouse_shock,
                                              slb_add_title = slb_add_title_round1,
                                              stg_subtitle=stg_subtitle, stg_caption=stg_caption,
                                              it_img_width=it_img_width, it_img_height=it_img_height,
                                              st_img_units='mm', it_img_res=it_img_res, it_img_pointsize=it_img_pointsize,
                                              st_img_suffix=st_img_suffix, bl_save_img=bl_save_img,
                                              spt_img_save=spt_img_save_use)

          # 1st and 2nd Stack Graph ------
          # Assume Round 1 actual duplicated in Round 2.
          # Labels
          allocate_type_levels <- c("Opti 2nd + Actual 1st" = "opti2nd_total", "Actual 1st" = "actual")
          slb_subtitle_stack1st2nd = ''
          # Update Consumption Frame
          df_alloc_i_long_covar_c_zs5_2nd_stack <- df_alloc_i_long_covar_c_zs5_2nd %>%
            group_by(rho, id_i) %>%
            pivot_wider(names_from = allocate_type,
                        values_from = checks) %>%
            mutate(total = optimal + actual) %>%
            rename(checks_optimal = optimal,
                   checks_actual = actual,
                   checks_opti2nd_total = total) %>%
            pivot_longer(cols = starts_with('checks'),
                         names_to = c('allocate_type'),
                         names_pattern = paste0("checks_(.*)"),
                         values_to = "checks") %>%
            filter(allocate_type != 'optimal') %>%
            mutate(allocate_type = fct_recode(allocate_type, !!!allocate_type_levels))

          # Update Value Frame
          df_alloc_i_long_covar_v_zs5_2nd_stack <- df_alloc_i_long_covar_v_zs5_2nd %>%
            group_by(rho, id_i) %>%
            pivot_wider(names_from = allocate_type,
                        values_from = checks) %>%
            mutate(total = optimal + actual) %>%
            rename(checks_optimal = optimal,
                   checks_actual = actual,
                   checks_opti2nd_total = total) %>%
            pivot_longer(cols = starts_with('checks'),
                         names_to = c('allocate_type'),
                         names_pattern = paste0("checks_(.*)"),
                         values_to = "checks") %>%
            filter(allocate_type != 'optimal') %>%
            mutate(allocate_type = fct_recode(allocate_type, !!!allocate_type_levels))

          # graph
          stg_subtitle_stack=paste0(slb_subtitle_stack_2nd, stg_subtitle)
          st_img_suffix <- paste0(st_file_type_withspouse_shock, '_', st_img_suf_age_ybin, paste0(st_img_suf, '_stack'))
          ls_pl_stack <- ffp_snw_graph_feasible(ar_rho=ar_rho,
                                                df_input_il_noninc_covar=df_input_il_noninc_covar_zs5_2nd,
                                                df_alloc_i_long_covar_c=df_alloc_i_long_covar_c_zs5_2nd_stack,
                                                df_alloc_i_long_covar_v=df_alloc_i_long_covar_v_zs5_2nd_stack,
                                                ls_st_gen_imgs = c('checks_c', 'checks_v'),
                                                ls_st_save_imgs = c('checks_c', 'checks_v'),
                                                bl_optimal_with_age_groups = bl_optimal_with_age_groups,
                                                st_file_type = st_file_type_withspouse_shock,
                                                slb_add_title = slb_add_title_round2,
                                                stg_subtitle=stg_subtitle_stack, stg_caption=stg_caption,
                                                it_img_width=it_img_width, it_img_height=it_img_height,
                                                st_img_units='mm', it_img_res=it_img_res, it_img_pointsize=it_img_pointsize,
                                                st_img_suffix=st_img_suffix, bl_save_img=bl_save_img,
                                                spt_img_save=spt_img_save_use)

          # 1st and 2nd Joint Graph ----
          allocate_type_levels <- c("1st round" = "1st_round", "2nd round" = "2nd_round")
          # Update Consumption Frame
          df_alloc_i_long_covar_c_zs5_2nd_joint <- df_alloc_i_long_covar_c_zs5_2nd %>%
            group_by(rho, id_i) %>%
            pivot_wider(names_from = allocate_type,
                        values_from = checks) %>%
            rename(optimal_2nd_round = optimal) %>%
            left_join(df_alloc_i_long_covar_c_zs5_1st %>%
                        filter( allocate_type == 'optimal') %>%
                        select(id_i, rho, checks) %>%
                        rename(optimal_1st_round = checks),
                      by=setNames(c('id_i', 'rho'), c('id_i', 'rho'))) %>%
            select(-actual) %>%
            pivot_longer(cols = starts_with('optimal'),
                         names_to = c('allocate_type'),
                         names_pattern = paste0("optimal_(.*)"),
                         values_to = "checks") %>%
            mutate(allocate_type = fct_recode(allocate_type, !!!allocate_type_levels))
          # Update Value Frame
          df_alloc_i_long_covar_v_zs5_2nd_joint <- df_alloc_i_long_covar_v_zs5_2nd %>%
            group_by(rho, id_i) %>%
            pivot_wider(names_from = allocate_type,
                        values_from = checks) %>%
            rename(optimal_2nd_round = optimal) %>%
            left_join(df_alloc_i_long_covar_v_zs5_1st %>%
                        filter( allocate_type == 'optimal') %>%
                        select(id_i, rho, checks) %>%
                        rename(optimal_1st_round = checks),
                      by=setNames(c('id_i', 'rho'), c('id_i', 'rho'))) %>%
            select(-actual) %>%
            pivot_longer(cols = starts_with('optimal'),
                         names_to = c('allocate_type'),
                         names_pattern = paste0("optimal_(.*)"),
                         values_to = "checks") %>%
            mutate(allocate_type = fct_recode(allocate_type, !!!allocate_type_levels))

          # Generate Graphs
          stg_subtitle_joint=paste0(slb_subtitle_joint_1st2nd, stg_subtitle)
          st_img_suffix <- paste0(st_file_type_withspouse_shock, '_', st_img_suf_age_ybin, paste0(st_img_suf, '_joint'))
          ls_pl_joint <- ffp_snw_graph_feasible(ar_rho=ar_rho,
                                                df_input_il_noninc_covar=df_input_il_noninc_covar_zs5_2nd,
                                                df_alloc_i_long_covar_c=df_alloc_i_long_covar_c_zs5_2nd_joint,
                                                df_alloc_i_long_covar_v=df_alloc_i_long_covar_v_zs5_2nd_joint,
                                                ls_st_gen_imgs = c('checks_c', 'checks_v'),
                                                ls_st_save_imgs = c('checks_c', 'checks_v'),
                                                bl_optimal_with_age_groups = bl_optimal_with_age_groups,
                                                st_file_type = st_file_type_withspouse_shock,
                                                slb_add_title = slb_add_title_round2,
                                                stg_subtitle=stg_subtitle_joint, stg_caption=stg_caption,
                                                it_img_width=it_img_width, it_img_height=it_img_height,
                                                st_img_units='mm', it_img_res=it_img_res, it_img_pointsize=it_img_pointsize,
                                                st_img_suffix=st_img_suffix, bl_save_img=bl_save_img,
                                                spt_img_save=spt_img_save_use)
        }
      }

      # End Timer
      tm_solu <- proc.time() - tm_start_solu
      print(paste0(st_solu_name, ', tm_solu:', tm_solu[['elapsed']], ' seconds'))

    }

    # CSV Save all Results first and second rounds --------------------
    df_alloc_all <- rbind(
      ls_prc_outputs_zs5_1st$df_alloc_i_long_covar_v %>%
        filter(rho_val == ar_rho[1]) %>%
        mutate(allocate_type = case_when(allocate_type == 'optimal' ~ 'optimal_v_1st',
                                         TRUE ~ allocate_type)),
      ls_prc_outputs_zs5_1st$df_alloc_i_long_covar_c %>%
        filter(rho_val == ar_rho[1]) %>%
        filter(allocate_type == 'optimal') %>%
        mutate(allocate_type = case_when(allocate_type == 'optimal' ~ 'optimal_c_1st',
                                         TRUE ~ allocate_type)),
      ls_prc_outputs_zs5_2nd$df_alloc_i_long_covar_v %>%
        filter(rho_val == ar_rho[1]) %>%
        filter(allocate_type == 'optimal') %>%
        mutate(allocate_type = case_when(allocate_type == 'optimal' ~ 'optimal_v_2nd',
                                         TRUE ~ allocate_type)),
      ls_prc_outputs_zs5_2nd$df_alloc_i_long_covar_c %>%
        filter(rho_val == ar_rho[1]) %>%
        filter(allocate_type == 'optimal') %>%
        mutate(allocate_type = case_when(allocate_type == 'optimal' ~ 'optimal_c_2nd',
                                         TRUE ~ allocate_type))) %>%
      arrange(id_i, allocate_type) %>%
      select(-rho, -F_star_i, -EH_star_i, -survive) %>%
      select(id_i, rho_val, everything())

    df_alloc_all <- df_alloc_all %>%
      pivot_wider(names_from = allocate_type,
                  values_from = checks)

    # export
    write.csv(df_alloc_all,
              paste0(srt_csv_path_use, paste0('df_alloc_all_',st_suffix_csvimg,'.csv'), row.names = TRUE))


    # CSV Export Grouped Average Statistics -------
    ls_svr_groups <- ar_svr_groups
    for (svr_group in ls_svr_groups) {

      # group mean
      df_alloc_combine_group_mean <- df_alloc_all %>%
        ungroup() %>% group_by(!!sym(svr_group)) %>%
        summarize(actual_mean = sum(actual*mass)/sum(mass),
                  optimal_c_1st_mean = sum(optimal_c_1st*mass)/sum(mass),
                  optimal_v_1st_mean = sum(optimal_v_1st*mass)/sum(mass),
                  optimal_c_2nd_mean = sum(optimal_c_2nd*mass)/sum(mass),
                  optimal_v_2nd_mean = sum(optimal_v_2nd*mass)/sum(mass),
                  mass_sum = sum(mass))

      # Export
      write.csv(df_alloc_combine_group_mean,
                paste0(srt_csv_path_use, "df_alloc_",svr_group, "_", st_suffix_csvimg, ".csv"),
                row.names = TRUE)

      # All but Group mean
      ls_svr_groups_oneless <- ls_svr_groups[ls_svr_groups != svr_group]
      df_alloc_combine_group_mean_oneless <- df_alloc_all %>%
        ungroup() %>% group_by(!!!syms(ls_svr_groups_oneless)) %>%
        summarize(actual_mean = sum(actual*mass)/sum(mass),
                  optimal_c_1st_mean = sum(optimal_c_1st*mass)/sum(mass),
                  optimal_v_1st_mean = sum(optimal_v_1st*mass)/sum(mass),
                  optimal_c_2nd_mean = sum(optimal_c_2nd*mass)/sum(mass),
                  optimal_v_2nd_mean = sum(optimal_v_2nd*mass)/sum(mass),
                  mass_sum = sum(mass))

      # Export
      write.csv(df_alloc_combine_group_mean_oneless,
                paste0(srt_csv_path_use, "df_alloc_",svr_group,"_without_",st_suffix_csvimg,".csv"),
                row.names = TRUE)

    }

    # CSV Average Statistics without Age -----------
    ls_svr_groups <- ar_svr_groups_noage
    for (svr_group in ls_svr_groups) {

      # All but Group mean
      ls_svr_groups_oneless <- ls_svr_groups[ls_svr_groups != svr_group]
      df_alloc_combine_group_mean_oneless <- df_alloc_all %>%
        ungroup() %>% group_by(!!!syms(ls_svr_groups_oneless)) %>%
        summarize(actual_mean = sum(actual*mass)/sum(mass),
                  optimal_c_1st_mean = sum(optimal_c_1st*mass)/sum(mass),
                  optimal_v_1st_mean = sum(optimal_v_1st*mass)/sum(mass),
                  optimal_c_2nd_mean = sum(optimal_c_2nd*mass)/sum(mass),
                  optimal_v_2nd_mean = sum(optimal_v_2nd*mass)/sum(mass),
                  mass_sum = sum(mass))

      # Export
      write.csv(df_alloc_combine_group_mean_oneless,
                paste0(srt_csv_path_use, "df_alloc_",svr_group,"_without_noage_", st_suffix_csvimg, ".csv"),
                row.names = TRUE)

    }

    ### CSV Print and save MASS REV results ---------------
    # Save REV to table, Stack them
    tb_rho_rev_mass_v1_tab <- ls_prc_outputs_zs5_1st$tb_rho_rev_v %>%
      mutate(objective = 'vlife',
             constraint = st_tfo_method,
             allocround = 'first',
             maxchecks = it_max_checks)
    tb_rho_rev_mass_c1_tab <- ls_prc_outputs_zs5_1st$tb_rho_rev_c %>%
      mutate(objective = 'c2020',
             constraint = st_tfo_method,
             allocround = 'first',
             maxchecks = it_max_checks)
    tb_rho_rev_mass_v2_tab <- ls_prc_outputs_zs5_2nd$tb_rho_rev_v %>%
      mutate(objective = 'vlife',
             constraint = st_tfo_method,
             allocround = 'second',
             maxchecks = it_max_checks)
    tb_rho_rev_mass_c2_tab <- ls_prc_outputs_zs5_2nd$tb_rho_rev_c %>%
      mutate(objective = 'c2020',
             constraint = st_tfo_method,
             allocround = 'second',
             maxchecks = it_max_checks)
    # Stack frames
    tb_rho_rev_mass_v_c_tab <- rbind(tb_rho_rev_mass_v1_tab, tb_rho_rev_mass_c1_tab,
                                     tb_rho_rev_mass_v2_tab, tb_rho_rev_mass_c2_tab)
    # export
    write.csv(tb_rho_rev_mass_v_c_tab,
              paste0(srt_csv_path_use, "rev_", st_suffix_csvimg, ".csv"),
              row.names = TRUE)

  }
}

# stopCluster(cl)
