# This file works out one particular optimal allocation problem, and draws allocation by full queue
# This is more time-consuming than normal, and generates potentially memory breaking dataframes, so should
# not be used when queue is too long or the number of candidate individuals is too large.

# the controlling loops here are identical to fs_opti.R, but the code here does not include graphing and csv
# generation files as there.

library(tidyverse)
library(REconTools)
# library(PrjOptiAlloc)
library(forcats)
#


# Show Queue Allocation i col queue row panel
bl_df_alloc_il = TRUE
# Return all queue positions even after resources
bl_return_allQ_V = TRUE


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
it_img_width=270
it_img_height=216
st_img_units='mm'
it_img_res=300
it_img_pointsize=7

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
ar_svr_groups_age_ymin <- c('age_group', 'ymin_group')
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

# Loop Over Soluiton Files
# it_solu_type <- 1
# it_age_type <- 1
# it_solu_round <- 1

ar_double_triple_alloc <- c(32)
# 22 is mulone feasible
# 32 is double kids feasible, 33 is double kids g4, 34 is double kids g47
# 42 is double adults feasible, 33 is double adults g4, 34 is double adults g47
for (it_solu_type in ar_double_triple_alloc) {

  for (it_age_type in c(2)) {
    if (it_age_type == 1) {
      it_max_age = 64
      it_min_age = 18
    }
    if (it_age_type == 2) {
      it_max_age = 99
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
      it_max_age = 100
      it_min_age = 22
    }
    if (it_age_type == 6) {
      it_max_age = 74
      it_min_age = 22
    }
    # Image Save Suffix
    st_img_suf_age_ybin <- paste0(it_min_age, 't', it_max_age)

    for (it_solu_round in c(2)) {

      it_solu_group <- it_solu_round*100 + it_solu_type
      print(paste0('solution group: ', it_solu_group, '=it_solu_group started'))

      # Start timer
      tm_start_solu <- proc.time()

      # Base settings
      bl_graph <- TRUE
      bl_threshold <- FALSE
      it_age_bins <- 1
      # ar_rho <- 1 - (10^(c(seq(-2,2, length.out=8))))
      # ar_rho <- unique(c(1,ar_rho))
      ar_rho <- c(1)
      ar_rho_g47 <- c(1)
      # ar_rho_g47 <- unique(c(1,ar_rho))
      bl_optimal_with_age_groups <- FALSE
      it_check_headorspouse <- 12
      it_check_perkids <- 5

      # Number of Checks to Consider
      if ((it_solu_group >= 101 & it_solu_group <=105) |
          (it_solu_group >= 201 & it_solu_group <=205)
      ) {
        # x02 to x04 is allocating up to 44 checks, common bounds for all, except for x01
        # x01 is threshold allocation, in 4400 folder
        it_max_checks <- 44
        spt_img_save_use <- paste0(spt_img_save, '_4400_', st_img_suf_age_ybin,'/')
        srt_csv_path_use <- paste0(srt_csv_path, '_4400_', st_img_suf_age_ybin,'/')
        st_solu_name_check_count <- '44 checks'
        st_image_name <- ''
      } else if ((it_solu_group >= 111 & it_solu_group <=115) |
                 (it_solu_group >= 211 & it_solu_group <=215)) {
        # x05 to x07 is allocating up to 200k, common bounds for all
        it_max_checks <- 200
        spt_img_save_use <- paste0(spt_img_save, '_20k_', st_img_suf_age_ybin,'/')
        srt_csv_path_use <- paste0(srt_csv_path, '_20k_', st_img_suf_age_ybin,'/')
        st_solu_name_check_count <- '200 checks'
        st_image_name <- '_20k'
      } else if ((it_solu_group >= 121 & it_solu_group <=125) |
                 (it_solu_group >= 221 & it_solu_group <=225)) {
        # same as actual allocation, so this is the same as threshold allocation but
        # now does threshold allocation with g4 and g47 as well.
        bl_threshold <- TRUE
        it_max_checks <- 44
        it_check_headorspouse <- 12
        it_check_perkids <- 5
        spt_img_save_use <- paste0(spt_img_save, '_mulone_', st_img_suf_age_ybin,'/')
        srt_csv_path_use <- paste0(srt_csv_path, '_mulone_', st_img_suf_age_ybin,'/')
        st_solu_name_check_count <- 'actual checks adults'
        st_image_name <- '_mulone'
      } else if ((it_solu_group >= 131 & it_solu_group <=135) |
                 (it_solu_group >= 231 & it_solu_group <=235)) {
        # double allocation per kid to 1000 vs 500, so 4 households = 4000 from kids + 2400
        # 64 checks for 4 kids households married
        bl_threshold <- TRUE
        it_max_checks <- 64
        it_check_headorspouse <- 12
        it_check_perkids <- 10
        spt_img_save_use <- paste0(spt_img_save, '_dblkid_', st_img_suf_age_ybin,'/')
        srt_csv_path_use <- paste0(srt_csv_path, '_dblkid_', st_img_suf_age_ybin,'/')
        st_solu_name_check_count <- 'double actual checks kids'
        st_image_name <- '_dblkid'
      } else if ((it_solu_group >= 141 & it_solu_group <=145) |
                 (it_solu_group >= 241 & it_solu_group <=245)) {
        # double adult allocation from 12 to 24
        # four kids married = 24*2 + 5*4 = 68
        bl_threshold <- TRUE
        it_max_checks <- 68
        it_check_headorspouse <- 24
        it_check_perkids <- 5
        spt_img_save_use <- paste0(spt_img_save, '_dbladt_', st_img_suf_age_ybin,'/')
        srt_csv_path_use <- paste0(srt_csv_path, '_dbladt_', st_img_suf_age_ybin,'/')
        st_solu_name_check_count <- 'double actual checks adults'
        st_image_name <- '_dbladt'
      } else if ((it_solu_group >= 151 & it_solu_group <=155) |
                 (it_solu_group >= 251 & it_solu_group <=255)) {
        # double adult allocation from 12 to 24, kids 5 to 10
        # four kids married = 24*2 + 10*4 = 88
        bl_threshold <- TRUE
        it_max_checks <- 88
        it_check_headorspouse <- 24
        it_check_perkids <- 10
        spt_img_save_use <- paste0(spt_img_save, '_dblbth_', st_img_suf_age_ybin,'/')
        srt_csv_path_use <- paste0(srt_csv_path, '_dblbth_', st_img_suf_age_ybin,'/')
        st_solu_name_check_count <- 'double actual checks adults + kids'
        st_image_name <- '_dblbth'
      } else if ((it_solu_group >= 161 & it_solu_group <=165) |
                 (it_solu_group >= 261 & it_solu_group <=265)) {
        # triple allocation per kid to 1500 vs 500, so 4 kids married = 1500*4 + 2400 = 8400
        # 64 checks for 4 kids households married
        bl_threshold <- TRUE
        it_max_checks <- 84
        it_check_headorspouse <- 12
        it_check_perkids <- 15
        spt_img_save_use <- paste0(spt_img_save, '_trbkid_', st_img_suf_age_ybin,'/')
        srt_csv_path_use <- paste0(srt_csv_path, '_trbkid_', st_img_suf_age_ybin,'/')
        st_solu_name_check_count <- 'triple actual checks kids'
        st_image_name <- '_trbkid'
      } else if ((it_solu_group >= 171 & it_solu_group <=175) |
                 (it_solu_group >= 271 & it_solu_group <=275)) {
        # double adult allocation from 12 to 36
        # four kids married = 36*2 + 5*4 = 92
        bl_threshold <- TRUE
        it_max_checks <- 92
        it_check_headorspouse <- 36
        it_check_perkids <- 5
        spt_img_save_use <- paste0(spt_img_save, '_trbadt_', st_img_suf_age_ybin,'/')
        srt_csv_path_use <- paste0(srt_csv_path, '_trbadt_', st_img_suf_age_ybin,'/')
        st_solu_name_check_count <- 'triple actual checks adults'
        st_image_name <- '_trbadt'
      } else if ((it_solu_group >= 181 & it_solu_group <=185) |
                 (it_solu_group >= 281 & it_solu_group <=285)) {
        # double adult allocation from 12 to 24, kids 5 to 10
        # four kids married = 36*2 + 15*4 = 132
        bl_threshold <- TRUE
        it_max_checks <- 132
        it_check_headorspouse <- 36
        it_check_perkids <- 15
        spt_img_save_use <- paste0(spt_img_save, '_trbbth_', st_img_suf_age_ybin,'/')
        srt_csv_path_use <- paste0(srt_csv_path, '_trbbth_', st_img_suf_age_ybin,'/')
        st_solu_name_check_count <- 'triple actual checks adults + kids'
        st_image_name <- '_trbbth'
      }

      # Generate Path if does not exist
      dir.create(file.path(spt_img_save_use), showWarnings = FALSE, recursive = TRUE)
      dir.create(file.path(srt_csv_path_use), showWarnings = FALSE, recursive = TRUE)

      # R split number
      ar_solu_group_digits <- as.numeric(strsplit(as.character(it_solu_group), "")[[1]])

      # Seven Types of Allocations
      if (it_solu_group == 101) {
        bl_threshold <- TRUE
        st_tfo_method <- 'thresold'
        st_suffix_csvimg <- paste0(st_tfo_method, st_image_name)
        st_img_suf <- paste0('_round1_', st_suffix_csvimg)
        st_solu_name <- paste('1st threshold', it_max_checks)
      }
      if (ar_solu_group_digits[1] == 1 & ar_solu_group_digits[3] == 2) {
        st_tfo_method <- 'feasible'
        st_suffix_csvimg <- paste0(st_tfo_method, st_image_name)
        st_img_suf <- paste0('_round1_', st_suffix_csvimg)
        st_solu_name <- paste('1st feasible', it_max_checks)
      }
      if (ar_solu_group_digits[1] == 1 & ar_solu_group_digits[3] == 3) {
        bl_optimal_with_age_groups <- TRUE
        it_age_bins <- 4
        st_tfo_method <- 'optimalg4'
        st_suffix_csvimg <- paste0(st_tfo_method, st_image_name)
        st_img_suf <- paste0('_round1_', st_suffix_csvimg)
        st_solu_name <- paste('1st optimal g4', it_max_checks)
      }
      if (ar_solu_group_digits[1] == 1 & ar_solu_group_digits[3] == 4) {
        bl_optimal_with_age_groups <- TRUE
        bl_graph <- FALSE
        it_age_bins <- it_max_age - it_min_age + 1
        st_tfo_method <- 'optimalg47'
        st_suffix_csvimg <- paste0(st_tfo_method, st_image_name)
        st_img_suf <- paste0('_round1_', st_suffix_csvimg)
        st_solu_name <- paste('1st optimal g47 1 planner', it_max_checks)
        ar_rho <- ar_rho_g47
      }
      if (ar_solu_group_digits[1] == 1 & ar_solu_group_digits[3] == 5) {
        bl_optimal_with_age_groups <- TRUE
        bl_graph <- FALSE
        it_age_bins <- it_max_age - it_min_age + 1
        st_tfo_method <- 'optimalg47'
        st_suffix_csvimg <- paste0(st_tfo_method, st_image_name)
        st_img_suf <- paste0('_round1_', st_suffix_csvimg)
        st_solu_name <- paste('1st optimal g47 all planners', it_max_checks)
      }

      # Seven Types of Allocations 2nd round
      if (it_solu_group == 201) {
        bl_threshold <- TRUE
        st_tfo_method <- 'thresold'
        it_check_headorspouse <- 12
        it_check_perkids <- 5
        st_suffix_csvimg <- paste0(st_tfo_method, st_image_name)
        st_img_suf <- paste0('_round2_', st_suffix_csvimg)
        st_solu_name <- paste('2nd threshold', it_max_checks)
      }
      if (ar_solu_group_digits[1] == 2 & ar_solu_group_digits[3] == 2) {
        st_tfo_method <- 'feasible'
        st_suffix_csvimg <- paste0(st_tfo_method, st_image_name)
        st_img_suf <- paste0('_round2_', st_suffix_csvimg)
        st_solu_name <- paste('2nd feasible', it_max_checks)
      }
      if (ar_solu_group_digits[1] == 2 & ar_solu_group_digits[3] == 3) {
        bl_optimal_with_age_groups <- TRUE
        it_age_bins <- 4
        st_tfo_method <- 'optimalg4'
        st_suffix_csvimg <- paste0(st_tfo_method, st_image_name)
        st_img_suf <- paste0('_round2_', st_suffix_csvimg)
        st_solu_name <- paste('2nd optimal g4', it_max_checks)
      }
      if (ar_solu_group_digits[1] == 2 & ar_solu_group_digits[3] == 4) {
        bl_optimal_with_age_groups <- TRUE
        bl_graph <- FALSE
        it_age_bins <- it_max_age - it_min_age + 1
        st_tfo_method <- 'optimalg47'
        st_suffix_csvimg <- paste0(st_tfo_method, st_image_name)
        st_img_suf <- paste0('_round2_', st_suffix_csvimg)
        st_solu_name <- paste('2nd optimal g47 1 planner', it_max_checks)
        ar_rho <- ar_rho_g47
      }
      if (ar_solu_group_digits[1] == 2 & ar_solu_group_digits[3] == 5) {
        bl_optimal_with_age_groups <- TRUE
        bl_graph <- FALSE
        it_age_bins <- it_max_age - it_min_age + 1
        st_tfo_method <- 'optimalg47'
        st_suffix_csvimg <- paste0(st_tfo_method, st_image_name)
        st_img_suf <- paste0('_round2_', st_suffix_csvimg)
        st_solu_name <- paste('2nd optimal g47 all planner', it_max_checks)
      }

      # First Round Allocation Solve ------------------------
      if (it_solu_group >= 101 & it_solu_group <= 199) {

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
          bl_df_alloc_il = bl_df_alloc_il,
          bl_return_allQ_V = bl_return_allQ_V,
          bl_threshold = bl_threshold,
          it_check_headorspouse = it_check_headorspouse,
          it_check_perkids = it_check_perkids,
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

        df_alloc_il_long_c_zs5_1st=ls_prc_outputs_zs5_1st$df_alloc_il_long_c
        df_alloc_il_long_v_zs5_1st=ls_prc_outputs_zs5_1st$df_alloc_il_long_v

        stg_subtitle=ls_prc_outputs_zs5_1st$stg_subtitle
        stg_caption=ls_prc_outputs_zs5_1st$stg_caption

        # Graphing
      }

      # Second Round Allocation Solve ------------------------
      if (it_solu_group >= 201 & it_solu_group <= 299) {

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
          bl_df_alloc_il = bl_df_alloc_il,
          bl_return_allQ_V = bl_return_allQ_V,
          bl_given_firstcheck = bl_given_firstcheck,
          bl_threshold = bl_threshold,
          it_check_headorspouse = it_check_headorspouse,
          it_check_perkids = it_check_perkids,
          bl_non_inc_adjust = TRUE,
          bl_print = FALSE,
          bl_print_verbose = FALSE)

        tb_rho_rev_c=ls_prc_outputs_zs5_2nd$tb_rho_rev_c
        tb_rho_rev_v=ls_prc_outputs_zs5_2nd$tb_rho_rev_v

        df_input_il_noninc_covar_zs5_2nd=ls_prc_outputs_zs5_2nd$df_input_il_noninc_covar
        df_alloc_i_long_covar_c_zs5_2nd=ls_prc_outputs_zs5_2nd$df_alloc_i_long_covar_c
        df_alloc_i_long_covar_v_zs5_2nd=ls_prc_outputs_zs5_2nd$df_alloc_i_long_covar_v

        df_alloc_il_long_c_zs5_2nd=ls_prc_outputs_zs5_2nd$df_alloc_il_long_c
        df_alloc_il_long_v_zs5_2nd=ls_prc_outputs_zs5_2nd$df_alloc_il_long_v

        print(paste0(it_solu_group, ': Second Round Allocation Finished'))

      }

      # End Timer
      tm_solu <- proc.time() - tm_start_solu
      print(paste0(st_solu_name, ', tm_solu:', tm_solu[['elapsed']], ' seconds'))

    }

  }
}

# stopCluster(cl)

# Queue household as columns to household id as category
df_alloc_il_long_longer <- df_alloc_il_long_c_zs5_2nd %>%
  pivot_longer(cols = starts_with('sid'),
               names_to = c('id_i'),
               names_pattern = paste0("sid_(.*)"),
               values_to = "checks") %>%
  mutate(id_i = as.numeric(id_i))

# Combine queue Rank info with household attributes and mass
df_alloc_il_long_longer_joined_id <- df_alloc_il_long_longer %>%
  left_join(df_alloc_i_long_covar_c_zs5_2nd %>%
              filter(allocate_type == 'optimal') %>%
              select(id_i, marital, kids, age_group, ymin_group, mass),
            by='id_i') %>%
  group_by(Q_il) %>%
  mutate(Q_il_mass = sum(mass*checks))

# Generate aggregate measure of checks at each queue point
df_group_graph <- df_alloc_il_long_longer_joined_id %>%
  group_by(Q_il, marital, kids) %>%
  summarize(checks_grp = sum(checks*mass),
            Q_il_mass = mean(Q_il_mass)) %>%
  mutate(percentage =  checks_grp/Q_il_mass)

# Graph group by kids and marital status
marry_levels <- c(Single = "0", Married = "1")
kids_levels <- c("no children" = "0", "one child" = "1",
                 "two children" = "2", "three children" = "3",
                 "four children" = "4")
df_group_graph <- df_group_graph %>%
  mutate(kids = as.factor(kids),
         marital = as.factor(marital)) %>%
  mutate(kids = fct_recode(kids, !!!kids_levels),
         marital = fct_recode(marital, !!!marry_levels))

# Graph
df_group_graph %>% ggplot(aes(x=Q_il_mass,
                              y=percentage,
                              colour=kids,
                              linetype=marital)) +
  geom_point() + ylim(0, 0.30) + xlim(0.2, 12)



