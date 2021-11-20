# .rs.restartR()
# rm(list = ls())
# rm(list = ls(all.names = TRUE))

# This is a R file, not RMD file, the loops through allocation problems to generate results.
# This is a R file, not RMD file, the loops through allocation problems to generate results.
library(tidyverse)
library(REconTools)

library(forcats)
library(foreach)
library(doParallel)

it_no_cores <- detectCores(logical = TRUE)
cl <- makeCluster(5)
registerDoParallel(cl)

# Various Lists
ls_st_file_suffix_temp <-
  c('snwx_bidenchk_moredense_a65zh266zs5_b1_xi0_manna_168')

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

ls_st_file_suffix_beta_heterogeneity <- c('snwx_bidenchk_docdense_e2hm2_b1_xi0_manna_88_bt99',
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

ls_st_file_suffix_test_multiple <- c('snwx_bidenchk_tiny_b1_xi0_manna_168',
                                     'snwx_bidenchk_tiny_b1_xi0_manna_168_unmarried',
                                     'snwx_bidenchk_tiny_b1_xi0_manna_168_married',
                                     'snwx_bidenchk_tiny_b1_xi0_manna_168_bt60',
                                     'snwx_bidenchk_tiny_b1_xi0_manna_168_bt95')
ls_st_file_suffix_test_multiple <- rev(ls_st_file_suffix_test_multiple)

ls_st_file_suffix_single_main <- c('snwx_bidenchk_moredense_a65zh266zs5_b1_xi0_manna_168')
ls_st_file_suffix_single_main <- c('snwx_bchklock_moredense_a65zh266zs5_b1_xi0_manna_168')
ls_st_file_suffix_test <- c('snwx_bidenchk_tiny_b1_xi0_manna_168')

# Run Collections-------------------------------------------------
it_run_collection <- 2

# OPTION 1 CSV SOLU to Run-------------------------------------------------
# list to run
ls_st_file_suffix <- c(ls_st_file_suffix_bidenchk,
                       ls_st_file_suffix_trumpchk,
                       ls_st_file_suffix_bidenchk_mixturealter,
                       ls_st_file_suffix_bchkbnoui,
                       ls_st_file_suffix_bchklock)
# ls_st_file_suffix <- c(ls_st_file_suffix_test_1)

# OPTION 2 Per Capita or Per Household-------------------------------------------------
# per capita or per household model, using what rho/elasticity values.
ls_bl_per_capita <- c(TRUE, FALSE)
ls_bl_per_capita <- c(TRUE)

# OPTION 3 Elasticity RHO Choices-------------------------------------------------
# Elasticities to try
# Set 1
# 10, 1, 0.1, 0.01; 5, 0.5. 0.05. I
# ar_elasticty <- c(10, 1, 0.1, 0.01, 5, 0.5, 0.05)
# ar_rho <- (ar_elasticty-1)/(ar_elasticty)
ar_rho_nine <- c(1, 0.9, 0.8,
                 0.05,
                 -1, -9, -19, -99)
# Set 2
ar_seq_equi <- seq(-2, 2, length.out=30)
ar_rho_30 <- 1 - (10^(c(ar_seq_equi)))
ar_rho_30 <- c(1, round(ar_rho_30, 4), -100)
# Set 3
ar_rho_segment <- 1 - (10^(c(seq(ar_seq_equi[23], ar_seq_equi[27], length.out=15))))
ar_rho_segment <- round(ar_rho_segment, 4)
# Set 4
ar_rho_intermediate <- c(-1)
# ar_rho all
ar_rho <- c(-1)
# ar_rho <- sort(unique(c(ar_rho_nine, ar_rho_30)))

# OPTION 4 Bounds to use-------------------------------------------------
# Here is a list of the 17 key results included in The Paper Tables
ar_double_triple_alloc_set1 <- c(22, 32, 42, 52, 62, 72, 82,
                                 23, 33, 43, 53, 63, 73, 83,
                                 12, 13,
                                 24)
ar_double_triple_alloc_set2 <- c(22, 52, 82, 23, 24)
ar_double_triple_alloc_set3 <- c(22, 52, 82)
ar_double_triple_alloc_set4 <- c(82)
ar_double_triple_alloc <- ar_double_triple_alloc_set3

# 12 is 170 checks, 13 is g4 for 170 checks
# 22 is mulone feasible, 23 is mulone g4, 24 is mulone g47
# 32 is double kids feasible, 33 is double kids g4, 34 is double kids g47
# 42 is double adults feasible, 33 is double adults g4, 34 is double adults g47
# 52 is double both feasible, 33 is double both g4, 34 is double both g47

# OPTION 5 Graph or not and other controls -------------------------------------------------
# graph or not
# bl_graph <- TRUE
bl_graph <- FALSE
ls_st_gen_agg_folders <- c('aggregator', 'aggcsvsolu', 'apcmpc', 'rev')
snm_main <- 'Results202111'

# OPTION 6 Perturb or Not ----------------------------------------------
fl_rand_adj_A_prop <- 0
it_rand_adj_A_rng_seed_start <- 0
ar_it_perturb_n <- -1

# Option 7 File saving options ---------------------------
bl_save_full_queue <- FALSE

# Run Collections ------------------------------
if (it_run_collection == 1) {

  # TEST RUN
  ls_st_file_suffix <- c(ls_st_file_suffix_test)
  ls_bl_per_capita <- c(TRUE)
  ar_rho <- c(1, -1, -99)
  ar_double_triple_alloc <- ar_double_triple_alloc_set3
  bl_graph <- TRUE
  bl_save_full_queue <- TRUE

} else if (it_run_collection == 2) {

  # Main run single filr
  ls_st_file_suffix <- c(ls_st_file_suffix_single_main)
  ls_bl_per_capita <- c(TRUE, FALSE)
  ar_rho <- ar_rho_nine
  ar_double_triple_alloc <- ar_double_triple_alloc_set1
  bl_graph <- TRUE

} else if (it_run_collection == 3) {
  # change fl_rho loop to paralle, use 5 workers, 6 max.
  # Main run many rhos
  ls_st_file_suffix <- c(ls_st_file_suffix_single_main)
  ls_bl_per_capita <- c(TRUE)
  ar_rho <- sort(unique(c(ar_rho_nine, ar_rho_30)))
  ar_double_triple_alloc <- ar_double_triple_alloc_set3
  bl_graph <- FALSE

} else if (it_run_collection == 4) {
  # run benchmark result with one single rho, 20 checks, generate also full queue
  ls_st_file_suffix <- c(ls_st_file_suffix_single_main)
  ls_bl_per_capita <- c(TRUE)
  ar_rho <- c(-1)
  ar_double_triple_alloc <- ar_double_triple_alloc_set4
  bl_graph <- FALSE
  bl_save_full_queue <- TRUE

}

# Run Collections Perturbation ------------------------------
if (it_run_collection >= 11 && it_run_collection <= 20) {
  ls_bl_per_capita <- c(TRUE)
  ls_st_gen_agg_folders <- c('')
  ar_rho <- ar_rho_intermediate
  snm_main <- 'Results202111_perturb'

  if (it_run_collection == 11 ) {

    # TEST RUN for perturbation exercise T
    ls_st_file_suffix <- c(ls_st_file_suffix_test)
    ar_double_triple_alloc <- ar_double_triple_alloc_set3
    # graph on to see if allocation results will look different
    bl_graph <- TRUE

    # random perturbation adjustments.
    fl_rand_adj_A_prop <- 0.1
    it_rand_adj_A_rng_seed_start <- 123
    ar_it_perturb_n <- seq(1, 5, by=1)

  } else if (it_run_collection == 12) {

    # Main run single file, 100 different tries
    ls_st_file_suffix <- c(ls_st_file_suffix_single_main)
    ar_double_triple_alloc <- ar_double_triple_alloc_set3

    # random perturbation adjustments.
    fl_rand_adj_A_prop <- 0.1

    # First 100 graph all
    it_rand_adj_A_rng_seed_start <- 123
    ar_it_perturb_n <- seq(1, 100, by=1)
    bl_graph <- TRUE

  } else if (it_run_collection == 13) {

    # Main run single file, 100 different tries
    ls_st_file_suffix <- c(ls_st_file_suffix_single_main)
    ar_double_triple_alloc <- ar_double_triple_alloc_set3

    # random perturbation adjustments.
    fl_rand_adj_A_prop <- 0.1

    # Run 400 more simulations, without graphs
    it_rand_adj_A_rng_seed_start <- 223
    ar_it_perturb_n <- seq(1, 400, by=1)
    bl_graph <- FALSE

  }
}

# Run -----------------
for (bl_per_capita in ls_bl_per_capita) {
  print(paste0('bl_per_capita=', bl_per_capita))

  # foreach (it_perturb_ctr=ar_it_perturb_n) %dopar% {
  for (it_perturb_ctr in ar_it_perturb_n) {
    it_rand_adj_A_rng_seed <- it_rand_adj_A_rng_seed_start + it_perturb_ctr

    foreach (fl_rho=ar_rho) %dopar% {
    # for (fl_rho in ar_rho) {
      print(paste0('bl_per_capita=', bl_per_capita, ', fl_rho=', fl_rho))
      try(rm(df_plan_v_tilde_full))

      # foreach (st_which_solu=ls_st_file_suffix) %dopar% {
      for (st_which_solu in ls_st_file_suffix) {
        print(paste0('bl_per_capita=', bl_per_capita, ', fl_rho=', fl_rho, ', st_which_solu=', st_which_solu))

        # File Names, Paths Etc -------
        if (it_perturb_ctr == -1) {
          it_rand_adj_A_rng_seed <- NULL
        }
        ls_output <- fs_opti_support_202111(st_which_solu,
                                            bl_per_capita=bl_per_capita,
                                            fl_rho=fl_rho,
                                            snm_main=snm_main,
                                            it_rand_draw_seed=it_rand_adj_A_rng_seed,
                                            ls_st_gen_agg_folders=ls_st_gen_agg_folders)
        st_file_type_withspouse_shock <- ls_output$st_file_type_withspouse_shock
        snm_simu_csv_withspouse_shock <- ls_output$snm_simu_csv_withspouse_shock
        srt_simu_path <- ls_output$srt_simu_path
        bl_save_img <- ls_output$bl_save_img
        spt_img_save <- ls_output$spt_img_save
        srt_csv_path <- ls_output$srt_csv_path
        ar_rho <- ls_output$ar_rho
        # ar_rho <- c(fl_rho)
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
        fl_max_phaseout = 200000
        it_bin_dollar_before_phaseout = 10000
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
            read.csv(file.path(srt_simu_path, snm_simu_csv_withspouse_shock), header=FALSE)) %>%
            rename_all(~c(ar_svr_csv))

          # df_plan_v_tilde_full <- vroom(paste0(srt_simu_path, snm_simu_csv_withspouse_shock), col_names=ar_svr_csv)
          tm_csv <- proc.time() - tm_start_csv
          print(paste0('df_plan_v_tilde_full loaded, tm_csv:', tm_csv))
        } else {
          print('already loaded df_plan_v_tilde_full previously')
        }

        # # Per capita allocation means shifting the C level per-capita level c
        # if (bl_per_capita) {
        #   # spt_img_save <- paste0(spt_img_save, '_rhoneg1')
        #   # srt_csv_path <- paste0(srt_csv_path, '_rhoneg1')
        #   # When cares about inequality aversion, cares about per capita
        #   df_plan_v_tilde_full <- df_plan_v_tilde_full %>%
        #     mutate(ctilde = ctilde/(kids+marital+1))
        # }

        for (it_solu_type in ar_double_triple_alloc) {
          # foreach (it_solu_type=ar_double_triple_alloc) %dopar% {
          print(paste0('bl_per_capita=', bl_per_capita, ', fl_rho=', fl_rho, ', it_solu_type=', it_solu_type))

          for (it_age_type in c(1)) {
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

            for (it_solu_round in c(1)) {

              it_solu_group <- it_solu_round*100 + it_solu_type
              print(paste0('solution group: ', it_solu_group, '=it_solu_group started'))

              # Start timer
              tm_start_solu <- proc.time()

              # Base settings
              bl_threshold <- FALSE
              it_age_bins <- 1
              # ar_rho <- 1 - (10^(c(seq(-2,2, length.out=8))))
              # ar_rho <- unique(c(1,ar_rho))
              # ar_rho <- c(1)

              ar_rho_g47 <- c(fl_rho)
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
                it_max_checks <- 84
                spt_img_save_use <- paste0(spt_img_save, '_4400_', st_img_suf_age_ybin,'/')
                srt_csv_path_use <- paste0(srt_csv_path, '_4400_', st_img_suf_age_ybin,'/')
                st_solu_name_check_count <- '84 checks'
                st_image_name <- ''
              } else if ((it_solu_group >= 111 & it_solu_group <=115) |
                         (it_solu_group >= 211 & it_solu_group <=215)) {
                # x05 to x07 is allocating up to 200k, common bounds for all
                it_max_checks <- 168
                spt_img_save_use <- paste0(spt_img_save, '_17k_', st_img_suf_age_ybin,'/')
                srt_csv_path_use <- paste0(srt_csv_path, '_17k_', st_img_suf_age_ybin,'/')
                st_solu_name_check_count <- '168 checks'
                st_image_name <- '_17k'
              } else if ((it_solu_group >= 121 & it_solu_group <=125) |
                         (it_solu_group >= 221 & it_solu_group <=225)) {
                # same as actual allocation, so this is the same as threshold allocation but
                # now does threshold allocation with g4 and g47 as well.
                # 14*4+14*2
                bl_threshold <- TRUE
                it_max_checks <- 84
                it_check_headorspouse <- 14
                it_check_perkids <- 14
                spt_img_save_use <- paste0(spt_img_save, '_14ca14ck_', st_img_suf_age_ybin,'/')
                srt_csv_path_use <- paste0(srt_csv_path, '_14ca14ck_', st_img_suf_age_ybin,'/')
                st_solu_name_check_count <- 'actual checks adults'
                st_image_name <- '_14ca14ck'
              } else if ((it_solu_group >= 131 & it_solu_group <=135) |
                         (it_solu_group >= 231 & it_solu_group <=235)) {
                # 17*4+14*2
                bl_threshold <- TRUE
                it_max_checks <- 96
                it_check_headorspouse <- 14
                it_check_perkids <- 17
                spt_img_save_use <- paste0(spt_img_save, '_17ckid_', st_img_suf_age_ybin,'/')
                srt_csv_path_use <- paste0(srt_csv_path, '_17ckid_', st_img_suf_age_ybin,'/')
                st_solu_name_check_count <- '14 checks adult 17 checks kids'
                st_image_name <- '_17ckid'
              } else if ((it_solu_group >= 141 & it_solu_group <=145) |
                         (it_solu_group >= 241 & it_solu_group <=245)) {
                # adult 18 checks
                # four kids married = 17*2 + 14*4
                bl_threshold <- TRUE
                it_max_checks <- 90
                it_check_headorspouse <- 17
                it_check_perkids <- 14
                spt_img_save_use <- paste0(spt_img_save, '_17cadt_', st_img_suf_age_ybin,'/')
                srt_csv_path_use <- paste0(srt_csv_path, '_17cadt_', st_img_suf_age_ybin,'/')
                st_solu_name_check_count <- '17 checks adult 14 checks kids'
                st_image_name <- '_17cadt'
              } else if ((it_solu_group >= 151 & it_solu_group <=155) |
                         (it_solu_group >= 251 & it_solu_group <=255)) {
                # double adult allocation from 12 to 24, kids 5 to 10
                # four kids married = 17*2 + 17*4
                bl_threshold <- TRUE
                it_max_checks <- 102
                it_check_headorspouse <- 17
                it_check_perkids <- 17
                spt_img_save_use <- paste0(spt_img_save, '_17ca17ck_', st_img_suf_age_ybin,'/')
                srt_csv_path_use <- paste0(srt_csv_path, '_17ca17ck_', st_img_suf_age_ybin,'/')
                st_solu_name_check_count <- '17 checks adult 17 checks kids'
                st_image_name <- '_17ca17ck'
              } else if ((it_solu_group >= 161 & it_solu_group <=165) |
                         (it_solu_group >= 261 & it_solu_group <=265)) {
                # 4*20 +14*2
                bl_threshold <- TRUE
                it_max_checks <- 108
                it_check_headorspouse <- 14
                it_check_perkids <- 20
                spt_img_save_use <- paste0(spt_img_save, '_20ckid_', st_img_suf_age_ybin,'/')
                srt_csv_path_use <- paste0(srt_csv_path, '_20ckid_', st_img_suf_age_ybin,'/')
                st_solu_name_check_count <- '14 checks adult 20 checks kids'
                st_image_name <- '_20ckid'
              } else if ((it_solu_group >= 171 & it_solu_group <=175) |
                         (it_solu_group >= 271 & it_solu_group <=275)) {
                # 4*14 +20*2
                bl_threshold <- TRUE
                it_max_checks <- 96
                it_check_headorspouse <- 20
                it_check_perkids <- 14
                spt_img_save_use <- paste0(spt_img_save, '_20cadt_', st_img_suf_age_ybin,'/')
                srt_csv_path_use <- paste0(srt_csv_path, '_20cadt_', st_img_suf_age_ybin,'/')
                st_solu_name_check_count <- '20 checks adult 14 checks kids'
                st_image_name <- '_20cadt'
              } else if ((it_solu_group >= 181 & it_solu_group <=185) |
                         (it_solu_group >= 281 & it_solu_group <=285)) {
                # triple adult allocation from 12 to 36, kids 5 to 15
                # four kids married = 20*2 + 20*4 = 132
                bl_threshold <- TRUE
                it_max_checks <- 120
                it_check_headorspouse <- 20
                it_check_perkids <- 20
                spt_img_save_use <- paste0(spt_img_save, '_20ca20ck_', st_img_suf_age_ybin,'/')
                srt_csv_path_use <- paste0(srt_csv_path, '_20ca20ck_', st_img_suf_age_ybin,'/')
                st_solu_name_check_count <- '20 checks adult 20 checks kids'
                st_image_name <- '_20ca20ck'
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
                it_check_headorspouse <- 14
                it_check_perkids <- 14
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
                  bl_per_capita = bl_per_capita,
                  it_max_age = it_max_age,
                  it_min_age = it_min_age,
                  it_age_bins = it_age_bins,
                  fl_rand_adj_A_prop = fl_rand_adj_A_prop,
                  it_rand_adj_A_rng_seed = it_rand_adj_A_rng_seed,
                  ar_svr_csv = ar_svr_csv,
                  ar_svr_groups = ar_svr_groups,
                  ar_svr_groups_stats = ar_svr_groups_stats,
                  svr_checks = svr_checks,
                  svr_v_value = svr_v_value,
                  svr_c_value = svr_c_value,
                  svr_mass = svr_mass,
                  ar_rho = ar_rho,
                  bl_threshold = bl_threshold,
                  it_check_headorspouse = it_check_headorspouse,
                  it_check_perkids = it_check_perkids,
                  bl_non_inc_adjust = TRUE,
                  bl_print = FALSE,
                  bl_print_verbose = FALSE)

                print(paste0(it_solu_group, ': First Round Allocation Finished'))

                # First Round Solution Results (Util Max)
                # The allocation Results here are based on utility maximization
                # How to allocate given resources available
                tb_rho_rev_c <- ls_prc_outputs_zs5_1st$tb_rho_rev_c
                tb_rho_rev_v <- ls_prc_outputs_zs5_1st$tb_rho_rev_v

                df_rho_gini_c <- ls_prc_outputs_zs5_1st$df_rho_gini_c
                df_rho_gini_v <- ls_prc_outputs_zs5_1st$df_rho_gini_v

                df_input_il_noninc_covar_zs5_1st=ls_prc_outputs_zs5_1st$df_input_il_noninc_covar
                df_alloc_i_long_covar_c_zs5_1st=ls_prc_outputs_zs5_1st$df_alloc_i_long_covar_c
                df_alloc_i_long_covar_v_zs5_1st=ls_prc_outputs_zs5_1st$df_alloc_i_long_covar_v

                df_queue_il_long_c <- ls_prc_outputs_zs5_1st$df_queue_il_long_c
                if (bl_save_full_queue) {
                  df_queue_il_long_c_full <- df_queue_il_long_c %>%
                    left_join(df_alloc_i_long_covar_c_zs5_1st %>%
                                select(id_i, marital, kids, ymin_group, mass), by='id_i')
                  write.csv(df_queue_il_long_c_full,
                            paste0(srt_csv_path, '/df_queue_il_long_c.csv'),
                            row.names = TRUE)
                }

                stg_subtitle=ls_prc_outputs_zs5_1st$stg_subtitle
                stg_caption=ls_prc_outputs_zs5_1st$stg_caption

                # First Round Solution Results (Exp Min)
                # The allocation results here are based on exp. minimization
                # The results rely on the same queue recorded in df_input_il_noninc_covar_zs5_1st
                # Once the queue has been solved for, can go back and forward to any point along
                # the queue. The REV file

                # EM1. Get Alternative Outcome Levels for C, scalar for utilitarian
                fl_utilitarian_AlterOutcome <- as.numeric(tb_rho_rev_c %>% filter(rho_val == ar_rho[1]) %>% pull(AlterOutcome))

                # EM2. Sbuset queue frame and count up to alteroutcome
                # df_alloc_id_expmin: id and checks columns, checks from expenditure minimization
                # df_alloc_id_expmin merged with allocation table to generate APC given exp. minimization
                df_alloc_id_expmin_1st <- df_queue_il_long_c %>% filter(rho_val == ar_rho[1]) %>%
                  mutate(D_Wbin_il_expmin :=
                           case_when(V_star_Q_il <= fl_utilitarian_AlterOutcome ~ 1,
                                     TRUE ~ 0)) %>%
                  group_by(id_i) %>%
                  summarize(checks := sum(D_Wbin_il_expmin)) %>%
                  ungroup() %>% arrange(id_i)

                # Graphing
                if (bl_graph) {

                  ls_it_kids_ct <- c(2, 4)
                  for (it_kids_ct in ls_it_kids_ct) {

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
                      it_kids_ct = it_kids_ct,
                      spt_img_save=spt_img_save_use)

                    # Generate Graphs
                    # if ((it_solu_group >= 101 & it_solu_group <=104) | (it_solu_group >= 201 & it_solu_group <=204)) {
                    ls_st_checks_cv <- c('checks_cjv', 'checks_c', 'checks_v')
                    # } else {
                    # ls_st_checks_cv <- c('checks_cjv', 'checks_c', 'checks_v')
                    # ls_st_checks_cv <- c('checks_cjv', 'checks_cjv_log', 'checks_c', 'checks_c_log', 'checks_v', 'checks_v_log')
                    # }
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
                      it_kids_ct = it_kids_ct,
                      spt_img_save=spt_img_save_use)
                  }
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
              ls_prc_outputs_zs5_1st$df_alloc_i_long_covar_c %>%
                filter(rho_val == ar_rho[1]) %>%
                filter(allocate_type == 'optimal') %>%
                mutate(allocate_type = case_when(allocate_type == 'optimal' ~ 'optiexpmin_c_1st',
                                                 TRUE ~ allocate_type)) %>%
                select(-checks) %>% left_join(df_alloc_id_expmin_1st, by='id_i')) %>%
              arrange(id_i, allocate_type) %>%
              select(-rho, -F_star_i, -EH_star_i, -survive) %>%
              select(id_i, rho_val, everything())

            df_alloc_all <- df_alloc_all %>%
              pivot_wider(names_from = allocate_type,
                          values_from = checks)

            # Stack three things together:
            # 1. select from df_input_il_noninc_covar_zs5_2nd: C without allocations
            # 2. merge df_alloc_i_long_covar_c_zs5_2nd if actual: C with actual allocations
            # 3. merge df_alloc_i_long_covar_c_zs5_2nd if optimal: C with optimal allocations
            # 4. select from (1),(2) and (3) only id and C column, merge with df_alloc_all

            # 1. select from df_input_il_noninc_covar_zs5_2nd: C without allocations
            df_c_no_allocation_1st <- df_input_il_noninc_covar_zs5_1st %>% ungroup() %>%
              filter(D_il==1) %>%
              mutate(optimal_c_1st_c_no_allocation = c_A_il) %>%
              select(id_i, optimal_c_1st_c_no_allocation)

            # 2. merge df_alloc_i_long_covar_c_zs5_2nd if actual: C with actual allocations
            # if actual = 0, use c_A_il when D_il == 1
            # if actual > 0, use c_A_il + c_alpha_il when D_il == 1
            df_c_allocation_actual_1st <- df_input_il_noninc_covar_zs5_1st %>% ungroup() %>%
              mutate(checks_il = checks) %>%
              select(id_i, checks_il, c_A_il, c_alpha_il) %>%
              left_join(df_alloc_i_long_covar_c_zs5_1st %>%
                          filter(rho_val == ar_rho[1]) %>%
                          filter(allocate_type == 'actual') %>%
                          select(id_i, checks), by = 'id_i') %>%
              mutate(optimal_c_1st_c_actual = case_when(
                checks==checks_il ~ c_A_il + c_alpha_il,
                checks==0 & checks_il == 1 ~ c_A_il)) %>%
              filter(!is.na(optimal_c_1st_c_actual)) %>%
              select(id_i, optimal_c_1st_c_actual)

            # 3. merge df_alloc_i_long_covar_c_zs5_2nd if optimal: C with optimal allocations
            # same precedure as above
            df_c_allocation_optimal_1st <- df_input_il_noninc_covar_zs5_1st %>% ungroup() %>%
              mutate(checks_il = checks) %>%
              select(id_i, checks_il, c_A_il, c_alpha_il) %>%
              left_join(df_alloc_i_long_covar_c_zs5_1st %>%
                          filter(rho_val == ar_rho[1]) %>%
                          filter(allocate_type == 'optimal') %>%
                          select(id_i, checks), by = 'id_i') %>%
              mutate(optimal_c_1st_c_optimal = case_when(
                checks==checks_il ~ c_A_il + c_alpha_il,
                checks==0 & checks_il == 1 ~ c_A_il)) %>%
              filter(!is.na(optimal_c_1st_c_optimal)) %>%
              select(id_i, optimal_c_1st_c_optimal)

            # optimal allocation c but for Expenditure Minimization 1st Round only
            # df_alloc_id_expmin_1st has ID and checks columns
            df_c_allocation_optiexpmin_1st <- df_input_il_noninc_covar_zs5_1st %>% ungroup() %>%
              mutate(checks_il = checks) %>%
              select(id_i, checks_il, c_A_il, c_alpha_il) %>%
              left_join(df_alloc_id_expmin_1st, by = 'id_i') %>%
              mutate(optimal_c_1st_c_optiexpmin = case_when(
                checks==checks_il ~ c_A_il + c_alpha_il,
                checks==0 & checks_il == 1 ~ c_A_il)) %>%
              filter(!is.na(optimal_c_1st_c_optiexpmin)) %>%
              select(id_i, optimal_c_1st_c_optiexpmin)

            # 4. Merge to oMain file
            df_alloc_all <- df_alloc_all %>%
              left_join(df_c_no_allocation_1st, by='id_i') %>%
              left_join(df_c_allocation_actual_1st, by='id_i') %>%
              left_join(df_c_allocation_optimal_1st, by='id_i') %>%
              left_join(df_c_allocation_optiexpmin_1st, by='id_i') %>%
              mutate(apc_actual_1st = (optimal_c_1st_c_actual-optimal_c_1st_c_no_allocation)/(actual*fl_percheck_dollar/fl_multiple),
                     apc_optimal_1st = (optimal_c_1st_c_optimal-optimal_c_1st_c_no_allocation)/(optimal_c_1st*fl_percheck_dollar/fl_multiple),
                     apc_optiexpmin_1st = (optimal_c_1st_c_optiexpmin-optimal_c_1st_c_no_allocation)/(optiexpmin_c_1st*fl_percheck_dollar/fl_multiple),
                     apc_gainratio_1st = apc_optimal_1st/apc_actual_1st)

            # Export Allocation for Every Cell -----------
            # For each cell under consideration, actual allocation and optimal Round1, Round 2 Allocations
            # Also store, EC without allocation, EC with allocation actual, EC with allocation optimal
            write.csv(df_alloc_all,
                      paste0(srt_csv_path_use, 'df_alloc_all_',st_suffix_csvimg,'.csv'),
                      row.names = TRUE)


            # CSV Export Grouped Average Statistics -------
            if (st_which_solu != "b1_manna_expmin") {
              ls_svr_groups <- ar_svr_groups
              for (svr_group in ls_svr_groups) {

                # group mean
                df_alloc_combine_group_mean <- df_alloc_all %>%
                  ungroup() %>% group_by(!!sym(svr_group)) %>%
                  summarize(actual_mean = sum(actual*mass)/sum(mass),
                            optimal_c_1st_mean = sum(optimal_c_1st*mass)/sum(mass),
                            optimal_v_1st_mean = sum(optimal_v_1st*mass)/sum(mass),
                            optimal_c_1st_c_no_allocation_mean = sum(optimal_c_1st_c_no_allocation*mass)/sum(mass),
                            optimal_c_1st_c_actual_mean = sum(optimal_c_1st_c_actual*mass)/sum(mass),
                            optimal_c_1st_c_optimal_mean = sum(optimal_c_1st_c_optimal*mass)/sum(mass),
                            apc_actual_1st = (optimal_c_1st_c_actual_mean-optimal_c_1st_c_no_allocation_mean)/(actual_mean*fl_percheck_dollar/fl_multiple),
                            apc_optimal_1st = (optimal_c_1st_c_optimal_mean-optimal_c_1st_c_no_allocation_mean)/(optimal_c_1st_mean*fl_percheck_dollar/fl_multiple),
                            apc_gainratio_1st = apc_optimal_1st/apc_actual_1st,
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
                            optimal_c_1st_c_no_allocation_mean = sum(optimal_c_1st_c_no_allocation*mass)/sum(mass),
                            optimal_c_1st_c_actual_mean = sum(optimal_c_1st_c_actual*mass)/sum(mass),
                            optimal_c_1st_c_optimal_mean = sum(optimal_c_1st_c_optimal*mass)/sum(mass),
                            apc_actual_1st = (optimal_c_1st_c_actual_mean-optimal_c_1st_c_no_allocation_mean)/(actual_mean*fl_percheck_dollar/fl_multiple),
                            apc_optimal_1st = (optimal_c_1st_c_optimal_mean-optimal_c_1st_c_no_allocation_mean)/(optimal_c_1st_mean*fl_percheck_dollar/fl_multiple),
                            apc_gainratio_1st = apc_optimal_1st/apc_actual_1st,
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
                            optimal_c_1st_c_no_allocation_mean = sum(optimal_c_1st_c_no_allocation*mass)/sum(mass),
                            optimal_c_1st_c_actual_mean = sum(optimal_c_1st_c_actual*mass)/sum(mass),
                            optimal_c_1st_c_optimal_mean = sum(optimal_c_1st_c_optimal*mass)/sum(mass),
                            apc_actual_1st = (optimal_c_1st_c_actual_mean-optimal_c_1st_c_no_allocation_mean)/(actual_mean*fl_percheck_dollar/fl_multiple),
                            apc_optimal_1st = (optimal_c_1st_c_optimal_mean-optimal_c_1st_c_no_allocation_mean)/(optimal_c_1st_mean*fl_percheck_dollar/fl_multiple),
                            apc_gainratio_1st = apc_optimal_1st/apc_actual_1st,
                            mass_sum = sum(mass))

                # Export
                write.csv(df_alloc_combine_group_mean_oneless,
                          paste0(srt_csv_path_use, "df_alloc_",svr_group,"_without_noage_", st_suffix_csvimg, ".csv"),
                          row.names = TRUE)

              }

              # for G4, also generate income and age grouping
              if (ar_solu_group_digits[3] == 3) {

                # All but Group mean
                df_alloc_combine_group_ymin_g4age <- df_alloc_all %>%
                  ungroup() %>% group_by(!!!syms(ar_svr_groups_age_ymin)) %>%
                  summarize(actual_mean = sum(actual*mass)/sum(mass),
                            optimal_c_1st_mean = sum(optimal_c_1st*mass)/sum(mass),
                            optimal_v_1st_mean = sum(optimal_v_1st*mass)/sum(mass),
                            optimal_c_1st_c_no_allocation_mean = sum(optimal_c_1st_c_no_allocation*mass)/sum(mass),
                            optimal_c_1st_c_actual_mean = sum(optimal_c_1st_c_actual*mass)/sum(mass),
                            optimal_c_1st_c_optimal_mean = sum(optimal_c_1st_c_optimal*mass)/sum(mass),
                            apc_actual_1st = (optimal_c_1st_c_actual_mean-optimal_c_1st_c_no_allocation_mean)/(actual_mean*fl_percheck_dollar/fl_multiple),
                            apc_optimal_1st = (optimal_c_1st_c_optimal_mean-optimal_c_1st_c_no_allocation_mean)/(optimal_c_1st_mean*fl_percheck_dollar/fl_multiple),
                            apc_gainratio_1st = apc_optimal_1st/apc_actual_1st,
                            mass_sum = sum(mass))

                # Export
                write.csv(df_alloc_combine_group_ymin_g4age,
                          paste0(srt_csv_path_use, "df_alloc_age_group_ymin_group_", st_suffix_csvimg, ".csv"),
                          row.names = TRUE)


              }

              ### CSV Print and save MASS REV results ---------------
              # Save REV to table, Stack them
              tb_rho_rev_mass_v1_tab <- ls_prc_outputs_zs5_1st$tb_rho_rev_v %>%
                mutate(objective = 'vlife',
                       constraint = st_tfo_method,
                       allocround = 'first',
                       maxchecks = it_max_checks) %>%
                left_join(ls_prc_outputs_zs5_1st$tb_rho_vstar_v %>%
                            select(rho_val, V_star_resource), by='rho_val') %>%
                left_join(ls_prc_outputs_zs5_1st$df_rho_gini_v %>%
                            select(rho_val,
                                   income_bound,
                                   gini_drv_EH_d0, gini_drv_EH_dact, gini_drv_EH_star,
                                   sd_log_EH_drv_d0,sd_log_EH_drv_dact,sd_log_EH_drv_star,
                                   atkinson_drv_EH_d0, atkinson_drv_EH_dact, atkinson_drv_EH_star,
                                   mean_EH_drv_d0, mean_EH_drv_dact, mean_drv_EH_star,
                                   min_EH_d0, min_EH_dact, min_EH_star,
                                   gini_drv_D_star,
                                   sd_D_star,
                                   mean_EH_star, sd_EH_star), by='rho_val')


              tb_rho_rev_mass_c1_tab <- ls_prc_outputs_zs5_1st$tb_rho_rev_c %>%
                mutate(objective = 'c2020',
                       constraint = st_tfo_method,
                       allocround = 'first',
                       maxchecks = it_max_checks) %>%
                left_join(ls_prc_outputs_zs5_1st$tb_rho_vstar_c %>%
                            select(rho_val, V_star_resource), by='rho_val') %>%
                left_join(ls_prc_outputs_zs5_1st$df_rho_gini_c %>%
                            select(rho_val,
                                   income_bound,
                                   gini_drv_EH_d0, gini_drv_EH_dact, gini_drv_EH_star,
                                   sd_log_EH_drv_d0,sd_log_EH_drv_dact,sd_log_EH_drv_star,
                                   atkinson_drv_EH_d0, atkinson_drv_EH_dact, atkinson_drv_EH_star,
                                   mean_EH_drv_d0, mean_EH_drv_dact, mean_drv_EH_star,
                                   min_EH_d0, min_EH_dact, min_EH_star,
                                   gini_drv_D_star,
                                   sd_D_star,
                                   mean_EH_star, sd_EH_star), by='rho_val')

              # Stack frames
              tb_rho_rev_mass_v_c_tab <- rbind(tb_rho_rev_mass_v1_tab, tb_rho_rev_mass_c1_tab)

              # Compute AAPC Aggregate Average Propensity to Consume -----------------
              # AAPC_Actual = (Increase in Aggregate Consumption)/(Total Checks Amount)
              # AAPC^{j}_Optimal = (Increase in Aggregate Consumption under Optimal Policy)/(Total Checks Amount)
              # Suppose AAPC_Actual = 0.10, that means overall 10 percent of all checks was spent, and 90 percent saved under actual allocation
              # Suppose AAPC_Optimal = 0.20, that means overall 20 percent of all checks were spent under optimal reallocation of the same total checks.
              # In this example, AAPC_Optimal/AAPC_Actual = 2, meaning optimal policy is able to double the consumption/stimulus effects with the same budget. The ratio AAPC_Optimal/AAPC_Actual is kind of the dual number of REV. REV is proportionally how much is saved with optimal allocation, AAPC_Optimal/AAPC_Actual is proportionally how much could be gained with optimal allocation.

              # Total Average Checks under Actual Policy
              fl_avg_checks_actual <- df_alloc_all %>% ungroup() %>%summarize(actual_mean = sum(actual*mass)/sum(mass))
              fl_checks_total_normalize <- as.numeric(fl_avg_checks_actual*fl_percheck_dollar/fl_multiple)

              # aggregation Consumption
              # c_aggregate and AlterOutcome are in the same units both need to be divided by mass
              # NoAllocOutcome is weighted
              fl_c_no_allocation_1st <- as.numeric(df_input_il_noninc_covar_zs5_1st %>% filter(D_il==1) %>% ungroup() %>% summarize(c_sum = sum(c_A_il*mass)/sum(mass)))

              # fl_c_no_allocation_unw <- as.numeric(df_input_il_noninc_covar_zs5_1st %>% filter(D_il==1) %>% ungroup() %>% summarize(c_sum = sum(c_A_il*mass)))
              fl_mass_sum <- as.numeric(df_input_il_noninc_covar_zs5_1st %>% filter(D_il==1)  %>% ungroup() %>% summarize(mass_sum = sum(mass)))
              tb_rho_rev_mass_v_c_tab <- tb_rho_rev_mass_v_c_tab %>%
                # mutate(c_aggregate = case_when( objective == 'c2020' & rho_val == 1 ~ V_star_resource)) %>%
                # mutate(NoAllocOutcome = fl_c_no_allocation_1st) %>%
                # mutate(GainRelaAlterOutcomeFrac = (c_aggregate-AlterOutcome)/AlterOutcome) %>%
                # mutate(GainRelaNoAllocOutcomeFrac = (c_aggregate/fl_mass_sum-NoAllocOutcome)/(NoAllocOutcome)) %>%
                mutate(AAPC_actual = case_when(
                  objective == 'c2020' & rho_val == ar_rho[1] & allocround == 'first' ~
                    (((AlterOutcome/fl_mass_sum-fl_c_no_allocation_1st)/fl_checks_total_normalize)))
                ) %>%
                mutate(AAPC_optimal = case_when(
                  objective == 'c2020' & rho_val == ar_rho[1] & allocround == 'first' ~
                    (((V_star_resource/fl_mass_sum-fl_c_no_allocation_1st)/fl_checks_total_normalize)))
                ) %>%
                mutate(AAPC_ratio = AAPC_optimal/AAPC_actual)

              # export
              spn_rev_csv_save <- paste0(srt_csv_path_use, "rev_", st_suffix_csvimg, ".csv")
              write.csv(tb_rho_rev_mass_v_c_tab, spn_rev_csv_save, row.names = TRUE)

              # CSV Expected APC Average Propensity to Consume Distribution -----------------
              # Expected consumption in 2020 (from 2019 perspective) given checks optimally allocated under optimal allocation policy j to individual i minus expected consumption without checks divided by total checks received. Only calculate this if i receives positive checks under j.
              # For AAPC, don't need to worry about expectation vs realization, since it is integrated aggregate. But perhaps important to keep in mind that our MPC and APC are all from 2019 perspective, they are expectations given COVID shocks etc. The checks are actually spent in 2020, but we don't get to condition check allocation by COVID realization or 2020 productivity realization.
              # We generate distribution from:  APC^{j}_{i} and MASS_{i}, save only these here
              # Store Results here that will be processed by Matlab MEconTools which can process Discrete Random Variables Better
              print(paste0('spn_rev_csv_save:', spn_rev_csv_save))
            }
          }
        }
      }
    }
  }
}
stopCluster(cl)
