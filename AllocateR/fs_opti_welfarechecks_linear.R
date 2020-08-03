ff_sup_clean_rmd <- function(srt_simu_path = 'C:/Users/fan/Documents/Dropbox (UH-ECON)/PrjNygaardSorensenWang/Output/',
                             st_file_type = 'moredense_a100zh266_e2m2',
                             fl_max_phaseout = 238000,
                             it_bin_dollar_before_phaseout = 2500,
                             fl_percheck_dollar = 100,
                             fl_multiple = 58056,
                             it_max_checks = 44,
                             fl_tax_hh = 128580000,
                             it_max_age = 64,
                             it_min_age = 18,
                             ar_svr_csv = c('age', 'marital', 'kids', 'checks',	'ymin', 'mass', 'survive', 'vtilde', 'ctilde'),
                             ar_svr_groups = c('marital', 'kids', 'age_group', 'ymin_group'),
                             ar_svr_groups_stats = c('mass', 'survive'),
                             svr_checks = 'checks',
                             svr_v_value = 'vtilde',
                             svr_c_value = 'ctilde',
                             bl_non_inc_adjust = FALSE,
                             bl_print = TRUE,
                             bl_print_verbose = FALSE) {
  #' This function cleans rmd: creates subfolder with pdf, html and R

  ## ----setup, include=FALSE------------------------------------------------------------------------------------------------------------------------------------
  # rm(list = ls())

  # knitr::opts_chunk$set(echo = TRUE)
  # knitr::opts_chunk$set(fig.width=12, fig.height=8)

  # library(tidyverse)
  # library(REconTools)

  # ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  # # Max Phase Out given 1200*2 + 500*4 = 4400
  # fl_max_phaseout = 238000
  # it_bin_dollar_before_phaseout = 2500
  # # Dollar Per Check
  # fl_percheck_dollar = 100
  # # Meaning of Ymin Ymax simulated interval of 1
  # fl_multiple = 58056
  # # Number of Max Checks
  # it_max_checks = 44
  # # Number of Tax Paying Households
  # fl_tax_hh = 128580000
  # # Number of Income Groups to Use: use 25 for 10,000 = 1
  # # Age Conditions
  # # it_max_age = 64
  # # it_min_age = 64
  # it_max_age = 64
  # it_min_age = 18


  # ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  # # File Path
  # srt_simu_path <- 'C:/Users/fan/Documents/Dropbox (UH-ECON)/PrjNygaardSorensenWang/Output/'
  # # File Name
  # # snm_simu_csv <- 'snwx_v_planner_small.csv'
  # # snm_simu_csv <- 'snwx_v_planner_small_dup.csv'
  # # snm_simu_csv <- 'snwx_v_planner_base.csv'
  # # snm_simu_csv <- 'snwx_v_planner_dense.csv'
  # # snm_simu_csv <- 'snwx_v_planner_densemore.csv'
  # # st_file_type <- 'small_dup'
  # # st_file_type <- 'dense'
  # # st_file_type <- 'moredense'
  # # st_file_type <- 'densemore_a55z133'
  # # st_file_type <- 'densemore_a55z266_e0m0'
  # # st_file_type <- 'moredense_a100z266_e0m0'
  # # st_file_type <- 'moredense_a55zh43zs11'
  # # st_file_type <- 'moredense_a75zh101zs5'
  #
  # # 1. e0m2 very dense a and zh test
  # st_file_type <- 'moredense_a100zh266_e1m1'
  #
  # # 2. e1m1 very dense a and zh test, both married and education, no spouse shock
  # st_file_type <- 'moredense_a100zh266_e2m2'
  #
  # # CSV Name
  # snm_simu_csv <- paste0('snwx_v_planner_',st_file_type,'.csv')
  #
  # # Column Names
  # ar_svr_csv <- c('age', 'marital', 'kids', 'checks',	'ymin', 'mass', 'survive', 'vtilde', 'ctilde')
  # # Variables That Identify Individual Types
  # ar_svr_groups <- c('marital', 'kids', 'age_group', 'ymin_group')
  # ar_svr_groups_stats <- c('mass', 'survive')
  # # Number of Checks and Planner Value
  # svr_checks <- 'checks'
  # svr_v_value <- 'vtilde'
  # svr_c_value <- 'ctilde'


  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  mt_plan_v_tilde <- read.csv(paste0(srt_simu_path, snm_simu_csv), header=FALSE)
  df_plan_v_tilde <- as_tibble(mt_plan_v_tilde) %>%
    rename_all(~c(ar_svr_csv)) %>%
    filter(vtilde != 0) %>%
    filter(checks <= it_max_checks) %>%
    filter(age <= it_max_age) %>%
    filter(age >= it_min_age)

  # Remove
  rm(mt_plan_v_tilde)

  # Column 1: Age (in year before COVID)
  # Column 2: Marital status (0 if not married; 1 if married)
  # Column 3: Nr of kids (0, 1, ..., 5) where 5 means 5 or more
  # Column 4: Number of welfare checks (here either equal to 0 or 1)
  # Column 5 and column 6 give income range
  # So the individual's income is at least as large as the value in column 5 but strictly less than the value in column 6
  # Column 7: Population weight Of that particular group (in the stationary distribution)
  # Column 8: Survival probability of that particular age (since the planner knows that some of the individuals will die before next period, so wasn't sure how you wanted me to include that. I did not already include it in V^tilde)
  # Column 9: Value of planner as in the slides (with the exception that I didn't multiply by the survival probability


  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  if (bl_print_verbose){
    REconTools::ff_summ_percentiles(df_plan_v_tilde)
  }

  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  if (bl_print){
    print(paste0('sum(df_plan_v_tilde %>% pull(mass)=', sum(df_plan_v_tilde %>% pull(mass))))
  }


  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  # Age Groups
  ar_agecut = seq(it_min_age-1, it_max_age, length.out=5)
  # Dimensions
  if (bl_print){
    print(paste0('dim(df_plan_v_tilde)=',dim(df_plan_v_tilde)))
  }

  # df_plan_v_tilde_yjm <- df_plan_v_tilde
  df_plan_v_tilde_yjm <- df_plan_v_tilde %>%
    mutate(age_group = (cut(age, ar_agecut))) %>%
    group_by(marital, kids, checks, ymin, age_group) %>%
    summarize(vtilde = sum(vtilde*mass)/sum(mass),
              ctilde = sum(ctilde*mass)/sum(mass),
              mass = sum(mass),
              survive = mean(survive)) %>%
    ungroup()

  # Remove
  rm(df_plan_v_tilde)

  # Summarize
  if (bl_print){
    print(paste0('dim(df_plan_v_tilde_yjm)=',dim(df_plan_v_tilde_yjm)))
  }
  if (bl_print_verbose){
    REconTools::ff_summ_percentiles(df_plan_v_tilde_yjm)
  }

  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  # Number of Cuts Before th Mas Phase Out Point.
  fl_thres = fl_max_phaseout/fl_multiple
  # Solution Grid Precision prior to max_phaseout point.
  inc_grid1 = seq(0, fl_thres, length.out=(fl_max_phaseout)/it_bin_dollar_before_phaseout)
  fl_grid_gap = inc_grid1[2] - inc_grid1[1]
  # Not all inc_grid1 points are valid, there are some elements that have no mass
  fl_min_ymin_posmass = min(df_plan_v_tilde_yjm %>% pull(ymin))
  # First group is min
  ar_ycut = c(0, inc_grid1[inc_grid1>fl_min_ymin_posmass+fl_grid_gap], 7)

  # alternative method
  # fl_max_full_phaseout = fl_max_phaseout/fl_multiple
  # ar_ycut2 = c(0, seq(
  #   min(df_plan_v_tilde_yjm %>% pull(ymin)),
  #   fl_max_full_phaseout,
  #   length.out=(it_inc_groups - 1))[2:(it_inc_groups - 1)], 7)

  # ar_ycut = c(0, seq(
  #   min(df_plan_v_tilde_yjm %>% pull(ymin)),
  #   7,
  #   length.out=(it_inc_groups - 1))[2:(it_inc_groups - 1)])

  if (bl_print){
    print(paste0('ar_ycut*fl_multiple:',ar_ycut*fl_multiple))
  }
  it_inc_groups = length(ar_ycut)
  fl_inc_gap = (ar_ycut[3]-ar_ycut[2])*fl_multiple
  fl_inc_min = min(df_plan_v_tilde_yjm %>% pull(ymin))*fl_multiple
  subtitle = paste0('1 unit along x-axis = $', round(fl_inc_gap),
                    ', x-axis min = $', round(fl_inc_min),
                    ', x-axis final group >= $', round(ar_ycut[length(ar_ycut)-1]*fl_multiple))
  if (bl_print){
    print(paste0('subtitle:',subtitle))
  }

  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  # df_plan_v_tilde_ygrpjm %>% filter(kids==0 & checks == 1)
  # table(df_plan_v_tilde_ygrpjm$ymin_group)
  # df_plan_v_tilde_ygrpjm <- df_plan_v_tilde_yjm %>% mutate(ymin_group = as.factor(ymin))
  df_plan_v_tilde_ygrpjm <- df_plan_v_tilde_yjm %>%
    mutate(ymin_group = (cut(ymin, ar_ycut))) %>%
    group_by(marital, kids, checks, age_group, ymin_group) %>%
    summarize(vtilde = sum(vtilde*mass)/sum(mass),
              ctilde = sum(ctilde*mass)/sum(mass),
              mass = sum(mass),
              survive = mean(survive)) %>%
    ungroup()

  # Remove
  rm(df_plan_v_tilde_yjm)

  # Summarize
  if (bl_print_verbose){
    REconTools::ff_summ_percentiles(df_plan_v_tilde_ygrpjm)
  }

  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  if (bl_print){
    print(paste0('sum(df_plan_v_tilde_ygrpjm %>% pull(mass))',sum(df_plan_v_tilde_ygrpjm %>% pull(mass))))
  }

  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  # group id
  svr_group_id <- 'group_id'
  # Define
  ls_svr_group_vars <- ar_svr_groups
  # panel dataframe following
  df_plan_v_tilde_id <- df_plan_v_tilde_ygrpjm %>%
    arrange(!!!syms(ls_svr_group_vars)) %>%
    group_by(!!!syms(ls_svr_group_vars)) %>%
    mutate(!!sym(svr_group_id) := (row_number()==1)*1) %>%
    ungroup() %>%
    rowid_to_column(var = "id") %>%
    mutate(!!sym(svr_group_id) := cumsum(!!sym(svr_group_id))) %>%
    select(one_of(svr_group_id, ls_svr_group_vars), everything())

  if (bl_print_verbose){
    REconTools::ff_summ_percentiles(df_plan_v_tilde_id)
  }

  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  if (bl_print_verbose){
    REconTools::ff_summ_count_unique_by_groups(df_plan_v_tilde_id,ls_svr_group_vars,svr_group_id)
  }
  if (bl_print_verbose){
    REconTools::ff_summ_percentiles(df_plan_v_tilde_id, bl_statsasrows = FALSE)
  }

  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  # Select Grouping by Variables
  df_id <- df_plan_v_tilde_id %>%
    select(one_of(svr_group_id, ls_svr_group_vars, ar_svr_groups_stats)) %>%
    group_by(!!!syms(svr_group_id)) %>%
    slice_head() %>% ungroup() %>%
    select(one_of(svr_group_id, ls_svr_group_vars, ar_svr_groups_stats)) %>%
    rename(id_i = !!sym(svr_group_id))
  ar_group_ids <- unique(df_id %>% pull(id_i))

  # Summarize
  if (bl_print_verbose){
    REconTools::ff_summ_percentiles(df_id)
  }

  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  if (bl_print){
    print(paste0('sum(df_id %>% pull(mass))', sum(df_id %>% pull(mass))))
  }

  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  # Select 4 variables
  df_value <- df_plan_v_tilde_id %>%
    select(one_of(svr_group_id, svr_checks, svr_v_value, svr_c_value))
  # remove
  rm(df_plan_v_tilde_id)
  # Summarize
  if (bl_print_verbose){
    REconTools::ff_summ_percentiles(df_value)
    REconTools::ff_summ_count_unique_by_groups(df_value, svr_group_id, svr_group_id)
  }

  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  # 1. id column and id_il
  df_il <- df_value %>% rename(id_i = !!sym(svr_group_id)) %>%
    mutate(id_il = row_number()) %>%
    select(id_i, id_il, everything())
  # 2. D_max_i and D_il
  df_il <- df_il %>%
    arrange(id_i, svr_checks) %>% group_by(id_i) %>%
    mutate(D_max_i = max(!!sym(svr_checks))) %>%
    rename(D_il = !!sym(svr_checks)) %>%
    mutate(beta_i = 1/n()) %>%
    select(id_i, id_il, D_max_i, D_il, everything())
  # Summarize
  if (bl_print_verbose){
    REconTools::ff_summ_percentiles(df_il)
  }

  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  # 3. A_il and alpha_il
  df_il_U <- df_il %>%
    mutate(c_alpha_il = lead(!!sym(svr_c_value)) - (!!sym(svr_c_value)),
           v_alpha_il = lead(!!sym(svr_v_value)) - (!!sym(svr_v_value))) %>%
    rename(c_A_il = !!sym(svr_c_value),
           v_A_il = !!sym(svr_v_value)) %>%
    ungroup()

  # 4. drop max check
  df_il_U <- df_il_U %>%
    filter(D_il != max(df_il$D_il)) %>%
    mutate(D_il = D_il + 1)

  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  # https://fanwangecon.github.io/PrjOptiAlloc/reference/df_opt_caschool_input_il.html
  if (bl_print_verbose){
    # id_i id_il D_max_i  D_il  A_il alpha_il  beta_i
    head(df_il_U, 50)
    tail(df_il_U, 50)
    # Summarize
    REconTools::ff_summ_percentiles(df_il_U)
  }

  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  # Rescale
  df_il_U <- df_il_U %>%
    mutate(v_A_il = v_A_il + 18) %>%
    mutate(v_A_il = case_when(v_A_il >= 1 ~ v_A_il,
                              v_A_il <  1 ~ 1 ))
  # Summarize
  if (bl_print){
    REconTools::ff_summ_percentiles(df_il_U)
  }

  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  ar_rho <- 1 - (10^(c(seq(-2,2, length.out=8))))
  ar_rho <- unique(ar_rho)
  ar_rho <- c(1)

  ## ---- fig.width=5, fig.height=5------------------------------------------------------------------------------------------------------------------------------
  # subset select
  set.seed(123)
  it_draw <- length(ar_group_ids)
  # it_draw <- 30
  ar_group_rand <- ar_group_ids[sample(length(ar_group_ids), it_draw, replace=FALSE)]
  df_input_il <- df_il_U %>%
    filter(id_i %in% ar_group_rand) %>%
    mutate(id_il = row_number())

  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  if (bl_non_inc_adjust){
    # Update c_alpha_il
    df_input_il_noninc <- df_input_il %>%
      group_by(id_i) %>%
      do(c_alpha_il_noninc = ffi_alpha_non_increasing_adj(.$c_alpha_il)) %>%
      unnest(c(c_alpha_il_noninc)) %>%
      group_by(id_i) %>%
      mutate(D_il = row_number()) %>%
      left_join(df_input_il, by=(c('id_i'='id_i', 'D_il'='D_il')))

    # Update v_alpha_il
    df_input_il_noninc <- df_input_il_noninc %>%
      group_by(id_i) %>%
      do(v_alpha_il_noninc = ffi_alpha_non_increasing_adj(.$v_alpha_il)) %>%
      unnest(c(v_alpha_il_noninc)) %>%
      group_by(id_i) %>%
      mutate(D_il = row_number()) %>%
      left_join(df_input_il_noninc, by=(c('id_i'='id_i', 'D_il'='D_il')))

    # replace
    df_input_il_noninc <- df_input_il_noninc %>%
      select(-c_alpha_il, -v_alpha_il) %>%
      rename(c_alpha_il = c_alpha_il_noninc) %>%
      rename(v_alpha_il = v_alpha_il_noninc)
  } else {
    df_input_il_noninc <- df_input_il
  }

  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  df_input_il_covid_actual <- df_input_il_noninc %>%
    left_join(df_id, by = "id_i") %>%
    mutate(ymin_group = as.numeric(ymin_group)) %>%
    ungroup() %>% rowwise() %>%
    mutate(actual_checks =
             ffi_vox_checks_ykm(ymin_group, marital, kids,
                                ar_ycut, fl_multiple, fl_percheck_dollar)) %>%
    ungroup()

  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  # Summarize:
  if (bl_print){
    REconTools::ff_summ_percentiles(df_input_il_covid_actual)
  }
  if (bl_print_verbose){
    # Group Summarize:
    ls_svr_group_vars_cact <- c('ymin_group', 'marital', 'kids')
    REconTools::ff_summ_bygroup(
      df_input_il_covid_actual, ls_svr_group_vars_cact, 'actual_checks')$df_table_grp_stats
    # Group Summarize:
    REconTools::ff_summ_bygroup(df_input_il_covid_actual, c('ymin_group'), 'actual_checks')$df_table_grp_stats
    REconTools::ff_summ_bygroup(df_input_il_covid_actual, c('marital'), 'actual_checks')$df_table_grp_stats
    REconTools::ff_summ_bygroup(df_input_il_covid_actual, c('kids'), 'actual_checks')$df_table_grp_stats
  }

  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  mass_sum <- df_input_il_covid_actual %>% summarize(mass_sum = sum(mass))
  fl_cost_max_checks_all <- mass_sum*fl_percheck_dollar*fl_tax_hh
  st_covid_check_cost <- paste0('Cost All Max Checks=$', round(fl_cost_max_checks_all/1000000000,2),
                                ' bil (assume ',round(fl_tax_hh/1000000,2),' mil tax households)')
  if (bl_print){
    print(st_covid_check_cost)
  }

  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  mass_sum_covid_vox_actual <- df_input_il_covid_actual %>%
    filter(actual_checks >= D_il) %>%
    summarize(mass_cumsum = sum(mass))
  fl_cost_actual <- mass_sum_covid_vox_actual*fl_percheck_dollar*fl_tax_hh
  st_covid_check_cost <- paste0('VOX Policy Cost=$', round(fl_cost_actual/1000000000,2),
                                ' bil (assume ',round(fl_tax_hh/1000000,2),
                                ' mil tax households, use SNW 2020 simulated P(Kids, Marry, INcome))')
  if (bl_print){
    print(st_covid_check_cost)
  }

  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  # 2020 consumption
  df_input_ib_c <- df_input_il_covid_actual %>%
    filter(case_when(actual_checks == 0 ~ D_il == 1, # if check = 0, filter D_il = 1
                     TRUE ~ D_il == actual_checks)) %>%
    rename(A_i_l0 = c_A_il) %>%
    mutate(alpha_o_i = case_when(actual_checks == 0 ~ 0,
                                 TRUE ~ c_alpha_il)) %>%
    select(id_i, A_i_l0, alpha_o_i, beta_i, actual_checks)

  # value
  df_input_ib_v <- df_input_il_covid_actual %>%
    filter(case_when(actual_checks == 0 ~ D_il == 1, # if check = 0, filter D_il = 1
                     TRUE ~ D_il == actual_checks)) %>%
    rename(A_i_l0 = v_A_il) %>%
    mutate(alpha_o_i = case_when(actual_checks == 0 ~ 0,
                                 TRUE ~ v_alpha_il)) %>%
    select(id_i, A_i_l0, alpha_o_i, beta_i, actual_checks)

  # summarize
  if (bl_print){
    REconTools::ff_summ_percentiles(df_input_ib_c)
    REconTools::ff_summ_percentiles(df_input_ib_v)
  }

  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  # Total checks
  it_total_checks <- df_input_il_covid_actual %>%
    filter(D_il == 1) %>%
    summarize(total_checks = sum(actual_checks))
  it_total_checks <- as.numeric(it_total_checks)
  # this is the measure of checks available given VOX allocation and simulated mass
  # summarize
  if (bl_print){
    print(paste0('mass_sum_covid_vox_actual=',mass_sum_covid_vox_actual))
  }
  # And this point, the number is not important
  fl_dis_w <- it_total_checks
  if (bl_print){
    print(paste0('fl_dis_w=',fl_dis_w))
  }

  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  # Inputs
  # 30 individuals 25 checks, half the amount is available
  df_input_il_c <- df_input_il_noninc %>%
    rename(A_il = c_A_il) %>%
    rename(alpha_il = c_alpha_il) %>%
    select(-v_A_il, -v_alpha_il)

  # Solve with Function
  ls_dis_solu_c <- suppressWarnings(suppressMessages(
    ffp_opt_anlyz_rhgin_dis(ar_rho,
                            fl_dis_w,
                            df_input_il_c,
                            bl_df_alloc_il = FALSE,
                            bl_return_V = TRUE,
                            bl_return_allQ_V = FALSE,
                            bl_return_inner_V = FALSE)))
  df_queue_il_long_c <-ls_dis_solu_c$df_queue_il_long
  df_alloc_i_long_c <- ls_dis_solu_c$df_alloc_i_long
  df_rho_gini_c <- ls_dis_solu_c$df_rho_gini


  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  # Inputs
  # 30 individuals 25 checks, half the amount is available
  df_input_il_v <- df_input_il_noninc %>%
    rename(A_il = v_A_il) %>%
    rename(alpha_il = v_alpha_il) %>%
    select(-c_A_il, -c_alpha_il)

  # Solve with Function
  ls_dis_solu_v <- suppressWarnings(suppressMessages(
    ffp_opt_anlyz_rhgin_dis(ar_rho,
                            fl_dis_w,
                            df_input_il_v,
                            bl_df_alloc_il = FALSE,
                            bl_return_V = TRUE,
                            bl_return_allQ_V = FALSE,
                            bl_return_inner_V = FALSE)))
  df_queue_il_long_v <-ls_dis_solu_v$df_queue_il_long
  df_alloc_i_long_v <- ls_dis_solu_v$df_alloc_i_long
  df_rho_gini_v <- ls_dis_solu_v$df_rho_gini


  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  REconTools::ff_summ_percentiles(df_queue_il_long_c, bl_statsasrows = FALSE)
  REconTools::ff_summ_percentiles(df_alloc_i_long_c, bl_statsasrows = FALSE)
  REconTools::ff_summ_percentiles(df_queue_il_long_v, bl_statsasrows = FALSE)
  REconTools::ff_summ_percentiles(df_alloc_i_long_v, bl_statsasrows = FALSE)


  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  tb_rho_rev_c <-
    PrjOptiAlloc::ffp_opt_anlyz_sodis_rev(ar_rho,
                                          fl_dis_w,
                                          df_input_ib = df_input_ib_c,
                                          df_queue_il_long_with_V = df_queue_il_long_c)


  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  # Display Results
  print(tb_rho_rev_c)


  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  tb_rho_rev_v <-
    PrjOptiAlloc::ffp_opt_anlyz_sodis_rev(ar_rho,
                                          fl_dis_w,
                                          df_input_ib = df_input_ib_v,
                                          df_queue_il_long_with_V = df_queue_il_long_v)

  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  # Display Results
  print(tb_rho_rev_v)


  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  # Consider only subsets of ages for graphing
  df_input_il_noninc_covar <- df_input_il_noninc %>%
    left_join(df_id, by = "id_i") %>%
    filter(kids <= 4) %>%
    mutate(ymin_group = as.numeric(ymin_group),
           kids = as.factor(kids),
           marital = as.factor(marital),
           age_group = as.factor(age_group),
           checks = D_il)
  # Summarize
  df_alloc_i_long_covar_v <- df_alloc_i_long_v %>%
    left_join(df_id, by = "id_i") %>%
    filter(kids <= 4)
  df_alloc_i_long_covar_c <- df_alloc_i_long_c %>%
    left_join(df_id, by = "id_i") %>%
    filter(kids <= 4)


  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  stg_source = paste0('SNW 2020 Simulation')
  stg_age = paste0('Age Between ', it_min_age, ' and ', it_max_age, '; Income Groups = ', it_inc_groups)
  stg_caption = paste0(stg_source, '\n', stg_age)


  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  df_alloc_i_long_covar_v <- df_alloc_i_long_covar_v %>%
    left_join(df_input_ib_v %>% select(id_i, actual_checks), by='id_i')
  df_alloc_i_long_covar_c <- df_alloc_i_long_covar_c %>%
    left_join(df_input_ib_c %>% select(id_i, actual_checks), by='id_i')


  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  df_alloc_i_long_covar_v %>%
    select(D_star_i, actual_checks, kids, marital, ymin_group) %>%
    arrange(ymin_group, marital, kids)


  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  # v table
  df_alloc_i_long_covar_v <- df_alloc_i_long_covar_v %>%
    rename(allocate_check_optimal = D_star_i,
           allocate_check_actual = actual_checks) %>%
    pivot_longer(cols = starts_with('allocate_check_'),
                 names_to = c('allocate_type'),
                 names_pattern = paste0("allocate_check_(.*)"),
                 values_to = "checks")
  # c table
  df_alloc_i_long_covar_c <- df_alloc_i_long_covar_c %>%
    rename(allocate_check_optimal = D_star_i,
           allocate_check_actual = actual_checks) %>%
    pivot_longer(cols = starts_with('allocate_check_'),
                 names_to = c('allocate_type'),
                 names_pattern = paste0("allocate_check_(.*)"),
                 values_to = "checks")


  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------

  # x-labels
  x.labels <- c('λ=0.99', 'λ=0.90', 'λ=0', 'λ=-10', 'λ=-100')
  x.breaks <- c(0.01, 0.10, 1, 10, 100)

  # title line 2
  # title_line1 <- sprintf("Percentage of Training Spots Misallocated, NSW Lalonde (AER, 1986)")
  # title_line2 <- sprintf("REV (Resource Equivalent Variation) Along Planner Spectrum")

  st_title <- sprintf(paste0('How much Fewer Resources are Needed (Shares) to Achieve the Same Welfare'))
  title_line1 <- sprintf("Covid Check: Compare Optimal vs Uniform Allocation (Andrew Yang)")

  # Graph Results--Draw
  line.plot <- tb_rho_rev_v %>%
    mutate(REV = 100*REV)  %>%
    mutate(one_minus_rho = 1 - rho_val)  %>%
    ggplot(aes(x=one_minus_rho, y=REV)) +
    geom_line() +
    geom_point() +
    # geom_vline(xintercept=c(1), linetype="dotted") +
    labs(title = st_title,
         subtitle = paste0(title_line1),
         x = 'log10 Rescale of λ, Log10(1-λ)\nλ=1 Utilitarian (Maximize Average), λ=-infty Rawlsian (Maximize Minimum)',
         y = paste0('100 x REV (Resource Equivalent Variations)'),
         caption = 'SNW 2020 Life Cycle Simulations.') +
    scale_x_continuous(trans='log10', labels = x.labels, breaks = x.breaks) +
    theme_bw(base_size=8) +
    ylim(0, 100)
  # +
  #   guides(colour=FALSE)


  # Print
  print(line.plot)

  spt_img_save <- 'img/'
  bl_save_img <- FALSE
  if (bl_save_img) {
    snm_cnts <- 'opti_vs_unif_rev.png'
    png(paste0(spt_img_save, snm_cnts),
        width = 135, height = 96, units='mm', res = 300, pointsize=7)
    print(line.plot)
    dev.off()
  }


  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  ls_svr_groups <- c('ymin_group', 'marital', 'kids')
  for (svr_group in ls_svr_groups) {

    # Group by variable
    print(paste0('current group = ', svr_group))

    # Summarize
    df <- df_alloc_i_long_covar_c
    vars.group <- c('rho_val', svr_group, 'allocate_type')
    var.numeric <- 'checks'
    str.stats.group <- 'allperc'
    ar.perc <- c(0.10, 0.25, 0.50, 0.75, 0.90)
    ls_summ_by_group <- REconTools::ff_summ_bygroup(
      df, vars.group, var.numeric, str.stats.group, ar.perc)
    print(ls_summ_by_group$df_table_grp_stats)

  }


  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  ls_svr_groups <- c('ymin_group', 'marital', 'kids')
  for (svr_group in ls_svr_groups) {

    # Group by variable
    print(paste0('current group = ', svr_group))

    # Summarize
    df <- df_alloc_i_long_covar_v
    vars.group <- c('rho_val', svr_group, 'allocate_type')
    var.numeric <- 'checks'
    str.stats.group <- 'allperc'
    ar.perc <- c(0.10, 0.25, 0.50, 0.75, 0.90)
    ls_summ_by_group <- REconTools::ff_summ_bygroup(
      df, vars.group, var.numeric, str.stats.group, ar.perc)
    print(ls_summ_by_group$df_table_grp_stats)

  }


  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  # graph mean check amount by income, marital status and kids counts
  lineplot_c <- df_alloc_i_long_covar_c %>%
    filter(rho_val == ar_rho[1]) %>% ungroup() %>%
    # filter(marital == 1) %>%
    mutate(ymin_group = as.numeric(ymin_group),
           kids = as.factor(kids),
           marital = as.factor(marital)) %>%
    ggplot(aes(x=ymin_group, y=mass ,
               colour=kids,
               shape=marital)) +
    facet_wrap( ~ marital + kids, ncol=5) +
    geom_point(size=3) +
    labs(title = paste0('Mass At Each Kids/Marry/Income Points, ',
                        st_file_type),
         subtitle = subtitle,
         x = 'Income Group',
         y = 'Mass for this Group',
         caption = stg_caption)
  print(lineplot_c)


  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  # df_input_il_noninc %>%
  #   left_join(df_id, by = "id_i")

  # graph mean check amount by income, marital status and kids counts
  lineplot_c <- df_input_il_noninc_covar %>%
    ggplot(aes(x=ymin_group, y=log(c_alpha_il),
               colour=age_group,
               shape=kids)) +
    facet_wrap( ~ marital + kids, ncol=5) +
    geom_point(size=1) +
    # scale_y_continuous(trans='log') +
    labs(title = paste0('Marginal Consumption Gain Each Check, Conditional on: Age+Marry+Kids+Income, ',
                        st_file_type, '\n', 'Top Row = Unmarried, Bottom Row = Married, Dark Blue = check 1, Light Blue = Final Check'),
         subtitle = subtitle,
         x = 'Income Group',
         y = 'Log(Marginal Consumption Gain Per Check)',
         caption = stg_caption)
  print(lineplot_c)


  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  # graph mean check amount by income, marital status and kids counts
  lineplot_c <- df_alloc_i_long_covar_c %>%
    filter(allocate_type == 'optimal') %>%
    filter(rho_val == ar_rho[1]) %>% ungroup() %>%
    mutate(ymin_group = as.numeric(ymin_group),
           kids = as.factor(kids),
           marital = as.factor(marital)) %>%
    ggplot(aes(x=ymin_group, y=checks,
               colour=age_group,
               shape=kids,
               linetype=marital)) +
    facet_wrap( ~ marital + kids, ncol=5) +
    geom_point(size=3) +
    geom_line() +
    labs(title = paste0(
      'Optimal First Check\n',
      'Optimize Expected 2020 Consumption, Conditional on: Age+Marry+Kids+Income, ',
      st_file_type, '\n', 'Colors Represent Actual Policy versus Optimal Policy'),
      subtitle = subtitle,
      x = 'Income Group',
      y = 'Number of Checks',
      caption = stg_caption)
  print(lineplot_c)


  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  # df_input_il_noninc %>%
  #   left_join(df_id, by = "id_i")

  # graph mean check amount by income, marital status and kids counts
  lineplot_c <- df_input_il_noninc_covar %>%
    ggplot(aes(x=ymin_group, y=log(v_alpha_il),
               colour=age_group,
               shape=kids)) +
    facet_wrap( ~ marital + kids, ncol=5) +
    geom_point(size=1) +
    labs(title = paste0('Marginal Value Gain Each Check, Conditional on: Age+Marry+Kids+Income, ',
                        st_file_type, '\n', 'Top Row = Unmarried, Bottom Row = Married, Dark Blue = check 1, Light Blue = Final Check'),
         subtitle = subtitle,
         x = 'Income Group',
         y = 'Log(Marginal Value Gain Per Check)',
         caption = stg_caption)
  print(lineplot_c)


  ## ------------------------------------------------------------------------------------------------------------------------------------------------------------
  # graph mean check amount by income, marital status and kids counts
  lineplot_v <- df_alloc_i_long_covar_v %>%
    filter(allocate_type == 'optimal') %>%
    filter(rho_val == ar_rho[1]) %>% ungroup() %>%
    mutate(ymin_group = as.numeric(ymin_group),
           kids = as.factor(kids),
           marital = as.factor(marital)) %>%
    ggplot(aes(x=ymin_group, y=checks,
               colour=age_group,
               shape=kids,
               linetype=marital)) +
    facet_wrap( ~ marital + kids, ncol=5) +
    geom_point(size=3) +
    geom_line() +
    labs(title =
           paste0('Optimize Expected 2020 Value, Conditional on: Age+Marry+Kids+Income, ',
                  st_file_type, '\n', 'Colors Represent Actual Policy versus Optimal Policy'),
         subtitle = subtitle,
         x = 'Income Group',
         y = 'Number of Checks',
         caption = stg_caption)
  print(lineplot_v)
}

ffi_alpha_non_increasing_adj <- function(ar_alpha, fl_min_inc_bd = 1e-20) {
  # alpha has tiny upticks sometimes due to approximation errors, need to be adjusted
  # Following theorem 1, alpha must be non-increasing.

  ar_cur <- ar_alpha
  ar_cur_diff <- diff(ar_cur)
  # New Array of NAs
  ar_cur_df_adj <- rep(NA, length(ar_cur_diff))
  # No changes needed if difference is negative, or zero, non-increasing
  ar_cur_df_adj[ar_cur_diff<=0] = 0
  # Record changes if positive
  ar_cur_df_adj[ar_cur_diff>0] = ar_cur_diff[ar_cur_diff>0]
  # Cumulative sum adjustment needed
  ar_cur_adj_cumsum <- cumsum(ar_cur_df_adj)
  ar_cur_adj_cumsum <- c(0, ar_cur_adj_cumsum);
  # Add adjustment to original array
  ar_cur_adj <- ar_cur - ar_cur_adj_cumsum
  # Make sure Adjusted changes are non-negative
  ar_cur_adj[ar_cur_adj < fl_min_inc_bd] = fl_min_inc_bd

  # return
  return(ar_cur_adj)
}


## ------------------------------------------------------------------------------------------------------------------------------------------------------------
ffi_vox_checks_ykm <- function(ymin_group, marital, kids,
                               ar_ycut,
                               fl_multiple = 58056, fl_percheck_dollar = 100,
                               it_inc_subgroups = 10) {
  # # Receive at Most
  # fl_max_checks <- 1200*2+4*500
  # # Max Phase Out Group
  # fl_mas_start_drop <- 150000
  # # Drop rate
  # fl_drop_rate <- 5/100
  # # max Phase
  # fl_phase_out_max <- fl_mas_start_drop + fl_max_checks/fl_drop_rate
  # print(fl_phase_out_max)

  # # Poorest married 4 kids
  # ffi_vox_checks_ykm(1,1,4, ar_ycut, fl_multiple, fl_percheck_dollar)
  # # Test Function
  # ffi_vox_checks_ykm(2,0,0, ar_ycut, fl_multiple, fl_percheck_dollar)
  # ffi_vox_checks_ykm(2,1,1, ar_ycut, fl_multiple, fl_percheck_dollar)
  # ffi_vox_checks_ykm(2,1,3, ar_ycut, fl_multiple, fl_percheck_dollar)
  # ffi_vox_checks_ykm(2,1,3, ar_ycut, fl_multiple, fl_percheck_dollar)
  # # Test Function
  # ffi_vox_checks_ykm(4,0,0, ar_ycut, fl_multiple, fl_percheck_dollar)
  # ffi_vox_checks_ykm(4,1,1, ar_ycut, fl_multiple, fl_percheck_dollar)
  # ffi_vox_checks_ykm(4,1,3, ar_ycut, fl_multiple, fl_percheck_dollar)
  # ffi_vox_checks_ykm(4,1,3, ar_ycut, fl_multiple, fl_percheck_dollar)
  # # Test Function
  # ffi_vox_checks_ykm(11,0,0, ar_ycut, fl_multiple, fl_percheck_dollar)
  # ffi_vox_checks_ykm(11,1,1, ar_ycut, fl_multiple, fl_percheck_dollar)
  # ffi_vox_checks_ykm(11,1,3, ar_ycut, fl_multiple, fl_percheck_dollar)
  # ffi_vox_checks_ykm(11,1,3, ar_ycut, fl_multiple, fl_percheck_dollar)

  # marital 0 or 1
  # kids 0,1,2,3
  # ymin_group index from 1 to N

  # income minimum and maximum,
  # by construction: length(ar_ycut_dollar) = max(ymin_group) + 1
  ar_ycut_dollar = ar_ycut*fl_multiple
  it_ygroups = length(ar_ycut_dollar)

  # Default no checks
  fl_check_dollor = 0

  # Only provide checks if not in the last income group, where all receive none
  if (ymin_group < it_ygroups - 1) {

    # start point household head
    fl_check_dollar = 1200

    # married household gets more, marital = 0, 1
    fl_check_dollar = fl_check_dollar + marital*1200

    # Households with kids: 0, 1,2,3,4
    fl_check_dollar = fl_check_dollar + kids*500


    # lower and upper bounds on income
    fl_inc_group_lower_bound = ar_ycut_dollar[ymin_group]
    fl_inc_group_upper_bound = ar_ycut_dollar[ymin_group+1]

    # A grid of income between these two points: 10 points
    ar_inc_grid = seq(fl_inc_group_lower_bound, fl_inc_group_upper_bound,
                      length.out=it_inc_subgroups)

    # What is the tax rate at each point of these incomes given marry and kids?
    ar_check_reduce_inc_grid = matrix(data=NA, nrow=length(ar_inc_grid), ncol=1)

    # as income increases, fl_check_dollar go down
    it_ctr = 0
    for (fl_inc in ar_inc_grid) {

      it_ctr = it_ctr + 1

      if (marital == 0 && kids == 0) {
        # The benefit would start decreasing at a rate of $5 for every additional $100 in income
        fl_check_reduce = ((max(fl_inc - 75000,0))/100)*5
      }

      # phaseout starts $112,500 for heads of household
      if (marital == 0 && kids != 0) {
        # The benefit would start decreasing at a rate of $5 for every additional $100 in income
        fl_check_reduce = ((max(fl_inc - 112500,0))/100)*5
      }

      # phaseout starts $150,000 for heads of household
      if (marital == 1 ) {
        # The benefit would start decreasing at a rate of $5 for every additional $100 in income
        fl_check_reduce = ((max(fl_inc - 150000,0))/100)*5
      }

      ar_check_reduce_inc_grid[it_ctr] = max(0, fl_check_dollar - fl_check_reduce)

    }

    fl_check_dollor = mean(ar_check_reduce_inc_grid)

  }

  # Check Numbers
  fl_avg_checks = round(fl_check_dollor/fl_percheck_dollar, 0)

  return(fl_avg_checks)
}
