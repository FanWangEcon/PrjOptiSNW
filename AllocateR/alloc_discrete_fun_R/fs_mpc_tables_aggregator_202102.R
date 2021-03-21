# This is an aggregation function that automates prior excel tasks
# in PrjNygaardSorensenWang\Results202101\ folder
# file apc_smooth_bykidsmarital20k_allchecks_18t64_comparison.csv
# Where we aggregate a file like apc_smooth_bykidsmarital20k_allchecks_18t64.csv
# in folder like \Results202101\20210106_snwx_v_planner_docdense_e2hm2_b1_xi0_manna_88_bt97\csv
# to generate aggregate APC or MPC at different levels of checks.
# This was done manually prior. This function is meant to automate the process.
#
# 1. Load a particular APC or MPC csv file
# 2. weighted average for all columns, representing different check amounts. single row
# 3. Gather the single rows together into a common file.
# 4. Reweighted rows for mixture models.
#

ls_st_file_collection <- c('bidenchk_docdense', 'trumpchk_docdense', 'mixture_bidenchk_trumpchk')
ls_st_file_collection <- c('mixture_bidenchk_trumpchk')

ls_snm_csv_raw_solu_mixturemodel <- c('snwx_trumpchk_moredense_a65zh266zs5_b1_xi0_manna_168_rho1',
                                      'snwx_trumpchk_moredense_a65zh266zs5_b1_xi0_manna_168_unmarried_rho1',
                                      'snwx_trumpchk_moredense_a65zh266zs5_b1_xi0_manna_168_married_rho1',
                                      'snwx_trumpchk_moredense_a65zh266zs5_b1_xi0_manna_168_bt60_rho1',
                                      'snwx_trumpchk_moredense_a65zh266zs5_b1_xi0_manna_168_bt95_rho1',
                                      'snwx_bidenchk_moredense_a65zh266zs5_b1_xi0_manna_168_rho1',
                                      'snwx_bidenchk_moredense_a65zh266zs5_b1_xi0_manna_168_unmarried_rho1',
                                      'snwx_bidenchk_moredense_a65zh266zs5_b1_xi0_manna_168_married_rho1',
                                      'snwx_bidenchk_moredense_a65zh266zs5_b1_xi0_manna_168_bt60_rho1',
                                      'snwx_bidenchk_moredense_a65zh266zs5_b1_xi0_manna_168_bt95_rho1')

st_solu_prefix <- 'bidenchk_docdense'
ls_snm_csv_raw_solu_bidenchk <- c(paste0('snwx_', st_solu_prefix, '_e2hm2_b1_xi0_manna_88_bt99'),
                                  paste0('snwx_', st_solu_prefix, '_e2hm2_b1_xi0_manna_88_bt97'),
                                  paste0('snwx_', st_solu_prefix, '_e2hm2_b1_xi0_manna_88_bt95'),
                                  paste0('snwx_', st_solu_prefix, '_e2hm2_b1_xi0_manna_88_bt90'),
                                  paste0('snwx_', st_solu_prefix, '_e2hm2_b1_xi0_manna_88_bt80'),
                                  paste0('snwx_', st_solu_prefix, '_e2hm2_b1_xi0_manna_88_bt70'),
                                  paste0('snwx_', st_solu_prefix, '_e2hm2_b1_xi0_manna_88_bt60'),
                                  paste0('snwx_', st_solu_prefix, '_e2hm2_b1_xi0_manna_88_bt50'),
                                  paste0('snwx_', st_solu_prefix, '_e1lm2_b1_xi0_manna_88_bt99'),
                                  paste0('snwx_', st_solu_prefix, '_e1lm2_b1_xi0_manna_88_bt97'),
                                  paste0('snwx_', st_solu_prefix, '_e1lm2_b1_xi0_manna_88_bt95'),
                                  paste0('snwx_', st_solu_prefix, '_e1lm2_b1_xi0_manna_88_bt90'),
                                  paste0('snwx_', st_solu_prefix, '_e1lm2_b1_xi0_manna_88_bt80'),
                                  paste0('snwx_', st_solu_prefix, '_e1lm2_b1_xi0_manna_88_bt70'),
                                  paste0('snwx_', st_solu_prefix, '_e1lm2_b1_xi0_manna_88_bt60'),
                                  paste0('snwx_', st_solu_prefix, '_e1lm2_b1_xi0_manna_88_bt50'))
st_solu_prefix <- 'trumpchk_docdense'
ls_snm_csv_raw_solu_trumpchk <- c(paste0('snwx_', st_solu_prefix, '_e2hm2_b1_xi0_manna_88_bt99'),
                                  paste0('snwx_', st_solu_prefix, '_e2hm2_b1_xi0_manna_88_bt97'),
                                  paste0('snwx_', st_solu_prefix, '_e2hm2_b1_xi0_manna_88_bt95'),
                                  paste0('snwx_', st_solu_prefix, '_e2hm2_b1_xi0_manna_88_bt90'),
                                  paste0('snwx_', st_solu_prefix, '_e2hm2_b1_xi0_manna_88_bt80'),
                                  paste0('snwx_', st_solu_prefix, '_e2hm2_b1_xi0_manna_88_bt70'),
                                  paste0('snwx_', st_solu_prefix, '_e2hm2_b1_xi0_manna_88_bt60'),
                                  paste0('snwx_', st_solu_prefix, '_e1lm2_b1_xi0_manna_88_bt99'),
                                  paste0('snwx_', st_solu_prefix, '_e1lm2_b1_xi0_manna_88_bt97'),
                                  paste0('snwx_', st_solu_prefix, '_e1lm2_b1_xi0_manna_88_bt95'),
                                  paste0('snwx_', st_solu_prefix, '_e1lm2_b1_xi0_manna_88_bt90'),
                                  paste0('snwx_', st_solu_prefix, '_e1lm2_b1_xi0_manna_88_bt80'),
                                  paste0('snwx_', st_solu_prefix, '_e1lm2_b1_xi0_manna_88_bt70'),
                                  paste0('snwx_', st_solu_prefix, '_e1lm2_b1_xi0_manna_88_bt60'),
                                  paste0('snwx_', st_solu_prefix, '_e1lm2_b1_xi0_manna_88_bt50'))

ls_snm_csv_mpcapc <- c('mpc_smooth_bykidsmarital20k_allchecks_18t64',
                       'mpc_smooth_bykidsmarital20k_allchecks_18t99',
                       'apc_smooth_bykidsmarital20k_allchecks_18t64',
                       'apc_smooth_bykidsmarital20k_allchecks_18t99')


bl_verbose <- TRUE
ar_checks_select <- c(0, 14, 28, 42, 56, 70, 84)

for (st_file_collection in ls_st_file_collection) {

  # Loop over MPC/APC file type
  for (snm_csv_mpcapc in ls_snm_csv_mpcapc) {

    it_df_ctr <- 0
    # Loop over beta/edu folders
    if (st_file_collection == 'bidenchk_docdense') {
      ls_snm_csv_raw_solu <- ls_snm_csv_raw_solu_bidenchk
      st_file_collection_prefix <- st_file_collection
      ar_checks_select <- seq(from=0, to=84, by=1)
    } else if (st_file_collection == 'trumpchk_docdense') {
      ls_snm_csv_raw_solu <- ls_snm_csv_raw_solu_trumpchk
      st_file_collection_prefix <- st_file_collection
      ar_checks_select <- seq(from=0, to=84, by=1)
    } else if (st_file_collection == 'mixture_bidenchk_trumpchk') {
      ls_snm_csv_raw_solu <- ls_snm_csv_raw_solu_mixturemodel
      st_file_collection_prefix <- st_file_collection
      ar_checks_select <- seq(from=0, to=168, by=1)
    } else {
      error('stop no st_file_collection specified')
    }

    # list to collect singles row df containing aggregate mpc/apc for beta/edu groups
    ls_df_aggregate_mpcapc <- vector(mode = "list", length = length(ls_snm_csv_mpcapc)*length(ls_snm_csv_raw_solu))

    for (snm_csv_raw_solu in ls_snm_csv_raw_solu) {

      # Get key strings
      st_which_solu <- snm_csv_raw_solu
      ls_output <- fs_opti_support_202102(st_which_solu)
      srt_results_root <- ls_output$srt_results_root
      srt_csv_path_root <- ls_output$srt_csv_path_root

      # Load file
      print(paste0('Starting to load file ', snm_csv_mpcapc, ' from folder ', srt_csv_path_root))
      df_mpc_apc_by_check <- as_tibble(
        read.csv(paste0(srt_csv_path_root, snm_csv_mpcapc, '.csv'), header=TRUE))

      # summary and str
      if (bl_verbose) {
        print(paste0('Starting to summarize ', snm_csv_mpcapc))
        REconTools::ff_summ_percentiles(df_mpc_apc_by_check)
      }

      # 1. summ over weight column, total weight from the mass column
      fl_total_mass <- sum(df_mpc_apc_by_check$mass)

      # 2. drop marital, kids, income group columns, which are not relevant
      df_mpc_apc_by_check <- df_mpc_apc_by_check %>% select(-marital, -kids, -ymin_group, -X)

      # 3. wide to long
      df_mpc_apc_by_check_long <- df_mpc_apc_by_check %>%
        pivot_longer(cols = starts_with('X'),
                     names_to = c('checks'),
                     names_pattern = paste0("X(.*)"),
                     values_to = "GroupMPCorAPC")

      # 4. group by checks, sort by id_i, weighted average betwen mass and groupMPCorAPC
      df_mpc_apc_by_check_long_summ <- df_mpc_apc_by_check_long %>%
        mutate(checks = as.numeric(checks)) %>%
        group_by(checks) %>%
        summarise(GroupMPCorAPCAggregate = sum(mass*GroupMPCorAPC)/sum(mass)) %>%
        filter(checks %in% ar_checks_select)

      # 5. reshape to wide
      df_mpc_apc_by_check_wide_agg <- df_mpc_apc_by_check_long_summ %>%
        pivot_wider(names_from = "checks",
                    values_from = "GroupMPCorAPCAggregate") %>%
        rename_at(vars(num_range('', ar_checks_select))
                  , list(~paste0('check_', . , ''))
        ) %>%
        mutate(snm_csv_mpcapc = snm_csv_mpcapc) %>%
        mutate(snm_csv_raw_solu = snm_csv_raw_solu)

      # 6. parse file name strings for more information and store
      if (grepl('trumpchk', snm_csv_raw_solu)) {
        st_trumpchkstatus <- 'WITHOUT trump chk (trumpchk file prefix)'
      } else if (grepl('bidenchk', snm_csv_raw_solu)) {
        st_trumpchkstatus <- 'GIVEN trump chk (bidenchk file prefix)'
      } else {
        st_trumpchkstatus <- NA
      }

      if (grepl('_manna', snm_csv_raw_solu)) {
        st_tax_type <- 'manna'
      } else {
        st_tax_type <- NA
      }
      if (grepl('_b1_', snm_csv_raw_solu)) {
        st_bbenefit_type <- 'b=1'
      } else {
        st_bbenefit_type <- NA
      }
      if (grepl('_xi0_', snm_csv_raw_solu)) {
        st_xishk_type <- 'xi=0'
      } else {
        st_xishk_type <- NA
      }

      if (grepl('_docdense', snm_csv_raw_solu)) {
        st_solu_type <- '81x5 shks'
      } else if (grepl('moredense_a65zh266zs5', snm_csv_raw_solu)) {
        st_solu_type <- '266x5 shks'
      } else {
        st_solu_type <- NA
      }

      if (grepl('_e1lm2', snm_csv_raw_solu)) {
        st_edu_type <- 'noncollege'
      } else if (grepl('_e2hm2', snm_csv_raw_solu)) {
        st_edu_type <- 'college'
      } else {
        st_edu_type <- 'alledu'
      }
      if (grepl('_bt99', snm_csv_raw_solu)) {
        st_beta_type <- '0.99'
      } else if (grepl('_bt97', snm_csv_raw_solu)) {
        st_beta_type <- '0.97'
      } else if (grepl('_bt95', snm_csv_raw_solu)) {
        st_beta_type <- '0.95'
      } else if (grepl('_bt90', snm_csv_raw_solu)) {
        st_beta_type <- '0.90'
      } else if (grepl('_bt80', snm_csv_raw_solu)) {
        st_beta_type <- '0.80'
      } else if (grepl('_bt70', snm_csv_raw_solu)) {
        st_beta_type <- '0.70'
      } else if (grepl('_bt60', snm_csv_raw_solu)) {
        st_beta_type <- '0.60'
      } else if (grepl('_bt50', snm_csv_raw_solu)) {
        st_beta_type <- '0.50'
      } else {
        st_beta_type <- NA
      }
      if (grepl('_married_', snm_csv_raw_solu)) {
        st_marriage_type <- 'married'
      } else if (grepl('_unmarried_', snm_csv_raw_solu)) {
        st_marriage_type <- 'unmarried'
      } else {
        st_marriage_type <- 'both'
      }

      df_mpc_apc_by_check_wide_agg <- df_mpc_apc_by_check_wide_agg %>%
        mutate(edu = st_edu_type, beta=st_beta_type, marital=st_marriage_type,
               tax = st_tax_type, unempbenefit=st_bbenefit_type, xishk=st_xishk_type,
               solution = st_solu_type,
               trumpbidenchk = st_trumpchkstatus) %>%
        select(snm_csv_mpcapc, snm_csv_raw_solu,
               solution,
               trumpbidenchk,
               tax, unempbenefit, xishk,
               edu, beta, marital,
               matches('check_'))

      # 7. Save to list
      if (it_df_ctr == 0) {
        df_aggregate_mpcapc <- df_mpc_apc_by_check_wide_agg
      } else {
        df_aggregate_mpcapc <- bind_rows(df_aggregate_mpcapc, df_mpc_apc_by_check_wide_agg)
      }
      it_df_ctr <- it_df_ctr + 1

    }

    # ls of dataframe to a ingle dataframe
    if (bl_verbose) {
      print(paste0('Starting to summarize df_aggregate_mpcapc'))
      REconTools::ff_summ_percentiles(df_aggregate_mpcapc)
    }


    # store aggregate MPC/APC file results
    snm_save_aggregate_mpcapc_csv <- paste0(st_file_collection_prefix, '_', snm_csv_mpcapc, '_agg.csv')
    write.csv(df_aggregate_mpcapc,
              file.path(srt_results_root, snm_save_aggregate_mpcapc_csv),
              row.names = TRUE)


  }
}
