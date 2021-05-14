# This file merges several MPC csv files together, that were solved for different beta
# and EDU combinations. Each file contains P(j,k,m,y|edu,beta) probability mass, and V and
# C weighted averages at each (j,k,m,y,edu,beta) point. Construct P(j,k,m,y) aggregate
# mass across edu and beta types. And generated weighted average C and V at each
# (j, k, m, y) point.
# - model solves conditional on beta and edu
# - merge to mixture with appropriate weights
# - generate marital type specific files
# - merge over edu by beta

# rm(list = ls())

library(tidyverse)
library(readr)

# library(foreach)
# library(doParallel)
# it_no_cores <- detectCores(logical = TRUE)
# cl <- makeCluster(3)
# registerDoParallel(cl)

# Parameters
# st_biden_or_trump <- 'bidenchk'
# st_biden_or_trump <- 'trumpchk'
# st_biden_or_trump <- 'bchklock'
st_biden_or_trump <- 'bchknoui'
n_welfchecksgrid = 169
# ls_it_betaetc_grp <- c(1,2,3,4,5)
ls_it_betaetc_grp <- c(1)
# ls_it_mixture_grp <- c(1,2,3,4,5)
ls_it_mixture_grp <- c(1)
fl_rho <- 1

st_rho <- ''
# if (fl_rho == 1) {
#   st_rho <- '_rho1'
# } else if (fl_rho == -1) {
#   st_rho <- '_rhoneg1'
# }

# Loop over, consider different beta groups
for (it_mixture_grp in ls_it_mixture_grp) {
  for (it_betaetc_grp in ls_it_betaetc_grp) {
  # foreach (it_betaetc_grp=ls_it_betaetc_grp) %dopar% {
    # Conditionally change beta lists to consider, with different names
    st_beta <- ''
    st_marital <- ''
    if (it_betaetc_grp <= 3) {
      ls_fl_beta_val=c(0.60, 0.95)
      if (it_betaetc_grp == 2) {
        st_marital <- '_unmarried'
      }
      if (it_betaetc_grp == 3) {
        st_marital <- '_married'
      }
    } else if (it_betaetc_grp == 4) {
      ls_fl_beta_val=c(0.60)
      st_beta <- '_bt60'
    } else if (it_betaetc_grp == 5) {
      ls_fl_beta_val=c(0.95)
      st_beta <- '_bt95'
    }

    # ls_fl_beta_val=c(0.95)
    ls_it_edu_simu_type=c(1, 2)

    # conditional probabilities
    if (it_mixture_grp == 1) {
      st_mixturegrp <- ''
      fl_wgt_college_frac <- 0.303
      fl_wgt_low_beta_college <- 0.10
      fl_wgt_low_beta_non_college <- 0.40
    }
    if (it_mixture_grp == 2) {
      st_mixturegrp <- '_3o6'
      fl_wgt_college_frac <- 0.303
      fl_wgt_low_beta_college <- 0.10*(3/6)
      fl_wgt_low_beta_non_college <- 0.40*(3/6)
    }
    if (it_mixture_grp == 3) {
      st_mixturegrp <- '_4o6'
      fl_wgt_college_frac <- 0.303
      fl_wgt_low_beta_college <- 0.10*(4/6)
      fl_wgt_low_beta_non_college <- 0.40*(4/6)
    }
    if (it_mixture_grp == 4) {
      st_mixturegrp <- '_8o6'
      fl_wgt_college_frac <- 0.303
      fl_wgt_low_beta_college <- 0.10*(8/6)
      fl_wgt_low_beta_non_college <- 0.40*(8/6)
    }
    if (it_mixture_grp == 5) {
      st_mixturegrp <- '_9o6'
      fl_wgt_college_frac <- 0.303
      fl_wgt_low_beta_college <- 0.10*(9/6)
      fl_wgt_low_beta_non_college <- 0.40*(9/6)
    }

    # File names
    ar_svr_grps <- c('age', 'marital', 'kids', 'ymin','checks')
    ar_svr_csv = c('age', 'marital', 'kids', 'checks',	'ymin',
                   'mass', 'survive', 'vtilde', 'ctilde')
    # File path and computesize
    srt_root <- 'D:/Dropbox (UH-ECON)/PrjNygaardSorensenWang/'
    # st_param_group_base <- 'default_tiny'
    # st_param_group_base <- 'default_small53'
    # st_param_group_base <- 'default_dense'
    # st_param_group_base <- 'default_moredense_a65zh133zs5'
    st_param_group_base <- 'default_moredense_a65zh266zs5'
    # save name

    snm_invoke_suffix = gsub(x = st_param_group_base,  pattern = "default_", replacement = "")
    snm_file_csvwgt_save = paste0('snwx_', st_biden_or_trump, '_',
                                  snm_invoke_suffix, '_b1_xi0_manna_', (n_welfchecksgrid-1),
                                  st_beta, st_marital, st_rho, st_mixturegrp);
    # Folders
    srt_results_root <- paste0(srt_root, 'Results202102')
    srt_simu_path <- paste0(srt_root, 'Output202102/')

    # A. 2d list of file names for beta (row) edu (column) specifi  c MPC results. -----
    # from PrjOptiSNW\invoke\202102\snw_res_b1_manna_mixture.m
    it_beta_cnt <- length(ls_fl_beta_val)
    it_edu_cnt <- length(ls_it_edu_simu_type)
    ls_betaedu_filenames <- vector(mode = "list", length = it_beta_cnt*it_edu_cnt)
    # consistent with m files, beta as rows, edu as columns
    dim(ls_betaedu_filenames) <- c(it_beta_cnt, it_edu_cnt)

    for (it_beta_ctr in seq(1,it_beta_cnt)) {
      fl_beta_val = ls_fl_beta_val[it_beta_ctr]
      for (it_edu_ctr in seq(1,it_edu_cnt)) {
        it_edu_simu_type = ls_it_edu_simu_type[it_edu_ctr]

        st_file_name_suffix = paste0('_bt', round(fl_beta_val*100, 0))
        if (it_edu_simu_type == 1) {
          st_param_group_suffix = '_e1lm2';
        }
        if (it_edu_simu_type == 2) {
          st_param_group_suffix = '_e2hm2';
        }

        # Construct file name
        st_param_group = paste0(st_param_group_base, st_param_group_suffix)
        snm_suffix = paste0('_b1_xi0_manna_', (n_welfchecksgrid-1), st_file_name_suffix);
        snm_invoke_suffix = gsub(x = st_param_group,  pattern = "default_", replacement = "")
        snm_file = paste0('snwx_', st_biden_or_trump, '_', snm_invoke_suffix, snm_suffix);

        # Save file name
        ls_betaedu_filenames[[it_beta_ctr, it_edu_ctr]] <- snm_file

      }
    }
    # name rows and columns
    dimnames(ls_betaedu_filenames)[[1]] <- paste0('beta', ls_fl_beta_val)
    dimnames(ls_betaedu_filenames)[[2]] <- paste0('edu', ls_it_edu_simu_type)

    # B. Probability weight beta and edu group -----
    # from PrjOptiSNW\fgov\snw_tax_steady_multitypes.m``
    fl_wgt_sum <- 0
    ls_betaedu_wgt <- vector(mode = "list", length = it_beta_cnt*it_edu_cnt)
    dim(ls_betaedu_wgt) <- c(it_beta_cnt, it_edu_cnt)
    for (it_beta_ctr in seq(1,it_beta_cnt)) {
      fl_beta_val = ls_fl_beta_val[it_beta_ctr]
      for (it_edu_ctr in seq(1,it_edu_cnt)) {
        it_edu_simu_type = ls_it_edu_simu_type[it_edu_ctr]

        if (it_edu_ctr == 1 && fl_beta_val == 0.6) {
          fl_wgt = (1-fl_wgt_college_frac)*fl_wgt_low_beta_non_college;
        } else if (it_edu_ctr == 1 && fl_beta_val == 0.95) {
          fl_wgt = (1-fl_wgt_college_frac)*(1-fl_wgt_low_beta_non_college);
        } else if (it_edu_ctr == 2 && fl_beta_val == 0.6) {
          fl_wgt = fl_wgt_college_frac*fl_wgt_low_beta_college;
        } else if (it_edu_ctr == 2 && fl_beta_val == 0.95) {
          fl_wgt = fl_wgt_college_frac*(1-fl_wgt_low_beta_college);
        }

        # Save weights
        ls_betaedu_wgt[[it_beta_ctr, it_edu_ctr]] <- fl_wgt
        fl_wgt_sum <- fl_wgt_sum + fl_wgt
      }
    }
    # re-weight
    for (it_beta_ctr in seq(1,it_beta_cnt)) {
      fl_beta_val = ls_fl_beta_val[it_beta_ctr]
      for (it_edu_ctr in seq(1,it_edu_cnt)) {
        fl_wgt <- ls_betaedu_wgt[[it_beta_ctr, it_edu_ctr]]
        ls_betaedu_wgt[[it_beta_ctr, it_edu_ctr]] <- fl_wgt/fl_wgt_sum
      }
    }

    # C. Loop through merge files together, each edu/beta MPC as own column ----
    # update v and c columns with mass*fl_wgt_betaedu
    for (it_betaedugrp in seq(1, length(ls_betaedu_filenames))) {

      # file name
      snm_simu_csv <- ls_betaedu_filenames[[it_betaedugrp]]
      fl_wgt_betedu <- ls_betaedu_wgt[[it_betaedugrp]]

      # load file as csv
      df_plan_v_tilde_full <- as_tibble(
        read.csv(file.path(srt_simu_path, paste0(snm_simu_csv, ".csv")), header=FALSE)) %>%
        rename_all(~c(ar_svr_csv))

      # sum of the Conditional Mass by Type P(kids,marry,inc|edu, beta)
      # value belwow is (Sum_{k,m,y} P(kids,marry,inc|edu, beta))*n_welfchecksgrid
      print(sum(df_plan_v_tilde_full$mass))
      # df_plan_v_tilde_full <- df_plan_v_tilde_full %>% mutate(mass = as.numeric(mass))

      # P(kids,marry,inc, edu, beta)) = P(kids,marry,inc|edu, beta))*P(edu, beta)
      df_plan_v_tilde_full <- df_plan_v_tilde_full %>%
        mutate(mass_joint = mass*fl_wgt_betedu) %>%
        rename(!!(paste0('mass_joint_g', it_betaedugrp)) := mass_joint,
               !!(paste0('survive_', it_betaedugrp)) := survive,
               !!(paste0('vtilde_', it_betaedugrp)) := vtilde,
               !!(paste0('ctilde_', it_betaedugrp)) := ctilde) %>%
        select(-mass)

      # File conditional processing
      if (it_betaedugrp == 1) {
        # This is the first file
        df_plan_v_tilde_full_betaedu <- df_plan_v_tilde_full
      } else {
        # This is other files, merge in
        df_plan_v_tilde_full_betaedu <- df_plan_v_tilde_full_betaedu %>%
          full_join(df_plan_v_tilde_full,
                    by=setNames(c('age', 'marital', 'kids', 'ymin', 'checks'),
                                c('age', 'marital', 'kids', 'ymin', 'checks')))
      }
    }

    # check and review individual mass
    if (it_betaetc_grp <= 3) {
      fl_total_mass <- sum(df_plan_v_tilde_full_betaedu$mass_joint_g1, na.rm=TRUE)+
        sum(df_plan_v_tilde_full_betaedu$mass_joint_g2, na.rm=TRUE)+
        sum(df_plan_v_tilde_full_betaedu$mass_joint_g3, na.rm=TRUE)+
        sum(df_plan_v_tilde_full_betaedu$mass_joint_g4, na.rm=TRUE)
    } else {
      fl_total_mass <- sum(df_plan_v_tilde_full_betaedu$mass_joint_g1, na.rm=TRUE)+
        sum(df_plan_v_tilde_full_betaedu$mass_joint_g2, na.rm=TRUE)
    }
    if (fl_total_mass != 1*n_welfchecksgrid) {
      warning(paste0('fl_total_mass=', fl_total_mass, ', not equal to 1.'))
    }

    # D. Generated conditional weights ----
    # Generate P(edu, beta|age, kids,marry,inc)) at each jkmy, weighted average over C and V
    # select data subset containing only keys and the joint mass
    # 1. Joint mass to long, summ within jkmy group total mass
    # 2. P(edu, beta|jkmy) through division
    # 3. long to wide, keep new conditional column

    # D1. file to long
    df_joint_mass_long <- df_plan_v_tilde_full_betaedu %>%
      select(!!!syms(ar_svr_grps), contains('mass_joint_g')) %>%
      pivot_longer(cols = starts_with('mass_joint_g'),
                   names_to = c('edubetagroup'),
                   names_pattern = paste0("mass_joint_g(.*)"),
                   values_to = "p_joint_jkmyeb") %>%
      group_by(!!!syms(ar_svr_grps)) %>%
      mutate(p_joint_jkmy = sum(p_joint_jkmyeb, na.rm = TRUE)) %>%
      mutate(p_betaedu_condi_jkmy = p_joint_jkmyeb/p_joint_jkmy)
    # p.joint.jkmy should not have nay NA values
    # REconTools::ff_summ_percentiles(df_joint_mass_long)

    # D2. jkmy joint mass sum: P(age, kids,marry,inc))
    # all long columns are the same, so keep 1
    ar_unique_ids <- sort(unique(df_joint_mass_long %>% pull(edubetagroup)))
    df_joint_mass_p_joint_jkmy <- df_joint_mass_long %>%
      select(!!!syms(ar_svr_grps), edubetagroup, p_joint_jkmy) %>%
      pivot_wider(names_from = edubetagroup,
                  values_from = p_joint_jkmy) %>%
      rename_at(vars(num_range('', ar_unique_ids))
                , list(~paste0('p_jkmy_', . , ''))) %>%
      select(!!!syms(ar_svr_grps), p_jkmy_1) %>%
      rename(p_jkmy = p_jkmy_1)
    # debug
    fl_summ_jkmy_mass <- sum(df_joint_mass_p_joint_jkmy$p_jkmy)
    if (fl_summ_jkmy_mass != 1*n_welfchecksgrid) {
      warning(paste0('fl_summ_jkmy_mass=', fl_summ_jkmy_mass, ', not equal to 1.'))
    }
    # REconTools::ff_summ_percentiles(df_joint_mass_p_joint_jkmy)

    # D3, conditional mass: P(edu, beta|age, kids,marry,inc))
    ar_unique_ids <- sort(unique(df_joint_mass_long %>% pull(edubetagroup)))
    df_joint_mass_p_joint_betaedu_condi_jkmy <- df_joint_mass_long %>%
      select(!!!syms(ar_svr_grps), edubetagroup, p_betaedu_condi_jkmy) %>%
      pivot_wider(names_from = edubetagroup,
                  values_from = p_betaedu_condi_jkmy) %>%
      rename_at(vars(num_range('', ar_unique_ids))
                , list(~paste0('p_betaedu_condi_jkmy_', . , ''))
      )

    # E. merge together, and weighted average ----
    # \sum_{edu,beta} \left( V_{edu, beta}*P(edu,beta|jkmy) \right)

    # E1. merge files together
    df_plan_v_tilde_full_betaedu_condiprob <- df_plan_v_tilde_full_betaedu %>%
      left_join(df_joint_mass_p_joint_jkmy,
                by=setNames(c('age', 'marital', 'kids', 'ymin', 'checks'),
                            c('age', 'marital', 'kids', 'ymin', 'checks'))) %>%
      left_join(df_joint_mass_p_joint_betaedu_condi_jkmy,
                by=setNames(c('age', 'marital', 'kids', 'ymin', 'checks'),
                            c('age', 'marital', 'kids', 'ymin', 'checks')))

    # df_plan_v_tilde_full_betaedu_states <-
    #   REconTools::ff_summ_percentiles(df_plan_v_tilde_full_betaedu_condiprob, FALSE)

    # E2. weighted average of V and U vectors
    # first convert all NA to zero
    if (it_betaetc_grp <= 3) {
      df_plan_v_tilde_full_wgted <- df_plan_v_tilde_full_betaedu_condiprob %>%
        mutate_at(vars(matches('p_betaedu_condi_jkmy_|vtilde_|ctilde_')),
                  list(~replace_na(., 0))) %>%
        mutate(vtilde_wgt =
                 (vtilde_1^fl_rho*p_betaedu_condi_jkmy_1 +
                    vtilde_2^fl_rho*p_betaedu_condi_jkmy_2 +
                    vtilde_3^fl_rho*p_betaedu_condi_jkmy_3 +
                    vtilde_4^fl_rho*p_betaedu_condi_jkmy_4)^(1/fl_rho)) %>%
        mutate(ctilde_wgt =
                 (ctilde_1^fl_rho*p_betaedu_condi_jkmy_1 +
                    ctilde_2^fl_rho*p_betaedu_condi_jkmy_2 +
                    ctilde_3^fl_rho*p_betaedu_condi_jkmy_3 +
                    ctilde_4^fl_rho*p_betaedu_condi_jkmy_4)^(1/fl_rho))
      # E3, Survive variable processing
      df_plan_v_tilde_full_wgted <- df_plan_v_tilde_full_wgted %>%
        mutate(survive_joint = case_when(survive_1 > 0 ~ survive_1,
                                         survive_1 <= 0 ~ survive_1,
                                         TRUE ~ survive_3))
    } else {
      df_plan_v_tilde_full_wgted <- df_plan_v_tilde_full_betaedu_condiprob %>%
        mutate_at(vars(matches('p_betaedu_condi_jkmy_|vtilde_|ctilde_')),
                  list(~replace_na(., 0))) %>%
        mutate(vtilde_wgt =
                 case_when(vtilde_1 !=0 && vtilde_2 !=0 ~
                             ((vtilde_1^fl_rho)*p_betaedu_condi_jkmy_1 + (vtilde_2^fl_rho)*p_betaedu_condi_jkmy_2
                             )^(1/fl_rho),
                           vtilde_1 !=0 && vtilde_2 ==0 ~ vtilde_1,
                           vtilde_1 ==0 && vtilde_2 !=0 ~ vtilde_2,
                           TRUE ~ 1)) %>%
        mutate(ctilde_wgt =
                 case_when(ctilde_1 !=0 && ctilde_2 !=0 ~
                             (ctilde_1^fl_rho)*p_betaedu_condi_jkmy_1,
                             # ((ctilde_1^fl_rho)*p_betaedu_condi_jkmy_1 + (ctilde_2^fl_rho)*p_betaedu_condi_jkmy_2
                             # )^(1/fl_rho),
                           ctilde_1 !=0 && ctilde_2 ==0 ~ ctilde_1,
                           ctilde_1 ==0 && ctilde_2 !=0 ~ ctilde_2,
                           TRUE ~ 1)) %>%
        mutate(survive_joint = survive_1)
    }
    # df_plan_v_tilde_full_betaedu_states <-
    #   REconTools::ff_summ_percentiles(df_plan_v_tilde_full_wgted, FALSE)
    # view(df_plan_v_tilde_full_betaedu_states)

    # E4. Column selection etc
    # c('age', 'marital', 'kids', 'checks',	'ymin', 'mass', 'survive', 'vtilde', 'ctilde')
    # column sequence must conform to the ar_svr_csv parameter in ffp_snw_process_inputs
    df_plan_v_tilde_full_wgted_final <- df_plan_v_tilde_full_wgted %>%
      select(age, marital, kids, checks, ymin, p_jkmy, survive_joint, vtilde_wgt, ctilde_wgt)

    # df_plan_v_tilde_full_wgted_final_stats <-
    #   REconTools::ff_summ_percentiles(df_plan_v_tilde_full_wgted_final, FALSE)
    #
    # df_plan_v_tilde_full_wgted_final_stats <-
    #   REconTools::ff_summ_percentiles(df_plan_v_tilde_full_wgted_final, FALSE)
    # view(df_plan_v_tilde_full_wgted_final_stats)

    # view(df_plan_v_tilde_full_wgted_final_stats)
    # ls_summ_by_group_ctilde_wgt <-
    #   REconTools::ff_summ_bygroup(df_plan_v_tilde_full_wgted_final,
    #                               c('ymin'), 'ctilde_wgt')
    # view(ls_summ_by_group_ctilde_wgt$df_table_grp_stats)
    # ls_summ_by_group_vtilde_wgt <-
    #   REconTools::ff_summ_bygroup(df_plan_v_tilde_full_wgted_final,
    #                               c('ymin'), 'vtilde_wgt')
    # view(ls_summ_by_group_vtilde_wgt$df_table_grp_stats)

    # select marital subset
    if (it_betaetc_grp == 2) {
      # unmarried
      print('select unmarried')
      df_plan_v_tilde_full_wgted_final <- df_plan_v_tilde_full_wgted_final %>% filter(marital == 0)
    } else if (it_betaetc_grp == 3) {
      # married
      print('select married')
      df_plan_v_tilde_full_wgted_final <- df_plan_v_tilde_full_wgted_final %>% filter(marital == 1)
    } else {
      # do nothing
      print('do not select marital subset')
    }

    df_plan_v_tilde_full <-
      write_csv(df_plan_v_tilde_full_wgted_final,
                file.path(srt_simu_path,
                          paste0(snm_file_csvwgt_save, ".csv")), col_names = FALSE)

  }
}
# stopCluster(cl)
