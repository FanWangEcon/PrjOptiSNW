# Feasible Allocation Graph, Compare C vs V Allocaitons.
# this is the 2021 04 version of the C vs V results

# Libraries
library(tidyverse)
library(REconTools)

ar_rho <- c(0.05, -1)
ar_rho_nine <- c(1, 0.9, 0.8,
                 0.05,
                 -1, -9, -19, -99)

for (fl_rho in ar_rho_nine) {
  # Parameter setting and controls -------
  st_which_solu <- c('snwx_bidenchk_moredense_a65zh266zs5_b1_xi0_manna_168')
  bl_per_capita <- TRUE
  fl_rho <- c(fl_rho)
  snm_main_noperturb <- 'Results202103'
  st_bound_files <- '20ca20ck'
  st_rho <- paste0(fl_rho)
  st_rho <- gsub(x = st_rho,  pattern = "-", replacement = "n")
  st_rho <- gsub(x = st_rho,  pattern = "\\.", replacement = "p")
  st_rho <- paste0(st_rho)
  st_rho <- paste0('_rh', st_rho)
  st_rho_nounderscore <- paste0('rh', st_rho)
  st_img_name_full <- paste0('Robustness_c_vs_v', st_rho, '_', st_bound_files , '_20210421.png')

  # Load file -----------------
  ls_output_noperturb <- fs_opti_support_202103(st_which_solu,
                                                bl_per_capita=bl_per_capita,
                                                fl_rho=fl_rho,
                                                snm_main=snm_main_noperturb,
                                                it_rand_draw_seed=NULL,
                                                ls_st_gen_agg_folders=c(''))
  srt_paper_appendix_textgraph <- ls_output$srt_paper_appendix_textgraph


  # Get folder name path
  spn_csv_full_path_noperturb <- file.path(ls_output_noperturb$srt_results_root,
                                           ls_output_noperturb$srt_folder, 'csv',
                                           paste0('b1_a64_', st_bound_files, '_18t64'),
                                           paste0('df_alloc_all_feasible_', st_bound_files, '.csv'))
  # rng = 99999 is the not perturbed result.
  df_alloc_cur <- as_tibble(read.csv(spn_csv_full_path_noperturb))

  # Modifying dataframe ------

  # Marital and Kids Level Labeling
  marry_levels <- c(Single = "0", Married = "1")
  kids_levels <- c("no children" = "0", "one child" = "1",
                   "two children" = "2", "three children" = "3",
                   "four children" = "4")
  kids_levels <- c("no children" = "0", "one child" = "1",
                   "two children" = "2")

  ## Select, as factor, and recode
  df_alloc_cur <- df_alloc_cur %>%
    mutate(kids = as.numeric(kids)) %>%
    filter(kids <= 2) %>%
    mutate(kids = as.factor(kids),
           marital = as.factor(marital)) %>%
    mutate(kids = fct_recode(kids, !!!kids_levels),
           marital = fct_recode(marital, !!!marry_levels))

  # Get value for minimum and maximum income levels for each bin
  df_alloc_cur <- df_alloc_cur %>%
    rowwise() %>%
    mutate(y_group_min = substring(strsplit(ymin_group, ",")[[1]][1], 2),
           y_group_max = gsub(strsplit(ymin_group, ",")[[1]][2],  pattern = "]", replacement = "")) %>%
    ungroup() %>%
    mutate(ymin_group = as.factor(ymin_group),
           kids = as.factor(kids),
           marital = as.factor(marital))

  # Get Optimal Allocation Result fro Graph
  # 1st period optimal c allocation
  df_alloc_cur_long <- df_alloc_cur %>%
    mutate(checks_c_optimal = optimal_c_1st,
           checks_v_optimal = optimal_v_1st) %>%
    select(-actual, -optimal_c_1st, -optimal_v_1st) %>%
    pivot_longer(cols = starts_with('checks'),
                 names_to = c('allocatetype'),
                 names_pattern = paste0("checks_(.*)"),
                 values_to = "checks")

  # Filtering data for allocation
  df_alloc_cur_long <- df_alloc_cur_long %>%
    mutate(y_group_min = as.numeric(y_group_min),
           checks = as.numeric(checks))

  # GGPLOT grahing ------
  # GGPLOT AES and FACET_WRAP
  plt_cur <- df_alloc_cur_long %>%
    ggplot(aes(x=y_group_min, y=checks*100,
               colour=allocatetype,
               shape=allocatetype)) +
    facet_wrap(~ marital + kids,
               ncol=3,
               scales = "free_x",
               labeller = label_wrap_gen(multi_line=FALSE))
  plt_cur <- plt_cur + geom_line() + geom_point(size=3)

  # GGPLOT legendings and coloring -------
  # X and Y labels and Legends etc.
  x.labels <- c('0', '50K', '100K', '150K', '200K')
  x.breaks <- c(0,
                50000/58056,
                100000/58056,
                150000/58056,
                200000/58056)

  # x-axis labeling
  plt_cur <- plt_cur +
    scale_x_continuous(labels = x.labels, breaks = x.breaks,
                       limits = c(0, 200000/58056))
  #
  # plt_cur <- plt_cur +
  #   scale_x_continuous(labels = x.labels, breaks = x.breaks,
  #                      limits = c(0, 175000/58056))
  # limits = c(0, 240000/58056)

  # Legend Labeling
  # \u2013 allows for en-dash
  ar_st_age_group_leg_labels <- c("C\u2013Optimal", "V\u2013Optimal")
  # color #F8766D is default red
  # color #F8766D is default bluish color

  if (fl_rho == 1) {
    st_color_c <- '#33cc33'
    it_shape <- 17
  } else if (fl_rho == -1) {
    st_color_c <- '#CD9600'
    it_shape <- 15
  } else if (fl_rho == -99) {
    st_color_c <- '#030e8c'
    it_shape <- 0
  } else {
    st_color_c <- '#808080'
    it_shape <- 3
  }
  ar_st_colours <- c(st_color_c, "#e75480")

  plt_cur <- plt_cur +
    theme(
      text = element_text(size = 12),
      legend.title = element_blank(),
      legend.position = c(0.14, 0.9),
      legend.background = element_rect(fill = "white", colour = "black", linetype='solid')) +
    scale_colour_manual(values = ar_st_colours, labels=ar_st_age_group_leg_labels) +
    scale_shape_manual(values = c(it_shape, 18), labels=ar_st_age_group_leg_labels)

  # X and Y Titling
  stg_x_title_hhinc <- 'Household income (thousands of 2012 USD)'
  stg_y_title_checks <- 'Stimulus check amount (2012 USD)'
  plt_cur <- plt_cur + labs(x = stg_x_title_hhinc, y = stg_y_title_checks)

  # Save Image Outputs -----
  # PNG and print
  if (bl_save_img) {
    png(paste0(srt_paper_appendix_textgraph, st_img_name_full),
        width = 270,
        height = 216, units='mm',
        res = 300, pointsize=7)
  }
  print(plt_cur)
  if (bl_save_img) {
    dev.off()
  }
}
