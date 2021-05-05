# Present graphically the relationship between per capita A and marginal alpha effects.
# This provides intuitions for the allocation problem across planners.
# Tested out graphing structure here: https://fanwangecon.github.io/R4Econ/tabgraph/ggscatter/htmlpdfr/fs_ggscatter_3cts_mdisc.html
# The code/loop structure follows AllocateR\alloc_discrete_fun_R\fs_mpc_tables_increments_202103.R
# fs_mpc_tables_increments_202103.R generates the A and alpha inputs needed here.

# try(dev.off(dev.list()["RStudioGD"]),silent=TRUE)
# try(dev.off(),silent=TRUE)

library(tidyverse)
library(REconTools)
# library(PrjOptiAlloc)

# library(forcats)
# library(foreach)
# library(doParallel)

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

bl_per_capita <- TRUE
fl_rho <- -1

# list to run
ls_st_file_suffix <- c(ls_st_file_suffix_bidenchk,
                       ls_st_file_suffix_trumpchk,
                       ls_st_file_suffix_bchklock)
# ls_st_file_suffix <- c(ls_st_file_suffix_trumpchk, ls_st_file_suffix_bidenchk,
#                        ls_st_file_suffix_bidenchk_betaedu, ls_st_file_suffix_trumpchk_betaedu)
# ls_st_file_suffix_test <- "snwx_bidenchk_moredense_a65zh266zs5_b1_xi0_manna_168"

ls_st_file_suffix <- c('snwx_bidenchk_moredense_a65zh266zs5_b1_xi0_manna_168')

# Solve iteratively
for (st_which_solu in ls_st_file_suffix) {
  # foreach (st_which_solu=ls_st_file_suffix) %dopar% {

  # Files:
  # source('fs_opti_support.R')
  # st_which_solu <- 'b1_manna'
  # st_which_solu <- paste0('b1_xi0_manna_88_', st_file_suffix)
  ls_output <- fs_opti_support_202103(st_which_solu, bl_per_capita=bl_per_capita, fl_rho=fl_rho)
  # Load folder names etc
  srt_results_root <- ls_output$srt_results_root
  bl_save_img <- ls_output$bl_save_img
  srt_folder <- ls_output$srt_folder
  bl_save_img <- TRUE
  srt_imgcsv_mpcapc_root <- ls_output$srt_imgcsv_mpcapc_root

  # file load
  srt_res_source_folder <- file.path(srt_folder, 'csv')
  srt_csv_allocate_subfolder <- paste0('b1_a64')
  srt_csv_file <- paste0('df_queue_il_long_c.csv')
  spn_csv_full_path <- file.path(srt_results_root, srt_res_source_folder, srt_csv_allocate_subfolder, srt_csv_file)

  # Load File
  df_queue_il_long_c <- as_tibble(read.csv(spn_csv_full_path, header=TRUE))


  # summarize
  # REconTools::ff_summ_percentiles(df_MPC_results, bl_statsasrows = FALSE)

  # # Gather the needed data
  # df_MPC_results_A_alpha <- df_MPC_results %>%
  #   select(X, marital, kids, ymin_group, mass, c_avg_chk0_usd, X0) %>%
  #   mutate(A = c_avg_chk0_usd/(1 + marital + kids),
  #          alpha = X0)

  # The continuous variables are A, alpha, mass. The categorical variables are kids and marital.
  # Modifying dataframe ------
  # Marital and Kids Level Labeling
  marry_levels <- c(Single = "0", Married = "1")
  kids_levels <- c("0" = "0", "1" = "1",
                   "2" = "2", "3" = "3",
                   "4" = "4")

  # ## Unique Kids Count
  # it_kids_marital_unique_n <- dim(as.matrix(df_alloc_combined
  #                                           %>% group_by(kids, marital) %>% summarize(freq=n())))[1]

  ## Select, as factor, and recode
  df_alloc_use <- df_queue_il_long_c %>%
    mutate(kids = as.numeric(kids)) %>%
    # filter(kids <= 2) %>%
    mutate(kids = as.factor(kids),
           marital = as.factor(marital)) %>%
    mutate(kids = fct_recode(kids, !!!kids_levels),
           marital = fct_recode(marital, !!!marry_levels))

  # Get value for minimum and maximum income levels for each bin
  df_alloc_use <- df_alloc_use %>%
    rowwise() %>%
    mutate(y_group_min = substring(strsplit(ymin_group, ",")[[1]][1], 2),
           y_group_max = gsub(strsplit(ymin_group, ",")[[1]][2],  pattern = "]", replacement = "")) %>%
    ungroup()

  # Generate factors
  df_alloc_use <- df_alloc_use %>%
    mutate(ymin_group = as.factor(ymin_group),
           y_group_min = as.numeric(y_group_min),
           kids = as.factor(kids),
           marital = as.factor(marital))

  # Transform scale
  df_alloc_use_mod <- df_alloc_use %>%
    group_by(id_i) %>%
    mutate(Q_il_start = min(Q_il))  %>% mutate(mass = mass*100)


  # Use R4Econ Scatter Code Template
  # Graphing
  plt_rank <-
    ggplot(df_alloc_use_mod, aes(x=y_group_min, y=Q_il)) +
    geom_point(aes(size=mass,
                    colour=kids)) +
    # geom_smooth(span = 0.50, se=FALSE) +
    theme(text = element_text(size = 16)) +
    facet_wrap(~ marital)

  # print(plt_rank)

  # Color controls
  # ar_st_colors <- c("#33cc33", "#F8766D")
  # ar_st_colors_label <- c("v-shaped", "straight")
  fl_legend_color_symbol_size <- 5
  st_leg_color_lab <- "Nr. of children"
  # Shape controls
  # ar_it_shapes <- c(9, 15)
  # ar_st_shapes_label <- c("auto", "manuel")
  fl_legend_shape_symbol_size <- 5
  st_leg_shape_lab <- "Marital status"
  # Control scatter point size
  fl_min_size <- 0.3125
  fl_max_size <- 3
  # fl_min_size <- 1.25
  # fl_max_size <- 12

  ar_size_range <- c(fl_min_size, fl_max_size)
  st_leg_size_lab <- "Pop. share (%)"

  # Labeling
  # st_title <- paste0('Distribution of HP and QSEC from mtcars')
  # st_subtitle <- paste0('https://fanwangecon.github.io/',
  #                       'R4Econ/tabgraph/ggscatter/htmlpdfr/fs_ggscatter_3cts_mdisc.html')
  # st_caption <- paste0('mtcars dataset, ',
  #                      'https://fanwangecon.github.io/R4Econ/')
  # if (bl_log_x_axis) {
    st_x_label <- 'Start queue position'
  # } else {
    st_x_label <- 'Minimum household income given income bin'
  # }
  st_y_label <- 'Queue positions'

  x.labels <- c('0', '50K', '100K', '150K', '200K')
  x.breaks <- c(0,
                50000/58056,
                100000/58056,
                150000/58056,
                200000/58056)

  # x-axis labeling
  plt_rank <- plt_rank +
    scale_x_continuous(labels = x.labels, breaks = x.breaks,
                       limits = c(0, 200000/58056))
  # limits = c(0, 240000/58056)

  # Add titles and labels
  plt_rank <- plt_rank +
    labs(x = st_x_label, y = st_y_label)

  # Color, shape and size controls
  plt_rank <- plt_rank +
    # scale_colour_manual(values=ar_st_colors, labels=ar_st_colors_label) +
    # scale_shape_manual(values=ar_it_shapes, labels=ar_st_shapes_label) +
    scale_size_continuous(range = ar_size_range)

  # replace the default labels for each legend segment
  plt_rank <- plt_rank +
    labs(colour = st_leg_color_lab,
         shape = st_leg_shape_lab,
         size = st_leg_size_lab) +
    theme(
      text = element_text(size = 14),
      legend.position = c(0.41, 0.21),
      legend.background = element_rect(fill = "white", colour = "black", linetype='solid'),
      legend.key.width = unit(2.5, "line"))


  # Control the order of legend display
  # Show color, show shape, then show size.
  plt_rank <- plt_rank + guides(
    # shape = guide_legend(order = 1, override.aes = list(size = fl_legend_shape_symbol_size)),
    colour = guide_legend(order = 2, override.aes = list(size = fl_legend_color_symbol_size)),
    size = guide_legend(order = 3))

  # x-axis
  # if (bl_log_x_axis) {
  #   print('do not label x-axis, not in level')
  # } else {
  #   x.labels <- c('0', '25K', '50k', '75k', '100K', '125K', '150K')
  #   x.breaks <- c(0,
  #                 25000,
  #                 50000,
  #                 75000,
  #                 100000,
  #                 125000,
  #                 150000)
  #   # x-axis labeling
  #   plt_rank <- plt_rank +
  #     scale_x_continuous(labels = x.labels, breaks = x.breaks)
  # }

  # Graph
  if (bl_save_img) {

    # add png
    snm_save_png <- paste0('queue_rank.png')
    ggsave(plt_rank,
           file=file.path(srt_imgcsv_mpcapc_root, snm_save_png),
           width = 270,
           height = 216, units='mm',
           dpi = 300)
  }

}
