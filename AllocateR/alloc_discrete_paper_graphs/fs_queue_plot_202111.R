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

# Types of allocation files to consider
bl_per_capita <- TRUE
# For bush stimulus problem fl_rho = 1
fl_rho <- 1

# list to run
ls_st_file_suffix <- c('snwx_bushchck_moredense_a65zh266zs5_b1_xi0_manna_96')

# Solve iteratively
for (st_which_solu in ls_st_file_suffix) {

  # Max Phase Out given 1200*2 + 500*4 = 4400
  fl_max_phaseout_trump_biden = 225000
  fl_max_phaseout_bush = 225000
  # Meaning of Ymin Ymax simulated interval of 1
  fl_multiple_trump_biden = 62502
  # fl_multiple_trump_biden = 58056
  fl_multiple_bush = 54831
  if (grepl('bushchck', st_which_solu)) {
    fl_multiple <- fl_multiple_bush
    fl_max_phaseout <- fl_max_phaseout_bush
    st_stimulus_or_rebate <- "Tax rebate"

  } else {
    fl_multiple <- fl_multiple_trump_biden
    fl_max_phaseout <- fl_max_phaseout_trump_biden

    if (st_v_or_c_opti == "vopti") {
      st_stimulus_or_rebate <- "Check"
    } else if (st_v_or_c_opti == "copti") {
      st_stimulus_or_rebate <- "Stimulus check"
    }

  }

  # foreach (st_which_solu=ls_st_file_suffix) %dopar% {

  # Files:
  # source('fs_opti_support.R')
  # st_which_solu <- 'b1_manna'
  # st_which_solu <- paste0('b1_xi0_manna_88_', st_file_suffix)
  ls_output <- fs_opti_support_202111(st_which_solu, bl_per_capita=bl_per_capita, fl_rho=fl_rho)
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

  bl_budget_lines_billions <- TRUE
  if (bl_budget_lines_billions) {
    # If to add in lines for budgets
    # hand calculated from cost column in df_queue_il_long_c in Results202111\capi_rh1\20210225_bushchck_a65zh266zs5_b1_xi0_manna_96\csv\b1_a64
    fl_100billion <- 1778
    fl_125billion <- 2470
    fl_150billion <- 3500
    plt_rank <- plt_rank +
      geom_hline(aes(yintercept=fl_100billion, linetype="100bil"), color = "black") +
      geom_hline(aes(yintercept=fl_125billion, linetype="125bil"), color = "black") +
      geom_hline(aes(yintercept=fl_150billion, linetype="150bil"), color = "black") +
      scale_linetype_manual(name = "Budget", values=c(1,2,3), guide = guide_legend(ooverride.aes = list(color = c("black", "black", "black"))))
  }


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
                50000/fl_multiple,
                100000/fl_multiple,
                150000/fl_multiple,
                200000/fl_multiple)

  # x-axis labeling
  plt_rank <- plt_rank +
    scale_x_continuous(labels = x.labels, breaks = x.breaks,
                       limits = c(0, fl_max_phaseout/fl_multiple))
  # limits = c(0, 240000/fl_multiple)

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
      legend.position = c(0.59, 0.70),
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
