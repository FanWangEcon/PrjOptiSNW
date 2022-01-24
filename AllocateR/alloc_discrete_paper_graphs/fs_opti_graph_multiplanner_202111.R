# Graphs where allocation lines from multiple planners are jointly compared.
# copied over from AllocateR\alloc_discrete_paper_graphs\fs_opti_graph_multiplanner.R

try(dev.off(dev.list()["RStudioGD"]),silent=TRUE)
try(dev.off(),silent=TRUE)

# Libraries
library(tidyverse)
library(REconTools)
library(scales)

# Parameters and Options --------

# Results from which planners should
ls_it_rho_combo_type <- c(1,2,3,4,5,6)

# Types of allocation files to consider
# Per capita or per household results
ls_bl_per_capita <- c(TRUE)

# Allocation bounds types for Trump/Biden Checks
ls_st_file_suffix <- c("snwx_bchklock_moredense_a65zh266zs5_b1_xi0_manna_168")
ls_st_bound_files <- c('14ca14ck', '17ca17ck', '20ca20ck', '20ca20ck', '17k',
                       '17cadt', '17ckid', '20cadt', '20ckid')
ls_st_bound_files <- c('20ca20ck', '20cadt', '20ckid', '14ca14ck')

# # Allocation bounds types for Bush checks
# ls_st_file_suffix <- c("snwx_bushchck_moredense_a65zh266zs5_b1_xi0_manna_96")
# ls_st_bound_files <- c('6ca3ck', '9ca6ck', '12ca9ck', '12cadt', '9cadt')

# Loop over file types1
for (st_v_or_c_opti in c('vopti','copti')) {
  for (bl_per_capita in ls_bl_per_capita) {
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

      for (st_bound_files in ls_st_bound_files) {
        for (it_rho_combo_type in ls_it_rho_combo_type) {

          # Four allocation graph display types
          # 1. Two lines: Rawlsian and Utilitarian.
          # 2. Three lines: Rawlsian, Utilitarian, and actual policy.
          # 3. Three lines: Rawlsian, Utilitarian, and benhcmark (intermediate=-1).
          # 4. Four lines: Rawlsian, Utilitarian, benchmark, and actual policy.
          fl_rho_val_select <- 1
          if (it_rho_combo_type == 1) {
            ar_rho <- c(1, -99)
            bl_include_actual <- FALSE

          } else if (it_rho_combo_type == 2) {
            ar_rho <- c(1, -99)
            bl_include_actual <- TRUE

          } else if (it_rho_combo_type == 3) {
            ar_rho <- c(1, -1, -99)
            bl_include_actual <- FALSE

          } else if (it_rho_combo_type == 4) {
            ar_rho <- c(1, -1, -99)
            bl_include_actual <- TRUE

          } else if (it_rho_combo_type == 5) {
            fl_rho_val_select <- -1
            ar_rho <- c(-1)
            bl_include_actual <- TRUE

          } else if (it_rho_combo_type == 6) {
            ar_rho <- c(1)
            bl_include_actual <- TRUE
          }

          # Generate files by aggregating over rho types
          # Results counter
          it_results_ctr <- 0
          # Loop over results
          for (fl_rho in ar_rho) {

            # Counter updating
            it_results_ctr <- it_results_ctr + 1

            # Load support file
            ls_output <- fs_opti_support_202111(st_which_solu, bl_per_capita=bl_per_capita, fl_rho=fl_rho)
            # Load folder names etc
            srt_results_root <- ls_output$srt_results_root
            bl_save_img <- ls_output$bl_save_img
            srt_folder <- ls_output$srt_folder
            bl_save_img <- TRUE
            srt_paper_appendix_textgraph <- ls_output$srt_paper_appendix_textgraph

            # folder that is not rho specific for storing aggregate folder
            srt_img_aggregator_save_root <- ls_output$srt_img_aggregator_save_root

            # Get Last folder which is per-capita/hhousehold and rho specific
            st_hhpercapita_rho <- tail(strsplit(srt_results_root, "/")[[1]],n=1)

            # Get folder name path
            srt_res_source_folder <- file.path(srt_folder, 'csv')
            srt_csv_allocate_subfolder <- paste0('b1_a64_', st_bound_files, '_18t64')
            srt_csv_file <- paste0('df_alloc_all_feasible_', st_bound_files, '.csv')
            spn_csv_full_path <- file.path(srt_results_root, srt_res_source_folder, srt_csv_allocate_subfolder, srt_csv_file)

            # Clean up each file
            svr_opti_alloc_col <- paste0('checks_optimal_', st_hhpercapita_rho)
            print(svr_opti_alloc_col)

            if (st_v_or_c_opti == "vopti") {
              df_alloc_cur <- as_tibble(read.csv(spn_csv_full_path)) %>%
                mutate(!!sym(svr_opti_alloc_col) := optimal_v_1st)
            } else if (st_v_or_c_opti == "copti") {
              df_alloc_cur <- as_tibble(read.csv(spn_csv_full_path)) %>%
                mutate(!!sym(svr_opti_alloc_col) := optimal_c_1st)
            }

            # Merging FIles Different Allocations are Additional Columns
            if (it_results_ctr == 1) {

              if (st_v_or_c_opti == "vopti") {
                df_alloc_combined <- df_alloc_cur %>% select(-X, -optimal_v_1st) %>% ungroup()
              } else if (st_v_or_c_opti == "copti") {
                df_alloc_combined <- df_alloc_cur %>% select(-X, -optimal_c_1st) %>% ungroup()
              }

              ls_st_opti_coll <- c(svr_opti_alloc_col)

            } else {
              df_alloc_combined <- df_alloc_combined %>%
                left_join(df_alloc_cur %>%
                            select(id_i, svr_opti_alloc_col) %>%
                            ungroup(), by='id_i')
              ls_st_opti_coll <- c(ls_st_opti_coll, svr_opti_alloc_col)
            }
          }

          # graphing
          # G4 Figure in Table, Figure Number 4 Main Paper, Double both Kids and Adult Max
          st_img_name_full <- paste0(srt_csv_allocate_subfolder,
            '_v', it_rho_combo_type, '_', st_v_or_c_opti, '.png')

          # Modifying dataframe ------
          # Marital and Kids Level Labeling
          marry_levels <- c(Single = "0", Married = "1")
          kids_levels <- c("no children" = "0", "one child" = "1",
                           "two children" = "2", "three children" = "3",
                           "four children" = "4")

          ## Unique Kids Count
          it_kids_marital_unique_n <- dim(as.matrix(df_alloc_combined
                                                    %>% group_by(kids, marital) %>% summarize(freq=n())))[1]

          ## Select, as factor, and recode
          df_alloc_use <- df_alloc_combined %>% mutate(kids = as.numeric(kids)) %>%
            filter(kids <= 2) %>%
            mutate(kids = as.factor(kids),
                   marital = as.factor(marital)) %>%
            mutate(kids = fct_recode(kids, !!!kids_levels),
                   marital = fct_recode(marital, !!!marry_levels))

          # Get value for minimum and maximum income levels for each bin
          df_alloc_use <- df_alloc_use %>%
            filter(rho_val == fl_rho_val_select) %>% rowwise() %>%
            mutate(y_group_min = substring(strsplit(ymin_group, ",")[[1]][1], 2),
                   y_group_max = gsub(strsplit(ymin_group, ",")[[1]][2],  pattern = "]", replacement = "")) %>%
            ungroup() %>%
            mutate(ymin_group = as.factor(ymin_group),
                   kids = as.factor(kids),
                   marital = as.factor(marital))

          # Get Optimal Allocation Result fro Graph
          # 1st period optimal c allocation
          if (bl_include_actual) {
            # include checks_actual, with checks prefix kept in transformation
            df_alloc_use <- df_alloc_use %>% mutate(checks_actual = actual)
          }

          if (st_v_or_c_opti == "vopti") {
            df_alloc_cur_long <- df_alloc_use %>% select(-actual, -optimal_c_1st)
          } else if (st_v_or_c_opti == "copti") {
            df_alloc_cur_long <- df_alloc_use %>% select(-actual, -optimal_v_1st)
          }

          df_alloc_cur_long <- df_alloc_cur_long %>%
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
                        50000/fl_multiple,
                        100000/fl_multiple,
                        150000/fl_multiple,
                        200000/fl_multiple)
          x.bounds <- c(0, fl_max_phaseout/fl_multiple)

          y.labels <- c('0', '2000', '4000', '6000', '8000')
          y.breaks <- c(0,
                        2000,
                        4000,
                        6000,
                        8000)
          y.bounds <- c(0, 8000)

          # x-axis labeling
          plt_cur <- plt_cur +
            scale_x_continuous(labels = x.labels, breaks = x.breaks,
                               limits = x.bounds)
          # limits = c(0, 240000/fl_multiple)

          # Legend Labeling
          # \u2013 allows for en-dash
          if (st_v_or_c_opti == "vopti") {
            st_efficient_planner <- 'lambda=1'
          } else if (st_v_or_c_opti == "copti") {
            st_efficient_planner <- 'Stimulus'
          }
          st_actual <- "Actual"
          st_actual_color <- "#F8766D"
          it_actual_shape <- 16

          st_util <- st_efficient_planner
          st_util_color <- "#33cc33"
          it_util_shape <- 17

          st_inter <- "Intermediate"
          st_inter_color <- "#CD9600"
          it_inter_shape <- 15

          st_rawl <- "Rawlsian"
          st_rawl_color <- "#030e8c"
          it_rawl_shape <- 0

          if (it_rho_combo_type == 1) {
            ar_st_age_group_leg_labels <- c(st_util, st_rawl)
            ar_st_colours <- c(st_util_color, st_rawl_color)
            ar_it_shape_val <- c(it_util_shape, it_rawl_shape)

          } else if (it_rho_combo_type == 2) {
            ar_st_age_group_leg_labels <- c(st_actual, st_util, st_rawl)
            ar_st_colours <- c(st_actual_color, st_util_color, st_rawl_color)
            ar_it_shape_val <- c(it_actual_shape, it_util_shape, it_rawl_shape)

          } else if (it_rho_combo_type == 3) {
            ar_st_age_group_leg_labels <- c(st_util, st_inter, st_rawl)
            ar_st_colours <- c(st_util_color, st_inter_color, st_rawl_color)
            ar_it_shape_val <- c(it_util_shape, it_inter_shape, it_rawl_shape)

          } else if (it_rho_combo_type == 4) {
            ar_st_age_group_leg_labels <- c(st_actual, st_util, st_inter, st_rawl)
            ar_st_colours <- c(st_actual_color, st_util_color, st_inter_color, st_rawl_color)
            ar_it_shape_val <- c(it_actual_shape, it_util_shape, it_inter_shape, it_rawl_shape)

          } else if (it_rho_combo_type == 5) {
            ar_st_age_group_leg_labels <- c(st_actual, st_inter)
            ar_st_colours <- c(st_actual_color, st_inter_color)
            ar_it_shape_val <- c(it_actual_shape, it_inter_shape)

          } else if (it_rho_combo_type == 6) {
            ar_st_age_group_leg_labels <- c(st_actual, st_util)
            ar_st_colours <- c(st_actual_color, st_util_color)
            ar_it_shape_val <- c(it_actual_shape, it_util_shape)

          }

          # custom theme
          theme_custom <- theme(
            text = element_text(size = 12),
            legend.title = element_blank(),
            legend.position = c(0.14, 0.9),
            legend.background = element_rect(fill = "white", colour = "black", linetype='solid'))
          # custom theme for single plot
          theme_custom_single <- theme(
            text = element_text(size = 12),
            legend.title = element_blank(),
            legend.position = c(0.35, 0.85),
            legend.background = element_rect(fill = "white", colour = "black", linetype='solid'))
          theme_custom_single_nolabels <- theme(
            text = element_text(size = 12),
            legend.title = element_blank(),
            legend.position = "none",
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())

          # color #F8766D is default red
          # color #F8766D is default bluish color
          plt_cur <- plt_cur + theme_custom +
            scale_colour_manual(values=ar_st_colours, labels=ar_st_age_group_leg_labels) +
            scale_shape_manual(values=ar_it_shape_val, labels=ar_st_age_group_leg_labels)

          # X and Y Titling
          stg_x_title_hhinc <- 'Household income (thousands of 2012 USD)'
          stg_y_title_checks <- paste0(st_stimulus_or_rebate, ' amount (2012 USD)')
          plt_cur <- plt_cur + labs(x = stg_x_title_hhinc, y = stg_y_title_checks)

          # Image height
          it_img_height <- 216
          if (it_kids_marital_unique_n <= 5) {
            # single line if only married or unmarried
            it_img_height <- it_img_height*0.55
          } else {
            # standard two lines solving problem jointly
            it_img_height <- it_img_height
          }

          # Save Image Outputs -----
          # PNG and print
          if (bl_save_img) {
            # png(paste0(srt_img_aggregator_save_root,
            #            st_img_name_full),
            #     width = 270,
            #     height = 216, units='mm',
            #     res = 300, pointsize=7)
            ggsave(plt_cur, file=paste0(srt_img_aggregator_save_root,
                          st_img_name_full),
                   width = 270,
                   height = it_img_height, units='mm',
                   dpi = 300, pointsize=7)
          }

  #         if (bl_save_img) {
  #           dev.off()
  #         }

          for (it_subplot_as_own_vsr in c(1,2)) {

            if (it_subplot_as_own_vsr == 1) {
              theme_custom_single_use <- theme_custom_single
              st_file_suffix <- '_haslegend'
              it_width <- 100
            } else if (it_subplot_as_own_vsr == 2) {
              theme_custom_single_use <- theme_custom_single_nolabels
              st_file_suffix <- '_nolegend'
              it_width <- 87
            }

            # Generate 1 by 1 figures ----------------------------
            df_alloc_cur_long_grps <- df_alloc_cur_long %>%
              group_by(marital, kids) %>% summarize(checks=mean(checks)) %>% ungroup()
            ls_plots <- apply(df_alloc_cur_long_grps, 1, function(ar_row_unique) {
              # 1. Graph main
              # plt_mtcars_scatter <-
              #   ggplot(tb_mtcars %>% filter(vs == st_vs_cate),
              #          aes(x=hp, y=qsec,
              #              colour=am, shape=am, linetype=am)) +
              #   geom_smooth(se = FALSE, lwd = 1.5) + # Lwd = line width
              #   geom_point(size = 5, stroke = 2)
              plt_cur <- df_alloc_cur_long %>%
                filter(marital == ar_row_unique[1], kids == ar_row_unique[2]) %>%
                ggplot(aes(x=y_group_min, y=checks*100,
                           colour=allocatetype,
                           shape=allocatetype)) +
                facet_wrap(~ marital + kids,
                           ncol=1,
                           scales = "free_x",
                           labeller = label_wrap_gen(multi_line=FALSE))
              plt_cur <- plt_cur + geom_line() + geom_point(size=3)

              # 2. Add titles and labels
              plt_cur <- plt_cur + labs(x = stg_x_title_hhinc, y = stg_y_title_checks)

              # 3. x and y ticks
              plt_cur <- plt_cur +
                scale_x_continuous(labels = x.labels, breaks = x.breaks,
                                   limits = x.bounds) +
                scale_y_continuous(labels = y.labels, breaks = y.breaks,
                                   limits = y.bounds)

              # 4. Color, shape and linetype controls
              plt_cur <- plt_cur +
                scale_colour_manual(values=ar_st_colours, labels=ar_st_age_group_leg_labels) +
                scale_shape_manual(values=ar_it_shape_val, labels=ar_st_age_group_leg_labels)

              # 5. replace the default labels for each legend segment)
              plt_cur <- plt_cur + theme_custom_single_use
            })

            it_fig_ctr <- 0
            for (plt_cur in ls_plots) {
              if (bl_save_img) {
                # png(paste0(srt_img_aggregator_save_root,
                #            st_img_name_full),
                #     width = 270,
                #     height = 216, units='mm',
                #     res = 300, pointsize=7)
                it_fig_ctr <- it_fig_ctr + 1

                srt_img_aggcsvsolu_save_subfolder <- file.path(srt_img_aggregator_save_root,
                                                               paste0(srt_csv_allocate_subfolder,
                                                                      '_v', it_rho_combo_type,
                                                                      '_', st_v_or_c_opti), '/')
                dir.create(file.path(srt_img_aggcsvsolu_save_subfolder), showWarnings = FALSE, recursive = TRUE)

                st_img_name_full_subfigure <- paste0(srt_csv_allocate_subfolder,
                                                     '_v', it_rho_combo_type,
                                                     '_', st_v_or_c_opti,
                                                     '_subfig', it_fig_ctr, st_file_suffix,
                                                     '.png')
                ggsave(plt_cur, file=paste0(srt_img_aggcsvsolu_save_subfolder,
                                            st_img_name_full_subfigure),
                       width = it_width,
                       height = 118.8, units='mm',
                       dpi = 150, pointsize=7)
              }
            }
          }

        }
      }
    }
  }
}
