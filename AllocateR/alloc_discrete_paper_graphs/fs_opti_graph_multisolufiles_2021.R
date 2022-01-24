# This is developed after fs_opti_graph_multiplanner.m, Rather than showing fl_multiple
# planners together given the same csv-solution file, now show for a single planner type
# but combine results frorm multiple csv-solution files. Specifically, files with Different
# mixture assumptions.

# try(dev.off(dev.list()["RStudioGD"]),silent=TRUE)
# try(dev.off(),silent=TRUE)

# Libraries
library(tidyverse)
library(REconTools)
library(scales)

# Parameters and Options --------
# Possible rho values to select from
# ar_rho <- c(1,
#             0.9, 0.8,
#             0.05, -1,
#             -9, -19, -99)
ar_rho <- c(-1)

# Results from which planners should
ls_it_csvsolu_combo_type <- c(1,2,3,4, 11,12, 21,22)
# ls_it_csvsolu_combo_type <- c(21,22)

# Per capita or per household results
ls_bl_per_capita <- c(TRUE)

# Allocation bounds types
# ls_st_bound_files <- c('14ca14ck', '17ca17ck', '20ca20ck', '20ca20ck', '17k',
#                        '17cadt', '17ckid', '20cadt', '20ckid')
ls_st_bound_files <- c('14ca14ck', '17ca17ck', '20ca20ck')


# Max Phase Out given 1200*2 + 500*4 = 4400
fl_max_phaseout_trump_biden = 225000
# Meaning of Ymin Ymax simulated interval of 1
fl_multiple_trump_biden = 62502
fl_max_phaseout <- fl_max_phaseout_trump_biden
fl_multiple <- fl_multiple_trump_biden

# Loop over file types1
for (st_v_or_c_opti in c('vopti','copti')) {
  for (bl_per_capita in ls_bl_per_capita) {
    for (st_bound_files in ls_st_bound_files) {
      for (fl_rho in ar_rho) {
        for (it_csvsolu_combo_type in ls_it_csvsolu_combo_type) {

          if (it_csvsolu_combo_type >= 1 & it_csvsolu_combo_type <= 10) {
            st_subfolder <- 'betaheter'
          } else if (it_csvsolu_combo_type >= 11 & it_csvsolu_combo_type <= 20) {
            st_subfolder <- 'uinoui'
          } else if (it_csvsolu_combo_type >= 21 & it_csvsolu_combo_type <= 30) {
            st_subfolder <- 'saveRlowhigh'
          }

          # A. beta heterogeneity -----------------------------------------
          # Four allocation graph display types
          # 1. Two lines: Rawlsian and Utilitarian.
          # 2. Three lines: Rawlsian, Utilitarian, and actual policy.
          # 3. Three lines: Rawlsian, Utilitarian, and benhcmark (intermediate=-1).
          # 4. Four lines: Rawlsian, Utilitarian, benchmark, and actual policy.
          if (it_csvsolu_combo_type == 1) {
            ls_st_file_name_short <- c('mix4o6', 'mix6o6', 'mix8o6')
            ls_st_file_suffix_bchklock_mixturealter <-
              c('snwx_bchklock_moredense_a65zh266zs5_b1_xi0_manna_168_4o6',
                'snwx_bchklock_moredense_a65zh266zs5_b1_xi0_manna_168',
                'snwx_bchklock_moredense_a65zh266zs5_b1_xi0_manna_168_8o6')
            bl_include_actual <- FALSE
            st_save_suffix <- '_4o6o8noactual'

          } else if (it_csvsolu_combo_type == 2) {
            ls_st_file_name_short <- c('mix4o6', 'mix6o6', 'mix8o6')
            ls_st_file_suffix_bchklock_mixturealter <-
              c('snwx_bchklock_moredense_a65zh266zs5_b1_xi0_manna_168_4o6',
                'snwx_bchklock_moredense_a65zh266zs5_b1_xi0_manna_168',
                'snwx_bchklock_moredense_a65zh266zs5_b1_xi0_manna_168_8o6')
            bl_include_actual <- TRUE
            st_save_suffix <- '_4o6o8'

          } else if (it_csvsolu_combo_type == 3) {
            ls_st_file_name_short <- c('mix3o6', 'mix6o6', 'mix9o6')
            ls_st_file_suffix_bchklock_mixturealter <-
              c('snwx_bchklock_moredense_a65zh266zs5_b1_xi0_manna_168_3o6',
                'snwx_bchklock_moredense_a65zh266zs5_b1_xi0_manna_168',
                'snwx_bchklock_moredense_a65zh266zs5_b1_xi0_manna_168_9o6')
            bl_include_actual <- FALSE
            st_save_suffix <- '_3o6o9noactual'

          } else if (it_csvsolu_combo_type == 4) {
            ls_st_file_name_short <- c('mix3o6', 'mix6o6', 'mix9o6')
            ls_st_file_suffix_bchklock_mixturealter <-
              c('snwx_bchklock_moredense_a65zh266zs5_b1_xi0_manna_168_3o6',
                'snwx_bchklock_moredense_a65zh266zs5_b1_xi0_manna_168',
                'snwx_bchklock_moredense_a65zh266zs5_b1_xi0_manna_168_9o6')
            bl_include_actual <- TRUE
            st_save_suffix <- '_3o6o9'
          }

          # B. With and without ui benefits -----------------------------------------
          if (it_csvsolu_combo_type == 11) {
            ls_st_file_name_short <- c('bcklknou', 'bchklock')
            ls_st_file_suffix_bchklock_mixturealter <-
              c('snwx_bcklknou_moredense_a65zh266zs5_b1_xi0_manna_168',
                'snwx_bchklock_moredense_a65zh266zs5_b1_xi0_manna_168')
            bl_include_actual <- FALSE
            st_save_suffix <- '_uinouiactual'

          } else if (it_csvsolu_combo_type == 12) {
            ls_st_file_name_short <- c('bcklknou', 'bchklock')
            ls_st_file_suffix_bchklock_mixturealter <-
              c('snwx_bcklknou_moredense_a65zh266zs5_b1_xi0_manna_168',
                'snwx_bchklock_moredense_a65zh266zs5_b1_xi0_manna_168')
            bl_include_actual <- TRUE
            st_save_suffix <- '_uinoui'
          }

          # C. With and without lockdown -----------------------------------------
          if (it_csvsolu_combo_type == 21) {
            ls_st_file_name_short <- c('bchklkr2', 'bchklock')
            ls_st_file_suffix_bchklock_mixturealter <-
              c('snwx_bchklkr2_moredense_a65zh266zs5_b1_xi0_manna_168',
                'snwx_bchklock_moredense_a65zh266zs5_b1_xi0_manna_168')
            bl_include_actual <- FALSE
            st_save_suffix <- '_lockactual'

          } else if (it_csvsolu_combo_type == 22) {
            ls_st_file_name_short <- c('bchklkr2', 'bchklock')
            ls_st_file_suffix_bchklock_mixturealter <-
              c('snwx_bchklkr2_moredense_a65zh266zs5_b1_xi0_manna_168',
                'snwx_bchklock_moredense_a65zh266zs5_b1_xi0_manna_168')
            bl_include_actual <- TRUE
            st_save_suffix <- '_lock'
          }

          ls_st_file_suffix <- ls_st_file_suffix_bchklock_mixturealter

          it_results_ctr <- 0
          for (st_which_solu in ls_st_file_suffix) {

            # Counter updating
            it_results_ctr <- it_results_ctr + 1
            st_file_name_short <- ls_st_file_name_short[it_results_ctr]

            # Generate files by aggregating over rho types
            # Results counter
            # Loop over results

            # Load support file
            ls_output <- fs_opti_support_202111(st_which_solu, bl_per_capita=bl_per_capita, fl_rho=fl_rho)
            # Load folder names etc
            srt_results_root <- ls_output$srt_results_root
            bl_save_img <- ls_output$bl_save_img
            srt_folder <- ls_output$srt_folder
            bl_save_img <- TRUE
            srt_paper_appendix_textgraph <- ls_output$srt_paper_appendix_textgraph

            # folder that is not rho specific for storing aggregate folder
            srt_img_aggcsvsolu_save_root <- ls_output$srt_img_aggcsvsolu_save_root

            # Get Last folder which is per-capita/hhousehold and rho specific
            # st_hhpercapita_rho <- tail(strsplit(srt_results_root, "/")[[1]],n=1)

            # Get folder name path
            srt_res_source_folder <- file.path(srt_folder, 'csv')
            srt_csv_allocate_subfolder <- paste0('b1_a64_', st_bound_files, '_18t64')
            srt_csv_file <- paste0('df_alloc_all_feasible_', st_bound_files, '.csv')
            spn_csv_full_path <- file.path(srt_results_root, srt_res_source_folder, srt_csv_allocate_subfolder, srt_csv_file)

            # Clean up each file
            svr_opti_alloc_col <- paste0('checks_optimal_', st_file_name_short)
            print(svr_opti_alloc_col)
            if (st_v_or_c_opti == "vopti") {
              st_stimulus_or_rebate <- "Check"
              df_alloc_cur <- as_tibble(read.csv(spn_csv_full_path)) %>%
                mutate(!!sym(svr_opti_alloc_col) := optimal_v_1st)
            } else if (st_v_or_c_opti == "copti") {
              st_stimulus_or_rebate <- "Stimulus check"
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

          # Generate policy specific subfolder, beta heter, uinoui or lockdown
          srt_img_aggcsvsolu_save_root <- file.path(srt_img_aggcsvsolu_save_root, st_subfolder, '/')
          dir.create(file.path(srt_img_aggcsvsolu_save_root), showWarnings = FALSE, recursive = TRUE)

          # graphing
          # G4 Figure in Table, Figure Number 4 Main Paper, Double both Kids and Adult Max
          st_img_name_full <- paste0(srt_csv_allocate_subfolder, st_save_suffix, '_', st_v_or_c_opti, '.png')

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
            rowwise() %>%
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

          # Legend Labeling
          # \u2013 allows for en-dash
          st_actual <- "Actual"
          st_actual_color <- "#F8766D"
          it_actual_shape <- 16

          st_util_color <- "#33cc33"
          it_util_shape <- 17

          st_inter <- "Benchmark model"
          st_inter_color <- "#CD9600"
          it_inter_shape <- 15

          st_rawl_color <- "#030e8c"
          it_rawl_shape <- 0


          # ZA, beta heterogeneity results --------------------------------
          st_util <- "Fewer impatient households"
          st_rawl <- "More impatient households"
          if (it_csvsolu_combo_type == 1) {
            # st_util <- "Lower low beta (2/3)"
            # st_rawl <- "Higher low beta (4/3)"
            ar_st_age_group_leg_labels <- c(st_util, st_inter, st_rawl)
            ar_st_colours <- c(st_util_color, st_inter_color, st_rawl_color)
            ar_it_shape_val <- c(it_util_shape, it_inter_shape, it_rawl_shape)

          } else if (it_csvsolu_combo_type == 2) {
            # st_util <- "Lower low beta (2/3)"
            # st_rawl <- "Higher low beta (4/3)"
            ar_st_age_group_leg_labels <- c(st_actual, st_util, st_inter, st_rawl)
            ar_st_colours <- c(st_actual_color, st_util_color, st_inter_color, st_rawl_color)
            ar_it_shape_val <- c(it_actual_shape, it_util_shape, it_inter_shape, it_rawl_shape)

          } else if (it_csvsolu_combo_type == 3) {
            # st_util <- "Lower low beta (1/2)"
            # st_rawl <- "Higher low beta (3/2)"
            ar_st_age_group_leg_labels <- c(st_util, st_inter, st_rawl)
            ar_st_colours <- c(st_util_color, st_inter_color, st_rawl_color)
            ar_it_shape_val <- c(it_util_shape, it_inter_shape, it_rawl_shape)

          } else if (it_csvsolu_combo_type == 4) {
            # st_util <- "Lower low beta (1/2)"
            # st_rawl <- "Higher low beta (3/2)"
            ar_st_age_group_leg_labels <- c(st_actual, st_util, st_inter, st_rawl)
            ar_st_colours <- c(st_actual_color, st_util_color, st_inter_color, st_rawl_color)
            ar_it_shape_val <- c(it_actual_shape, it_util_shape, it_inter_shape, it_rawl_shape)

          }

          # ZB, UI or NO UI results --------------------------------
          st_noui <- "Without UI (b=0)"
          st_withui <- "With UI (b=1)"
          if (it_csvsolu_combo_type == 11) {
            ar_st_age_group_leg_labels <- c(st_noui, st_withui)
            ar_st_colours <- c(st_util_color, st_inter_color)
            ar_it_shape_val <- c(it_util_shape, it_inter_shape)

          } else if (it_csvsolu_combo_type == 12) {
            ar_st_age_group_leg_labels <- c(st_actual, st_noui, st_withui)
            ar_st_colours <- c(st_actual_color, st_util_color, st_inter_color)
            ar_it_shape_val <- c(it_actual_shape, it_util_shape, it_inter_shape)

          }

          # ZC, Lockdown or no lockdown results --------------------------------
          st_low_r <- "With r=0.02"
          st_high_r <- "With r=0.04"
          if (it_csvsolu_combo_type == 21) {
            ar_st_age_group_leg_labels <- c(st_low_r, st_high_r)
            ar_st_colours <- c(st_util_color, st_inter_color)
            ar_it_shape_val <- c(it_util_shape, it_inter_shape)

          } else if (it_csvsolu_combo_type == 22) {
            ar_st_age_group_leg_labels <- c(st_actual, st_low_r, st_high_r)
            ar_st_colours <- c(st_actual_color, st_util_color, st_inter_color)
            ar_it_shape_val <- c(it_actual_shape, it_util_shape, it_inter_shape)

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
            ggsave(plt_cur, file=paste0(srt_img_aggcsvsolu_save_root,
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

                srt_img_aggcsvsolu_save_subfolder <- file.path(srt_img_aggcsvsolu_save_root,
                                                               paste0(srt_csv_allocate_subfolder,
                                                                      st_save_suffix,
                                                                      '_', st_v_or_c_opti),
                                                               '/')
                dir.create(file.path(srt_img_aggcsvsolu_save_subfolder), showWarnings = FALSE, recursive = TRUE)

                st_img_name_full_subfigure <- paste0(srt_csv_allocate_subfolder,
                                                     st_save_suffix,
                                                     '_', st_v_or_c_opti,
                                                     '_subfig', it_fig_ctr, st_file_suffix, '.png')
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
