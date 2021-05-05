# Developed based on *fs_opti_graph_multiplanner.R*, aggregate across subfolders
# generated from the same data file, same planner, but different perturbed, with
# varying shock draws.
try(dev.off(dev.list()["RStudioGD"]),silent=TRUE)
try(dev.off(),silent=TRUE)

# Libraries
library(tidyverse)
library(REconTools)
library(scales)

# Parameters and Options --------
# graph types
ls_it_graph_types <- c(1, 2, 3, 4, 5)
# snm_main folder
snm_main <- 'Results202104'
# Possible rho values to select from
ar_rho <- c(-1)
# Types of allocation files to consider
ls_st_file_suffix <- c('snwx_bidenchk_moredense_a65zh266zs5_b1_xi0_manna_168')
# Per capita or per household results
ls_bl_per_capita <- c(TRUE)
# Allocation bounds types
ls_st_bound_files <- c('14ca14ck', '17ca17ck', '20ca20ck')
# starting seed and pertrub sequence
it_rand_adj_A_rng_seed_start <- 123
ar_it_perturb_n <- seq(1, 500, by=1)
bl_gen_agg_folders <- TRUE
# aggregation folder for store
ls_st_gen_agg_folders <- c('perturb')

# Parameters and Options not perturbed original --------------------
# Also allow loading in the original non-perturbed file for easier graphing and joint plotting
snm_main_noperturb <- 'Results202103'
it_disturb <- 99999

# Loop over file types1
for (bl_per_capita in ls_bl_per_capita) {
  for (st_which_solu in ls_st_file_suffix) {
    for (st_bound_files in ls_st_bound_files) {
      for (fl_rho in ar_rho) {

        # Load in NOT perturbed file -------------------
        ls_output_noperturb <- fs_opti_support_202103(st_which_solu,
                                                      bl_per_capita=bl_per_capita,
                                                      fl_rho=fl_rho,
                                                      snm_main=snm_main_noperturb,
                                                      it_rand_draw_seed=NULL,
                                                      ls_st_gen_agg_folders=c(''))
        # Get folder name path
        spn_csv_full_path_noperturb <- file.path(ls_output_noperturb$srt_results_root,
                                                 ls_output_noperturb$srt_folder, 'csv',
                                                 paste0('b1_a64_', st_bound_files, '_18t64'),
                                                 paste0('df_alloc_all_feasible_', st_bound_files, '.csv'))
        # rng = 99999 is the not perturbed result.
        svr_opti_alloc_col <- paste0('checks_optimal_capi_rhn1_perturb_rng', it_disturb)
        df_alloc_cur_optimal_nopertrub <- as_tibble(read.csv(spn_csv_full_path_noperturb)) %>%
          mutate(!!sym(svr_opti_alloc_col) := optimal_c_1st)

        # Start aggregate file
        df_alloc_combined <- df_alloc_cur_optimal_nopertrub %>%
          select(-X, -optimal_c_1st) %>% ungroup()
        ls_st_opti_coll <- c(svr_opti_alloc_col)

        # Load in and stack perturbed files ---------------------------------
        # Results counter
        it_results_ctr <- 0
        # Loop over results
        for (it_perturb_ctr in ar_it_perturb_n) {
          it_rand_adj_A_rng_seed <- it_rand_adj_A_rng_seed_start + it_perturb_ctr

          # Counter updating
          it_results_ctr <- it_results_ctr + 1

          # Load support file
          ls_output <- fs_opti_support_202103(st_which_solu,
                                              bl_per_capita=bl_per_capita,
                                              fl_rho=fl_rho,
                                              snm_main=snm_main,
                                              it_rand_draw_seed=it_rand_adj_A_rng_seed,
                                              ls_st_gen_agg_folders=ls_st_gen_agg_folders)

          # Load folder names etc
          srt_results_root <- ls_output$srt_results_root
          bl_save_img <- ls_output$bl_save_img
          srt_folder <- ls_output$srt_folder
          bl_save_img <- TRUE
          srt_paper_appendix_textgraph <- ls_output$srt_paper_appendix_textgraph

          # folder that is not rho specific for storing aggregate folder
          srt_imgcsv_perturb_root <- ls_output$srt_imgcsv_perturb_root

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
          df_alloc_cur <- as_tibble(read.csv(spn_csv_full_path)) %>%
            mutate(!!sym(svr_opti_alloc_col) := optimal_c_1st)

          # Merging FIles Different Allocations are Additional Columns
          # if (it_results_ctr == 1) {
          #   df_alloc_combined <- df_alloc_cur %>%
          #     select(-X, -optimal_c_1st) %>% ungroup()
          #   ls_st_opti_coll <- c(svr_opti_alloc_col)
          # } else {
          df_alloc_combined <- df_alloc_combined %>%
            left_join(df_alloc_cur %>%
                        select(id_i, svr_opti_alloc_col) %>%
                        ungroup(), by='id_i')
          ls_st_opti_coll <- c(ls_st_opti_coll, svr_opti_alloc_col)
          # }
        }

        # Column Selection
        df_alloc_combined <- df_alloc_combined %>%
          select(id_i, rho_val, marital, kids, age_group, ymin_group,
                 mass, hhsize, actual,
                 contains('perturb'))

        # Transform Wide to Long
        df_alloc_combined_long <- df_alloc_combined %>%
          select(-actual) %>%
          pivot_longer(cols = starts_with('checks_optimal_capi_rhn1_perturb_rng'),
                       names_to = c('shockdraw'),
                       names_pattern = paste0("checks_optimal_capi_rhn1_perturb_rng(.*)"),
                       values_to = "checks") %>%
          rowwise() %>%
          mutate(y_group_min = substring(strsplit(ymin_group, ",")[[1]][1], 2),
                 y_group_max = gsub(strsplit(ymin_group, ",")[[1]][2],  pattern = "]", replacement = "")) %>%
          ungroup() %>%
          mutate(ymin_group = as.factor(ymin_group),
                 y_group_min = as.numeric(y_group_min))

        # Keep disturbed and non-disturbed results differently
        df_alloc_combined_long_disturbed <- df_alloc_combined_long %>% filter(shockdraw != it_disturb)
        df_alloc_combined_long_undisturbed <- df_alloc_combined_long %>% filter(shockdraw == it_disturb)

        # Distributional review
        # Group Summarize:
        str.stats.group <- 'allperc'
        ar.perc <- c(0.01, 0.05, 0.10, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99)
        REconTools::ff_summ_bygroup(
          df_alloc_combined_long_disturbed, c('marital'), 'checks', str.stats.group, ar.perc)$df_table_grp_stats
        REconTools::ff_summ_bygroup(
          df_alloc_combined_long_disturbed, c('kids'), 'checks', str.stats.group, ar.perc)$df_table_grp_stats
        REconTools::ff_summ_bygroup(
          df_alloc_combined_long_disturbed, c('age_group'), 'checks', str.stats.group, ar.perc)$df_table_grp_stats
        REconTools::ff_summ_bygroup(
          df_alloc_combined_long_disturbed, c('ymin_group'), 'checks', str.stats.group, ar.perc)$df_table_grp_stats

        # Joint group key stats, simply average due to random simple sampling of draws
        df_table_grp_stats_across_draws <- REconTools::ff_summ_bygroup(
          df_alloc_combined_long_disturbed, c('id_i', 'marital', 'kids', 'ymin_group'),
          'checks', str.stats.group, ar.perc)$df_table_grp_stats

        # Save all the allocation results jointly
        # the raw file
        snm_save_csv <- paste0('alloc_perturb_', st_bound_files, '.csv')
        write.csv(df_alloc_combined,
                  paste0(srt_imgcsv_perturb_root, snm_save_csv),
                  row.names = TRUE)
        # the distributional file
        snm_save_csv <- paste0('alloc_perturb_', st_bound_files, '_perc_by_group.csv')
        write.csv(df_table_grp_stats_across_draws,
                  paste0(srt_imgcsv_perturb_root, snm_save_csv),
                  row.names = TRUE)

        # graphing
        # G4 Figure in Table, Figure Number 4 Main Paper, Double both Kids and Adult Max
        # st_img_name_full <- paste0('alloc_perturb_', st_bound_files, '_perc_by_group.png')

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
        df_alloc_use <- df_alloc_combined_long_disturbed %>% mutate(kids = as.numeric(kids)) %>%
          filter(kids <= 2) %>%
          mutate(kids = as.factor(kids),
                 marital = as.factor(marital)) %>%
          mutate(kids = fct_recode(kids, !!!kids_levels),
                 marital = fct_recode(marital, !!!marry_levels))

        # Get value for minimum and maximum income levels for each bin
        df_alloc_use <- df_alloc_use %>%
          # rowwise() %>%
          # mutate(y_group_min = substring(strsplit(ymin_group, ",")[[1]][1], 2),
          #        y_group_max = gsub(strsplit(ymin_group, ",")[[1]][2],  pattern = "]", replacement = "")) %>%
          # ungroup() %>%
          mutate(kids = as.factor(kids),
                 marital = as.factor(marital))

        # factor to numeric for plotting
        df_alloc_use <- df_alloc_use %>%
          mutate(checks = as.numeric(checks))

        # Group by and generate aggregate statistics
        df_alloc_use_summ <- df_alloc_use %>%
          group_by(id_i, marital, kids, y_group_min) %>%
          summarize(stats_P05 = quantile(checks, probs = .05),
                    stats_P10 = quantile(checks, probs = .10),
                    stats_P25 = quantile(checks, probs = .25),
                    stats_P50 = quantile(checks, probs = .50),
                    stats_P75 = quantile(checks, probs = .75),
                    stats_P90 = quantile(checks, probs = .90),
                    stats_P95 = quantile(checks, probs = .95),
                    stats_Mean = mean(checks))

        # Group joining
        df_alloc_use_summ_all_wide <- df_alloc_use_summ %>%
          left_join(df_alloc_combined_long_undisturbed %>%
                      select(id_i, checks) %>%
                      rename(stats_nodisturb = checks), by='id_i')

        # Stats wide to long
        df_alloc_use_summ_all <- df_alloc_use_summ_all_wide %>%
          pivot_longer(cols = starts_with('stats_'),
                       names_to = c('stats'),
                       names_pattern = paste0("stats_(.*)"),
                       values_to = "checks") %>%
          mutate(checks = checks*100) %>%
          mutate(stats = case_when(stats == 'nodisturb' ~ 'Benchmark model',
                                   TRUE ~ stats))

        for (it_graph_types in ls_it_graph_types) {

          fl_norm_size <- 0.5
          fl_thck_size <- 0.75
          if (it_graph_types == 1) {
            df_alloc_use_graph <- df_alloc_use_summ_all %>%
              filter(stats=='P05' | stats=='P25' | stats=='P50' |
                       stats=='P75' | stats=='P95' | stats=='Mean')

            ls_scale_colour_manual <- c("red", "black", "blue", "blue", "blue", "black")
            scale_linetype_manual <- c("twodash", "dotted", "dashed", "solid", "dashed", "dotted")
            ls_scale_size_manual <- c(fl_norm_size, fl_thck_size, fl_norm_size, fl_norm_size, fl_norm_size, fl_thck_size)

          } else if (it_graph_types == 2) {

            df_alloc_use_graph <- df_alloc_use_summ_all %>%
              filter(stats=='P05' | stats=='P95' | stats=='Benchmark model')

            ls_scale_colour_manual <- c("#CD9600", "black", "black")
            scale_linetype_manual <- c("solid", "dotted", "dotted")
            ls_scale_size_manual <- c(fl_norm_size, fl_thck_size, fl_thck_size)

          } else if (it_graph_types == 3) {

            df_alloc_use_graph <- df_alloc_use_summ_all %>%
              filter(stats=='P05' | stats=='P95'  | stats=='Mean' | stats=='Benchmark model')
            ls_scale_colour_manual <- c("#CD9600", "red", "black", "black")
            scale_linetype_manual <- c("solid", "dashed", "dotted", "dotted")
            ls_scale_size_manual <- c(fl_norm_size, fl_norm_size, fl_thck_size, fl_thck_size)

          } else if (it_graph_types == 4) {

            df_alloc_use_graph <- df_alloc_use_summ_all %>%
              filter(stats=='P05' | stats=='P25' | stats=='P75' | stats=='P95' |
                       stats=='Mean' | stats=='Benchmark model')

            ls_scale_colour_manual <- c("#CD9600", "red", "black", "blue", "blue", "black")
            scale_linetype_manual <- c("solid", "twodash", "dotted", "dashed", "dashed", "dotted")
            ls_scale_size_manual <- c(fl_norm_size, fl_norm_size, fl_thck_size, fl_norm_size, fl_norm_size, fl_thck_size)

          } else if (it_graph_types == 5) {

            df_alloc_use_graph <- df_alloc_use_summ_all %>%
              filter(stats=='P05' | stats=='P25' | stats=='P75' | stats=='P95' |
                       stats=='Benchmark model')

            ls_scale_colour_manual <- c("#CD9600", "black", "blue", "blue", "black")
            scale_linetype_manual <- c("solid", "dotted", "dashed", "dashed", "dotted")
            ls_scale_size_manual <- c(fl_norm_size, fl_thck_size, fl_norm_size, fl_norm_size, fl_thck_size)
          }

          plt_cur <- df_alloc_use_graph %>%
            ggplot(aes(x=y_group_min, y=checks,
                       colour=stats,
                       linetype=stats,
                       size=stats)) +
            facet_wrap(~ marital + kids,
                       ncol=3,
                       scales = "free_x",
                       labeller = label_wrap_gen(multi_line=FALSE))
          plt_cur <- plt_cur + geom_line()

          if (it_graph_types != 1) {
            plt_cur <- plt_cur + geom_point(data=df_alloc_use_graph %>%
                                              filter(stats=='Benchmark model'),
                                            aes(x=y_group_min, y=checks),
                                            size=3, shape=15,
                                            show.legend = FALSE)
          }

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
          # limits = c(0, 240000/58056)


          # color #F8766D is default red
          # color #F8766D is default bluish color
          plt_cur <- plt_cur +
            theme(
              text = element_text(size = 12),
              legend.title = element_blank(),
              legend.position = c(0.14, 0.9),
              legend.background = element_rect(fill = "white", colour = "black", linetype='solid'),
              legend.key.width = unit(3, "line")) +
            scale_colour_manual(values=ls_scale_colour_manual) +
            scale_linetype_manual(values=scale_linetype_manual) +
            scale_size_manual(values=ls_scale_size_manual, guide=FALSE)

          # X and Y Titling
          stg_x_title_hhinc <- 'Household income (thousands of 2012 USD)'
          stg_y_title_checks <- 'Stimulus check amount (2012 USD)'
          plt_cur <- plt_cur + labs(x = stg_x_title_hhinc, y = stg_y_title_checks,
                                    colour = "abc", linetype = "abc")

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
            st_img_name_full <- paste0('alloc_perturb_', st_bound_files, '_perc_by_group_v', it_graph_types, '.png')
            ggsave(plt_cur, file=paste0(srt_imgcsv_perturb_root,
                                        st_img_name_full),
                   width = 270,
                   height = it_img_height, units='mm',
                   dpi = 300, pointsize=7)
          }

        }
      }
    }
  }
}
