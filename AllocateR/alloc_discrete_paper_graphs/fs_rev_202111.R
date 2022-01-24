#  Graphs where allocation lines from multiple planners are jointly compared.
# "Paper\Appendix_graphs_and_tables\Robustness_Tax_model_Actual_vs_feasible_c_allocation_round1_18-64_year-olds_by_income.png"
try(dev.off(dev.list()["RStudioGD"]),silent=TRUE)
try(dev.off(),silent=TRUE)

# Libraries
library(tidyverse)
library(REconTools)
library(scales)

# Parameters and Options --------

# Possible rho values to select from
ar_seq_equi <- seq(-2, 2, length.out=30)
# initial points
ar_rho_init <- 1 - (10^(c(seq(-2, 2, length.out=30))))
# ar_rho_init <- ar_rho_init[seq(1, 30, by=2)]
ar_rho_init <- c(1, round(ar_rho_init, 4), -100)
# extra points
ar_rho_extra <- 1 - (10^(c(seq(ar_seq_equi[23], ar_seq_equi[27], length.out=15))))
# ar_rho_extra <- round(ar_rho_extra, 4)
ar_rho_five <- c(1, 0.9, -1, -19, -99)
# combine together
ar_rho <- c(ar_rho_five)
ar_rho <- sort(ar_rho)

# log10((-1)*(ar_rho - 1))
# 1/(1-ar_rho)

# Results from which planners should
ls_it_rho_combo_type <- c(1)

# Per capita or per household results
ls_bl_per_capita <- c(TRUE, FALSE)
ls_bl_per_capita <- c(TRUE)


# # Types of allocation files to consider
# ls_st_file_suffix <- c('snwx_bchklock_moredense_a65zh266zs5_b1_xi0_manna_168')
# # Allocation bounds types
# ls_st_bound_files <- c('14ca14ck', '17ca17ck', '20ca20ck')
# Types of allocation files to consider
ls_st_file_suffix <- c('snwx_bushchck_moredense_a65zh266zs5_b1_xi0_manna_96')
# Allocation bounds types
ls_st_bound_files <- c('6ca3ck', '9ca6ck', '12ca9ck')
ls_st_bound_files <- c('6ca3ck', '9ca6ck', '12ca9ck',
                       '6ckid', '9cadt', '9ckid', '12cadt')

# show results for all age as one group or optimally for each of the four age groups
ls_st_all_or_g4 <- c('feasible', 'optimalg4')

# number of y-ticks
it_y_nbreaks <- 5
# Which graph to draw
bl_rev_graph <- FALSE
bl_atkinson_graph <- FALSE
bl_gini_outcome_graph <- FALSE
bl_gini_checks_graph <- FALSE
bl_sd_log_c_graph <- FALSE
bl_outcome_mean_graph <- FALSE
bl_outcome_min_graph <- FALSE

# Loop over file types1
for (bl_per_capita in ls_bl_per_capita) {
  for (st_which_solu in ls_st_file_suffix) {
    for (st_all_or_g4 in ls_st_all_or_g4) {
      for (bl_approximate in c(TRUE)) {
        for (bl_zoomout0t1 in c(FALSE)) {
          for (bl_nonone in c(FALSE, TRUE)) {
            for (it_income_bounds in c(1,2,3)) {

              if (it_income_bounds == 1) {
                fl_income_bound <- 1000000
                st_income_bound <- '_yleqAll'
              } else if (it_income_bounds == 2) {
                fl_income_bound <- 200000
                st_income_bound <- '_yleq200k'
              } else if (it_income_bounds == 3) {
                fl_income_bound <- 160000
                st_income_bound <- '_yleq160k'
              }

              # Generate files by aggregating over rho types
              # Results counter
              it_ls_results_ctr <- 0
              ls_of_ls_rev_results <- vector(mode = "list", length = length(ls_st_bound_files))

              for (st_bound_files in ls_st_bound_files) {
                it_ls_results_ctr <- it_ls_results_ctr + 1

                # Generate files by aggregating over rho types
                # Results counter
                it_results_ctr <- 0
                ls_rev_results <- vector(mode = "list", length = length(ar_rho))
                # Loop over results
                for (fl_rho in ar_rho) {

                  # Counter updating
                  it_results_ctr <- it_results_ctr + 1

                  # Load support file
                  ls_output <- fs_opti_support_202111(st_which_solu, bl_per_capita=bl_per_capita, fl_rho=fl_rho)
                  # Load folder names etc
                  srt_results_root <- ls_output$srt_results_root
                  bl_save_img <- ls_output$bl_save_img
                  srt_imgcsv_rev_root <- ls_output$srt_imgcsv_rev_root
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
                  srt_csv_file <- paste0('rev_', st_all_or_g4, '_', st_bound_files, '.csv')
                  spn_csv_full_path <- file.path(srt_results_root, srt_res_source_folder,
                                                 srt_csv_allocate_subfolder, srt_csv_file)

                  # Clean up each file
                  df_alloc_cur <- as_tibble(read.csv(spn_csv_full_path))
                  ls_rev_results[[it_results_ctr]] <- df_alloc_cur
                }

                # bind files together
                mt_rev_across_rhos <- do.call(rbind, ls_rev_results)
                ls_of_ls_rev_results[[it_ls_results_ctr]] <- mt_rev_across_rhos %>%
                  mutate(bounds = srt_csv_allocate_subfolder)

                # CSV Save
                snm_save_csv <- paste0('rev_', st_all_or_g4, '_', st_bound_files, '_rhomultiple.csv')
                write.csv(mt_rev_across_rhos,
                          paste0(srt_imgcsv_rev_root, snm_save_csv),
                          row.names = TRUE)

              }

              # cAll results from different bounds to the same file
              mt_rev_across_rhos_bounds <- do.call(rbind, ls_of_ls_rev_results)

              # Filter by income thresholds
              mt_rev_across_rhos_bounds <- mt_rev_across_rhos_bounds %>%
                filter(income_bound == fl_income_bound)

              # FIGURE Shared -----------------------------------
              # graph, copied from https://fanwangecon.github.io/PrjOptiAlloc/articles/ffv_opt_sobin_rkone_allrw_training.html
              # x-labels
              x.labels <- c('\u03bb \U2248 1.00', '\u03bb = 0.90', '\u03bb = 0', '\u03bb = -10', '\u03bb = -100')
              x.breaks <- c(0.01, 0.10, 1, 10, 100)

              # FIGURE 1 REV -----------------------------------
              if (bl_rev_graph && it_income_bounds == 1) {
                # collect only relevant columns
                mt_rev_across_rhos_bounds_selected <- mt_rev_across_rhos_bounds %>%
                  filter(objective == "c2020") %>%
                  mutate(REV = REV*100) %>%
                  mutate(one_minus_rho = 1 - rho_val) %>%
                  select(one_minus_rho, REV, bounds) %>%
                  rename(Outcome = bounds, rev = REV)

                # Relabel Variable
                Outcome_levels <- c("$1,400" = "b1_a64_14ca14ck_18t64",
                                    "$2,000" = "b1_a64_20ca20ck_18t64")
                mt_rev_across_rhos_bounds_selected_fig <- mt_rev_across_rhos_bounds_selected %>%
                  filter(Outcome == 'b1_a64_14ca14ck_18t64'
                         | Outcome == 'b1_a64_20ca20ck_18t64') %>%
                  mutate(Outcome = as_factor(Outcome)) %>%
                  mutate(Outcome = fct_recode(Outcome, !!!Outcome_levels))

                if (bl_approximate) {
                  # Graph Results--Draw
                  line.plot <- mt_rev_across_rhos_bounds_selected_fig %>%
                    ggplot(aes(x=one_minus_rho, y=rev,
                               group=Outcome,
                               colour=Outcome,
                               linetype=Outcome,
                               shape=Outcome)) +
                    # geom_point() +
                    # geom_line() +
                    geom_smooth(span = 0.35, se=FALSE, size=2)
                } else {
                  # Graph Results--Draw
                  line.plot <- mt_rev_across_rhos_bounds_selected_fig %>%
                    ggplot(aes(x=one_minus_rho, y=rev,
                               group=Outcome,
                               colour=Outcome,
                               linetype=Outcome,
                               shape=Outcome)) +
                    geom_point(size=6) +
                    geom_line(size=2)
                  # geom_smooth(span = 0.35, se=FALSE, size=2)
                  # geom_vline(xintercept=c(1), linetype="dotted") +
                }

                line.plot <- line.plot +
                  labs(x = 'Planner inequality aversion, \u03bb',
                       y = paste0('Resource equivalent variation (percent)')) +
                  scale_x_continuous(trans='log10', labels = x.labels, breaks = x.breaks) +
                  scale_y_continuous(n.breaks = it_y_nbreaks) +
                  theme_bw(base_size=8) +
                  ylim(0, 100)

                # legend area
                line.plot <- line.plot +
                  theme(
                    text = element_text(size = 16),
                    legend.position = c(0.22, 0.75),
                    legend.background = element_rect(fill = "white", colour = "black",
                                                     linetype='solid'))

                line.plot <- line.plot +
                  labs(colour = "Check limit per adult and child",
                       shape = "Check limit per adult and child",
                       linetype = "Check limit per adult and child")

                # # Labeling
                # line.plot$labels$linetype <- "REV\nOptimal\nvs.\nAlternatives"
                # line.plot$labels$colour <- "REV\nOptimal\nvs.\nAlternatives"
                # line.plot$labels$shape <- "REV\nOptimal\nvs.\nAlternatives"

                # Print
                print(line.plot)

                # add png
                if (bl_approximate) {
                  snm_save_png <- paste0('rev_feasible_rhomultiple_smooth', st_income_bound, '.png')
                } else {
                  snm_save_png <- paste0('rev_feasible_rhomultiple_pointline', st_income_bound, '.png')
                }

                # Save
                ggsave(line.plot,
                       file=file.path(srt_imgcsv_rev_root, snm_save_png),
                       width = 270,
                       height = 216, units='mm',
                       dpi = 300)

              }


              # FIGURE 2 Atkinson -----------------------------------
              if (bl_atkinson_graph && it_income_bounds == 1) {

                # collect only relevant columns
                mt_atkinson_across_rhos_bounds_selected <- mt_rev_across_rhos_bounds %>%
                  filter(objective == "c2020") %>%
                  mutate(one_minus_rho = 1 - rho_val) %>%
                  select(one_minus_rho, atkinson_drv_EH_star, bounds) %>%
                  rename(Outcome = bounds)
                # Actual outcome
                mt_atkinson_across_rhos_bounds_actual <- mt_rev_across_rhos_bounds %>%
                  filter(objective == "c2020") %>%
                  mutate(one_minus_rho = 1 - rho_val) %>%
                  select(one_minus_rho, atkinson_drv_EH_dact, bounds) %>%
                  rename(Outcome = bounds) %>%
                  filter(Outcome == 'b1_a64_14ca14ck_18t64') %>%
                  select(one_minus_rho, atkinson_drv_EH_dact) %>%
                  rename(atkinson_drv_EH_star = atkinson_drv_EH_dact) %>%
                  mutate(Outcome = 'Actual')
                # None outcome
                mt_atkinson_across_rhos_bounds_none <- mt_rev_across_rhos_bounds %>%
                  filter(objective == "c2020") %>%
                  mutate(one_minus_rho = 1 - rho_val) %>%
                  select(one_minus_rho, atkinson_drv_EH_d0, bounds) %>%
                  rename(Outcome = bounds) %>%
                  filter(Outcome == 'b1_a64_14ca14ck_18t64') %>%
                  select(one_minus_rho, atkinson_drv_EH_d0) %>%
                  rename(atkinson_drv_EH_star = atkinson_drv_EH_d0) %>%
                  mutate(Outcome = 'None')
                # Stack Actual and Optimal
                mt_atkinson_across_rhos_bounds_selected <-
                  rbind(mt_atkinson_across_rhos_bounds_none,
                        mt_atkinson_across_rhos_bounds_actual,
                        mt_atkinson_across_rhos_bounds_selected)

                # Relabel Variable
                Outcome_levels <- c("No checks" = "None",
                                    "Actual allocations" = "Actual",
                                    "$1,400 limit optimal" = "b1_a64_14ca14ck_18t64",
                                    "$2,000 limit optimal" = "b1_a64_20ca20ck_18t64")
                mt_atkinson_across_rhos_bounds_selected_fig <- mt_atkinson_across_rhos_bounds_selected %>%
                  filter(Outcome == 'b1_a64_14ca14ck_18t64'
                         | Outcome == 'b1_a64_20ca20ck_18t64'
                         | Outcome == 'Actual'
                         | Outcome == 'None') %>%
                  mutate(Outcome = as_factor(Outcome)) %>%
                  mutate(Outcome = fct_recode(Outcome, !!!Outcome_levels))

                if (bl_approximate) {
                  # Graph Results--Draw
                  line.plot <- mt_atkinson_across_rhos_bounds_selected_fig %>%
                    ggplot(aes(x=one_minus_rho, y=atkinson_drv_EH_star,
                               group=Outcome,
                               colour=Outcome,
                               linetype=Outcome,
                               shape=Outcome)) +
                    # geom_point() +
                    # geom_line() +
                    geom_smooth(span = 0.35, se=FALSE, size=2)
                } else {
                  # Graph Results--Draw
                  line.plot <- mt_atkinson_across_rhos_bounds_selected_fig %>%
                    ggplot(aes(x=one_minus_rho, y=atkinson_drv_EH_star,
                               group=Outcome,
                               colour=Outcome,
                               linetype=Outcome,
                               shape=Outcome)) +
                    geom_point(size=6) +
                    geom_line(size=2)
                  # geom_smooth(span = 0.35, se=FALSE, size=2)
                  # geom_vline(xintercept=c(1), linetype="dotted") +
                }

                line.plot <- line.plot +
                  labs(x = 'Planner inequality aversion, \u03bb',
                       y = paste0('Atkinson Index')) +
                  scale_x_continuous(trans='log10', labels = x.labels, breaks = x.breaks) +
                  scale_y_continuous(n.breaks = it_y_nbreaks) +
                  theme_bw(base_size=8) +
                  ylim(0, 1)

                # legend area
                line.plot <- line.plot +
                  theme(
                    text = element_text(size = 16),
                    legend.position = c(0.22, 0.75),
                    legend.background = element_rect(fill = "white", colour = "black",
                                                     linetype='solid'))

                line.plot <- line.plot +
                  labs(colour = "Check limit per adult and child",
                       shape = "Check limit per adult and child",
                       linetype = "Check limit per adult and child")

                # # Labeling
                # line.plot$labels$linetype <- "REV\nOptimal\nvs.\nAlternatives"
                # line.plot$labels$colour <- "REV\nOptimal\nvs.\nAlternatives"
                # line.plot$labels$shape <- "REV\nOptimal\nvs.\nAlternatives"

                # Print
                print(line.plot)

                # add png
                if (bl_approximate) {
                  snm_save_png <- paste0('atkinson_feasible_rhomultiple_smooth', st_income_bound, '.png')
                } else {
                  snm_save_png <- paste0('atkinson_feasible_rhomultiple_pointline', st_income_bound, '.png')
                }

                # Save
                ggsave(line.plot,
                       file=file.path(srt_imgcsv_rev_root, snm_save_png),
                       width = 270,
                       height = 216, units='mm',
                       dpi = 300)

              }

              # FIGURE 3 Outcome Gini -----------------------------------
              if (bl_gini_outcome_graph) {

                # collect only relevant columns
                mt_gini_across_rhos_bounds_selected <- mt_rev_across_rhos_bounds %>%
                  filter(objective == "c2020") %>%
                  mutate(one_minus_rho = 1 - rho_val) %>%
                  select(one_minus_rho, gini_drv_EH_star, bounds) %>%
                  rename(Outcome = bounds)
                # Zero outcome
                mt_gini_across_rhos_bounds_none <- mt_rev_across_rhos_bounds %>%
                  filter(objective == "c2020") %>%
                  mutate(one_minus_rho = 1 - rho_val) %>%
                  select(one_minus_rho, gini_drv_EH_d0, bounds) %>%
                  rename(Outcome = bounds) %>%
                  filter(Outcome == 'b1_a64_14ca14ck_18t64') %>%
                  select(one_minus_rho, gini_drv_EH_d0) %>%
                  rename(gini_drv_EH_star = gini_drv_EH_d0) %>%
                  mutate(Outcome = 'None')
                # Actual outcome
                mt_gini_across_rhos_bounds_actual <- mt_rev_across_rhos_bounds %>%
                  filter(objective == "c2020") %>%
                  mutate(one_minus_rho = 1 - rho_val) %>%
                  select(one_minus_rho, gini_drv_EH_dact, bounds) %>%
                  rename(Outcome = bounds) %>%
                  filter(Outcome == 'b1_a64_14ca14ck_18t64') %>%
                  select(one_minus_rho, gini_drv_EH_dact) %>%
                  rename(gini_drv_EH_star = gini_drv_EH_dact) %>%
                  mutate(Outcome = 'Actual')

                # Stack Actual and Optimal
                mt_gini_across_rhos_bounds_selected <-
                  rbind(mt_gini_across_rhos_bounds_none,
                        mt_gini_across_rhos_bounds_actual,
                        mt_gini_across_rhos_bounds_selected)

                # Relabel Variable
                if (bl_nonone) {
                  Outcome_levels <- c("Actual allocations" = "Actual",
                                      "$1,400 limit optimal" = "b1_a64_14ca14ck_18t64",
                                      "$1,700 limit optimal" = "b1_a64_17ca17ck_18t64",
                                      "$2,000 limit optimal" = "b1_a64_20ca20ck_18t64")
                  mt_gini_across_rhos_bounds_selected_fig <- mt_gini_across_rhos_bounds_selected %>%
                    filter(Outcome == 'b1_a64_14ca14ck_18t64'
                           | Outcome == 'b1_a64_20ca20ck_18t64'
                           | Outcome == 'Actual') %>%
                    mutate(Outcome = as_factor(Outcome)) %>%
                    mutate(Outcome = fct_recode(Outcome, !!!Outcome_levels))
                } else {
                  Outcome_levels <- c("No checks" = "None",
                                      "Actual allocations" = "Actual",
                                      "$1,400 limit optimal" = "b1_a64_14ca14ck_18t64",
                                      "$1,700 limit optimal" = "b1_a64_17ca17ck_18t64",
                                      "$2,000 limit optimal" = "b1_a64_20ca20ck_18t64")
                  mt_gini_across_rhos_bounds_selected_fig <- mt_gini_across_rhos_bounds_selected %>%
                    filter(Outcome == 'b1_a64_14ca14ck_18t64'
                           | Outcome == 'b1_a64_17ca17ck_18t64'
                           | Outcome == 'b1_a64_20ca20ck_18t64'
                           | Outcome == 'Actual'
                           | Outcome == 'None') %>%
                    mutate(Outcome = as_factor(Outcome)) %>%
                    mutate(Outcome = fct_recode(Outcome, !!!Outcome_levels))
                }

                if (bl_approximate) {
                  # Graph Results--Draw
                  line.plot <- mt_gini_across_rhos_bounds_selected_fig %>%
                    ggplot(aes(x=one_minus_rho, y=gini_drv_EH_star,
                               group=Outcome,
                               colour=Outcome,
                               linetype=Outcome,
                               shape=Outcome)) +
                    # geom_point() +
                    # geom_line() +
                    geom_smooth(span = 0.35, se=FALSE, size=2)
                } else {
                  # Graph Results--Draw
                  line.plot <- mt_gini_across_rhos_bounds_selected_fig %>%
                    ggplot(aes(x=one_minus_rho, y=gini_drv_EH_star,
                               group=Outcome,
                               colour=Outcome,
                               linetype=Outcome,
                               shape=Outcome)) +
                    geom_point(size=6) +
                    geom_line(size=2)
                  # geom_smooth(span = 0.35, se=FALSE, size=2)
                  # geom_vline(xintercept=c(1), linetype="dotted") +
                }

                fl_gini_min <- min(mt_gini_across_rhos_bounds_selected_fig$gini_drv_EH_star)
                fl_gini_max <- max(mt_gini_across_rhos_bounds_selected_fig$gini_drv_EH_star)

                line.plot <- line.plot +
                  labs(x = 'Planner inequality aversion, \u03bb',
                       y = paste0('GINI coefficient for consumption given stimulus checks')) +
                  scale_x_continuous(trans='log10', labels = x.labels, breaks = x.breaks) +
                  scale_y_continuous(n.breaks = it_y_nbreaks) +
                  theme_bw(base_size=8)

                if (bl_zoomout0t1) {
                  line.plot <- line.plot + ylim(0, 1)
                } else {
                  line.plot <- line.plot +
                    ylim(fl_gini_min - abs(fl_gini_max-fl_gini_min)/20, fl_gini_max + abs(fl_gini_max-fl_gini_min)/20)
                }

                # legend area
                line.plot <- line.plot +
                  theme(
                    text = element_text(size = 16),
                    legend.position = c(0.22, 0.75),
                    legend.background = element_rect(fill = "white", colour = "black",
                                                     linetype='solid'))

                line.plot <- line.plot +
                  labs(colour = "Check limit per adult and child",
                       shape = "Check limit per adult and child",
                       linetype = "Check limit per adult and child")

                # # Labeling
                # line.plot$labels$linetype <- "REV\nOptimal\nvs.\nAlternatives"
                # line.plot$labels$colour <- "REV\nOptimal\nvs.\nAlternatives"
                # line.plot$labels$shape <- "REV\nOptimal\nvs.\nAlternatives"

                # Print
                print(line.plot)

                # add png
                if (bl_zoomout0t1) {
                  st_zoom <- '_0t1'
                } else {
                  st_zoom <- '_zoomed'
                }
                if (bl_nonone) {
                  st_actual <- '_nonone'
                } else {
                  st_actual <- '_withnone'
                }
                if (bl_approximate) {
                  snm_save_png <- paste0('gini_EC_feasible_rhomultiple_smooth', st_zoom, st_actual, st_income_bound, '.png')
                } else {
                  snm_save_png <- paste0('gini_EC_feasible_rhomultiple_pointline', st_zoom, st_actual, st_income_bound, '.png')
                }

                # Save
                ggsave(line.plot,
                       file=file.path(srt_imgcsv_rev_root, snm_save_png),
                       width = 270,
                       height = 216, units='mm',
                       dpi = 300)

              }

              # FIGURE 4 Allocation Gini -----------------------------------
              if (bl_gini_checks_graph && it_income_bounds == 1) {

                # collect only relevant columns
                mt_gini_across_rhos_bounds_selected <- mt_rev_across_rhos_bounds %>%
                  filter(objective == "c2020") %>%
                  mutate(one_minus_rho = 1 - rho_val) %>%
                  select(one_minus_rho, gini_drv_D_star, bounds) %>%
                  rename(Outcome = bounds)
                # Actual outcome
                # see G:\repos\PrjOptiSNW\AllocateR\alloc_discrete_paper_graphs\fs_actual_checks_gini.R for 0.2886926
                mt_gini_across_rhos_bounds_actual <- mt_rev_across_rhos_bounds %>%
                  filter(objective == "c2020") %>%
                  mutate(one_minus_rho = 1 - rho_val) %>%
                  select(one_minus_rho, bounds) %>%
                  rename(Outcome = bounds) %>%
                  filter(Outcome == 'b1_a64_14ca14ck_18t64') %>%
                  select(one_minus_rho) %>%
                  mutate(gini_drv_D_star = 0.346999380993645, Outcome = 'Actual')
                # Stack Actual and Optimal
                mt_gini_across_rhos_bounds_selected <-
                  rbind(mt_gini_across_rhos_bounds_selected, mt_gini_across_rhos_bounds_actual)

                # Relabel Variable
                Outcome_levels <- c("$1,400 limit optimal" = "b1_a64_14ca14ck_18t64",
                                    "$2,000 limit optimal" = "b1_a64_20ca20ck_18t64",
                                    "Actual (all pop) allocations" = "Actual")
                mt_gini_across_rhos_bounds_selected_fig <- mt_gini_across_rhos_bounds_selected %>%
                  filter(Outcome == 'b1_a64_14ca14ck_18t64'
                         | Outcome == 'b1_a64_20ca20ck_18t64'
                         | Outcome == 'Actual') %>%
                  mutate(Outcome = as_factor(Outcome)) %>%
                  mutate(Outcome = fct_recode(Outcome, !!!Outcome_levels))

                if (bl_approximate) {
                  # Graph Results--Draw
                  line.plot <- mt_gini_across_rhos_bounds_selected_fig %>%
                    ggplot(aes(x=one_minus_rho, y=gini_drv_D_star,
                               group=Outcome,
                               colour=Outcome,
                               linetype=Outcome,
                               shape=Outcome)) +
                    # geom_point() +
                    # geom_line() +
                    geom_smooth(span = 0.35, se=FALSE, size=2)
                } else {
                  # Graph Results--Draw
                  line.plot <- mt_gini_across_rhos_bounds_selected_fig %>%
                    ggplot(aes(x=one_minus_rho, y=gini_drv_D_star,
                               group=Outcome,
                               colour=Outcome,
                               linetype=Outcome,
                               shape=Outcome)) +
                    geom_point(size=6) +
                    geom_line(size=2)
                  # geom_smooth(span = 0.35, se=FALSE, size=2)
                  # geom_vline(xintercept=c(1), linetype="dotted") +
                }

                fl_gini_min <- min(mt_gini_across_rhos_bounds_selected_fig$gini_drv_EH_star)
                fl_gini_max <- max(mt_gini_across_rhos_bounds_selected_fig$gini_drv_EH_star)

                line.plot <- line.plot +
                  labs(x = 'Planner inequality aversion, \u03bb',
                       y = paste0('Gini coefficient for the allocation of stimulus checks')) +
                  scale_x_continuous(trans='log10', labels = x.labels, breaks = x.breaks) +
                  scale_y_continuous(n.breaks = it_y_nbreaks) +
                  theme_bw(base_size=8)

                if (bl_zoomout0t1) {
                  line.plot <- line.plot + ylim(0, 1)
                } else {
                  line.plot <- line.plot +
                    ylim(fl_gini_min - abs(fl_gini_max-fl_gini_min)/20, fl_gini_max + abs(fl_gini_max-fl_gini_min)/20)
                }

                # legend area
                line.plot <- line.plot +
                  theme(
                    text = element_text(size = 16),
                    legend.position = c(0.22, 0.75),
                    legend.background = element_rect(fill = "white", colour = "black",
                                                     linetype='solid'))

                line.plot <- line.plot +
                  labs(colour = "Check limit per adult and child",
                       shape = "Check limit per adult and child",
                       linetype = "Check limit per adult and child")

                # # Labeling
                # line.plot$labels$linetype <- "REV\nOptimal\nvs.\nAlternatives"
                # line.plot$labels$colour <- "REV\nOptimal\nvs.\nAlternatives"
                # line.plot$labels$shape <- "REV\nOptimal\nvs.\nAlternatives"

                # Print
                print(line.plot)

                # add png
                if (bl_zoomout0t1) {
                  st_zoom <- '_0t1'
                } else {
                  st_zoom <- '_zoomed'
                }
                if (bl_approximate) {
                  snm_save_png <- paste0('gini_checks_feasible_rhomultiple_smooth', st_zoom, st_income_bound, '.png')
                } else {
                  snm_save_png <- paste0('gini_checks_feasible_rhomultiple_pointline', st_zoom, st_income_bound,  '.png')
                }

                # Save
                ggsave(line.plot,
                       file=file.path(srt_imgcsv_rev_root, snm_save_png),
                       width = 270,
                       height = 216, units='mm',
                       dpi = 300)

              }

              # FIGURE 5 s.d. of log of c -----------------------------------
              if (bl_sd_log_c_graph) {
                # collect only relevant columns
                mt_sd_log_c_across_rhos_bounds_selected <- mt_rev_across_rhos_bounds %>%
                  filter(objective == "c2020") %>%
                  mutate(one_minus_rho = 1 - rho_val) %>%
                  select(one_minus_rho, sd_log_EH_drv_star, bounds) %>%
                  rename(Outcome = bounds)
                # Zero outcome
                mt_sd_log_c_across_rhos_bounds_none <- mt_rev_across_rhos_bounds %>%
                  filter(objective == "c2020") %>%
                  mutate(one_minus_rho = 1 - rho_val) %>%
                  select(one_minus_rho, sd_log_EH_drv_d0, bounds) %>%
                  rename(Outcome = bounds) %>%
                  filter(Outcome == 'b1_a64_14ca14ck_18t64') %>%
                  select(one_minus_rho, sd_log_EH_drv_d0) %>%
                  rename(sd_log_EH_drv_star = sd_log_EH_drv_d0) %>%
                  mutate(Outcome = 'None')
                # Actual outcome
                mt_sd_log_c_across_rhos_bounds_actual <- mt_rev_across_rhos_bounds %>%
                  filter(objective == "c2020") %>%
                  mutate(one_minus_rho = 1 - rho_val) %>%
                  select(one_minus_rho, sd_log_EH_drv_dact, bounds) %>%
                  rename(Outcome = bounds) %>%
                  filter(Outcome == 'b1_a64_14ca14ck_18t64') %>%
                  select(one_minus_rho, sd_log_EH_drv_dact) %>%
                  rename(sd_log_EH_drv_star = sd_log_EH_drv_dact) %>%
                  mutate(Outcome = 'Actual')

                # Stack Actual and Optimal
                mt_sd_log_c_across_rhos_bounds_selected <-
                  rbind(mt_sd_log_c_across_rhos_bounds_none,
                        mt_sd_log_c_across_rhos_bounds_actual,
                        mt_sd_log_c_across_rhos_bounds_selected)

                # Relabel Variable
                if (bl_nonone) {
                  Outcome_levels <- c("Actual allocations" = "Actual",
                                      "$1,400 limit optimal" = "b1_a64_14ca14ck_18t64",
                                      "$1,700 limit optimal" = "b1_a64_17ca17ck_18t64",
                                      "$2,000 limit optimal" = "b1_a64_20ca20ck_18t64")
                  mt_sd_log_c_across_rhos_bounds_selected_fig <- mt_sd_log_c_across_rhos_bounds_selected %>%
                    filter(Outcome == 'b1_a64_14ca14ck_18t64'
                           | Outcome == 'b1_a64_17ca17ck_18t64'
                           | Outcome == 'b1_a64_20ca20ck_18t64'
                           | Outcome == 'Actual') %>%
                    mutate(Outcome = as_factor(Outcome)) %>%
                    mutate(Outcome = fct_recode(Outcome, !!!Outcome_levels))
                } else {
                  Outcome_levels <- c("No checks" = "None",
                                      "Actual allocations" = "Actual",
                                      "$1,400 limit optimal" = "b1_a64_14ca14ck_18t64",
                                      "$1,700 limit optimal" = "b1_a64_17ca17ck_18t64",
                                      "$2,000 limit optimal" = "b1_a64_20ca20ck_18t64")
                  mt_sd_log_c_across_rhos_bounds_selected_fig <- mt_sd_log_c_across_rhos_bounds_selected %>%
                    filter(Outcome == 'b1_a64_14ca14ck_18t64'
                           | Outcome == 'b1_a64_17ca17ck_18t64'
                           | Outcome == 'b1_a64_20ca20ck_18t64'
                           | Outcome == 'Actual'
                           | Outcome == 'None') %>%
                    mutate(Outcome = as_factor(Outcome)) %>%
                    mutate(Outcome = fct_recode(Outcome, !!!Outcome_levels))
                }

                if (bl_approximate) {
                  # Graph Results--Draw
                  line.plot <- mt_sd_log_c_across_rhos_bounds_selected_fig %>%
                    ggplot(aes(x=one_minus_rho, y=sd_log_EH_drv_star,
                               group=Outcome,
                               colour=Outcome,
                               linetype=Outcome,
                               shape=Outcome)) +
                    # geom_point() +
                    # geom_line() +
                    geom_smooth(span = 0.35, se=FALSE, size=2)
                } else {
                  # Graph Results--Draw
                  line.plot <- mt_sd_log_c_across_rhos_bounds_selected_fig %>%
                    ggplot(aes(x=one_minus_rho, y=sd_log_EH_drv_star,
                               group=Outcome,
                               colour=Outcome,
                               linetype=Outcome,
                               shape=Outcome)) +
                    geom_point(size=6) +
                    geom_line(size=2)
                  # geom_smooth(span = 0.35, se=FALSE, size=2)
                  # geom_vline(xintercept=c(1), linetype="dotted") +
                }

                fl_gini_min <- min(mt_sd_log_c_across_rhos_bounds_selected_fig$gini_drv_EH_star)
                fl_gini_max <- max(mt_sd_log_c_across_rhos_bounds_selected_fig$gini_drv_EH_star)

                line.plot <- line.plot +
                  labs(x = 'Planner inequality aversion, \u03bb',
                       y = paste0('Standard deviation of log of consumption given stimulus checks')) +
                  scale_x_continuous(trans='log10', labels = x.labels, breaks = x.breaks) +
                  scale_y_continuous(n.breaks = it_y_nbreaks) +
                  theme_bw(base_size=8)

                if (bl_zoomout0t1) {
                  line.plot <- line.plot + ylim(0, 1)
                } else {
                  line.plot <- line.plot +
                    ylim(fl_gini_min - abs(fl_gini_max-fl_gini_min)/20, fl_gini_max + abs(fl_gini_max-fl_gini_min)/20)
                }

                # legend area
                line.plot <- line.plot +
                  theme(
                    text = element_text(size = 16),
                    legend.position = c(0.22, 0.75),
                    legend.background = element_rect(fill = "white", colour = "black",
                                                     linetype='solid'))

                line.plot <- line.plot +
                  labs(colour = "Check limit per adult and child",
                       shape = "Check limit per adult and child",
                       linetype = "Check limit per adult and child")

                # # Labeling
                # line.plot$labels$linetype <- "REV\nOptimal\nvs.\nAlternatives"
                # line.plot$labels$colour <- "REV\nOptimal\nvs.\nAlternatives"
                # line.plot$labels$shape <- "REV\nOptimal\nvs.\nAlternatives"

                # Print
                print(line.plot)

                # add png
                if (bl_zoomout0t1) {
                  st_zoom <- '_0t1'
                } else {
                  st_zoom <- '_zoomed'
                }
                if (bl_nonone) {
                  st_actual <- '_nonone'
                } else {
                  st_actual <- '_withnone'
                }

                if (bl_approximate) {
                  snm_save_png <- paste0('sd_log_c_feasible_rhomultiple_smooth', st_zoom, st_actual, st_income_bound, '.png')
                } else {
                  snm_save_png <- paste0('sd_log_c_feasible_rhomultiple_pointline', st_zoom, st_actual, st_income_bound, '.png')
                }

                # Save
                ggsave(line.plot,
                       file=file.path(srt_imgcsv_rev_root, snm_save_png),
                       width = 270,
                       height = 216, units='mm',
                       dpi = 300)

              }

              # FIGURE 6 Outcome Mean -----------------------------------
              if (bl_outcome_mean_graph) {

                # collect only relevant columns
                mt_mean_EH_across_rhos_bounds_selected <- mt_rev_across_rhos_bounds %>%
                  filter(objective == "c2020") %>%
                  mutate(one_minus_rho = 1 - rho_val) %>%
                  select(one_minus_rho, mean_drv_EH_star, bounds) %>%
                  rename(Outcome = bounds)
                # Actual outcome
                mt_mean_EH_across_rhos_bounds_actual <- mt_rev_across_rhos_bounds %>%
                  filter(objective == "c2020") %>%
                  mutate(one_minus_rho = 1 - rho_val) %>%
                  select(one_minus_rho, mean_EH_drv_dact, bounds) %>%
                  rename(Outcome = bounds) %>%
                  filter(Outcome == 'b1_a64_14ca14ck_18t64') %>%
                  select(one_minus_rho, mean_EH_drv_dact) %>%
                  rename(mean_drv_EH_star = mean_EH_drv_dact) %>%
                  mutate(Outcome = 'Actual')
                # Zero outcome
                mt_mean_EH_across_rhos_bounds_none <- mt_rev_across_rhos_bounds %>%
                  filter(objective == "c2020") %>%
                  mutate(one_minus_rho = 1 - rho_val) %>%
                  select(one_minus_rho, mean_EH_drv_d0, bounds) %>%
                  rename(Outcome = bounds) %>%
                  filter(Outcome == 'b1_a64_14ca14ck_18t64') %>%
                  select(one_minus_rho, mean_EH_drv_d0) %>%
                  rename(mean_drv_EH_star = mean_EH_drv_d0) %>%
                  mutate(Outcome = 'None')
                # Stack Actual and Optimal
                mt_mean_EH_across_rhos_bounds_selected <-
                  rbind(mt_mean_EH_across_rhos_bounds_none,
                        mt_mean_EH_across_rhos_bounds_actual,
                        mt_mean_EH_across_rhos_bounds_selected)

                # Relabel Variable
                if (bl_nonone) {
                  Outcome_levels <- c("Actual allocations" = "Actual",
                                      "$1,400 limit optimal" = "b1_a64_14ca14ck_18t64",
                                      "$2,000 limit optimal" = "b1_a64_20ca20ck_18t64")
                  mt_mean_EH_across_rhos_bounds_selected_fig <- mt_mean_EH_across_rhos_bounds_selected %>%
                    filter(Outcome == 'b1_a64_14ca14ck_18t64'
                           | Outcome == 'b1_a64_20ca20ck_18t64'
                           | Outcome == 'Actual') %>%
                    mutate(Outcome = as_factor(Outcome)) %>%
                    mutate(Outcome = fct_recode(Outcome, !!!Outcome_levels)) %>%
                    mutate(mean_drv_EH_star = mean_drv_EH_star*58056)
                } else {
                  Outcome_levels <- c("No checks" = "None",
                                      "Actual allocations" = "Actual",
                                      "$1,400 limit optimal" = "b1_a64_14ca14ck_18t64",
                                      "$2,000 limit optimal" = "b1_a64_20ca20ck_18t64")
                  mt_mean_EH_across_rhos_bounds_selected_fig <- mt_mean_EH_across_rhos_bounds_selected %>%
                    filter(Outcome == 'b1_a64_14ca14ck_18t64'
                           | Outcome == 'b1_a64_20ca20ck_18t64'
                           | Outcome == 'Actual'
                           | Outcome == 'None') %>%
                    mutate(Outcome = as_factor(Outcome)) %>%
                    mutate(Outcome = fct_recode(Outcome, !!!Outcome_levels)) %>%
                    mutate(mean_drv_EH_star = mean_drv_EH_star*58056)
                }


                if (bl_approximate) {
                  # Graph Results--Draw
                  line.plot <- mt_mean_EH_across_rhos_bounds_selected_fig %>%
                    ggplot(aes(x=one_minus_rho, y=mean_drv_EH_star,
                               group=Outcome,
                               colour=Outcome,
                               linetype=Outcome,
                               shape=Outcome)) +
                    # geom_point() +
                    # geom_line() +
                    geom_smooth(span = 0.35, se=FALSE, size=2)
                } else {
                  # Graph Results--Draw
                  line.plot <- mt_mean_EH_across_rhos_bounds_selected_fig %>%
                    ggplot(aes(x=one_minus_rho, y=mean_drv_EH_star,
                               group=Outcome,
                               colour=Outcome,
                               linetype=Outcome,
                               shape=Outcome)) +
                    geom_point(size=6) +
                    geom_line(size=2)
                  # geom_smooth(span = 0.35, se=FALSE, size=2)
                  # geom_vline(xintercept=c(1), linetype="dotted") +
                }

                fl_mean_min <- min(mt_mean_EH_across_rhos_bounds_selected$mean_drv_EH_star)
                fl_mean_max <- max(mt_mean_EH_across_rhos_bounds_selected$mean_drv_EH_star)
                fl_mean_min <- fl_mean_min*58056
                fl_mean_max <- fl_mean_max*58056

                line.plot <- line.plot +
                  labs(x = 'Planner inequality aversion, \u03bb',
                       y = paste0('Mean consumption per-capita given stimulus checks (USD)')) +
                  scale_x_continuous(trans='log10', labels = x.labels, breaks = x.breaks) +
                  scale_y_continuous(n.breaks = it_y_nbreaks) +
                  theme_bw(base_size=8)

                line.plot <- line.plot +
                  ylim(fl_mean_min - abs(fl_mean_max-fl_mean_min)/20, fl_mean_max + abs(fl_mean_max-fl_mean_min)/20)


                # legend area
                line.plot <- line.plot +
                  theme(
                    text = element_text(size = 16),
                    legend.position = c(0.22, 0.75),
                    legend.background = element_rect(fill = "white", colour = "black",
                                                     linetype='solid'))

                line.plot <- line.plot +
                  labs(colour = "Check limit per adult and child",
                       shape = "Check limit per adult and child",
                       linetype = "Check limit per adult and child")

                # # Labeling
                # line.plot$labels$linetype <- "REV\nOptimal\nvs.\nAlternatives"
                # line.plot$labels$colour <- "REV\nOptimal\nvs.\nAlternatives"
                # line.plot$labels$shape <- "REV\nOptimal\nvs.\nAlternatives"

                # Print
                print(line.plot)

                if (bl_nonone) {
                  st_actual <- '_nonone'
                } else {
                  st_actual <- '_withnone'
                }
                # add png
                if (bl_approximate) {
                  snm_save_png <- paste0('mean_EH_feasible_rhomultiple_smooth', st_actual, st_income_bound, '.png')
                } else {
                  snm_save_png <- paste0('mean_EH_feasible_rhomultiple_pointline', st_actual, st_income_bound, '.png')
                }

                # Save
                ggsave(line.plot,
                       file=file.path(srt_imgcsv_rev_root, snm_save_png),
                       width = 270,
                       height = 216, units='mm',
                       dpi = 300)

              }

              # FIGURE 7 Outcome Min -----------------------------------
              if (bl_outcome_min_graph && it_income_bounds == 1) {

                # collect only relevant columns
                mt_min_EH_across_rhos_bounds_selected <- mt_rev_across_rhos_bounds %>%
                  filter(objective == "c2020") %>%
                  mutate(one_minus_rho = 1 - rho_val) %>%
                  select(one_minus_rho, min_EH_star, bounds) %>%
                  rename(Outcome = bounds)
                # Actual outcome
                mt_min_EH_across_rhos_bounds_actual <- mt_rev_across_rhos_bounds %>%
                  filter(objective == "c2020") %>%
                  mutate(one_minus_rho = 1 - rho_val) %>%
                  select(one_minus_rho, min_EH_dact, bounds) %>%
                  rename(Outcome = bounds) %>%
                  filter(Outcome == 'b1_a64_14ca14ck_18t64') %>%
                  select(one_minus_rho, min_EH_dact) %>%
                  rename(min_EH_star = min_EH_dact) %>%
                  mutate(Outcome = 'Actual')
                # Zero outcome
                mt_min_EH_across_rhos_bounds_none <- mt_rev_across_rhos_bounds %>%
                  filter(objective == "c2020") %>%
                  mutate(one_minus_rho = 1 - rho_val) %>%
                  select(one_minus_rho, min_EH_d0, bounds) %>%
                  rename(Outcome = bounds) %>%
                  filter(Outcome == 'b1_a64_14ca14ck_18t64') %>%
                  select(one_minus_rho, min_EH_d0) %>%
                  rename(min_EH_star = min_EH_d0) %>%
                  mutate(Outcome = 'None')
                # Stack Actual and Optimal
                mt_min_EH_across_rhos_bounds_selected <-
                  rbind(mt_min_EH_across_rhos_bounds_none,
                        mt_min_EH_across_rhos_bounds_actual,
                        mt_min_EH_across_rhos_bounds_selected)

                # Relabel Variable
                if (bl_nonone) {
                  Outcome_levels <- c("Actual allocations" = "Actual",
                                      "$1,400 limit optimal" = "b1_a64_14ca14ck_18t64",
                                      "$2,000 limit optimal" = "b1_a64_20ca20ck_18t64")
                  mt_min_EH_across_rhos_bounds_selected_fig <- mt_min_EH_across_rhos_bounds_selected %>%
                    filter(Outcome == 'b1_a64_14ca14ck_18t64'
                           | Outcome == 'b1_a64_20ca20ck_18t64'
                           | Outcome == 'Actual') %>%
                    mutate(Outcome = as_factor(Outcome)) %>%
                    mutate(Outcome = fct_recode(Outcome, !!!Outcome_levels)) %>%
                    mutate(min_EH_star = min_EH_star*58056)
                } else {
                  Outcome_levels <- c("No checks" = "None",
                                      "Actual allocations" = "Actual",
                                      "$1,400 limit optimal" = "b1_a64_14ca14ck_18t64",
                                      "$2,000 limit optimal" = "b1_a64_20ca20ck_18t64")
                  mt_min_EH_across_rhos_bounds_selected_fig <- mt_min_EH_across_rhos_bounds_selected %>%
                    filter(Outcome == 'b1_a64_14ca14ck_18t64'
                           | Outcome == 'b1_a64_20ca20ck_18t64'
                           | Outcome == 'Actual'
                           | Outcome == 'None') %>%
                    mutate(Outcome = as_factor(Outcome)) %>%
                    mutate(Outcome = fct_recode(Outcome, !!!Outcome_levels)) %>%
                    mutate(min_EH_star = min_EH_star*58056)
                }

                if (bl_approximate) {
                  # Graph Results--Draw
                  line.plot <- mt_min_EH_across_rhos_bounds_selected_fig %>%
                    ggplot(aes(x=one_minus_rho, y=min_EH_star,
                               group=Outcome,
                               colour=Outcome,
                               linetype=Outcome,
                               shape=Outcome)) +
                    # geom_point() +
                    # geom_line() +
                    geom_smooth(span = 0.35, se=FALSE, size=2)
                } else {
                  # Graph Results--Draw
                  line.plot <- mt_min_EH_across_rhos_bounds_selected_fig %>%
                    ggplot(aes(x=one_minus_rho, y=min_EH_star,
                               group=Outcome,
                               colour=Outcome,
                               linetype=Outcome,
                               shape=Outcome)) +
                    geom_point(size=6) +
                    geom_line(size=2)
                  # geom_smooth(span = 0.35, se=FALSE, size=2)
                  # geom_vline(xintercept=c(1), linetype="dotted") +
                }

                fl_gini_min <- min(mt_gini_across_rhos_bounds_selected_fig$gini_drv_EH_star)
                fl_gini_max <- max(mt_gini_across_rhos_bounds_selected_fig$gini_drv_EH_star)

                line.plot <- line.plot +
                  labs(x = 'Planner inequality aversion, \u03bb',
                       y = paste0('Min consumption per-capita given stimulus checks (USD)')) +
                  scale_x_continuous(trans='log10', labels = x.labels, breaks = x.breaks) +
                  scale_y_continuous(n.breaks = it_y_nbreaks) +
                  theme_bw(base_size=8)

                line.plot <- line.plot +
                  ylim(fl_gini_min - abs(fl_gini_max-fl_gini_min)/20, fl_gini_max + abs(fl_gini_max-fl_gini_min)/20)

                # legend area
                line.plot <- line.plot +
                  theme(
                    text = element_text(size = 16),
                    legend.position = c(0.22, 0.75),
                    legend.background = element_rect(fill = "white", colour = "black",
                                                     linetype='solid'))

                line.plot <- line.plot +
                  labs(colour = "Check limit per adult and child",
                       shape = "Check limit per adult and child",
                       linetype = "Check limit per adult and child")

                # # Labeling
                # line.plot$labels$linetype <- "REV\nOptimal\nvs.\nAlternatives"
                # line.plot$labels$colour <- "REV\nOptimal\nvs.\nAlternatives"
                # line.plot$labels$shape <- "REV\nOptimal\nvs.\nAlternatives"

                # Print
                print(line.plot)

                if (bl_nonone) {
                  st_actual <- '_nonone'
                } else {
                  st_actual <- '_withnone'
                }
                # add png
                if (bl_approximate) {
                  snm_save_png <- paste0('min_EH_feasible_rhomultiple_smooth', st_actual, st_income_bound, '.png')
                } else {
                  snm_save_png <- paste0('min_EH_feasible_rhomultiple_pointline', st_actual, st_income_bound, '.png')
                }

                # Save
                ggsave(line.plot,
                       file=file.path(srt_imgcsv_rev_root, snm_save_png),
                       width = 270,
                       height = 216, units='mm',
                       dpi = 300)

              }
            }
          }
        }
      }
    }
  }
}
