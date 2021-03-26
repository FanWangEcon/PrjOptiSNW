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
# combine together
ar_rho <- c(ar_rho_init)
ar_rho <- sort(ar_rho)

# log10((-1)*(ar_rho - 1))
# 1/(1-ar_rho)

# Results from which planners should
ls_it_rho_combo_type <- c(1)

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

ls_st_file_suffix_bchklock <-
  c('snwx_bchklock_moredense_a65zh266zs5_b1_xi0_manna_168_bt95',
    'snwx_bchklock_moredense_a65zh266zs5_b1_xi0_manna_168_bt60',
    'snwx_bchklock_moredense_a65zh266zs5_b1_xi0_manna_168_married',
    'snwx_bchklock_moredense_a65zh266zs5_b1_xi0_manna_168_unmarried',
    'snwx_bchklock_moredense_a65zh266zs5_b1_xi0_manna_168')
ls_st_file_suffix_bchklock <- rev(ls_st_file_suffix_bchklock)

# list to run
ls_st_file_suffix <- c(ls_st_file_suffix_trumpchk,
                       ls_st_file_suffix_bidenchk,
                       ls_st_file_suffix_bchklock)
ls_st_file_suffix <- c('snwx_bidenchk_tiny_b1_xi0_manna_168')
ls_st_file_suffix <- c('snwx_bidenchk_moredense_a65zh266zs5_b1_xi0_manna_168')
# Per capita or per household results
ls_bl_per_capita <- c(TRUE)

# Allocation bounds types
ls_st_bound_files <- c('14ca14ck', '17ca17ck', '20ca20ck')

# Loop over file types1
for (bl_per_capita in ls_bl_per_capita) {
  for (st_which_solu in ls_st_file_suffix) {

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
        ls_output <- fs_opti_support_202103(st_which_solu, bl_per_capita=bl_per_capita, fl_rho=fl_rho)
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
        srt_csv_file <- paste0('rev_feasible_', st_bound_files, '.csv')
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
      snm_save_csv <- paste0('rev_feasible_', st_bound_files, '_rhomultiple.csv')
      write.csv(mt_rev_across_rhos,
                paste0(srt_imgcsv_rev_root, snm_save_csv),
                row.names = TRUE)

      # Graph
      mt_rev_across_rhos %>% filter(objective == "c2020")

    }

    # cAll results from different bounds to the same file
    mt_rev_across_rhos_bounds <- do.call(rbind, ls_of_ls_rev_results)

    # collect only relevant columns
    mt_rev_across_rhos_bounds_selected <- mt_rev_across_rhos_bounds %>%
      filter(objective == "c2020") %>%
      mutate(REV = REV*100) %>%
      mutate(one_minus_rho = 1 - rho_val) %>%
      select(one_minus_rho, REV, bounds) %>%
      rename(Outcome = bounds, rev = REV)

    # graph, copied from https://fanwangecon.github.io/PrjOptiAlloc/articles/ffv_opt_sobin_rkone_allrw_training.html
    # x-labels
    x.labels <- c('\u03bb \U2248 1.00', '\u03bb = 0.90', '\u03bb = 0', '\u03bb = -10', '\u03bb = -100')
    x.breaks <- c(0.01, 0.10, 1, 10, 100)

    # title line 2
    # title_line1 <- sprintf("Percentage of Training Spots Misallocated, NSW Lalonde (AER, 1986)")
    # title_line2 <- sprintf("REV (Resource Equivalent Variation) Along Planner Spectrum")

    st_title <- sprintf(paste0('How much Fewer Resources are Needed (Shares) to Achieve the Same Welfare'))
    title_line1 <- sprintf("Compare alternative allocations to optimal allocations given observables and estimates")
    title_line2 <- sprintf("Solid Red Line: train 297 random NSW treatment individuals vs. optimally allocating 297 spots")
    title_line3 <- sprintf("Dashed Blue Line: train 297 lowest baseline wage individuals vs. optimally allocating 297 spots")

    # Relabel Variable
    Outcome_levels <- c("$1,400" = "b1_a64_14ca14ck_18t64",
                        "$2,000" = "b1_a64_20ca20ck_18t64")
    mt_rev_across_rhos_bounds_selected_fig <- mt_rev_across_rhos_bounds_selected %>%
      filter(Outcome == 'b1_a64_14ca14ck_18t64'
             | Outcome == 'b1_a64_20ca20ck_18t64') %>%
      mutate(Outcome = as_factor(Outcome)) %>%
      mutate(Outcome = fct_recode(Outcome, !!!Outcome_levels))

    # Graph Results--Draw
    line.plot <- mt_rev_across_rhos_bounds_selected_fig %>%
      ggplot(aes(x=one_minus_rho, y=rev,
                 group=Outcome,
                 colour=Outcome,
                 linetype=Outcome,
                 shape=Outcome)) +
      # geom_point() +
      # geom_line() +
      geom_smooth(span = 0.35, se=FALSE, size=2) +
      # geom_vline(xintercept=c(1), linetype="dotted") +
      labs(x = 'Planner inequality aversion, \u03bb',
           y = paste0('Resource equivalent variation (percent)')) +
      scale_x_continuous(trans='log10', labels = x.labels, breaks = x.breaks) +
      theme_bw(base_size=8) +
      ylim(0, 100)

    # legend area
    line.plot <- line.plot +
      theme(
        text = element_text(size = 16),
        legend.position = c(0.18, 0.85),
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
    snm_save_png <- paste0('rev_feasible_rhomultiple.png')

    # Save
    ggsave(line.plot,
           file=file.path(srt_imgcsv_rev_root, snm_save_png),
           width = 270,
           height = 216, units='mm',
           dpi = 300)
  }
}

