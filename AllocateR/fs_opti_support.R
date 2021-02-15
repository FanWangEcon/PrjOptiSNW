fs_opti_support <- function(st_which_solu='b1_manna') {
  #' Gateway function to access string objects including paths
  #'
  #' @description
  #' Get file path, depending on simulation. There are three key ones: 1,
  #' allocation without tax changes b=1, xi=anything, meaning there is full
  #' unemployment benefits under covid, and no tax change; 2, allocation
  #' with tax change b=1, xi=0, max impact of covid shock if shock no income
  #' but fully pay unemployment benefits; 3, allocation without tax change
  #' and b=0 and xi=0.25, meaning that no unemployment benefit, very bad shock.
  #'
  #' @params st_which_solu string which one of the three problems been solved.
  #'

  # Path
  if (dir.exists(file.path('D:/Dropbox (UH-ECON)/PrjNygaardSorensenWang/'))) {
    srt_root <- 'D:/Dropbox (UH-ECON)/PrjNygaardSorensenWang/'
  } else {
    srt_root <- 'C:/Users/fan/Documents/Dropbox (UH-ECON)/PrjNygaardSorensenWang/'
  }

  bl_save_img <- TRUE

  # Preferences
  ar_rho <- 1 - (10^(c(seq(-2,2, length.out=8))))
  ar_rho <- unique(c(1,ar_rho))
  # ar_rho <- c(1)

  srn_csv <- 'csv/'
  # Allocation Type specific folder etc
  ls_st_submit2020 <- c('b1_manna', 'b1_manna_expmin', 'b1_xi0_tax', 'b0_xi0p25_manna')

  if (grepl('snwx_v_planner_', st_which_solu)) {
    # supplying full string names
    # started doing this when docdense files are been used.

    if (grepl('b1_', st_which_solu)) {
      st_b0b1 <- 'b1'
    }
    if (grepl('docdense_e1lm2_b1_xi0_manna', st_which_solu) ||
        grepl('docdense_e2hm2_b1_xi0_manna', st_which_solu)) {
      srt_folder <- paste0('20210106_', st_which_solu)
    }

    st_file_type_withspouse_shock <- gsub(x = st_which_solu,pattern = "snwx_v_planner_", replacement = "")

  } else {
    # use string prefix

    if (st_which_solu %in% ls_st_submit2020) {
        if (st_which_solu == 'b1_manna') {
          st_b0b1 <- 'b1'
          srt_folder <- '2020-08-23-b1_manna'
          st_file_type_withspouse_shock <- paste0('moredense_a65zh266zs5_e2m2_', st_which_solu, '_244')
        } else if (st_which_solu == 'b1_manna_expmin') {
          st_b0b1 <- 'b1'
          srt_folder <- '2020-08-23-b1_manna'
          st_file_type_withspouse_shock <- paste0('moredense_a65zh266zs5_e2m2_', st_which_solu, '_244')
          srn_csv <- 'csv_with_expmin_apc/'
          bl_save_img <- FALSE
          ar_rho <- c(1)
        } else if (st_which_solu == 'b1_xi0_tax') {
          st_b0b1 <- 'b1'
          srt_folder <- '2020-08-23-b1_xi0_tax'
          st_file_type_withspouse_shock <- paste0('moredense_a65zh266zs5_e2m2_', st_which_solu, '_244')
        } else if (st_which_solu == 'b0_xi0p25_manna') {
          st_b0b1 <- 'b0'
          srt_folder <- '2020-08-23-b0_xi0p25_manna'
          st_file_type_withspouse_shock <- paste0('moredense_a65zh266zs5_e2m2_', st_which_solu, '_244')
        }
    } else {
        st_b0b1 <- 'b1'
        srt_folder <- paste0('20210106_', st_which_solu)
        st_file_type_withspouse_shock <- paste0('moredense_a65zh266zs5_e2m2_', st_which_solu)
    }
  }

  # st_file_type_withspouse_shock <- paste0('dense_ybin2500')
  snm_simu_csv_withspouse_shock <- paste0('snwx_v_planner_',st_file_type_withspouse_shock,'.csv')

  # Results and Output folder
  if (st_which_solu %in% ls_st_submit2020) {
      srt_results_root <- paste0(srt_root, 'Results')
      srt_img_save_root <- paste0(srt_root, 'Results/', srt_folder, '/Graphs/')
      srt_csv_path_root <- paste0(srt_root, 'Results/', srt_folder, '/', srn_csv)
      srt_simu_path <- paste0(srt_root, 'Output/')
  } else {
      srt_results_root <- paste0(srt_root, 'Results202101')
      srt_img_save_root <- paste0(srt_root, 'Results202101/', srt_folder, '/Graphs/')
      srt_csv_path_root <- paste0(srt_root, 'Results202101/', srt_folder, '/', srn_csv)
      srt_simu_path <- paste0(srt_root, 'Output202101/')
  }

  srt_paper_main_textgraph <- paste0(srt_root, 'Paper/Main_text_graphs_and_tables/')
  srt_paper_appendix_textgraph <- paste0(srt_root, 'Paper/Appendix_graphs_and_tables/')

  # Ages
  # ar_it_max_age <- c(69, 74, 79, 84)
  # ar_it_max_age <- c(69, 74)
  # for (it_max_age in ar_it_max_age) {

  # Change max age
  it_max_age <- 64
  st_file_folder <- paste0(st_b0b1, '_a', it_max_age, '')

  # CSV and Image Paths
  spt_img_save = paste0(srt_img_save_root, st_file_folder)
  srt_csv_path = paste0(srt_csv_path_root, st_file_folder)
  dir.create(file.path(spt_img_save), showWarnings = FALSE, recursive = TRUE)
  dir.create(file.path(srt_csv_path), showWarnings = FALSE, recursive = TRUE)

  ls_output <- list(st_file_type_withspouse_shock=st_file_type_withspouse_shock,
                    snm_simu_csv_withspouse_shock=snm_simu_csv_withspouse_shock,
                    st_b0b1=st_b0b1,
                    srt_simu_path=srt_simu_path,
                    bl_save_img=bl_save_img,
                    srt_results_root=srt_results_root,
                    srt_img_save_root=srt_img_save_root,
                    srt_csv_path_root=srt_csv_path_root,
                    srt_paper_main_textgraph=srt_paper_main_textgraph,
                    srt_paper_appendix_textgraph=srt_paper_appendix_textgraph,
                    spt_img_save=spt_img_save,
                    srt_csv_path=srt_csv_path,
                    ar_rho=ar_rho,
                    it_max_age=it_max_age)

  return(ls_output)
}
