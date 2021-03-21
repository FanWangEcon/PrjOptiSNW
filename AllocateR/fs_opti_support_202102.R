fs_opti_support_202102 <- function(st_which_solu='snwx_trumpchk_moredense_a65zh266zs5_b1_xi0_manna_168_married') {
  #' Same as fs_opti_support, but shorter for just the beta/edu weighted results
  #'
  #' @description
  #' see fs_opti_support
  #'
  #' @params st_which_solu string name of the csv file.
  #'

  # Path
  if (dir.exists(file.path('D:/Dropbox (UH-ECON)/PrjNygaardSorensenWang/'))) {
    srt_root <- 'D:/Dropbox (UH-ECON)/PrjNygaardSorensenWang/'
  } else {
    srt_root <- 'C:/Users/fan/Documents/Dropbox (UH-ECON)/PrjNygaardSorensenWang/'
  }

  bl_save_img <- TRUE

  # Preferences
  # ar_rho <- 1 - (10^(c(seq(-2,2, length.out=8))))
  # ar_rho <- unique(c(1,ar_rho))
  ar_rho <- c(1)

  srn_csv <- 'csv/'
  # Allocation Type specific folder etc
  st_b0b1 <- 'b1'

  # cleaning up csv file name, to achieve two goals:
  # 1. shorten file names slightly so that folder path does not get too long
  # 2. add on date prefix in case there are date-specific things
  snm_which_solu_folder <- gsub(x = st_which_solu,  pattern = "moredense_a65zh266zs5", replacement = "a65zh266zs5")
  snm_which_solu_folder <- gsub(x = snm_which_solu_folder,  pattern = "snwx_", replacement = "")
  srt_folder <- paste0('20210225_', snm_which_solu_folder)
  # st_file_type_withspouse_shock <- paste0('moredense_', st_which_solu)
  # st_file_type_withspouse_shock <- paste0(st_which_solu)
  st_file_type_withspouse_shock <- snm_which_solu_folder

  # CSV file to load
  # snm_simu_csv_withspouse_shock <- paste0('snwx_bidenchk_',st_file_type_withspouse_shock,'.csv')
  snm_simu_csv_withspouse_shock <- paste0(st_which_solu,'.csv')

  # Results and Output folder
  srt_results_root <- paste0(srt_root, 'Results202102')
  srt_img_save_root <- paste0(srt_root, 'Results202102/', srt_folder, '/Graphs/')
  srt_csv_path_root <- paste0(srt_root, 'Results202102/', srt_folder, '/', srn_csv)
  srt_simu_path <- paste0(srt_root, 'Output202102/')

  srt_paper_main_textgraph <- paste0(srt_root, 'Paper/Main_text_graphs_and_tables/')
  srt_paper_appendix_textgraph <- paste0(srt_root, 'Paper/Appendix_graphs_and_tables/')

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
