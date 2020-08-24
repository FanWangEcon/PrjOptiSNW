fs_opti_support <- function() {
# other functions might overrid what is specified here.

# ls_output <- fs_opti_support()
# st_file_type_withspouse_shock <- ls_output$st_file_type_withspouse_shock
# snm_simu_csv_withspouse_shock <- ls_output$snm_simu_csv_withspouse_shock
# srt_simu_path <- ls_output$srt_simu_path
# bl_save_img <- ls_output$bl_save_img
# spt_img_save <- ls_output$spt_img_save
# srt_csv_path <- ls_output$srt_csv_path

# This function store a number of support strings. for folder output path,
# file name, etc that are shared across feasiable, optimal and threshold allocations

srt_root <- 'D:/Dropbox (UH-ECON)/PrjNygaardSorensenWang/'
# srt_root <- 'C:/Users/fan/Documents/Dropbox (UH-ECON)/PrjNygaardSorensenWang/'
# st_b0b1 <- 'b1_spouseinc'
st_b0b1 <- 'b1'
st_file_type_withspouse_shock <- paste0('moredense_a65zh266zs5_e2m2_',st_b0b1)
# st_file_type_withspouse_shock <- paste0('dense_ybin2500')
snm_simu_csv_withspouse_shock <- paste0('snwx_v_planner_',st_file_type_withspouse_shock,'_spouseinc.csv')

srt_simu_path <- paste0(srt_root, 'Output/')

bl_save_img <- TRUE
srt_img_save_root <- paste0(srt_root, 'Results/2020-08-23/Graphs/')
srt_csv_path_root <- paste0(srt_root, 'Results/2020-08-23/csv/')
srt_paper_main_textgraph <- paste0(srt_root, 'Paper/Main_text_graphs_and_tables/')
srt_paper_appendix_textgraph <- paste0(srt_root, 'Paper/Appendix_graphs_and_tables/')

# Preferences
ar_rho <- 1 - (10^(c(seq(-2,2, length.out=8))))
ar_rho <- unique(c(1,ar_rho))
# ar_rho <- c(1)

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
