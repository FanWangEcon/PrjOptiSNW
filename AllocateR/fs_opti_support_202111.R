fs_opti_support_202111 <- function(st_which_solu='snwx_trumpchk_moredense_a65zh266zs5_b1_xi0_manna_168_married',
                                   bl_per_capita=FALSE, fl_rho=1,
                                   snm_main='Results202111',
                                   it_rand_draw_seed=NULL,
                                   ls_st_gen_agg_folders=c('aggregator', 'aggcsvsolu', 'apcmpc', 'rev')) {
  #' Similar to fs_opti_support_202103, but for Nov 2021 revisions.
  #'
  #' Random perturbation changes the A values for each allocation group.
  #'
  #' @description
  #' see fs_opti_support
  #'
  #' @params st_which_solu string name of the csv file.
  #' @params bl_per_capita boolean if true per-capita folder name, if false, do not append.
  #' @params fl_rho foat the CEs parameter, default is 1, no changes to folder name, but otherwise change 1
  #' is for the Utilitarian planner.
  #' @params it_rand_draw_seed is the random seed used for perturbation draws
  #' @params ls_st_gen_agg_folders list of strings of types of aggregation folders to generate
  #' in the current results folder. Possible values include: c('aggregator', 'aggcsvsolu', 'apcmpc', 'rev')
  #'

  # Path
  if (dir.exists(file.path('D:/Dropbox (UH-ECON)/PrjNygaardSorensenWang/'))) {
    srt_root <- 'D:/Dropbox (UH-ECON)/PrjNygaardSorensenWang/'
  } else {
    srt_root <- 'C:/Users/fan/Documents/Dropbox (UH-ECON)/PrjNygaardSorensenWang/'
  }

  # Folder suffix for per capita
  if (bl_per_capita) {
    st_per_capita <- 'capi'
  } else {
    st_per_capita <- 'hhds'
  }

  # Seed
  if(!is.null(it_rand_draw_seed)) {
    st_rng_seed <- paste0('_perturb_rng', it_rand_draw_seed)
  } else {
    st_rng_seed <- ''
  }


  # Preferences
  # Folder suffix CES parameter
  st_rho <- paste0(fl_rho)
  st_rho <- gsub(x = st_rho,  pattern = "-", replacement = "n")
  st_rho <- gsub(x = st_rho,  pattern = "\\.", replacement = "p")
  st_rho <- paste0(st_rho)
  st_rho <- paste0('_rh', st_rho)
  st_rho_nounderscore <- paste0('rh', st_rho)

  # ar_rho <- 1 - (10^(c(seq(-2,2, length.out=8))))
  # ar_rho <- unique(c(1,ar_rho))
  # ar_rho <- c(1)
  ar_rho <- c(fl_rho)

  # Folder name, main folder below root storing results.

  # rho specific folder
  snm_full_main_folder <- file.path(snm_main, paste0(st_per_capita, st_rho, st_rng_seed))
  dir.create(file.path(snm_full_main_folder), showWarnings = FALSE, recursive = TRUE)

  # these are unlikely to be needed for perturb exercises.
  # bl_gen_agg_folders <- FALSE
  # snm_image_aggregator_folder: folder where various rho results aggregated together
  snm_image_aggregator_folder <- file.path(snm_main, paste0(st_per_capita, '_rhmultiple'))
  # snm_image_aggcsvsolu_folder: folder where results over multiple csv solu stored together
  snm_image_aggcsvsolu_folder <- file.path(snm_main, paste0(st_per_capita, '_csvsolumultiple'))
  # snm_apcmpc_folder: folder where various rho results aggregated together
  snm_apcmpc_folder <- file.path(snm_main, paste0('apc_mpc'))
  # snm_apcmpc_folder: folder where various rho results aggregated together
  snm_rev_folder <- file.path(snm_main, paste0(st_per_capita, '_rhomultiple_rev'))
  # Folder for perturbed results aggregated
  snm_perturb_folder <- file.path(snm_main, paste0(st_per_capita, '_perturb'))

  if ('aggregator' %in% ls_st_gen_agg_folders){
    dir.create(file.path(snm_image_aggregator_folder), showWarnings = FALSE, recursive = TRUE)
  }
  if ('aggcsvsolu' %in% ls_st_gen_agg_folders){
    dir.create(file.path(snm_image_aggcsvsolu_folder), showWarnings = FALSE, recursive = TRUE)
  }
  if ('apcmpc' %in% ls_st_gen_agg_folders){
    dir.create(file.path(snm_apcmpc_folder), showWarnings = FALSE, recursive = TRUE)
  }
  if ('rev' %in% ls_st_gen_agg_folders){
    dir.create(file.path(snm_rev_folder), showWarnings = FALSE, recursive = TRUE)
  }
  if ('perturb' %in% ls_st_gen_agg_folders){
    dir.create(file.path(snm_perturb_folder), showWarnings = FALSE, recursive = TRUE)
  }

  bl_save_img <- TRUE

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
  srt_results_root <- file.path(srt_root, snm_full_main_folder, '/')
  srt_img_save_root <- file.path(srt_root, snm_full_main_folder, srt_folder, 'Graphs', '/')
  srt_csv_path_root <- file.path(srt_root, snm_full_main_folder, srt_folder, 'csv', '/')
  srt_simu_path <- file.path(srt_root, 'Output202111', '/')

  # Aggregate image save folder, joint over rho
  srt_img_aggregator_save_root <- file.path(srt_root, snm_image_aggregator_folder, srt_folder, '/')
  # Aggregate image save folder, but now joint over csvsolu
  srt_img_aggcsvsolu_save_root <- file.path(srt_root, snm_image_aggcsvsolu_folder, st_rho_nounderscore, '/')
  # Aggregate MPC and Related Scatter plot folder
  srt_imgcsv_mpcapc_root <- file.path(srt_root, snm_apcmpc_folder, srt_folder, '/')
  # REV folder
  srt_imgcsv_rev_root <- file.path(srt_root, snm_rev_folder, srt_folder, '/')
  # Perturb folder
  srt_imgcsv_perturb_root <- file.path(srt_root, snm_perturb_folder, srt_folder, '/')
  if ('aggregator' %in% ls_st_gen_agg_folders){
    dir.create(file.path(srt_img_aggregator_save_root), showWarnings = FALSE, recursive = TRUE)
  }
  if ('aggcsvsolu' %in% ls_st_gen_agg_folders){
    dir.create(file.path(srt_img_aggcsvsolu_save_root), showWarnings = FALSE, recursive = TRUE)
  }
  if ('apcmpc' %in% ls_st_gen_agg_folders){
    dir.create(file.path(srt_imgcsv_mpcapc_root), showWarnings = FALSE, recursive = TRUE)
  }
  if ('rev' %in% ls_st_gen_agg_folders){
    dir.create(file.path(srt_imgcsv_rev_root), showWarnings = FALSE, recursive = TRUE)
  }
  if ('perturb' %in% ls_st_gen_agg_folders){
    dir.create(file.path(srt_imgcsv_perturb_root), showWarnings = FALSE, recursive = TRUE)
  }


  srt_paper_main_textgraph <- file.path(srt_root, 'Paper', 'Main_text_graphs_and_tables', '/')
  srt_paper_appendix_textgraph <- file.path(srt_root, 'Paper', 'Appendix_graphs_and_tables', '/')

  # Change max age
  it_max_age <- 64
  st_file_folder <- paste0(st_b0b1, '_a', it_max_age, '')

  # CSV and Image Paths
  srt_csv_path = file.path(srt_csv_path_root, st_file_folder)
  dir.create(file.path(srt_csv_path), showWarnings = FALSE, recursive = TRUE)
  spt_img_save = file.path(srt_img_save_root, st_file_folder)
  if (bl_save_img) {
    dir.create(file.path(spt_img_save), showWarnings = FALSE, recursive = TRUE)
  }

  ls_output <- list(st_file_type_withspouse_shock=st_file_type_withspouse_shock,
                    snm_simu_csv_withspouse_shock=snm_simu_csv_withspouse_shock,
                    st_b0b1=st_b0b1,
                    srt_simu_path=srt_simu_path,
                    bl_save_img=bl_save_img,
                    srt_results_root=srt_results_root,
                    srt_img_save_root=srt_img_save_root,
                    srt_img_aggregator_save_root=srt_img_aggregator_save_root,
                    srt_img_aggcsvsolu_save_root=srt_img_aggcsvsolu_save_root,
                    srt_mpcapc_baseroot=file.path(srt_root, snm_apcmpc_folder, '/'),
                    srt_imgcsv_mpcapc_root=srt_imgcsv_mpcapc_root,
                    srt_imgcsv_rev_root=srt_imgcsv_rev_root,
                    srt_imgcsv_perturb_root=srt_imgcsv_perturb_root,
                    srt_csv_path_root=srt_csv_path_root,
                    srt_folder=srt_folder,
                    srt_paper_main_textgraph=srt_paper_main_textgraph,
                    srt_paper_appendix_textgraph=srt_paper_appendix_textgraph,
                    spt_img_save=spt_img_save,
                    srt_csv_path=srt_csv_path,
                    ar_rho=ar_rho,
                    it_max_age=it_max_age)

  return(ls_output)
}
