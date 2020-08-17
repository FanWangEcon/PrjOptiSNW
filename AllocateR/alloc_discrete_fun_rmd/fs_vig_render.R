rm(list = ls())
knitr::opts_chunk$set(echo = TRUE)
knitr::opts_chunk$set(fig.width=12, fig.height=8)
library(tidyverse)
library(REconTools)
library(PrjOptiAlloc)

srt_root <- 'D:/Dropbox (UH-ECON)/PrjNygaardSorensenWang/'
# srt_root <- 'C:/Users/fan/Documents/Dropbox (UH-ECON)/PrjNygaardSorensenWang/'

source(paste0(srt_root, 'PrjOptiSNW/AllocateR/alloc_discrete_fun/fs_opti_support.R'))

ar_spt_root <- c(paste0(srt_root, 'PrjOptiSNW/AllocateR/alloc_discrete_fun/'))
ar_spn_skip <- c('_main', 'mpc', 'g47')
# file directory for HTML output
# Update Directory for Image output inside each of the RMd files.
# spt_out_directory <- 'C:/Users/fan/Documents/Dropbox (UH-ECON)/PrjNygaardSorensenWang/Results/2020-07-29/HTML_files/'
spt_out_directory <- paste0(srt_root, 'Results/2020-08-05/HTML_files/',
                            'b0_a64/')
st_save_add_suffix <- ''
ls_bool_convert <- list(bl_pdf=FALSE, bl_html=TRUE, bl_R=FALSE)
bl_verbose <- TRUE
svr_checks <- 'checks'
REconTools::ff_sup_rmd2htmlpdfr(ar_spt_root=ar_spt_root, ar_spn_skip=ar_spn_skip,
                                spt_out_directory=spt_out_directory,
                                st_save_add_suffix=st_save_add_suffix,
                                ls_bool_convert=ls_bool_convert,
                                bl_verbose=bl_verbose)
