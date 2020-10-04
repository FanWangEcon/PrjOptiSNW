# Threshold Allocation Graph, Figure 1 in Main Text
# This is to generate the figure stored with path below:
# "Paper\Main_text_graphs_and_tables\Actual_vs_threshold_C_allocation_by_income_kids_marital_status.png"

# Libraries
library(tidyverse)
library(REconTools)

# File Names, Paths Etc -------
ls_output <- fs_opti_support('b1_manna')
srt_csv_path_root <- ls_output$srt_csv_path_root
bl_save_img <- ls_output$bl_save_img
srt_paper_main_textgraph <- ls_output$srt_paper_main_textgraph

# G4 Figure in Table, Figure Number 4 Main Paper, Double both Kids and Adult Max
srt_main_folder <- ''
for (MPC_type in c(1,2)) {
  
  if (MPC_type == 1) {
    snm_sub_folder <- 'mpc_smooth_bykidsmarital20k_allchecks.csv'
    st_img_name_full <- 'mpc_smooth_bykidsmarital20k_allchecks.png'
    
  } else {
    snm_sub_folder <- 'mpc_raw_bykidsmarital20k_allchecks.csv'
    st_img_name_full <- 'mpc_raw_bykidsmarital20k_allchecks.png'
  }

  # Load in CSV
  spn_image_csv_path <- file.path(srt_csv_path_root, srt_main_folder, snm_sub_folder)
  df_MPC_results <- as_tibble(read.csv(spn_image_csv_path)) %>% select(-X)
  
  # Reshape Wide to Long
  df_MPC_graph <- df_MPC_results %>%
    rename_at(vars(num_range('', 0:243)), list(~paste0('', . , '')))
  df_MPC_graph <- df_MPC_graph %>% 
    pivot_longer(cols = starts_with('X'),
                 names_to = c('X'),
                 names_pattern = paste0("X(.*)"),
                 values_to = "MPC") %>%
    mutate(X = as.numeric(X)) %>%
    mutate(X = X + 1) %>% rename(checks = X)
  
  # Extract Income Group Info
  df_MPC_graph <- df_MPC_graph %>% rowwise() %>%
    mutate(y_lower = 58056*as.numeric(substring(strsplit(ymin_group, ",")[[1]][1], 2))) %>% 
    mutate(y_upper = 58056*as.numeric(gsub(strsplit(ymin_group, ",")[[1]][2],  pattern = "]", replacement = ""))) %>% 
    mutate(y_lower = round(y_lower/10000)*10000, y_upper = round(y_upper/10000)*10000) %>% 
    ungroup %>%
    mutate(y_group = paste(y_lower, ' to ', y_upper)) %>%
    mutate(y_group = case_when(y_group == "0  to  20000" ~ "0 to 20k",
                               y_group == "20000  to  40000" ~ "20k to 40k",
                               y_group == "40000  to  60000" ~ "40k to 60k",
                               y_group == "60000  to  80000" ~ "60k to 80k",
                               y_group == "80000  to  1e+05" ~ "80k to 100k",
                               TRUE ~ ">100k")) %>%
    filter(y_group != ">100k")
  
  # Marital and Kids Level Labeling
  marry_levels <- c(Single = "0", Married = "1")
  kids_levels <- c("no children" = "0", "one child" = "1",
                   "two children" = "2", "three children" = "3",
                   "four children" = "4")
  
  ## Select, as factor, and recode
  df_MPC_graph <- df_MPC_graph %>%
    mutate(kids = as.factor(kids),
           marital = as.factor(marital)) %>%
    mutate(kids = fct_recode(kids, !!!kids_levels),
           marital = fct_recode(marital, !!!marry_levels))
  
  # Get value for minimum and maximum income levels for each bin
  df_MPC_graph <- df_MPC_graph %>%
    mutate(kids = as.factor(kids),
           marital = as.factor(marital))
  
  # GGPLOT grahing ------
  # GGPLOT AES and FACET_WRAP
  plt_cur <- df_MPC_graph %>%
    ggplot(aes(x=checks, y=MPC,
               colour=y_group,
               shape=y_group)) +
    facet_wrap(~ marital + kids,
               ncol=5,
               scales = "free_x",
               labeller = label_wrap_gen(multi_line=FALSE))
  plt_cur <- plt_cur + geom_line() + geom_point(size=2)
  
  
  # GGPLOT legendings etc -------
  # X and Y labels and Legends etc.
  x.labels <- c('0k', '5k', '10k', '20k', '25K')
  x.breaks <- c(0, 50, 100, 200, 250)
  
  # x-axis labeling
  plt_cur <- plt_cur +
    scale_x_continuous(labels = x.labels, breaks = x.breaks)
  # limits = c(0, 240000/58056)
  
  # # Legend Labeling
  # ar_st_age_group_leg_labels <- c("Actual", "Optimal")
  # # color #F8766D is default red
  # # color #F8766D is default bluish color
  # ar_st_colours <- c("#F8766D", "#33cc33")
  
  plt_cur <- plt_cur +
    theme(
      text = element_text(size = 12),
      legend.title = element_blank(),
      legend.position = c(0.14, 0.9),
      legend.background = element_rect(fill = "white", colour = "black", linetype='solid'))
  
  # +
  #   scale_colour_manual(values=ar_st_colours, labels=ar_st_age_group_leg_labels) +
  #   scale_shape_discrete(labels=ar_st_age_group_leg_labels)
  
  # X and Y Titling
  stg_x_title_hhinc <- 'Check Dollar Amounts'
  stg_y_title_checks <- 'MPC of the Next $100 Check'
  plt_cur <- plt_cur + labs(x = stg_x_title_hhinc, y = stg_y_title_checks)
  
  # Save Image Outputs -----
  # PNG and print
  if (bl_save_img) {
    png(paste0(srt_paper_main_textgraph, st_img_name_full),
        width = 270,
        height = 216, units='mm',
        res = 300, pointsize=7)
  }
  print(plt_cur)
  if (bl_save_img) {
    dev.off()
  }
  
}