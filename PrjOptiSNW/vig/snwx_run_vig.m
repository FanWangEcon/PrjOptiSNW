%% File Paths
clc
% git_repo = 'C:\Users\fan\PrjOptiSNW\PrjOptiSNW';
git_repo = 'C:\Users\fan\Documents\Dropbox (UH-ECON)\PrjNygaardSorensenWang\PrjOptiSNW\PrjOptiSNW';
git_repo_vig_main = fullfile(git_repo, 'vig');

dropbox_htmlout = 'C:\Users\fan\Documents\Dropbox (UH-ECON)\PrjNygaardSorensenWang\';
dorpbox_htmlout_vig_main = 'MatlabVig\07-07-2020';

% st_mlx_search_name = '*_densemore_*.mlx';
% st_mlx_search_name = 'snwx_v_planner_densemore.mlx';
st_mlx_search_name = 'snw_evuvw19_jmky*.mlx';
st_pub_format = 'html';
bl_run_mlx = true;
bl_run_mlx_only = false;
bl_verbose = true;

% cl_st_subfolder = {'sdist', 'splanner','svalpol'};
cl_st_subfolder = {'splanneralt'};

%% Run amin Vignette
st_proj_folder = git_repo_vig_main;
st_out_folder = fullfile(dropbox_htmlout, dorpbox_htmlout_vig_main);
ff_mlx2htmlpdf_runandexport(...
    st_proj_folder, cl_st_subfolder, ...
    st_mlx_search_name, st_out_folder, st_pub_format, ...
    bl_run_mlx, bl_run_mlx_only, ...
    bl_verbose);