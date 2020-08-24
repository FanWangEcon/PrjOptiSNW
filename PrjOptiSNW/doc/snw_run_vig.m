%% File Paths
git_repo = 'C:\Users\fan\Documents\Dropbox (UH-ECON)\PrjNygaardSorensenWang\PrjOptiSNW\PrjOptiSNW';
% git_repo = 'D:\Dropbox (UH-ECON)\PrjNygaardSorensenWang\PrjOptiSNW\PrjOptiSNW';
git_repo_vig_main = fullfile(git_repo, 'doc');

dropbox_htmlout = 'C:\Users\fan\Documents\Dropbox (UH-ECON)\PrjNygaardSorensenWang';
% dropbox_htmlout = 'D:\Dropbox (UH-ECON)\PrjNygaardSorensenWang';
dorpbox_htmlout_vig_main = 'MatlabVig\08-19-2020';
project = 'PrjOptiSNW\doc';

st_mlx_search_name = '*.mlx';
st_pub_format = 'html';
bl_run_mlx = true;
bl_run_mlx_only = false;
bl_verbose = true;

%% Run Vignettes
st_proj_folder = git_repo_vig_main;
cl_st_subfolder = {'calibrate', 'params', 'sdist', 'splannercheckval', ...
    'splannerjaeemk', 'splannerjmky', 'svalpol', 'svalpolsmall', 'svalpolunemploy'};
st_out_folder = fullfile(dropbox_htmlout, dorpbox_htmlout_vig_main);
ff_mlx2htmlpdf_runandexport(...
    st_proj_folder, cl_st_subfolder, ...
    st_mlx_search_name, st_out_folder, st_pub_format, ...
    bl_run_mlx, bl_run_mlx_only, ...
    bl_verbose);

