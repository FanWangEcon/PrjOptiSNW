%% File Paths

git_repo = 'C:\Users\fan\PrjOptiSNW\PrjOptiSNW';
git_repo_vig_main = fullfile(git_repo, 'vig');

dropbox_htmlout = 'C:\Users\fan\Documents\Dropbox (UH-ECON)\PrjNygaardSorensenWang\';
dorpbox_htmlout_vig_main = 'MatlabVig\07-07-2020';

% st_mlx_search_name = '*_densemore_*.mlx';
st_mlx_search_name = 'snwx_v_planner_densemore.mlx';
st_pub_format = 'html';
bl_run_mlx = true;
bl_run_mlx_only = false;
bl_verbose = true;

%% Run amin Vignette
% st_proj_folder = git_repo_vig_main;
% cl_st_subfolder = {'sdist'};
% st_out_folder = fullfile(dropbox_htmlout, dorpbox_htmlout_vig_main);
% ff_mlx2htmlpdf_runandexport(...
%     st_proj_folder, cl_st_subfolder, ...
%     st_mlx_search_name, st_out_folder, st_pub_format, ...
%     bl_run_mlx, bl_run_mlx_only, ...
%     bl_verbose);

%% Run ds_savingsunkc Vignette
st_proj_folder = git_repo_vig_main;
cl_st_subfolder = {'splanner'};
st_out_folder = fullfile(dropbox_htmlout, dorpbox_htmlout_vig_main);
ff_mlx2htmlpdf_runandexport(...
    st_proj_folder, cl_st_subfolder, ...
    st_mlx_search_name, st_out_folder, st_pub_format, ...
    bl_run_mlx, bl_run_mlx_only, ...
    bl_verbose);

%% Run savingfc Vignette
% st_proj_folder = git_repo_vig_main;
% cl_st_subfolder = {'svalpol'};
% st_out_folder = fullfile(dropbox_htmlout, dorpbox_htmlout_vig_main);
% ff_mlx2htmlpdf_runandexport(...
%     st_proj_folder, cl_st_subfolder, ...
%     st_mlx_search_name, st_out_folder, st_pub_format, ...
%     bl_run_mlx, bl_run_mlx_only, ...
%     bl_verbose);