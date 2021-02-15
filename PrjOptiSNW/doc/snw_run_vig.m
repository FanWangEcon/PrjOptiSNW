%% File Paths
mp_path = snw_mp_path('fan');
git_repo_vig_main = mp_path('spt_simu_codem_doc');

st_mlx_search_name = '*.mlx';
st_pub_format = 'html';
bl_run_mlx = true;
bl_run_mlx_only = false;
bl_verbose = true;

%% Run Vignettes
st_proj_folder = git_repo_vig_main;
% cl_st_subfolder = {'calibrate', 'params', 'sdist', 'splannercheckval', ...
%     'splannerjaeemk', 'splannerjmky', 'svalpol', 'svalpolsmall', 'svalpolunemploy'};
st_out_folder = mp_path('spt_simu_outputs_vig');
cl_st_subfolder = {'params', 'svalpol', 'svalpolsmall', 'svalpolunemploy', 'sdist'};
cl_st_subfolder = {'svalpolsmall', 'svalpolunemploy', 'sdist'};
cl_st_subfolder = {'splannercheckval', 'splannerjaeemk', 'splannerjmky', 'calibrate'};
cl_st_subfolder = {'svalpolsmall'};
ff_mlx2htmlpdf_runandexport(...
    st_proj_folder, cl_st_subfolder, ...
    st_mlx_search_name, st_out_folder, st_pub_format, ...
    bl_run_mlx, bl_run_mlx_only, ...
    bl_verbose);