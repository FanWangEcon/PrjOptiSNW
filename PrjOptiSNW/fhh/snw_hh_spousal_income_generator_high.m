% see file:///G:/repos/M4Econ/function/anonymous/htmlpdfm/fs_dyna_generate_func.html
cl_code = {'function F=snw_hh_spousal_income(j,educ,kids,earn,SS_inc,jret)', ...
     '% both', ...
     '%it_edu_grid_type = 1;', ...
     '% non-college only', ...
     '%it_edu_grid_type = 2;', ...
     '% college only', ...
     'it_edu_grid_type = 3;', ...
     'if (it_edu_grid_type == 1 || it_edu_grid_type == 2)', ...
     '    % no changes needed when both, or low education type', ...
     '    % low type = 1, educ = 1', ...
     '    F = snw_hh_spousal_income_func(j,educ,kids,earn,SS_inc,jret);', ...
     'elseif (it_edu_grid_type == 3)', ...
     '    % problem high type, educ should equal 2, but solving for high edu only, educ=1', ...
     '    F = snw_hh_spousal_income_func(j,2,kids,earn,SS_inc,jret);', ...
     'end', ...
     'end'};

% Check folder to use
spn_current_mlx_path = matlab.desktop.editor.getActiveFilename;
[spt_filepath, snm_name, st_ext] = fileparts(spn_current_mlx_path);
[spt_filepath] = fileparts(spt_filepath);
% disp(['spt_filepath=' spt_filepath]);

spt_filepath = fullfile(spt_filepath, 'fhh');
if ~exist(spt_filepath, 'dir')
    mkdir(spt_filepath);
end
% disp(['spt_filepath:' spt_filepath]);

spn_path_file = fullfile(spt_filepath, 'snw_hh_spousal_income.m');
f = fopen(spn_path_file, 'w');
fprintf(f, '%s\n', cl_code{:});
fclose(f);

addpath(spt_filepath);
