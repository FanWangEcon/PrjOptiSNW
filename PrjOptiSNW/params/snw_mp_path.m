%% SNW_MP_PATH Controls and Keeps Track of Paths
%    SNW_MP_PATH controls and keeps track of paths. Keeps track of external
%    path for storing larger data and simulation files. 
%
%    ST_COMPUTER, whose computer this is on: 'fan'
%
%    mp_paths = snw_mp_path('fan');
%    spt_simu_codem = mp_paths('spt_simu_codem');
%    spt_simu_outputs = mp_paths('spt_simu_outputs');
%    spt_simu_outputs_log = mp_paths('spt_simu_outputs_log');
%    spt_simu_results_csv = mp_paths('spt_simu_results_csv');
%

%%
function varargout = snw_mp_path(varargin)
%% Parse Main Inputs and Set Defaults
if (~isempty(varargin))

    if (length(varargin)==1)
        st_computer = varargin{:};
        bl_print_mp_path = false;
    elseif (length(varargin)==2)
        [st_computer, bl_print_mp_path] = varargin{:};
    end    

else

    st_computer = 'fan';
    bl_print_mp_path = true;
    
end

%% Dropbox Root Auto Detect by Computer
if (strcmp(st_computer, 'fan'))
     if (exist('D:/Dropbox (UH-ECON)', 'dir')>0)
         spt_dropbox_root = fullfile('D:/Dropbox (UH-ECON)/PrjNygaardSorensenWang/');
     elseif (exist('C:/Users/fan/Documents/Dropbox (UH-ECON)/', 'dir')>0)
         spt_dropbox_root = fullfile('C:/Users/fan/Documents/Dropbox (UH-ECON)/PrjNygaardSorensenWang/');
     end
end

%% Parametesr Grid Points
% store output mat files, this is to save speed during check calculation,
% re-calculation

spt_simu_codem = fullfile(spt_dropbox_root, 'PrjOptiSNW', filesep);
spt_simu_codem_doc = fullfile(spt_dropbox_root, 'PrjOptiSNW', 'PrjOptiSNW', 'doc', filesep);
spt_simu_outputs = fullfile(spt_dropbox_root, 'Output202102', filesep);
spt_simu_outputs_vig = fullfile(spt_dropbox_root, 'MatlabVig', '02-21-2021', filesep);
spt_simu_outputs_log = fullfile(spt_simu_outputs, 'log', filesep);
spt_simu_outputs_mat = fullfile(spt_simu_outputs, 'mat', filesep);
% spt_simu_results_csv = fullfile(spt_dropbox_root, 'Results', '2020-08-23-b1_manna', 'csv_with_expmin_apc');
% spt_simu_results_csv = fullfile(spt_dropbox_root, 'Results', '2020-08-23-b1_manna', 'csv', filesep);
% spt_simu_results_csv = fullfile(spt_dropbox_root, 'Results', '2021-01-02-b1_manna', 'csv', filesep);
spt_simu_results_csv = fullfile(spt_dropbox_root, 'Results', '2021-02-21-b1_manna', 'csv', filesep);

%% Set Parameter Maps
mp_path_external = containers.Map('KeyType', 'char', 'ValueType', 'any');
mp_path_external('spt_dropbox_root') = spt_dropbox_root;
mp_path_external('spt_simu_codem') = spt_simu_codem;
mp_path_external('spt_simu_codem_doc') = spt_simu_codem_doc;

mp_path_external('spt_simu_outputs') = spt_simu_outputs;
mp_path_external('spt_simu_outputs_vig') = spt_simu_outputs_vig;
mp_path_external('spt_simu_outputs_log') = spt_simu_outputs_log;
mp_path_external('spt_simu_outputs_mat') = spt_simu_outputs_mat;
mp_path_external('spt_simu_results_csv') = spt_simu_results_csv;

%% Combine Maps
mp_path = [mp_path_external];
mp_path('st_computer') = string(st_computer);

% MP_PARAMS = [MP_PARAMS_PREFTECHPRICE; MP_PARAMS_STATESGRID ; ...
%     MP_PARAMS_EXOTRANS; MP_PARAMS_TYPELIFE; MP_PARAMS_INTLEN];

%% Print
if (bl_print_mp_path)
    ff_container_map_display(mp_path_external);
end

%% Return
if (nargout==1)
    varargout = cell(nargout,0);
    varargout{1} = mp_path;
elseif (nargout==2)
    varargout = cell(nargout,0);
    varargout{1} = mp_path;
    varargout{2} = mp_path_external;
end

end
