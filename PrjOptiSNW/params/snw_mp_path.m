%% SNW_MP_PATH Controls and Keeps Track of Paths
%    SNW_MP_PATH controls and keeps track of paths. Keeps track of external
%    path for storing larger data and simulation files. 
%
%    ST_COMPUTER, whose computer this is on: 'fan'
%
%    MP_PATHS = SNW_MP_PATH(ST_COMPUTER) returns paths 
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
         spt_dropbox_root = 'D:/Dropbox (UH-ECON)/PrjNygaardSorensenWang/';
     elseif (exist('C:/Users/fan/Documents/Dropbox (UH-ECON)/', 'dir')>0)
         spt_dropbox_root = 'C:/Users/fan/Documents/Dropbox (UH-ECON)/PrjNygaardSorensenWang/';
     end
end

%% Parametesr Grid Points
spt_simu_outputs = [spt_dropbox_root 'output/'];

%% Set Parameter Maps
mp_path_external = containers.Map('KeyType', 'char', 'ValueType', 'any');
mp_path_external('spt_simu_outputs') = spt_simu_outputs;

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
