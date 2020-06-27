%% SNW_MP_CONTROL Organizes and Sets Various Solution Simu Control Parameters
%    SNW_MP_CONTROL opitmizer control, graph, print, and other controls
%
%    MP_CONTROLS = SNW_MP_CONTROL() get default parameters all in the
%    same container map
%
%    MP_CONTROLS = SNW_MP_CONTROL(ST_PARAM_GROUP)
%    generates default parameters for the type ST_PARAM_GROUP. 
%
%    MP_CONTROLS = SNW_MP_CONTROL(ST_PARAM_GROUP, bl_print_mp_controls) generates
%    default parameters for the type ST_PARAM_GROUP, display parameter map
%    details if bl_print_mp_controls is true.
%    
%    See also SNWX_MP_CONTROLS
%

%%
function varargout = snw_mp_control(varargin)
%% Parse Main Inputs and Set Defaults
if (~isempty(varargin))
    
    if (length(varargin)==1)
        st_param_group = varargin{:};
        bl_print_mp_controls = false;        
    elseif (length(varargin)==2)
        [st_param_group, bl_print_mp_controls] = varargin{:};
    end
    
else
    
    st_param_group = 'default_base';   
    st_param_group = 'default_test';   
    bl_print_mp_controls = true;
    
end

%% Control Optimization 
%amin=0;
%amax=agrid(end);
A_aux=[];
B_aux=[];
Aeq=[];
Beq=[];

nonlcon=[];
options=optimoptions('fmincon','Display', 'off');
options2=optimoptions('fsolve','Display','off');

%% Control Calibration 
err=1;
tol=0.005;

%% Control What to Calculate 
bl_compute_drv_stats = true;

%% Control Profiling and Display
bl_timer = true;

%% Controls Print
bl_print_vfi = true;
bl_print_ds = true;

if (strcmp(st_param_group, 'default_test')) 
    bl_print_ds_verbose = true;
    bl_print_vfi_verbose = true;
else
    bl_print_ds_verbose = false;
    bl_print_vfi_verbose = false;
end

%% Control Storage
if (strcmp(st_param_group, 'default_test')) 
    bl_ds_store_all = true;
    bl_vfi_store_all = true;
else
    bl_ds_store_all = false;
    bl_vfi_store_all = false;
end

%% Control Optimization 
mp_minimizer_controls = containers.Map('KeyType', 'char', 'ValueType', 'any');
mp_minimizer_controls('A_aux') = A_aux;
mp_minimizer_controls('B_aux') = B_aux;
mp_minimizer_controls('Aeq') = Aeq;
mp_minimizer_controls('Beq') = Beq;
mp_minimizer_controls('nonlcon') = nonlcon;
mp_minimizer_controls('options') = options;
mp_minimizer_controls('options2') = options2;

mp_calibrate = containers.Map('KeyType', 'char', 'ValueType', 'any');
mp_calibrate('err') = err;
mp_calibrate('tol') = tol;

mp_compute_stats = containers.Map('KeyType', 'char', 'ValueType', 'any');
mp_compute_stats('bl_compute_drv_stats') = bl_compute_drv_stats;

mp_profile = containers.Map('KeyType', 'char', 'ValueType', 'any');
mp_profile('bl_timer') = bl_timer;

mp_display = containers.Map('KeyType', 'char', 'ValueType', 'any');
mp_display('bl_print_vfi') = bl_print_vfi;
mp_display('bl_print_ds') = bl_print_vfi;
mp_display('bl_print_ds_verbose') = bl_print_ds_verbose;
mp_display('bl_print_vfi_verbose') = bl_print_vfi_verbose;

mp_store = containers.Map('KeyType', 'char', 'ValueType', 'any');
mp_store('bl_ds_store_all') = bl_ds_store_all;
mp_store('bl_vfi_store_all') = bl_vfi_store_all;

%% Combine Maps
mp_controls = [mp_minimizer_controls; ...
    mp_calibrate; mp_compute_stats; mp_profile; mp_display; mp_store];
mp_controls('mp_params_name') = string(st_param_group);

%% Print 
if (bl_print_mp_controls)
    ff_container_map_display(mp_controls);
end

%% Return
if (nargout==1)
    varargout = cell(nargout,0);
    varargout{1} = mp_controls;
elseif (nargout==3)
    varargout = cell(nargout,0);
    varargout{1} = mp_controls;
    varargout{2} = mp_minimizer_controls;
    varargout{3} = mp_profile;
    varargout{4} = mp_display;    
end

end
