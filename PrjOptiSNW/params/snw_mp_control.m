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
bl_print_a4chk = true;
bl_print_vu_vw = true;
bl_print_v_planner = true;
bl_print_precompute = true;
bl_print_evuvw20_jaeemk = true;
bl_print_evuvw19_jaeemk = true;
bl_print_evuvw19_jymk = true;
bl_print_evuvw19_jmky_mass = true;
bl_print_evuvw19_jmky_allchecks = true;

if (strcmp(st_param_group, 'default_test')) 
    bl_print_ds_verbose = true;
    bl_print_vfi_verbose = true;
    bl_print_a4chk_verbose = true;
    bl_print_vu_vw_verbose = true;
    bl_print_v_planner_verbose = true;    
    bl_print_precompute_verbose = true;
    bl_print_evuvw20_jaeemk_verbose = true;
    bl_print_evuvw19_jaeemk_verbose = true;
    bl_print_evuvw19_jmky_verbose = true;
    bl_print_evuvw19_jmky_mass_verbose = true;
    bl_print_evuvw19_jmky_allchecks_verbose = true;
else
    bl_print_ds_verbose = false;
    bl_print_vfi_verbose = false;
    bl_print_a4chk_verbose = false;
    bl_print_vu_vw_verbose = false;
    bl_print_v_planner_verbose = false;
    bl_print_precompute_verbose = false;
    bl_print_evuvw20_jaeemk_verbose = false;
    bl_print_evuvw19_jaeemk_verbose = false;
    bl_print_evuvw19_jmky_verbose = false;
    bl_print_evuvw19_jmky_mass_verbose = false;
    bl_print_evuvw19_jmky_allchecks_verbose = false;
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

mp_m4check_controls = containers.Map('KeyType', 'char', 'ValueType', 'any');
mp_m4check_controls('fl_max_trchk_perc_increase') = 1.50;

mp_calibrate = containers.Map('KeyType', 'char', 'ValueType', 'any');
mp_calibrate('err') = err;
mp_calibrate('tol') = tol;

%% What to Print Out, Store, Etc
mp_compute_stats = containers.Map('KeyType', 'char', 'ValueType', 'any');
mp_compute_stats('bl_compute_drv_stats') = bl_compute_drv_stats;

mp_profile = containers.Map('KeyType', 'char', 'ValueType', 'any');
mp_profile('bl_timer') = bl_timer;

mp_display = containers.Map('KeyType', 'char', 'ValueType', 'any');
mp_display('bl_print_vfi') = bl_print_vfi;
mp_display('bl_print_vfi_verbose') = bl_print_vfi_verbose;

mp_display('bl_print_ds') = bl_print_ds;
mp_display('bl_print_ds_verbose') = bl_print_ds_verbose;

mp_display('bl_print_a4chk') = bl_print_a4chk;
mp_display('bl_print_a4chk_verbose') = bl_print_a4chk_verbose;

mp_display('bl_print_vu_vw') = bl_print_vu_vw;
mp_display('bl_print_vu_vw_verbose') = bl_print_vu_vw_verbose;

mp_display('bl_print_v_planner') = bl_print_v_planner;
mp_display('bl_print_v_planner_verbose') = bl_print_v_planner_verbose;
 
mp_display('bl_print_precompute') = bl_print_precompute;
mp_display('bl_print_precompute_verbose') = bl_print_precompute_verbose;

mp_display('bl_print_evuvw20_jaeemk') = bl_print_evuvw20_jaeemk;
mp_display('bl_print_evuvw20_jaeemk_verbose') = bl_print_evuvw20_jaeemk_verbose;

mp_display('bl_print_evuvw19_jaeemk') = bl_print_evuvw19_jaeemk;
mp_display('bl_print_evuvw19_jaeemk_verbose') = bl_print_evuvw19_jaeemk_verbose;

mp_display('bl_print_evuvw19_jmky') = bl_print_evuvw19_jymk;
mp_display('bl_print_evuvw19_jmky_verbose') = bl_print_evuvw19_jmky_verbose;

mp_display('bl_print_evuvw19_jmky_mass') = bl_print_evuvw19_jmky_mass;
mp_display('bl_print_evuvw19_jmky_mass_verbose') = bl_print_evuvw19_jmky_mass_verbose;

mp_display('bl_print_evuvw19_jmky_allchecks') = bl_print_evuvw19_jmky_allchecks;
mp_display('bl_print_evuvw19_jmky_allchecks_verbose') = bl_print_evuvw19_jmky_allchecks_verbose;

%% Combine Maps
mp_controls = [mp_minimizer_controls; mp_m4check_controls; ...
    mp_calibrate; mp_compute_stats; mp_profile; mp_display];
mp_controls('mp_params_name') = string(st_param_group);

%% Print 
if (bl_print_mp_controls)
    ff_container_map_display(mp_controls);
end

%% Return
varargout = cell(nargout,0);
for it_k = 1:nargout
    if (it_k==1)
        ob_out_cur = mp_controls;
    elseif (it_k==2)
        ob_out_cur = mp_minimizer_controls;
    elseif (it_k==3)
        ob_out_cur = mp_m4check_controls;
    elseif (it_k==4)
        ob_out_cur = mp_calibrate;
    elseif (it_k==5)
        ob_out_cur = mp_compute_stats;        
    elseif (it_k==6)
        ob_out_cur = mp_profile;        
    elseif (it_k==7)
        ob_out_cur = mp_display;        
    elseif (it_k==8)
        ob_out_cur = mp_store;        
    end
    varargout{it_k} = ob_out_cur;
end


end
