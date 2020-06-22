%% SNW_MP_CONTROLS Organizes and Sets Various Solution Simu Control Parameters
%    SNW_MP_CONTROLS opitmizer control, graph, print, and other controls
%
%    MP_CONTROLS = SNW_MP_CONTROLS() get default parameters all in the
%    same container map
%
%    MP_CONTROLS = SNW_MP_CONTROLS(ST_PARAM_GROUP)
%    generates default parameters for the type ST_PARAM_GROUP. 
%
%    MP_CONTROLS = SNW_MP_CONTROLS(ST_PARAM_GROUP, bl_print_mp_controls) generates
%    default parameters for the type ST_PARAM_GROUP, display parameter map
%    details if bl_print_mp_controls is true.
%    
%    See also SNWX_MP_CONTROLS
%

%%
function varargout = snw_mp_controls(varargin)
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
    bl_print_mp_controls = true;
    
end

%% For minimization problem
%amin=0;
%amax=agrid(end);
A_aux=[];
B_aux=[];
Aeq=[];
Beq=[];

nonlcon=[];
options=optimoptions('fmincon','Display', 'off');
options2=optimoptions('fsolve','Display','off');

%% Timer and Profiler
bl_timer = true;

%% Print Controls
bl_print_iter = true;

%% Set Parameter Maps
mp_minimizer_controls = containers.Map('KeyType', 'char', 'ValueType', 'any');
mp_minimizer_controls('A_aux') = A_aux;
mp_minimizer_controls('B_aux') = B_aux;
mp_minimizer_controls('Aeq') = Aeq;
mp_minimizer_controls('Beq') = Beq;
mp_minimizer_controls('nonlcon') = nonlcon;
mp_minimizer_controls('options') = options;
mp_minimizer_controls('options2') = options2;

mp_profile = containers.Map('KeyType', 'char', 'ValueType', 'any');
mp_profile('bl_timer') = bl_timer;

mp_display = containers.Map('KeyType', 'char', 'ValueType', 'any');
mp_display('bl_print_iter') = bl_print_iter;

%% Combine Maps
mp_controls = [mp_minimizer_controls; mp_profile; mp_display];

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
