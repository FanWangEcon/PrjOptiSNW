%% SNW_A4CHK_WRK_BISEC (loop) solves for Asset Position Corresponding to Check Level
%    What is the value of a check? From the perspective of the value
%    function? We have Asset as a state variable, in a cash-on-hand sense,
%    how much must the asset (or think cash-on-hand) increase by, so that
%    it is equivalent to providing the household with a check? This is not
%    the same as the check amount because of tax as well as interest rates.
%    Interest rates means that you might need to offer a smaller a than the
%    check amount. The tax rate means that we might need to shift a by
%    larger than the check amount.
%
%    This function mainly exists to make sure to the bisection solution
%    returns the same answer as the fzero solution. Given that these are
%    the same, the vectorized version of this code should be used:
%    SNW_A4CHK_WRK_BISEC_VEC. Note that here, FF_OPTIM_BISEC_SAVEZRONE or
%    fzero speed differences is not substantial it seems.
%
%    This file solves the problem using bisection with the function
%    FF_OPTIM_BISEC_SAVEZRONE from MEconTools. Also solves via fsolve. Can
%    solve using both, and check if the answers are the same. This function
%    requires a common zero and 1 starting point. For the problem here, 0
%    bound is for a change for check to be zero, the 1 bound is for a
%    change for check to be 130 percent of the check value, which considers
%    the possible upper bound of taxes.
%
%    * WELF_CHECKS integer the number of checks
%    * TR float the value of each check
%    * V_SS ndarray the value matrix along standard state-space dimensions:
%    (n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid)
%    * MP_PARAMS map with model parameters
%    * MP_CONTROLS map with control parameters
%
%    % bl_fzero is true if solve via fzero
%    mp_controls('bl_fzero') = false;
%    $ bl_ff_bisec is true is solve via ff_optim_bisec_savezrone, if both
%    true solve both and return error message if answers do not match up.
%    mp_controls('bl_ff_bisec') = true;
%
%    [V_W, EXITFLAG_FSOLVE] = SNW_A4CHK_WRK_BISEC(WELF_CHECKS, TR,
%    V_SS, MP_PARAMS, MP_CONTROLS) solves for working value given V_SS
%    value function results, for number of check WELF_CHECKS, and given the
%    value of each check equal to TR.
%
%    See also SNW_A4CHK_WRK_BISEC_VEC, SNW_A4CHK_WRK, FIND_A_WORKING
%

%%
function [V_W, exitflag_fsolve]=snw_a4chk_wrk_bisec(varargin)

%% Default and Parse
if (~isempty(varargin))
    
    if (length(varargin)==3)
        [welf_checks, TR, V_ss] = varargin{:};
        mp_controls_ext = snw_mp_control('default_base');
    elseif (length(varargin)==5)
        [welf_checks, TR, V_ss, mp_params, mp_controls_ext] = varargin{:};
    end
    
else
    close all;
    
    % Solve the VFI Problem and get Value Function
    mp_params = snw_mp_param('default_tiny');
    mp_controls_ext = snw_mp_control('default_test');
    [V_ss,~,~,~] = snw_vfi_main_bisec_vec(mp_params, mp_controls_ext);
    welf_checks = 2;
    TR = 100/58056;
    
    % run fzero
    mp_controls_ext('bl_fzero') = true;    
    % run ff_optim_bisec_savezrone bisect as well to compare results
    mp_controls_ext('bl_ff_bisec') = true;
    
end

%% Reset All globals
% globals = who('global');
% clear(globals{:});
% Parameters used in this code directly
global agrid n_jgrid n_agrid n_etagrid n_educgrid n_marriedgrid n_kidsgrid
% Used in find_a_working function
global theta r agrid epsilon eta_H_grid eta_S_grid SS Bequests bequests_option throw_in_ocean

%% Parse Model Parameters
params_group = values(mp_params, {'theta', 'r'});
[theta,  r] = params_group{:};

params_group = values(mp_params, {'Bequests', 'bequests_option', 'throw_in_ocean'});
[Bequests, bequests_option, throw_in_ocean] = params_group{:};

params_group = values(mp_params, {'agrid', 'eta_H_grid', 'eta_S_grid'});
[agrid, eta_H_grid, eta_S_grid] = params_group{:};

params_group = values(mp_params, {'epsilon', 'SS'});
[epsilon, SS] = params_group{:};

params_group = values(mp_params, ...
    {'n_jgrid', 'n_agrid', 'n_etagrid', 'n_educgrid', 'n_marriedgrid', 'n_kidsgrid'});
[n_jgrid, n_agrid, n_etagrid, n_educgrid, n_marriedgrid, n_kidsgrid] = params_group{:};

%% Control Map Function Specific Local Defaults
mp_controls = containers.Map('KeyType', 'char', 'ValueType', 'any');
mp_controls('bl_fzero') = false;
mp_controls('bl_ff_bisec') = true;

if (length(varargin)>=2 || isempty(varargin))
    mp_controls = [mp_controls; mp_controls_ext];
end

%% Parse Model Controls
% Which to Run Controls
params_group = values(mp_controls, {'bl_fzero', 'bl_ff_bisec'});
[bl_fzero, bl_ff_bisec] = params_group{:};

% Minimizer Controls
params_group = values(mp_controls, {'fl_max_trchk_perc_increase', 'options2'});
[fl_max_trchk_perc_increase, options2] = params_group{:};

% Profiling Controls
params_group = values(mp_controls, {'bl_timer'});
[bl_timer] = params_group{:};

% Display Controls
params_group = values(mp_controls, {'bl_print_a4chk','bl_print_a4chk_verbose'});
[bl_print_a4chk, bl_print_a4chk_verbose] = params_group{:};

%% Timing and Profiling Start

if (bl_timer)
    tic
end

%% Loop over states and find A state for a Particular Check Level

V_W=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);
exitflag_fsolve=NaN(n_jgrid,n_agrid,n_etagrid,n_educgrid,n_marriedgrid,n_kidsgrid);

for j=1:n_jgrid % Age
    for a=1:n_agrid % Assets
        for eta=1:n_etagrid % Productivity
            for educ=1:n_educgrid % Educational level
                for married=1:n_marriedgrid % Marital status
                    for kids=1:n_kidsgrid % Number of kids
                        
                        % A. solve using ff_optim_bisec_savezrone
                        if (bl_ff_bisec)
                            % A1. Construct Bisect Percentage Function Handle 
                            fc_ffi_frac0t1_find_a_working = @(x) ffi_frac0t1_find_a_working(...
                                x, agrid(a), ...
                                j,a,eta,educ,married,kids,...
                                welf_checks, TR, fl_max_trchk_perc_increase);

                            % A2. Solve via Bisection
                            [a_aux_bisec_frac, a_aux_bisec] = ...
                                ff_optim_bisec_savezrone(fc_ffi_frac0t1_find_a_working);
                        end
                        
                        % B. Solve via fzero
                        if (bl_fzero)
                            x0=agrid(a)+TR*welf_checks; % Initial guess for a
                            
                            [a_aux_fzero,~,it_flag_fzero]=...
                                fsolve(@(x)find_a_working(x,j,a,eta,educ,married,kids,TR,welf_checks),x0,options2);
                            exitflag_fsolve(j,a,eta,educ,married,kids) = it_flag_fzero;
                        end
                        
                        % C. Check if Solution Answers are the Same
                        if (bl_ff_bisec && bl_fzero)
                            if (abs(a_aux_fzero - a_aux_bisec) >= 10e-5)
                                st_message = ['SNW_A4CHK_WRK_BISEC: a_aux=' num2str(a_aux_fzero) ...
                                    ', but, a_aux_bisec=' num2str(a_aux_bisec)];
                                error(st_message)
                            end
                            a_aux = a_aux_bisec;
                        elseif (bl_ff_bisec)
                            a_aux = a_aux_bisec;
                        else
                            a_aux = a_aux_fzero;
                        end
                        
                        % D. Error Check
                        if a_aux<0
                            disp(a_aux)
                            error('Check code! Should not allow for negative welfare checks')
                        elseif a_aux>agrid(n_agrid)
                            a_aux=agrid(n_agrid);
                        end
                        
                        % Linear interpolation
                        ind_aux=find(agrid<=a_aux,1,'last');
                        
                        if a_aux==0
                            inds(1)=1;
                            inds(2)=1;
                            vals(1)=1;
                            vals(2)=0;
                            
                        elseif a_aux==agrid(n_agrid)
                            inds(1)=n_agrid;
                            inds(2)=n_agrid;
                            vals(1)=1;
                            vals(2)=0;
                            
                        else
                            inds(1)=ind_aux;
                            inds(2)=ind_aux+1;
                            vals(1)=1-((a_aux-agrid(inds(1)))/(agrid(inds(2))-agrid(inds(1))));
                            vals(2)=1-vals(1);
                            
                        end
                        
                        V_W(j,a,eta,educ,married,kids)=vals(1)*V_ss(j,inds(1),eta,educ,married,kids)+vals(2)*V_ss(j,inds(2),eta,educ,married,kids);
                        
                    end
                end
            end
        end
    end
    
    if (bl_print_a4chk)
        disp(strcat(['SNW_A4CHK_WRK: Finished Age Group:' num2str(j) ' of ' num2str(n_jgrid)]));
    end
    
end

%% Timing and Profiling End
if (bl_timer)
    toc;
    st_complete_a4chk = strjoin(...
        ["Completed SNW_A4CHK_WRK", ...
         ['welf_checks=' num2str(welf_checks)], ...
         ['TR=' num2str(TR)], ...
         ['SNW_MP_PARAM=' char(mp_params('mp_params_name'))], ...
         ['SNW_MP_CONTROL=' char(mp_controls('mp_params_name'))] ...
        ], ";");
    disp(st_complete_a4chk);
end


%% Compare Difference between V_ss and V_W

if (bl_print_a4chk_verbose)
    mn_V_gain_check = V_W - V_ss;
    mn_V_gain_frac_check = (V_W - V_ss)./V_ss;
    mp_container_map = containers.Map('KeyType','char', 'ValueType','any');
    mp_container_map('V_W') = V_W;
    mp_container_map('V_ss') = V_ss;
    mp_container_map('V_W_minus_V_ss') = mn_V_gain_check;
    mp_container_map('V_W_minus_V_ss_divide_V_ss') = mn_V_gain_frac_check;
    ff_container_map_display(mp_container_map);
end

end

function [ar_root_zero, ar_a_aux_amz] = ...
    ffi_frac0t1_find_a_working(...
    ar_aux_change_frac_amz, ar_a_state_level_amz, ...
    j_amz,a_amz,eta_amz,educ_amz,married_amz,kids_amz, ...
    welf_checks, TR, fl_max_trchk_perc_increase)
    
    fl_a_aux_max = TR*welf_checks*fl_max_trchk_perc_increase;    
    ar_a_aux_amz = ar_a_state_level_amz + ar_aux_change_frac_amz.*fl_a_aux_max;
    
    ar_root_zero = find_a_working(...
        ar_a_aux_amz,...
        j_amz,a_amz,eta_amz,educ_amz,married_amz,kids_amz,...
        TR,welf_checks);    
end