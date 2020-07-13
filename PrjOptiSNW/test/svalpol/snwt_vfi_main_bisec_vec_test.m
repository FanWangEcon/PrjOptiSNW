% Test Policy Function Monotonicity
% the value function, and the policy function for c and ap should be
% monotonically increasing in savings state as well as shock state.

% Call function
mp_params = snw_mp_param('default_small');
mp_controls = snw_mp_control('default_test');
[V_VFI, ap_VFI, cons_VFI, mp_valpol_more] = snw_vfi_main_bisec_vec(mp_params, mp_controls);

jret = mp_params('jret');
n_eta_S_grid = mp_params('n_eta_S_grid');


%% Test 1: Value is increasing in asset, dim2

assert(min(diff(V_VFI,1,2), [], 'all') > 0, 'V is strictly increasing in a');

%% Test 2: Value is increasing in shock, dim3, before retirement
% but this is only the case when we only have one shock
if (n_eta_S_grid==1)
    assert(min(diff(V_VFI(1:(jret-1),:,:,:,:,:),1,3), [], 'all') > 0, ...
        'V is strictly increasing in husband incomoe before retirement');
end

%% Test 3: Value is Invariant to shocks, after retirement
% after retirement, when spouse has shock, shock still has impact post
% retirement. 

if (n_eta_S_grid==1)
    assert(sum(diff(V_VFI(jret:end,:,:,:,:,:),1,3), 'all') <= 1e-10, ...
        'V is strictly increasing in z before retirement');
end

%% Test 5: Ap is increasing in ap + c
% more asset, higher resource availability, must lead to either higher sum
% of consumption and savings choices. For the actual underlying model,
% rising a should lead to monotonically strictly increasing c and ap both,
% however, due to the discretization of the asset grid and the step-wise EV
% derivative dEV/dap, dc/dz could be zero. 

assert(min(diff(ap_VFI,1,2) + diff(cons_VFI,1,2), [], 'all') > 0, ...
    'ap+c is strictly increasing in a');

%% Test 6: Consumption and AP are non-decreasing in asset 
% should be strictly increasing, but due to step-function given linear
% approximation of EV, consumption might not increase when a increases.
% Savings choices should also be non-decreasing

assert(min(diff(cons_VFI,1,2), [], 'all') <= 1e-10, 'C is non-decreasing in a');
assert(min(diff(ap_VFI,1,2), [], 'all') <= 1e-10, 'AP is non-decreasing in a');

