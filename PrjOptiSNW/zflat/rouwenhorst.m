function [y_grid,p] = rouwenhorst(rho,sig_u,ngrid)
%%
% PURPOSE   Discretize an AR(1) of the form z = rho*z(-1)+u
% USAGE     [y_grid,p] = rouwenhorst(rho,sig_u,ngrid)
% INPUTS    rho     : first-order serial autocorrelation (1x1)           
%           sig_u   : standard deviation of u (1x1)
%           ngrid   : number of states (1x1)    
% OUTPUTS   y_grid  : discretized values
%           p       : Markov transition matrix
% REFERENCES
%           1.Rouwenhorst, K.G (1995): "Asset Pricing Implications of Equilibrium 
%           Business Cycle Models." In Frontiers of Business Cycle
%           Research, Chapter 10.
% USES      mkv_recursive.m
%
sig_y = sqrt(sig_u^2/(1-rho^2));
% contructing the grid
y_end = sqrt(ngrid-1)*sig_y;
y_beg = -y_end;
step = (y_end-y_beg)/(ngrid-1); % Could have used the command linspace

if ngrid == 1
    y_grid = 0;
else
    y_grid = (y_beg:step:y_end)';
end

p0 = 0.5*(1+rho);                                                           % from calibration rho = p+q-1 together wit p=q
q0 = p0;
p = mkv_recursive(p0,q0,ngrid);

end

%% Own notes
% This procedure depends on the no. of nodes (i.e., no. of grid points) we
% want to use. Hence, m in Tauchen is equal to sqrt(N-1) in Rouwenhorst