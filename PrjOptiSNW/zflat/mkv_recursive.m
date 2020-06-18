function [m] = mkv_recursive(p,q,h)
%%
% PURPOSE   Computes the Marvov transition matrix due to Rouwenhorst(1995)
%           recursively for any given number of states h
% USAGE     m = mkv_recursive(p,q,h)
% INPUTS    p,q     : p+q-1 is the first-order serial autocorrelation.                     
%                     The matrix for h=2 looks like [p (1-p);(1-q) q]
%           h       : number of states        
% OUTPUTS   m       : Markov transition matrix
%           s_grid  : standard deviation of the discretized approximation of the
%           process
% REFERENCES
%           Rouwenhorst, K.G (1995): "Asset Pricing Implications of Equilibrium 
%           Business Cycle Models." Frontiers of Business Cycle Research,
%           Chapter 10.
%
if h == 1;
    m = 1;
else
    m = 1;
    for i = 1:h-1;
        zv = zeros(1,i)';
        m = p*[[m zv];[zv' 0]] + (1-p)*[[zv m];[0 zv']] + (1-q)*[[zv' 0];[m zv]] + q*[[0 zv'];[zv m]]; 
        m = m/2;
        m(1,:) = 2*m(1,:);
        m(length(m),:) = 2*m(length(m),:);
    end    
end

end