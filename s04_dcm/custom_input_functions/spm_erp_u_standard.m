function [u] = spm_erp_u(t,P,M)
% Custom input function for standard input (at 564ms)

% This function is the original SPM spm_erp_u.m.

% returns the [scalar] input for EEG models (Gaussian function)
% FORMAT [u] = spm_erp_u(t,P,M)
% t      - PST (seconds)
% P      - parameter structure
% P.R  - scaling of [Gaussian] parameters
%
% u   - stimulus-related (subcortical) input
%
%--------------------------------------------------------------------------
% Preliminaries - check durations (ms)
%--------------------------------------------------------------------------
try
    if length(M.dur) ~= length(M.ons)
        M.dur = M.dur(1) + M.ons - M.ons;
    end
catch
    M.dur = 32 + M.ons - M.ons;
end

%--------------------------------------------------------------------------
% stimulus - gaussian (subcortical) impulse 
%--------------------------------------------------------------------------
nu = length(M.ons);           
u  = sparse(length(t), nu);   
t  = t * 1000; 
               
for i = 1:nu
    
    %----------------------------------------------------------------------
    % Gaussian bump function 
    %----------------------------------------------------------------------
    delay = M.ons(i) + 128 * P.R(i,1);
    scale = M.dur(i) * exp(P.R(i,2));

    % gaussian
    U = exp(-(t - delay).^2 / (2 * scale^2));

    u(:,i) = 32 * U; 
    
end

end
