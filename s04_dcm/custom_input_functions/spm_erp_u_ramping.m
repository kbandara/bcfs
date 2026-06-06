function [u] = spm_erp_u(t,P,M)
% Custom input function for ramping input
% based on spm_erp_u.m 

%{

 used for gradual face emergence in bCFS.
 ramp starts at 0ms (fixed), duration is a free parameter.

 this models the scalar input (u) in dcm as a gradual slope (rather than the standard gaussian bump) 
 here we use a simple model based on Stevens' power law (1957):

  u(t) = (t / D)^P,  constrained to [0, 1]
  where,
     D = duration to reach full contrast
     P = power/curvature exponent

  t - time vector (seconds)
  P - parameter structure with P.R (learned parameters)
         P.R(1) - Exponent/curvature 
         P.R(2) - Duration (time to reach full visibility)
  M - model structure with:
       M.ons - onset (NOT USED - fixed at 0)
       M.dur - base duration (ms) for ramping
  u - input matrix 
%}

% preliminaries - set base duration
%--------------------------------------------------------------------------
base_duration = 2125;  % hard code the default duration (not taken from DCM.options.Dur) this is the midpoint bw 1250ms and 3000ms
ramp_onset = 0;        % hard code onset to 0 

t_ms = t(:) * 1000;  

% extract parameters
if isvector(P.R)
    pr = P.R(:);
else
    pr = P.R(1,:)';
end

%--------------------------------------------------------------------------
% Parameter 1 - exponent
%--------------------------------------------------------------------------

try
    P_exp = exp(0.5 * pr(1));
catch
    P_exp = 1;  % default is linear
end

% Soft constraints on exponent
P_exp = max(0.2, min(P_exp, 5));  % 0.2 to 5

%--------------------------------------------------------------------------
% Parameter 2 - duration/slope 
%--------------------------------------------------------------------------

try
    D = base_duration * exp(0.5*pr(2));
catch
    D = base_duration;
end

% Soft constraints on duration
D = max(500, min(D, 3000));

%--------------------------------------------------------------------------
% create ramp
%--------------------------------------------------------------------------

% normalized time 
t_norm = (t_ms - ramp_onset) / D;

% clip to valid range
t_norm(t_norm < 0) = 0;
t_norm(t_norm > 1) = 1;

% power law
u = t_norm .^ P_exp;

u = 32 * u; %scale
u = u(:); %formatting

end
