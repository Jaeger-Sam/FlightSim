% This function computes acceleration due to gravity given an altitude H in
% ft for earth.
%
% g = grav(H)
% 
% INPUTS:
%   H: scalar of altitude in ft
% 
% OUTPUTS:
%   g: acceleration due to gravity in ft/s^2
%
% Sam Jaeger
% jaege246@umn.edu
% 1/4/2024

function g = grav(H)
    % constants for earth in US units
    G = 1.068846e-9; % gravitational constant (ft^3 lb^-1 s^-2) 
    M = 1.317e25; % mass of earth (lb)
    r_e = 20856299.2; % radius of earth (ft)
    r= r_e + H; 
    
    g = G*M/(r^2);
end