% stab_to_body.fcn calculates the body fixed velocities: u,v,w from the
% angle of attack, angle of sideslip, and freestream velocity. Input units
% must be consistent with output units. See Equation 2.3-6a in Stevens and
% Lewis Flight Simulation and Control.
%
% x = stab_to_body(alpha,beta,V_inf)
%
% INPUTS:
%   alpha: angle of attack (deg)
%   beta: angle of sideslip (deg)
%   V_inf: freestream velocity
%
% OUTPUTS:
%   x(1): u
%   x(2): v
%   x(3): w
%
% Sam Jaeger
% jaege246@umn.edu
% 3/6/2024

function x = stab_to_body(alpha,beta,V_inf)
    x(1,1) = V_inf*cosd(alpha)*cosd(beta);
    x(2,1) = V_inf*sind(beta);
    x(3,1) = V_inf*sind(alpha)*cosd(beta);
end