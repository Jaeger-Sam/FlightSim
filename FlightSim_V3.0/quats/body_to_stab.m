% body_to_stab.fcn calculates the angle of attack, angle of side sideslip, 
% freestream velocity as well as the flank angle from the body fixed 
% velocities. Input units must be consistent with outputs.
%
% [alpha,beta,V_inf, beta_f] = body_to_stab(x)
%
%
% INPUTS:
%   x(1): u
%   x(2): v
%   x(3): w
%
% OUTPUTS:
%   alpha: angle of attack (deg)
%   beta: angle of sideslip (deg)
%   V_inf: freestream velocity
%   beta_f: flank angle (deg)
%
% Sam Jaeger
% jaege246@umn.edu
% 3/6/2024

function [alpha,beta,V_inf, beta_f] = body_to_stab(x)
    V_inf = norm(x(1:3,1));
    alpha = atan2( x(3,1), x(1,1))*180/pi;
    beta_f = atan2( x(2,1), x(1,1))*180/pi; % flank angle
    beta = atan(cos( alpha*pi/180).*tan( beta_f*pi/180))*180/pi; % true sideslip
end