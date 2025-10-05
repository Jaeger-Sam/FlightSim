% force_to_stab.fcn converts body fixed forces (Fx, Fy, Fz) to standard
% aerodynamic forces (L,D,Y) using alpha and beta. See Equations 2.3-3a in
% Stevens and Lewis Flight Simulation and Control.
%
% F_s = force_to_stab(F_b,alpha,beta)
%
% INPUTS:
%   F_b: vector of body fixed forces [Fx; Fy; Fz]
%   alpha: angle of attack (deg)
%   beta: angle of sideslip (deg)
%
% OUTPUTS:
%   Fs: vector of forces in stability axis
%       Fs(1) = -D
%       Fs(2) = -Y
%       Fs(3) = -L
%
% Sam Jaeger
% jaege246@umn.edu
% 3/6/2024

function F_s = force_to_stab(F_b,alpha,beta)
    rot_mat = [cosd(alpha)*cosd(beta), sind(beta), sind(alpha)*cosd(beta);
              -cosd(alpha)*sind(beta), cosd(beta),-sind(alpha)*sind(beta);
              -sind(alpha),                     0,             cosd(alpha)];
    F_s = rot_mat*F_b;
end
