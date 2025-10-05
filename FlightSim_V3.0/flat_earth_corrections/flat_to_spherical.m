% flat_to_spherical.fcn computes the latitude (PHI), longitude (PSI), and
% altitude (H) from flat earth velocities (x_dot, y_dot, z_dot). This
% function integrates equation 11.12.22 from Phillips Mechanics of Flight.
% Lat/long are all expressed in radians.
%
% PHI_PSI_H = flat_to_spherical(t_out,x_f_dot,PHI_PSI_H_0)
%
% INPUTS:
%   t_out: Vector of time corresponding to x_f_dot
%   x_f_dot: Matrix of x_dot, y_dot, z_dot of flight path
%   PHI_PSI_H_0: Vector of inital latitude, longitude, and altitudes.
%
% OUTPUTS:
%   PHI_PSI_H: Matrix of integrated latitude, longitude, and altitude.
%
% Sam Jaeger
% jaege246@umn.edu
% 1/9/2024

function PHI_PSI_H = flat_to_spherical(t_out,x_f_dot,PHI_PSI_H_0)
    
    [~,PHI_PSI_H] = ode45(@(t,x)spherical_EOMS(t,x,x_f_dot),t_out,PHI_PSI_H_0);

    function x_dot = spherical_EOMS(t,x,x_f_dot)
        PHI = x(1);
        %PSI = x(2);
        H = x(3);

        R_e = 6366.707*3280.84; % Equitorial radius of earth (km to ft)

        x_dot(1) = x_f_dot(1)./(R_e + H);
        x_dot(2) = x_f_dot(2)./((R_e + H)*cos(PHI));
        x_dot(3) = - x_f_dot(3);

        x_dot = x_dot'; % Return column vector
    end
end
