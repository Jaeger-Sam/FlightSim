% flat_to_ellipsodial.fcn computes the latitude (PHI), longitude (PSI), and
% altitude (H) from flat earth velocities (x_dot, y_dot, z_dot). This
% function integrates equation 11.12.20 from Phillips Mechanics of Flight.
% Lat/long are all expressed in radians.
%
% PHI_PSI_H = flat_to_ellipsodial(t_out,x_f_dot,PHI_PSI_H_0)
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

function PHI_PSI_H = flat_to_ellipsodial(t_out,x_f_dot,PHI_PSI_H_0)
    
    [~,PHI_PSI_H] = ode45(@(t,x)ellipsoidal_EOMS(t,x,x_f_dot),t_out,PHI_PSI_H_0);

    function x_dot = ellipsoidal_EOMS(t,x,x_f_dot)
        PHI = x(1);
        %PSI = x(2);
        H = x(3);

        R_e = 6378.1363*3280.84; % Equitorial radius of earth (km to ft)
        eps2 = 0.006694385; % eccentricity of earth

        R_x = R_e*(1-eps2)./(1 - (eps2*(sin(PHI).^2))).^(3/2);
        R_y = R_e./(1 - (eps2*(sin(PHI).^2))).^(1/2);

        x_dot(1) = x_f_dot(1)./(R_x + H);
        x_dot(2) = x_f_dot(2)./((R_y + H)*cos(PHI));
        x_dot(3) = - x_f_dot(3);

        x_dot = x_dot'; % Return column vector
    end
end
