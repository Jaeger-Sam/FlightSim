% PropForces.fcn computes the propulsion forces and moments given the state
% of the aircraft, control_input, and propulsion propertes. This script
% assumes that the propulsion system is placed in the centerline of the
% aircraft.
%
% [F_P, M_P] = PropForces(x, rho, delta_T, propulsion)
%
% INPUTS:
%   x: state vector of aircraft
%   rho: atmospheric density (slugs/ft^3)
%   delta_T: input of throttle from 0 to 1
%   propulsion: data structure of propulsion properties from initalize_sim
%       propulsion.Kmotor: mapping from throttle input to rev/s of prop
%       propulsion.d: diameter of propeller (ft)
%       propulsion.p_CT: polynomial of fitted thrust coeffients with
%           respect to advance ratio J
%       propulsion.p_f: polynomial of fitted function coeffients for normal
%           force computation (see Etkin B.7, Riber (1944))
%       propulsion.x_prop: x location of propeller with respect to the CG
%       propulsion.z_prop: z location of propeller with respect to the CG
%       
%       
%
% OUTPUTS:
%   F_P: vector of propulsive forces in body fixed x,y,z coords
%   M_P: vector of moments in body fixed x,y,z coords
%
% Written By:
% Sam Jaeger
% jaege246@umn.edu
% 10/20/2023
%   Revised: 1/5/2024
%   Revised: 1/15/2024, added state vector input
%   Revised: 9/18/2025, changed to propeller
%
% TO DO: add propeller dynamics / gyroscopic effects

function [F_P, M_P] = PropForces(x, rho, delta_T, propulsion)
     n = propulsion.Kmotor*delta_T; % rev/s
     d = propulsion.d; % diameter of propeller
     
     u = x(1);
     v = x(2);
     w = x(3);
     %altitude = -x(9);

    % Compute relevant quantities
    V = norm([u,v,w]); % Total Freestream Velocity
    alpha = atan2(w,u);
    beta = atan2(v,V);
    % if altitude <= 3000
    %     rho = 0.0023770;
    % else
    %     rho = ATMOS_1976(altitude,'US',0);
    % end
    qbar = 0.5*rho*V^2;

    % thrust
    J = V/n/d;
    if isnan(J)
        J = 0;
    elseif n==0
        J=0;
    end
    CT = polyval(propulsion.p_CT,J);
    if isnan(CT)
        CT=0;
    end
    % if CT <= 0 % no negative thrust
    %     CT = 0;
    % end

    if propulsion.prop_moments == true
        % prop normal force at an angle of attack
        %   Etkin, Appendix B
        C_Np = polyval(propulsion.p_f,CT/J^2)*0.08; % normal force coefficient
        Np = C_Np*qbar*pi*(d^2)/4;
        if isnan(Np)
            Np = 0;
        elseif n==0
            Np = 0;
        end
    else
        Np=0;
    end

    % forces in xyz
    F_P(1) = CT*rho*(n^2)*(d^4);
    F_P(2) = Np*sin(beta)*cos(alpha); % approximate alpha, beta direction
    F_P(3) = -Np*cos(beta)*sin(alpha); % approximate alpha, beta direction

    % moments
    M_P(1) = 0; % rolling moment
    M_P(2) = - F_P(1)*propulsion.z_prop + propulsion.x_prop*F_P(3); % pitching moment
    M_P(3) = propulsion.x_prop*F_P(2); % yawing moment
    F_P=F_P';
    M_P=M_P';
end

