% AeroForces.fcn return the aerodynamic forces and moments for a particular
% aircraft at a particular state. Gives in body fixed coordinates.
%
% [F_b, M_b] = AeroForces(x, rho, alphabeta_dot, control_vec, aircraft)
%
% INPUTS:
%   x: state vector
%       x(1) = u
%       x(2) = v
%       x(3) = w
%       x(4) = p
%       x(5) = q
%       x(6) = r
%       x(7) = x_f
%       x(8) = y_f
%       x(9) = z_f
%       x(10) = e_0
%       x(11) = e_x
%       x(12) = e_y
%       x(13) = e_z
%   rho: atmospheric density (slugs/ft^3)
%   alphabeta_dot: column vector of alphabeta_dot
%   control_vec: u input forces
%   aircraft: data structure of aircraft properties
% 
% Outputs: 
%   F_b: 3x1 vector of forces corresponding to body fixed coordinates x,y,z
%   M_b: 3x1 vector of moments corresponding to body fixed coords x,y,z
%   
% Sam Jaeger
% jaege246@umn.edu
% 9/26/2023
%   Revised 10/20/23
%   Revised 1/5/23
%   Revised 1/15/2024, added state vector
%   Revised: 3/6/2024, correct function call for NewtAero
%   Revised: 9/18/2025, made rho as an input, changed
%       AeroForces,AeroCeofs,PropForces,


function [F_b, M_b] = AeroForces(x, rho, alphabeta_dot, control_vec, aircraft)
    if aircraft.aero.Newtonian_Aerodynamics == true
        [F_b, M_b] = NewtAero(x,control_vec,aircraft);
    else
        u = x(1);
        v = x(2);
        w = x(3);
        p = x(4);
        q = x(5);
        r = x(6);
        %altitude = -x(9);
    
        % Compute relevant quantities
        V = norm([u,v,w]); % Total Freestream Velocity
        alpha = atan2(w,u);
        beta = atan2(v,V);
        % if altitude <= sim_options.const_dens_alt
        %     rho = 0.0023770; % in slugs/ft^3
        % else
        %     rho = ATMOS_1976(altitude,'US',0);
        % end
    
        % alphabeta_dot
        alpha_dot = alphabeta_dot(1);
        beta_dot = alphabeta_dot(2);

        % non dim rates
        p_b = p*aircraft.geom.b_w/2/V; % roll
        q_b = q*aircraft.geom.c_b_w/2/V; % pitch
        r_b = r*aircraft.geom.b_w/2/V; % yaw
    
        % Propulsive Forces
        [F_P, M_P] = PropForces(x, rho, control_vec(1), aircraft.propulsion);
        n = aircraft.propulsion.Kmotor*control_vec(1); % rev/s
        J = V/n/aircraft.propulsion.d; 
       
        % Aero Coefficients of Vehicle
        [C_L, C_D, C_Y, C_l, C_m, C_n]= AeroCoefs(alpha, beta, alpha_dot, beta_dot, p_b, q_b, r_b, control_vec, aircraft);
        Delta_CL = change_lift(alpha,J,C_L,aircraft);
        C_L = C_L + Delta_CL; % lift addition due propulsion effects

        % body 
        %C_S = C_Y + C_D*sin(beta); % approximate for small beta
%         C_x = (C_L*sin(alpha) - C_Y*cos(alpha)*sin(beta) - C_D*cos(alpha)*cos(beta));
%         C_z = ( -C_L*cos(alpha) - C_Y*sin(alpha)*sin(beta) - C_D*sin(alpha)*cos(beta));
        
        %C_x = (C_L*sin(alpha) - C_Y*cos(alpha)*sin(beta) - C_D*cos(alpha)*cos(beta));
        %C_z = ( -C_L*cos(alpha) - C_Y*sin(alpha)*sin(beta) - C_D*sin(alpha)*cos(beta));
        
        %C_Y = -cos(beta)*C_Y - sin(beta)*C_D;

        C_x = cos(beta)*(C_L*cos(beta)*sin(alpha) - cos(alpha)*(C_D + C_Y*sin(beta)))/((sin(beta)^2)*(1 - 2*sin(alpha)^2) + 1);
        C_z = - (C_L*cos(alpha)*(sin(beta)^2 + 1) + cos(beta)*sin(alpha)*(C_Y*sin(beta) + C_D))/((sin(beta)^2)*(1 - 2*sin(alpha)^2) + 1);

        % body
        %C_z = -C_L*cos(alpha) - C_D*sin(alpha)*cos(beta);
        %C_x = C_L*sin(alpha) - C_D*cos(alpha)*cos(beta);

        % (negative) axial force
        F_x = F_P(1) + 0.5*rho*(V^2)*aircraft.geom.S_w*C_x;
    
        % side force
        F_y = F_P(2) + 0.5*rho*(V^2)*aircraft.geom.S_w*C_Y;
    
        % (negative) normal force
        F_z = F_P(3) + 0.5*rho*(V^2)*aircraft.geom.S_w*C_z;
    
    
        % Rolling Moment
        M_x = M_P(1) + 0.5*rho*(V^2)*aircraft.geom.S_w*aircraft.geom.b_w*C_l;
    
        % Pitching Moment
        M_y = M_P(2) + 0.5*rho*(V^2)*aircraft.geom.S_w*aircraft.geom.c_b_w*C_m;
    
        % Yawing Moment
        M_z = M_P(3) + 0.5*rho*(V^2)*aircraft.geom.S_w*aircraft.geom.b_w*C_n;
    
        % body fixed forces and moments
        F_b = [F_x; F_y; F_z];
        M_b = [M_x; M_y; M_z];
    end
end