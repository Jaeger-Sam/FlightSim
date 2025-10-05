% terminate_flight.fcn terminates the simulation if the ground is hit or
% the aircraft departs from controlled flight. Called from RK4 function.
%
% [terminate] = terminate_flight(t,altitude, alpha, beta, phi_theta_psi,sim_options)
%
% INPUTS:
%   t: time
%   altitude: altitude (ft)
%   alpha: angle of attack (deg)
%   beta: angle of sideslip (deg)
%   phi_theta_psi: vector of bank, pitch, and yaw angles (rad)
%   sim_options: structure that defines the departure values
%
% OUTPUTS:
%   terminate: true or false
%
% Sam Jaeger
% jaege246@umn.edu
% 2/18/2024
%   Revised 9/21/25: removed load factor

function [terminate] = terminate_flight(t,altitude, alpha, beta, phi_theta_psi,sim_options)
    
    %[phi_theta_psi] = attitude_from_quats(x(10:13)');
    departure_deg = sim_options.departure_deg;
    
    %%%%%% GROUND HIT %%%%%%%%%%%%%%%%%%%%%
    if altitude<0 
        terminate = true;
        warning(append('Ground hit at t = ', num2str(t),' (s)'))
    else
        terminate = false;
    end
    if terminate == true
        return
    end
    %%%%%%%%%%%%%%%%% ALPHA %%%%%%%%%%%%%%%
    if alpha > departure_deg.alpha_depart_deg(2) % ALPHA positive
        terminate = true; % yes, halt integration
        warning(append('Alpha + departure at t = ', num2str(t),' (s),  alpha = ', num2str(alpha),' (deg)'))
    elseif alpha < departure_deg.alpha_depart_deg(1) % ALPHA negative
        terminate = true; % yes, halt integration
        warning(append('Alpha - departure at t = ', num2str(t),' (s),  alpha = ', num2str(alpha),' (deg)'))
    else
        terminate = false;
    end
    if terminate == true
        return
    end

    %%%%%%%%%%%%% BETA %%%%%%%%%%%%%%%%%%%%
    if beta > departure_deg.beta_depart_deg(2) % beta positive
        terminate = 1; % yes, halt integration
        warning(append('Beta + departure at t = ', num2str(t),' (s),  beta = ', num2str(beta),' (deg)'))
    elseif beta < departure_deg.beta_depart_deg(1) % beta negative
        terminate = 1; % yes, halt integration
        warning(append('Alpha - departure at t = ', num2str(t),' (s),  beta = ', num2str(beta),' (deg)'))
    else
        terminate = false;
    end
    if terminate == true
        return
    end

    %%%%%%%%%%%%% THETA %%%%%%%%%%%%%%%%%%%
    if phi_theta_psi(2)*180/pi > departure_deg.theta_depart_deg(2)
        terminate = true;
        warning(append('Theta + departure at t = ', num2str(t),' (s),  theta = ', num2str(phi_theta_psi(2)*180/pi),' (deg)'))
    elseif phi_theta_psi(2)*180/pi < departure_deg.theta_depart_deg(1)
        terminate = true;
        warning(append('Theta - departure at t = ', num2str(t),' (s),  theta = ', num2str(phi_theta_psi(2)*180/pi),' (deg)'))
    else
        terminate = false;
    end
    if terminate == true
        return
    end

    %%%%%%%%%%%%% PHI %%%%%%%%%%%%%%%%%%%%%
    if phi_theta_psi(1)*180/pi > departure_deg.phi_depart_deg(2)
        terminate = true;
         warning(append('Phi + departure at t = ', num2str(t),' (s),  phi = ', num2str(phi_theta_psi(1)*180/pi),' (deg)'))
    elseif phi_theta_psi(1)*180/pi < departure_deg.phi_depart_deg(1)
        terminate = true;
        warning(append('Phi - departure at t = ', num2str(t),' (s),  phi = ', num2str(phi_theta_psi(1)*180/pi),' (deg)'))
    else
        terminate = false;
    end
    if terminate == true
        return
    end

    %%%%%%%%%%%%% n %%%%%%%%%%%%%%%%%%%%%
    % if n > 15
    %     terminate = true;
    %     warning(append('load factor critcial at t = ', num2str(t),' (s),  n = ', num2str(n),' (g)'))
    % else
    %     terminate = false;
    % end
    % if terminate == true
    %     return
    % end

end