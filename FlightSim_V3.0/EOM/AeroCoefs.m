% AeroCoefs computes the linear aerodynamic coefficients given alpha,beta,
% alpha_dot,beta_dot, p_b,q_b,r_b, control_vec. This can be either in a
% matrix form or a lookup table with an iterpolation scheme. Matrix form is
% significantly faster. Lookup table only has longitudinal coefficients 
% supported. The matrix formulation includes a 
% longitudinal poststall model and transforms steady dependence to a CG 
% location. The unsteady dervatives to be evaluated at the CG. 
% Propulsion corrections are not applied.
%
% [C_L, C_D, C_Y, C_l, C_m, C_n] = AeroCoefs(alpha, beta, alpha_dot, beta_dot, p_b, q_b, r_b, control_vec, aircraft)
%
% INPUTS:
%   alpha: angle of attack in radians
%   beta: angle of sideslip in radians
%   alpha_dot: time derivative of angle of attack in rad/s
%   beta_dot: time derative of angle of sideslip in rad/s
%   p_b: non dimensional roll rate
%   q_b: non dimensional pitch rate
%   r_b: non dimensional yaw rate
%   control_vec: u input
%       control_vec(1): delta_T (not used)
%       control_vec(2): delta_e
%       control_vec(3): delta_a
%       control_vec(4): delta_r
%       control_vec(5): delta_f
%   aircraft: data structure with the information
%       mass.xcg: x location of cg (ft)
%       mass.zcg: z location of cg (ft)
%       geom.c_b_w: mean aero chord used to nondimensionalize (ft)
%       geom.b_w: wingspan used to nondimensionalize (ft)
%       aero.unsteady: true/false to use unsteady terms
%       aero.nonlin: true/false to use nonlinear terms
%       aero.aero_data_mat: true/false if aerocoefs are in the form of a
%           lookup table, currently only impelemnted for longitudinal
%       aero.S: matrix of linear derivatives
%       aero.Su: matrix of unsteady derivatives
%       aero.Sn: matrix of nonlinear derivatives of alpha*beta, (alpha^2)*beta
%       aero.CD_mat: diagonal matrix of drag dependence
%       aero.C_D_0: parasitic drag coefficient 
%       aero.C_m_0: static pitching moment coefficient
%       aero.C_L_0: static lift coefficient
%       aero.alpha_stall: stall angle of attack (deg)
%       aero.alpha_b: stall blending angle 
%
% OUTPUTS:
%   C_L: lift coefficient
%   C_D: Drag Coefficient
%   C_Y: Side Force Coefficient
%   C_l: Rolling moment Coefficient
%   C_m: Pitching moment Coefficient
%   C_n: Yawing moment Coefficient
%
% Written By:
% Sam Jaeger
% jaege246@umn.edu
% 10/20/2023
%   Revised: 1/5/2024, Added input for aircraft data structure
%   Revised: 1/6/2024, Added input for Mach and Reynolds dependence
%   Revised: 1/15/2024, Added alpha_dot and beta_dot inputs
%   Revised: 11/3/2024, Changed interp2 to have a method input
%   Revised: 9/18/2025, Changed to matrix structure, removed Mach/Re dep,
%       added stall model to CL,Cm,CD

function [C_L, C_D, C_Y, C_l, C_m, C_n] = AeroCoefs(alpha, beta, alpha_dot, beta_dot, p_b, q_b, r_b, control_vec, aircraft)
    aero = aircraft.aero;
    
    if aero.aero_data_mat == true
        delta_e = control_vec(2);
        % some conditioning so matlab doesn't break itself
        if alpha > aero.alpha_table(end)
            alpha = aero.alpha_table(end);
        elseif alpha < aero.alpha_table(1)
            alpha = aero.alpha_table(1);
        end
        if delta_e > aero.delta_table(end)
            delta_e = aero.delta_table(end);
        elseif delta_e < aero.delta_table(1)
            delta_e = aero.delta_table(1);
        end

        C_L = interp2( aero.delta_table, aero.alpha_table, aero.CL, delta_e, alpha, aero.method);
        C_D = interp2( aero.delta_table, aero.alpha_table, aero.CD, delta_e, alpha, aero.method);
        C_m = interp2( aero.delta_table, aero.alpha_table, aero.Cm, delta_e, alpha, aero.method);
        C_Y = 0;
        C_l = 0;
        C_n = 0;

    else
        xaero = [alpha;beta;control_vec(2:end)]; % steady dependence
        C_lin = aero.S*xaero;
    
        if aero.unsteady== true % unsteady dependence
            xaerou = [alpha_dot;beta_dot;p_b;q_b;r_b]; 
            C_unsteady = aero.Su*xaerou;
        end
        if aero.nonlin == true % nonlinear dependence
            xaeron = [alpha*beta; (alpha^2)*beta]; 
            C_nonlin =  aero.Sn*xaeron;
            C_steady = C_lin + C_nonlin;
        else
            C_steady = C_lin ;
        end
    
        % steady force coefficients
        C_Ls = C_steady(1) + aero.C_L_0;
        C_Ys = C_steady(2);
        C_D = xaero'*aero.CD_mat*xaero + aero.C_D_0;
    
        % calulate moments at CG
        C_ls = C_steady(3) + (aircraft.mass.zcg/aircraft.geom.b_w)*C_Ys;
        C_ms = (C_steady(4) + aero.C_m_0) + ...
            (aircraft.mass.xcg/aircraft.geom.c_b_w)*(-C_Ls*cos(alpha) - C_D*sin(alpha)*cos(beta)) ...
            - (aircraft.mass.zcg/aircraft.geom.c_b_w)*(C_Ls*sin(alpha) - C_D*cos(alpha)*cos(beta));
        C_ns = C_steady(5) - (aircraft.mass.xcg/aircraft.geom.b_w)*C_Ys;
        
        % add in unsteady components
        if aero.unsteady == true
            C_L = C_Ls + C_unsteady(1);
            C_Y = C_Ys + C_unsteady(2);
            C_l = C_ls + C_unsteady(3);
            C_m = C_ms + C_unsteady(4);
            C_n = C_ns + C_unsteady(5);
        else
            C_L = C_Ls;
            C_Y = C_Ys;
            C_l = C_ls;
            C_m = C_ms;
            C_n = C_ns;
        end
    
        % stall model for longitudinal coefs
        if alpha > aero.alpha_stall % stall model
            M = aero.M;
            alpha_b = aero.alpha_b; % blending function parameter
            sigma = (1 + exp(-M*(alpha - alpha_b)) + exp(M*(alpha+alpha_b)) )/(1 + exp(-M*(alpha-alpha_b)) )/(1 + exp(M*(alpha+alpha_b)) );
            
            C_L_plate = 2*sign(alpha)*(sin(alpha)^2)*cos(alpha);
            C_D_plate = 2*sin(abs(alpha))^(3/2);
            C_m_plate = -0.8*sin(alpha);
        
            C_L = (1 - sigma)*C_L + sigma*C_L_plate;
            C_D = (1 - sigma)*C_D + sigma*C_D_plate;
            C_m = (1 - sigma)*C_m + sigma*C_m_plate;
        end
    end
end