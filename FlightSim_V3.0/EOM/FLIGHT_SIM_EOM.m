% Gives derivative of state given the current state x, control vector u and 
% disturbance vector w. Formulation has alpha_dot, beta_dot as an input.
% Also needs data structure of aircraft, sim_options
%
% x_dot = FLIGHT_SIM_EOM(t, x, control_vec, distrub_vec, alphabeta_dot, aircraft, sim_options)
%
% Flat Eart Quaternion state space EOM. See Phillips Mechanics of Flight
% Sections 11.11, 11.8.  There are different constraint EOMs. Additionally,
% see Chapters 1 & 2 in Stevens & Lewis Flight Control and Simulation.
% 
% Everything must be in US units (ft-lb-s).
%
% INPUTS:
%   t: time
%   x: state vector
%       x(1) = u
%       x(2) = v
%       x(3) = w
%       x(4) = p
%       x(5) = q
%       x(6) = r
%       x(7) = x_f (north)
%       x(8) = y_f (east)
%       x(9) = z_f (down)
%       x(10) = e_0
%       x(11) = e_x
%       x(12) = e_y
%       x(13) = e_z
%       x(14) = du1 (ONLY FOR GUSTS)
%       x(15) = dv1 (ONLY FOR GUSTS)
%       x(16) = dv2 (ONLY FOR GUSTS)
%       x(17) = dw1 (ONLY FOR GUSTS)
%       x(18) = dw2 (ONLY FOR GUSTS)
%       x(19) = dp1 (ONLY FOR GUSTS)
%       x(20) = dq1 (ONLY FOR GUSTS)
%       x(21) = dr1 (ONLY FOR GUSTS)
%   u: control_vec
%       u(1): delta_T
%       u(2): delta_e
%       u(3): delta_a
%       u(4): delta_r
%   disturb_vec:
%       disturb_vec(1): delta_u noise in x
%       disturb_vec(2): delta_v noise in y
%       disturb_vec(3): delta_w noise in z
%       disturb_vec(4): delta_p noise in roll
%   alphabeta_dot: vector of alpha_dot and beta_dot
%   aircraft: data structure of aircraft
%   sim_options: data structure of simulation options
%
% OUTPUTS:
%   x_dot: vector of derivative of state space
%       
% Sam Jaeger
% jaege246@umn.edu
% 2/18/2024
%   Revised: 3/6/2024
%   Revised: 11/3/2024
%       Added WGS84 gravity model, Coriolis accel
%   Revised: 9/21/2025
%       Added disturbances, alpha_dot/beta_dot as inputs
%   Revised: 10/5/2025
%       Added Von Karman turbulence model

function x_dot = FLIGHT_SIM_EOM(t, x, control_vec, disturb_vec, alphabeta_dot, aircraft, sim_options)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define variables for ease of use
    if sim_options.disturb.gusts == true
        % disturbances -> Von Karman Turbulence Model
        Wu = disturb_vec(1);
        Wv = disturb_vec(2);
        Ww = disturb_vec(3);
        Wp = disturb_vec(4);
    
        V = norm(x(1:3));
        Lfa = 1750; % length scale of atmosphere (ft)
        if -x(9) < Lfa % within the planetary boundary layer
            Lu = (Lfa^(2/3))*((-x(9))^(1/3));
            Lv = Lu;
            Lw = -x(9);
        else
            Lu = Lfa;
            Lv = Lfa;
            Lw = Lfa;
        end
    
        u1 = x(14);
        v1 = x(15);
        v2 = x(16);
        w1 = x(17);
        w2 = x(18);
        p1 = x(19);
        q1 = x(20);
        r1 = x(21);
    
        % u1dot
        gust_dot(1) = (-V/Lu)*u1 + Wu;
        
        % v1,v2 dot
        gust_dot(2) = v2;
        gust_dot(3) = (-(V/Lv)^2)*v1 + (-2*V/Lv)*v2 + Wv;
        
        % w1,w2 dot
        gust_dot(4) = w2;
        gust_dot(5) = (-(V/Lw)^2)*w1 + (-2*V/Lw)*w2 + Ww;
        
        % p1 dot
        gust_dot(6) = -pi*V/4/aircraft.geom.b_w*p1 + Wp;
    
        du = sqrt(2*V/pi/Lu)*sim_options.disturb.sig_u;
        dv = sqrt(3*V/pi/Lv)*sim_options.disturb.sig_v*[V/sqrt(3)/Lv, 1]*[v1; v2];
        dw = sqrt(4*V/pi/Lw)*sim_options.disturb.sig_w*[V/sqrt(3)/Lw, 1]*[w1; w2];
        dp = ((pi/4/aircraft.geom.b_w)^(7/6))*sqrt(4*V/5)/(Lw^(1/3))*sim_options.disturb.sig_w*p1;
        dq = q1;
        dr = r1;
    
        % qdot,rdot
        gust_dot(7) = -pi*V/4/aircraft.geom.b_w*q1 + pi/4/aircraft.geom.b_w*dw;
        gust_dot(8) = -pi*V/3/aircraft.geom.b_w*r1 + pi/3/aircraft.geom.b_w*dv;

    else
        du = 0;
        dv = 0;
        dw = 0;
        dp = 0;
        dq = 0;
        dr = 0;
    end

    u = x(1) + du;
    v = x(2) + dv;
    w = x(3) + dw;
    
    p = x(4) + dp;
    q = x(5) + dq;
    r = x(6) + dr;
    %omega = [p; q; r]; % rotation rate
    
    %x_f = x(7);
    %y_f = x(8);
    z_f = x(9);
    
    E = [x(10);x(11);x(12);x(13)]; % vector of quaterions
    E = E/(norm(E)); % renomalize for drift error, should have length 1.
    e_0 = E(1);
    e_x = E(2);
    e_y = E(3);
    e_z = E(4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Aircraft & Enviroment Properties %%%%
    mass = aircraft.mass; % Aircraft Mass
    g0 = 32.17404855643; % at sea level ft/s^2
    V_w = sim_options.disturb.V_w; % Wind in earth fixed coords
    
    if sim_options.WGS84_gravity == true
        % Acceleration due to gravity (and centripedal accel)
        g = grav_WGS84(-z_f, sim_options.earth.PHI_PSI_H_0(1)*180/pi,'US'); 
        omega_E = 7.2921150e-5; %rad/s sidereal rate of rotation of Earth
        % rotation of Earth in North-East-Down Earth CS
        omega_E = [cos(sim_options.earth.PHI_PSI_H_0(1)); 0; -sin(sim_options.earth.PHI_PSI_H_0(1))]*omega_E;
    else
        if -z_f <= sim_options.const_grav_alt % don't use variable gravity below threshold
            g = g0;
        else
            g = grav(-z_f);
        end
    end

    % Atmospheric density
    if -z_f <= sim_options.const_dens_alt
        rho = 0.0023770; % in slugs/ft^3
    else
        rho = ATMOS_1976(altitude,'US',0);
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Forces from Aerodynamics & Propulsion %%%%%
    [F_b, M_b] = AeroForces(x, rho, alphabeta_dot, control_vec, aircraft);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% !Derivatives! %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if sim_options.constrain_roll == true
        p=0;
        omega = [p; q; r]; % rotation rate

        % In inertial (flat earth fixed) reference frame!
        DCM_uvw = [ (e_x^2 + e_0^2 - e_y^2 - e_z^2), 2*(e_x*e_y - e_z*e_0),           2*(e_x*e_z + e_y*e_0);
                    2*(e_x*e_y + e_z*e_0),           (e_y^2 + e_0^2 - e_x^2 - e_z^2), 2*(e_y*e_z - e_x*e_0);
                    2*(e_x*e_z - e_y*e_0),            2*(e_y*e_z + e_x*e_0),        (e_z^2 + e_0^2 - e_x^2 - e_y^2)];
        
        % ground track velocity
        xyz_dot = DCM_uvw*[u;v;w] + V_w; % 11.11.3
        if sim_options.WGS84_gravity == true
            cor_accel = - 2*cross(omega_E,xyz_dot);
            uvw_dot = g0/mass.W*F_b + DCM_uvw*g + [((r*v)-(q*w)); (-(r*u)); ((q*u))] + cor_accel;
        else
        % Eq 11.11.1
            uvw_dot = g0/mass.W*F_b + g*[2*((e_x*e_z) - (e_y*e_0)); 2*((e_y*e_z) + (e_x*e_0)); (e_z^2 + e_0^2 - e_x^2 - e_y^2)] + [((r*v)-(q*w)); (-(r*u)); ((q*u))];
        end
        
        % Euler's full equations {omega} x [I]*{omega}
        euler_terms(1) = 0;
        euler_terms(2) = (mass.Inertia(3,3) - mass.Inertia(1,1))*p*r + mass.Inertia(1,3)*(r^2 - p^2) + mass.Inertia(1,2)*q*r - mass.Inertia(2,3)*p*q;
        euler_terms(3) = (mass.Inertia(1,1) - mass.Inertia(2,2))*p*q + mass.Inertia(1,2)*(p^2 - q^2) + mass.Inertia(2,3)*p*r - mass.Inertia(1,3)*q*r;
        
        omega_dot = mass.inv_Inertia*(mass.gyro*omega + euler_terms' + M_b); % Eq 11.11.2
        omega_dot(1) = 0;

        % DCM formulation seems to be slightly faster. Too many function calls?
        % xyz_dot = body_to_earth(E,[u,v,w])' + V_w; % Try Eq 11.6.8 for computational speed

        % Eq. 11.11.4
        E_dot = 0.5*[-e_x, -e_y, -e_z; e_0, -e_z, e_y; e_z, e_0, -e_x; -e_y, e_x, e_0]*omega;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif sim_options.constrain_pitch == true
        q=0;
        omega = [p; q; r]; % rotation rate
    
        % Flat Earth to body
        %   In inertial (flat earth fixed) reference frame!
        DCM_uvw = [ (e_x^2 + e_0^2 - e_y^2 - e_z^2), 2*(e_x*e_y - e_z*e_0),           2*(e_x*e_z + e_y*e_0);
                    2*(e_x*e_y + e_z*e_0),           (e_y^2 + e_0^2 - e_x^2 - e_z^2), 2*(e_y*e_z - e_x*e_0);
                    2*(e_x*e_z - e_y*e_0),            2*(e_y*e_z + e_x*e_0),        (e_z^2 + e_0^2 - e_x^2 - e_y^2)];
    
        % ground track velocity
        xyz_dot = DCM_uvw*[u;v;w] + V_w; % 11.11.3

        if sim_options.WGS84_gravity == true
            cor_accel = - 2*cross(omega_E,xyz_dot);
            uvw_dot = g0/mass.W*F_b + DCM_uvw*g + [((r*v)); ((p*w)-(r*u)); (-(p*v))] + cor_accel;
        else
        % Eq 11.11.1
            uvw_dot = g0/mass.W*F_b + g*[2*((e_x*e_z) - (e_y*e_0)); 2*((e_y*e_z) + (e_x*e_0)); (e_z^2 + e_0^2 - e_x^2 - e_y^2)] + [((r*v)); ((p*w)-(r*u)); (-(p*v))];
        end

        % Euler's full equations {omega} x [I]*{omega}
        euler_terms(1) = mass.Inertia(2,3)*( - r^2) - mass.Inertia(1,2)*p*r;
        euler_terms(2) = 0;
        euler_terms(3) = mass.Inertia(1,2)*(p^2) + mass.Inertia(2,3)*p*r;
        
        omega_dot = mass.inv_Inertia*(mass.gyro*omega + euler_terms' + M_b); % Eq 11.11.2
        omega_dot(2) = 0;
        
        % DCM formulation seems to be slightly faster. Too many function calls?
        % xyz_dot = body_to_earth(E,[u,v,w])' + V_w; % Try Eq 11.6.8 for computational speed

        % Eq. 11.11.4
        E_dot = 0.5*[-e_x, -e_y, -e_z; e_0, -e_z, e_y; e_z, e_0, -e_x; -e_y, e_x, e_0]*omega;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif sim_options.constrain_yaw == true
        r=0;
        omega = [p; q; r]; % rotation rate

        % Flat Earth to body
        %   In inertial (flat earth fixed) reference frame!
        DCM_uvw = [ (e_x^2 + e_0^2 - e_y^2 - e_z^2), 2*(e_x*e_y - e_z*e_0),           2*(e_x*e_z + e_y*e_0);
                    2*(e_x*e_y + e_z*e_0),           (e_y^2 + e_0^2 - e_x^2 - e_z^2), 2*(e_y*e_z - e_x*e_0);
                    2*(e_x*e_z - e_y*e_0),            2*(e_y*e_z + e_x*e_0),        (e_z^2 + e_0^2 - e_x^2 - e_y^2)];

        % ground track velocity
        xyz_dot = DCM_uvw*[u;v;w] + V_w; % 11.11.3
        
        if sim_options.WGS84_gravity == true
            cor_accel = - 2*cross(omega_E,xyz_dot);
            uvw_dot = g0/mass.W*F_b + DCM_uvw*g + [(-(q*w)); ((p*w)); ((q*u)-(p*v))] + cor_accel;
        else
        % Eq 11.11.1
            uvw_dot = g0/mass.W*F_b + g*[2*((e_x*e_z) - (e_y*e_0)); 2*((e_y*e_z) + (e_x*e_0)); (e_z^2 + e_0^2 - e_x^2 - e_y^2)] + [(-(q*w)); ((p*w)); ((q*u)-(p*v))];
        end

        % Euler's full equations {omega} x [I]*{omega}
        euler_terms(1) = mass.Inertia(2,3)*(q^2) + mass.Inertia(1,3)*p*q;
        euler_terms(2) = mass.Inertia(1,3)*(-p^2) - mass.Inertia(2,3)*p*q;
        euler_terms(3) = 0;
        
        omega_dot = mass.inv_Inertia*(mass.gyro*omega + euler_terms' + M_b); % Eq 11.11.2
        omega_dot(3) = 0;

        % DCM formulation seems to be slightly faster. Too many function calls?
        % xyz_dot = body_to_earth(E,[u,v,w])' + V_w; % Try Eq 11.6.8 for computational speed

        % Eq. 11.11.4
        E_dot = 0.5*[-e_x, -e_y, -e_z; e_0, -e_z, e_y; e_z, e_0, -e_x; -e_y, e_x, e_0]*omega;

        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif sim_options.pure_rolling == true
        q=0;
        r=0;
        omega = [p; q; r]; % rotation rate
        %e_y=0; e_z=0; E = [e_0;e_x;e_y;e_z];

        % flat Earth to body
        % In inertial (flat earth fixed) reference frame!
        DCM_uvw = [ (e_x^2 + e_0^2),              0,                0;
                    0,              (e_0^2 - e_x^2),     2*(-e_x*e_0);
                    0,                 2*( e_x*e_0), (e_0^2 - e_x^2)];

        % ground track velocity
        xyz_dot = DCM_uvw*[u;v;w] + V_w; % 11.11.3
        if sim_options.WGS84_gravity == true
            cor_accel = - 2*cross(omega_E,xyz_dot);
            uvw_dot = g0/mass.W*F_b + DCM_uvw*g + [0; ((p*w)); (-(p*v))] + [0; cor_accel(2); cor_accel(3)];
        else
        % Eq 11.11.1
            uvw_dot = g0/mass.W*F_b + g*[0; 2*( (e_x*e_0)); ( e_0^2 - e_x^2 )] + [0; ((p*w)); (-(p*v))];
        end
        
        % Euler's equations - constraints for 2 rotations makes everything zero
        euler_terms = [0,0,0];
        omega_dot = mass.inv_Inertia*(mass.gyro*omega + euler_terms' + M_b); % Eq 11.11.2
        omega_dot(2)=0;
        omega_dot(3)=0;

        % DCM formulation seems to be slightly faster. Too many function calls?
        % xyz_dot = body_to_earth(E,[u,v,w])' + V_w; % Try Eq 11.6.8 for computational speed
        
        % Eq. 11.11.4
        E_dot = 0.5*[-e_x, 0, 0; e_0, 0, 0; 0, e_0, -e_x; 0, e_x, e_0]*omega;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif sim_options.pure_pitching == true
        p=0;
        r=0;
        omega = [p; q; r]; % rotation rate
        %e_x=0; e_z=0; E = [e_0;e_x;e_y;e_z];

        % flat Earth to body
        %   In inertial (flat earth fixed) reference frame!
        DCM_uvw = [ (e_0^2 - e_y^2), 0,           2*( e_y*e_0);
                    0,           (e_y^2 + e_0^2), 0;
                    2*( - e_y*e_0),            0,        ( e_0^2 - e_y^2)];

        % ground velocity
        xyz_dot = DCM_uvw*[u;v;w] + V_w; % 11.11.3

        if sim_options.WGS84_gravity == true
            cor_accel = - 2*cross(omega_E,xyz_dot);
            uvw_dot = g0/mass.W*F_b + DCM_uvw*g + [(-(q*w)); 0; ((q*u))] + [cor_accel(1); 0; cor_accel(3)];
        else
        % Eq 11.11.1
            uvw_dot = g0/mass.W*F_b + g*[2*(- (e_y*e_0)); 0; ( e_0^2 - e_y^2)] + [(-(q*w)); 0; ((q*u))];
        end
        
        % Euler's equations - constraints for 2 rotations makes everything zero
        euler_terms = [0,0,0];
        omega_dot = mass.inv_Inertia*(mass.gyro*omega + euler_terms' + M_b); % Eq 11.11.2
        omega_dot(1)=0;
        omega_dot(3)=0;

        % DCM formulation seems to be slightly faster. Too many function calls?
        % xyz_dot = body_to_earth(E,[u,v,w])' + V_w; % Try Eq 11.6.8 for computational speed
        
        % Eq. 11.11.4
        E_dot = 0.5*[0, -e_y, 0; e_0, 0, e_y; 0, e_0, 0; -e_y, 0, e_0]*omega;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif sim_options.pure_yawing == true
        p=0;
        q=0;
        omega = [p; q; r]; % rotation rate
        %e_x=0; e_y=0; E = [e_0;e_x;e_y;e_z];

        % flat Earth to body
        %   In inertial (flat earth fixed) reference frame!
        DCM_uvw = [ (e_0^2 - e_z^2), 2*( - e_z*e_0),           0;
                    2*( e_z*e_0),           ( e_0^2 - e_z^2),  0;
                    0,            0,        (e_z^2 + e_0^2)];
    
        % ground track velocity 
        xyz_dot = DCM_uvw*[u;v;w] + V_w; % 11.11.3

        if sim_options.WGS84_gravity == true
            cor_accel = - 2*cross(omega_E,xyz_dot);
            uvw_dot = g0/mass.W*F_b + DCM_uvw*g + [((r*v)); (-(r*u)); 0] + [cor_accel(1); cor_accel(2); 0];
        else
        % Eq 11.11.1
            uvw_dot = g0/mass.W*F_b + g*[0; 0; (e_z^2 + e_0^2)] + [((r*v)); (-(r*u)); 0];
        end
        
        % Euler's equations - constraints for 2 rotations makes everything zero
        euler_terms = [0,0,0];
        omega_dot = mass.inv_Inertia*(mass.gyro*omega + euler_terms' + M_b); % Eq 11.11.2
        omega_dot(1)=0;
        omega_dot(2)=0;

        % DCM formulation seems to be slightly faster. Too many function calls?
        % xyz_dot = body_to_earth(E,[u,v,w])' + V_w; % Try Eq 11.6.8 for computational speed
        
        % Eq. 11.11.4
        E_dot = 0.5*[0, 0, -e_z; e_0, -e_z, 0; e_z, e_0, 0; 0, 0, e_0]*omega;
%%%%%%%%%%%%%%%%%%%%%%%%%
    else % unconstrained
        omega = [p; q; r]; % rotation rate

        % flat Earth to body fixed, (11.5.8) (11.7.1)
        %   In inertial (flat earth fixed) reference frame!
        DCM_uvw = [ (e_x^2 + e_0^2 - e_y^2 - e_z^2), 2*(e_x*e_y - e_z*e_0),           2*(e_x*e_z + e_y*e_0);
                    2*(e_x*e_y + e_z*e_0),           (e_y^2 + e_0^2 - e_x^2 - e_z^2), 2*(e_y*e_z - e_x*e_0);
                    2*(e_x*e_z - e_y*e_0),            2*(e_y*e_z + e_x*e_0),        (e_z^2 + e_0^2 - e_x^2 - e_y^2)];

        % ground track velocity 
        xyz_dot = DCM_uvw*[u;v;w] + V_w; % 11.11.3

        % DCM formulation seems to be slightly faster. Too many function calls?
        % xyz_dot = body_to_earth(E,[u,v,w])' + V_w; % Try Eq 11.6.8 for computational speed

        if sim_options.WGS84_gravity == true
            uvw_dot = g0/mass.W*F_b + DCM_uvw*g + [((r*v)-(q*w)); ((p*w)-(r*u)); ((q*u)-(p*v))] - 2*cross(omega_E,xyz_dot);
        else
        % Eq 11.11.1
            uvw_dot = g0/mass.W*F_b + g*[2*((e_x*e_z) - (e_y*e_0)); 2*((e_y*e_z) + (e_x*e_0)); (e_z^2 + e_0^2 - e_x^2 - e_y^2)] + [((r*v)-(q*w)); ((p*w)-(r*u)); ((q*u)-(p*v))];
        end
        
        % Euler's full equations {omega} x [I]*{omega}
        euler_terms(1) = (mass.Inertia(2,2) - mass.Inertia(3,3))*q*r + mass.Inertia(2,3)*(q^2 - r^2) + mass.Inertia(1,3)*p*q - mass.Inertia(1,2)*p*r;
        euler_terms(2) = (mass.Inertia(3,3) - mass.Inertia(1,1))*p*r + mass.Inertia(1,3)*(r^2 - p^2) + mass.Inertia(1,2)*q*r - mass.Inertia(2,3)*p*q;
        euler_terms(3) = (mass.Inertia(1,1) - mass.Inertia(2,2))*p*q + mass.Inertia(1,2)*(p^2 - q^2) + mass.Inertia(2,3)*p*r - mass.Inertia(1,3)*q*r;
        
        omega_dot = mass.inv_Inertia*(mass.gyro*omega + euler_terms' + M_b); % Eq 11.11.2
        
        
        
        % Eq. 11.11.4
        E_dot = 0.5*[-e_x, -e_y, -e_z; e_0, -e_z, e_y; e_z, e_0, -e_x; -e_y, e_x, e_0]*omega;
    end
  
%     % calculate alphabeta_dot
%     Vt = norm([u,v,w]);
%     Vtdot = (u*uvw_dot(1) + v*uvw_dot(2) + w*uvw_dot(3))/Vt; % Stevens 2.3-10c
%     alphabeta_dot(1,1) = (u*uvw_dot(3) - w*uvw_dot(1))/(u^2 + w^2); % Stevens 2.3-10a
%     alphabeta_dot(2,1) = (uvw_dot(2)*norm([u,v,w]) - v*Vtdot)/Vt/sqrt(u^2 + w^2); % Stevens 2.3-10b;

    % State vector
    %x_dot = [uvw_dot; omega_dot; xyz_dot; E_dot; alphabeta_dot];
    if sim_options.disturb.gusts == true
    x_dot = [uvw_dot; omega_dot; xyz_dot; E_dot; gust_dot'];
    else
    x_dot = [uvw_dot; omega_dot; xyz_dot; E_dot];
    end
end