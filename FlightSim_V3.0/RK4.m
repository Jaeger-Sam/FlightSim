% RK4 is a self written Runge-Kutta algorithm to integrate the equations of
% motion of FLIGHT_SIM_EOM.fcn. Described in Section 3.6 of Stengel Flight 
% Dynamics. See also Section 3.4 in Stevens and Lewis Flight Simulation and
% Control. 
%
% [sim_data_out] = RK4(aircraft,sim_options)
%
% =========================================================================
% =========================================================================
% INPUTS:
%   aircraft: data sturcture from INITALIZE_SIMULATION
%   sim_options: data structure from INITALIZE_SIMULATION
%
% =========================================================================
% OUTPUTS:
%   sim_data_out: data structure with the included variables
%           altitude: height CG is above ground (ft)
%           control_vec: Matrix of control deflections time history
%           Velocity: true airspeed (ft/s)
%           alpha: angle of attack (deg)
%           beta_f: flank angle (deg)
%           beta: side slip angle (deg)
%           alpha_dot: AoA rate (deg/s)
%           beta_f_dot: flank angle rate (deg/s)
%           beta_dot: AoS rate (deg/s)
%           t_out: vector of time (s)
%           X_out: matrix of integrated state variables
%           {u,v,w}: matrix of integrated air velocities (ft/s)
%           {p,q,r}: matrix of integrated angular velocities (rad/s)
%           {X,Y,Z}: matrix of integrated positions on earth (ft)
%           es_out: {e0,ex,ey,ez} quaternions
%           {phi,theta,psi}, matrix of recovered attitude (deg)
%           phi: bank angle (deg)
%           theta: pitch attitude (deg)
%           psi: heading (deg)
%           gamma: flight path angle (deg)
%           delta_T: Time history of thrust deflections [0,1]
%           delta_e: Time history of elevator deflections (rad)
%           delta_a: Time history of aileron deflections (rad)
%           delta_r: Time history of rudder deflections (rad)
%           delta_f: Time history of flap deflections (rad)
%           PHI_PSI_H: Latitude, longitude, altitude, for flight path
%               sim_options.earth.ellipsoidal_earth == true OR
%               sim_options.earth.spherical_earth == true OTHERWISE
%               don't calculate 
%       OPTIONAL OUTPUTS:==================================================
%        sim_options.save_atmos == true------------------------------------
%           rho: freestream density (slugs/ft^3)
%           T: freestream temperature (f)
%           p: freestream pressure (psi)
%           a: freestream speed of sound (ft/s)
%           nu: freestream Kinematic viscosity.
%           g: gravity (ft/s^2)
%        sim_options.save_Mach == true-------------------------------------
%           Mach: Mach number
%        sim_options.save_Re == true---------------------------------------
%           Reynolds: Reynolds number
%        sim_options.save_xdot == true ------------------------------------
%           x_dot: Matrix of state derivatives
%           ax: x acceleration at CG
%           ay: y acceleration at CG
%           az: z acceleration at CG
%        sim_options.save_forces == true-----------------------------------
%           F_b: Matrix of body fixed forces in x,y,z directions
%           M_b: Matrix of body fixed moments in x,y,z directions
%        sim_options.save_prop_forces == true------------------------------
%           F_p: Matrix of body fixed propulsion forces in x,y,z
%           M_p: Matrix of body fixed propulsion moments in x,y,z
%        sim_options.save_coefs == true -----------------------------------
%           CL: lift coefficient (time history
%           CD: drag coefficient (time history)
%           CY: side force coefficient (time history)
%           Cl: rolling moment coefficient (time history)
%           Cm: pitching moment coefficient (time history)
%           Cn: yawing moment coefficient (time history)
%           p_b: non-dim roll rate time hisotry
%           q_b: non-dim pitch rate time history
%           r_b: non-dim yaw rate time histtory
%           n: propeller rev/s
%           J: advance ratio
%           CT: thtrust coefficient
%           DCL: change in lift coefficient due to propeller
%        sim_options.save_load == true-------------------------------------
%           nx_g: g load in body fixed x direction
%           ny_g: g load in body fixed y direction
%           nz_g: g load in body fixed z direction
%
% Sam Jaeger
% jaege246@umn.edu
% 2/18/2024
%   Revised: 3/4/2024
%   Revised: 3/6/2024
%   Revised: 9/21/2025
%   Revised: 10/5/2025: added VK turb. & documentation

function [sim_data_out] = RK4(aircraft,sim_options)
    ICS = sim_options.ICS;
    t = sim_options.t;
    if size(ICS,1) ~= 13 && sim_options.disturb.gusts == false
        error('ICS vector must be 13 x 1!')
    elseif size(ICS,1) ~= 21 && sim_options.disturb.gusts == true
        error('ICS vector must be 21 x 1!')
    elseif size(t,1) > 1
        error('t must be a vector!')
    end
   
    n_step = length(t);
    %%%%%%%%%%%%%%%%%% INITALIZE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if sim_options.disturb.gusts == true
    x = NaN(21,n_step);
    else
    x = NaN(13,n_step);
    end
    x(:,1) = ICS;

    altitude = NaN(n_step,1); altitude(1) = - ICS(9);
    if sim_options.save_atmos == true
        rho = NaN(n_step,1); 
        T = NaN(n_step,1); 
        p = NaN(n_step,1); 
        a = NaN(n_step,1);
        nu = NaN(n_step,1); 
        g = NaN(n_step,1);
        gb = NaN(3,n_step);
        if altitude(1) <= sim_options.const_dens_alt % fix atmospheric variables below desired alt
            [rho(1), T(1), p(1), a(1), nu(1), ~, ~] = ATMOS_1976(0,'US',sim_options.disp_alt_warn);
            g(1) = 32.17404855643;
        else
            [rho(1), T(1), p(1), a(1), nu(1), g(1), ~] = ATMOS_1976(altitude(1),'US',sim_options.disp_alt_warn);
        end
    end
    if sim_options.save_grav == true
        gb(:,1) = earth_to_body(x(10:13,1)',[0;0;g(1)]');
    end

    control_vec = NaN(5,n_step); control_vec(:,1) = zeros(5,1);
    if sim_options.save_xdot == true
        if sim_options.disturb.gusts == true
        x_dot = NaN(21,n_step);
        else
        x_dot = NaN(13,n_step); 
        end
        x_dot(:,1) =  FLIGHT_SIM_EOM(t(1), x(:,1), control_vec(:,1), zeros(1,4), [0,0], aircraft, sim_options);
        ax = NaN(1,n_step); ax(1) = 0;
        ay = NaN(1,n_step); ay(1) = 0;
        az = NaN(1,n_step); az(1) = 0;
    end
    

    phi_theta_psi = NaN(n_step,3);
    phi_theta_psi(1,:) = attitude_from_quats(ICS(10:13)')';
    
    Velocity = NaN(n_step,1); 
    alpha  = NaN(n_step,1);
    beta_f = NaN(n_step,1);
    beta = NaN(n_step,1);
    alpha_dot = NaN(n_step,1);  alpha_dot(1) = 0;
    beta_f_dot = NaN(n_step,1); beta_f_dot(1) = 0;
    beta_dot = NaN(n_step,1); beta_dot(1) = 0;

    [alpha(1),beta(1),Velocity(1),beta_f(1)] = body_to_stab(x(:,1));

    if sim_options.save_forces == true
        F_b = NaN(3,n_step); 
        M_b  = NaN(3,n_step);
        [F_b(:,1), M_b(:,1)] = AeroForces(x(:,1), rho(1), [0,0], control_vec(:,1), aircraft);
    end
    if sim_options.save_coefs == true
        CL = NaN(n_step,1);CL(1)=0;
        CD = NaN(n_step,1);CD(1)=0;
        CY = NaN(n_step,1);CY(1)=0;
        Cl = NaN(n_step,1);Cl(1)=0;
        Cm = NaN(n_step,1);Cm(1)=0;
        Cn = NaN(n_step,1);Cn(1)=0;
        p_b = NaN(n_step,1);p_b(1)=0;
        q_b = NaN(n_step,1);q_b(1)=0;
        r_b = NaN(n_step,1);r_b(1)=0;
        n = NaN(n_step,1);n(1)=0;
        J = NaN(n_step,1);J(1)=0;
        CT = NaN(n_step,1);CT(1)=0;
        DCL = NaN(n_step,1);DCL(1)=0;
    end
    if sim_options.save_prop_forces == true
        F_p = NaN(3,n_step); 
        M_p  = NaN(3,n_step);
        [F_p(:,1), M_p(:,1)] = PropForces(x(:,1), rho(1), control_vec(1,1), aircraft.propulsion);
    end
    if sim_options.save_Mach == true
        Mach = NaN(n_step,1); Mach(1) = Velocity(1)/a(1);
    end
    if sim_options.save_Re == true
        Reynolds = NaN(n_step,1); Reynolds(1) = aircraft.geom.c_b_w*Velocity(1)/nu(1);
    end
    
    %%%%%%%%%%%%%%%%%% TIME MARCH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if sim_options.time_integration == true
        tic
    end
    for kk=2:n_step 
        %%%%%%%%%%%%%%%%%%%%%%%% RK4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        dt = t(kk) - t(kk-1);
        t_12 = t(kk-1) + dt/2; % midpoint of the time interval
        alphabeta_dot = [alpha_dot(kk-1);beta_dot(kk-1)]*pi/180; % use previous time derivative (rad/s)
        w = sim_options.disturb.dist_vec(kk-1,:);

        % dx1
        u = control_input(t(:,kk-1),x(:,kk-1),aircraft.control);
        dx1 = FLIGHT_SIM_EOM(t(:,kk-1), x(:,kk-1), u, w, alphabeta_dot, aircraft, sim_options)*dt;

        % dx2
        u = control_input(t_12,x(:,kk-1) + dx1/2,aircraft.control);
        %w  = sim_options.disturb.dist_vec(kk-1,:);
        dx2 = FLIGHT_SIM_EOM(t_12, x(:,kk-1) + dx1/2, u, w, alphabeta_dot, aircraft, sim_options)*dt;
        
        % dx3
        u = control_input(t_12,x(:,kk-1) + dx2/2,aircraft.control);
        %w  = sim_options.disturb.dist_vec(kk-1,:);
        dx3 = FLIGHT_SIM_EOM(t_12, x(:,kk-1) + dx2/2, u, w, alphabeta_dot, aircraft, sim_options)*dt;
        
        % dx4
        u = control_input(t(kk),x(:,kk-1) + dx3,aircraft.control);
        %w  = sim_options.disturb.dist_vec(kk-1,:);
        dx4 = FLIGHT_SIM_EOM(t(kk), x(:,kk-1) + dx3, u, w, alphabeta_dot, aircraft, sim_options)*dt;

        % Integrate -> Eq 3.6-16 Stengel
        x(:,kk) = x(:,kk-1) + (dx1 + 2*dx2 + 2*dx3 + dx4)/6;

        %%%%%%%%%%%%%%%%%%%% STORAGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % compute quantities for storage
        altitude(kk) = -x(9,kk);
        control_vec(:,kk) = control_input(t(kk),x(:,kk),aircraft.control);

        % aerodynamic angles (deg)
        [alpha(kk),beta(kk),Velocity(kk),beta_f(kk)] = body_to_stab(x(:,kk));

         % aero angle derivatives (in deg)
         alpha_dot(kk) = (alpha(kk) - alpha(kk-1))/dt;
         beta_f_dot(kk) = (beta_f(kk) - beta_f(kk-1))/dt;
         beta_dot(kk) = (beta(kk) - beta(kk-1))/dt;

        % attitude
        phi_theta_psi(kk,:) = attitude_from_quats(x(10:13,kk)')';

        % additional optional quantities to be saved
        if sim_options.save_atmos == true
            if altitude(kk) <= sim_options.const_dens_alt
                rho(kk) = rho(1);
                T(kk) = T(1);
                p(kk) = p(1);
                a(kk) = a(1);
                nu(kk) = nu(1);
                g(kk) = g(1);
            else
                [rho(kk), T(kk), p(kk), a(kk), nu(kk), g(kk), ~] = ATMOS_1976(altitude(kk),'US',sim_options.disp_alt_warn);
            end
        end
        if sim_options.save_xdot == true
            x_dot(:,kk) = FLIGHT_SIM_EOM(t(kk), x(:,kk), control_vec(:,kk), w, alphabeta_dot, aircraft, sim_options);

            % accelerations -> morelli 3.52
            ax(kk) = x_dot(1,kk) - x(6,kk)*x(2,kk) + x(5,kk)*x(3,kk) + g(kk)*sin(phi_theta_psi(kk,2));
            ay(kk) = x_dot(2,kk) - x(4,kk)*x(3,kk) + x(6,kk)*x(1,kk) - g(kk)*cos(phi_theta_psi(kk,2))*sin(phi_theta_psi(kk,1));
            az(kk) = x_dot(3,kk) - x(5,kk)*x(1,kk) + x(4,kk)*x(2,kk) - g(kk)*cos(phi_theta_psi(kk,2))*cos(phi_theta_psi(kk,1));
        end
        if sim_options.save_forces == true
            [F_b(:,kk), M_b(:,kk)] = AeroForces(x(:,kk), rho(kk), alphabeta_dot, control_vec(:,kk), aircraft);
        end
        if sim_options.save_coefs == true
            % non dim rates
            p_b(kk) = x(4,kk)*aircraft.geom.b_w/2/Velocity(kk); % roll
            q_b(kk) = x(5,kk)*aircraft.geom.c_b_w/2/Velocity(kk); % pitch
            r_b(kk) = x(6,kk)*aircraft.geom.b_w/2/Velocity(kk); % yaw
            [CL(kk), CD(kk), CY(kk), Cl(kk), Cm(kk), Cn(kk)] = AeroCoefs(alpha(kk)*pi/180, beta(kk)*pi/180, alpha_dot(kk)*pi/180, beta_dot(kk)*pi/180, p_b(kk), q_b(kk), r_b(kk), control_vec(:,kk), aircraft);
            n(kk) = aircraft.propulsion.Kmotor*control_vec(1,kk);
            J(kk) = Velocity(kk)/n(kk)/aircraft.propulsion.d;
            CT(kk) = polyval(aircraft.propulsion.p_CT,J(kk));
            DCL(kk) = change_lift(alpha(kk)*pi/180,J(kk), CL(kk),aircraft);
            CL(kk) = CL(kk) + DCL(kk);
        end
        if sim_options.save_prop_forces == true
            [F_p(:,kk), M_p(:,kk)] = PropForces(x(:,kk), rho(kk), control_vec(1,kk), aircraft.propulsion);
        end

        % non dim numbers
        if sim_options.save_Mach == true
            Mach(kk) = Velocity(kk)/a(kk);
        end
        if sim_options.save_Re == true
            Reynolds(kk) = aircraft.geom.c_b_w*Velocity(kk)/nu(kk);
        end
        if sim_options.save_grav == true
            % body fixed accel due to gravity
            gb(:,kk) = earth_to_body(x(10:13,kk)',[0;0;g(kk)]');
        end


        %%%%%%%%%%%%%%%%%%%% TERMINATE FLIGHT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Check state to terminate integration
        terminate = terminate_flight(t(kk),altitude(kk),alpha(kk),beta(kk),phi_theta_psi(kk,:),sim_options);
        if terminate == true
            break
        end

        %%%%%%%%%%%%%%%%%%%% PLOT INTEGRATED QUANTITIES %%%%%%%%%%%%%%%%%%%
        if floor(kk/sim_options.n_tstep_disp) == kk/sim_options.n_tstep_disp && sim_options.display_integration_time == true
            disp(append('time = ',num2str(t(kk))))
            % plot relevant quantities realtime
            if sim_options.plot_real_time == true
                figure(6969); plot3(x(7,kk),x(8,kk),-x(9,kk),'.b'); grid on; xlabel('$x_f$ (ft)','Interpreter','latex','FontSize',20); ylabel('$y_f$ (ft)','Interpreter','latex','FontSize',20); zlabel('$-z_f$ (ft)','Interpreter','latex','FontSize',20); title('Integrated Positon Realtime','FontSize',20); axis equal; hold on; drawnow
            end
        end
    end
    if sim_options.plot_real_time == true
        hold off
    end
    if sim_options.time_integration == true
        toc
    end
    %%%%%%%%%%%%%%% ASSEMBLE RESULTS INTO DATA STRUCTURE %%%%%%%%%%%%%%%%%%
    % REMOVE NaN VALUES 
    sim_data_out.altitude = rmmissing(altitude)';
    sim_data_out.control_vec = control_vec(:,1:kk)';
    sim_data_out.Velocity = rmmissing(Velocity)';
    sim_data_out.alpha  = rmmissing(alpha)';
    sim_data_out.beta_f = rmmissing(beta_f)';
    sim_data_out.beta = rmmissing(beta)';
    sim_data_out.alpha_dot = rmmissing(alpha_dot)';
    sim_data_out.beta_f_dot = rmmissing(beta_f_dot)';
    sim_data_out.beta_dot = rmmissing(beta_dot)';

    % state vector
    sim_data_out.t_out = t(1:kk);
    sim_data_out.X_out = x(:,1:kk);

    sim_data_out.uvw = [x(1,1:kk); x(2,1:kk); x(3,1:kk)]'; %ft/s (kts 1.688)
    sim_data_out.pqr = [x(4,1:kk); x(5,1:kk); x(6,1:kk)]'*180/pi; % deg/s
    sim_data_out.XYZ = [x(7,1:kk); x(8,1:kk); x(9,1:kk)]';
    sim_data_out.es_out = [x(10,1:kk); x(11,1:kk); x(12,1:kk); x(13,1:kk)]';

    % compute attitude
    sim_data_out.phi_theta_psi = (phi_theta_psi(1:kk,:))*180/pi;  % convert to deg
    sim_data_out.phi = sim_data_out.phi_theta_psi(1:end,1)';
    sim_data_out.theta = sim_data_out.phi_theta_psi(1:end,2)';
    sim_data_out.psi = sim_data_out.phi_theta_psi(1:end,3)';

    % flight path angle
    sim_data_out.gamma = sim_data_out.theta - sim_data_out.alpha;

    % control vector
    sim_data_out.delta_T = sim_data_out.control_vec(:,1); 
    sim_data_out.delta_e = sim_data_out.control_vec(:,2)*180/pi;
    sim_data_out.delta_a = sim_data_out.control_vec(:,3)*180/pi;
    sim_data_out.delta_r = sim_data_out.control_vec(:,4)*180/pi;
    sim_data_out.delta_f = sim_data_out.control_vec(:,5)*180/pi;


    % optional outputs
    if sim_options.save_atmos == true
        sim_data_out.rho = rmmissing(rho)';
        sim_data_out.T = rmmissing(T)';
        sim_data_out.p = rmmissing(p)';
        sim_data_out.a = rmmissing(a)';
        sim_data_out.nu = rmmissing(nu)';
        sim_data_out.g = rmmissing(g)';
    end
    if sim_options.save_grav == true
        sim_data_out.gb = gb(:,1:kk);
    end
    if sim_options.save_xdot == true
        sim_data_out.x_dot = x_dot(:,1:kk)';
        sim_data_out.ax = ax(1:kk)';
        sim_data_out.ay = ay(1:kk)';
        sim_data_out.az = az(1:kk)';
    end
    
    if sim_options.save_forces == true
        sim_data_out.F_b = rmmissing(F_b)';
        sim_data_out.M_b  = rmmissing(M_b)';
    end
    if sim_options.save_prop_forces == true
        sim_data_out.F_p = rmmissing(F_p');
        sim_data_out.M_p = rmmissing(M_p');
    end
    
    if sim_options.save_Mach == true
        sim_data_out.mach = rmmissing(Mach)';
    else
        sim_data_out.mach = zeros(length(t),1);
    end
    if sim_options.save_Re == true
        sim_data_out.Reynolds = rmmissing(Reynolds)';
    end
    if sim_options.save_coefs == true
        sim_data_out.CL = rmmissing(CL);
        sim_data_out.CY = rmmissing(CY);
        sim_data_out.CD = rmmissing(CD);
        sim_data_out.Cl = rmmissing(Cl);
        sim_data_out.Cm = rmmissing(Cm);
        sim_data_out.Cn = rmmissing(Cn);
        sim_data_out.p_b = rmmissing(p_b);
        sim_data_out.q_b = rmmissing(q_b);
        sim_data_out.r_b = rmmissing(r_b);
        sim_data_out.J = rmmissing(J);
        sim_data_out.n = rmmissing(n);
        sim_data_out.CT = rmmissing(CT);
        sim_data_out.DCL = rmmissing(DCL);
    end

    % load factors
    if sim_options.save_load == true && sim_options.save_xdot
        sim_data_out.nx_g = ( ax(1:kk))/aircraft.mass.W;
        sim_data_out.ny_g = ( ay(1:kk))/aircraft.mass.W;
        sim_data_out.nz_g = ( -az(1:kk))/aircraft.mass.W;
    end

    % earth information
    if sim_options.earth.ellipsoidal_earth == true && sim_options.save_xdot == true
        sim_data_out.PHI_PSI_H = flat_to_ellipsodial(sim_data_out.t_out, sim_data_out.x_dot(7:9,:), sim_options.earth.PHI_PSI_H_0);
    elseif sim_options.earth.spherical_earth == true && sim_options.save_xdot == true
        sim_data_out.PHI_PSI_H = flat_to_spherical(sim_data_out.t_out, sim_data_out.x_dot(7:9,:), sim_options.earth.PHI_PSI_H_0);
    else
        sim_data_out.PHI_PSI_H = sim_options.earth.PHI_PSI_H_0;
    end
end