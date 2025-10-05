% SIMULATION_PLOTS.fcn plots output data from SIMULATE_FLIGHT.fcn.
% Options to plot are specified in sim_options. Plots include 3d flight
% path, state vector plot, mach vs time, flight path angle vs time, 
% alpha/beta vs time, control input vs time, and 3d animation. There are no
% outputs.
%
% SIMULATION_PLOTS(sim_data_out,aircraft,sim_options)
%
% INPUTS:
%   sim_dat_out: Data structure from RK4.fcn
%   sim_options: Data structure from initalize_sim.fcn
%   aircraft: Data structure from initalize_sim.fcn
%
% OUTPUTS:
%
% Sam Jaeger
% jaege246@umn.edu
% 1/5/2024
%   Revised: added alpha, beta plots as well as control input plots
%   Revised: 1/25/2024 added aircraft data structure
%   Revised: 3/4/2024
%   Revised: 10/20/2024
%   Revised: 10/1/2025, Updated 3d animation section

function SIMULATION_PLOTS(sim_data_out,aircraft,sim_options)
    
    if sim_options.trajectory_3d_plot == true
        figure(100)
        plot3(sim_data_out.XYZ(:,1),sim_data_out.XYZ(:,2),-sim_data_out.XYZ(:,3))
        title('3D flight path','FontSize',20,'Interpreter','latex')
        xlabel('$X$ (ft)','FontSize',20,'Interpreter','latex')
        ylabel('$Y$ (ft)','FontSize',20,'Interpreter','latex')
        zlabel('$-Z$ (ft)','FontSize',20,'Interpreter','latex')
        axis equal
        grid on
    end

    if sim_options.state_vec_plot == true
        figure(10)
        subplot(4,1,1)
        plot(sim_data_out.t_out,sim_data_out.uvw(:,1) ,sim_data_out.t_out,sim_data_out.uvw(:,2), sim_data_out.t_out, sim_data_out.uvw(:,3)); 
        grid on; 
        xlabel('time (s)','FontSize',15,'Interpreter','latex'); 
        ylabel('velocity (ft/s)','FontSize',15,'Interpreter','latex'); 
        legend('u','v','w','location','eastoutside','Interpreter','latex')
        
        subplot(4,1,2)
        plot(sim_data_out.t_out, sim_data_out.pqr(:,1), sim_data_out.t_out, sim_data_out.pqr(:,2), sim_data_out.t_out, sim_data_out.pqr(:,3)); 
        grid on; 
        xlabel('time (s)','FontSize',15,'Interpreter','latex'); 
        ylabel('angular velocity (deg/s)','FontSize',15,'Interpreter','latex'); 
        legend('$p$','$q$','$r$','location','eastoutside','Interpreter','latex')
        
        subplot(4,1,3)
        plot(sim_data_out.t_out,sim_data_out.XYZ(:,1) , sim_data_out.t_out,sim_data_out.XYZ(:,2), sim_data_out.t_out,-sim_data_out.XYZ(:,3)); 
        grid on; 
        xlabel('time (s)','FontSize',15,'Interpreter','latex'); 
        ylabel('Location (ft)','FontSize',15,'Interpreter','latex'); 
        legend('$x_f$','$y_f$','$-z_f$','location','eastoutside','Interpreter','latex')
        
        subplot(4,1,4)
        plot(sim_data_out.t_out,sim_data_out.phi ,   sim_data_out.t_out,sim_data_out.theta, sim_data_out.t_out,sim_data_out.psi); 
        grid on; 
        xlabel('time (s)','FontSize',15,'Interpreter','latex'); 
        ylabel('angle (deg) ','FontSize',15,'Interpreter','latex'); 
        ylim([-180 180])
        legend('$\phi$','$\theta$', '$\psi$', 'location','eastoutside','Interpreter','latex')
    end

    if sim_options.mach_time_plot == true
        figure(11)
        plot(sim_data_out.t_out, sim_data_out.mach)
        xlabel('time (s)','FontSize',20,'Interpreter','latex')
        ylabel(' Mach ','FontSize',20,'Interpreter','latex')
        grid on
        title(' Mach Number vs Time','FontSize',20,'Interpreter','latex')
    end

    if sim_options.Reynolds_time_plot == true
        figure(11)
        plot(sim_data_out.t_out, sim_data_out.Reynolds)
        xlabel('time (s)','FontSize',20,'Interpreter','latex')
        ylabel(' Re ','FontSize',20,'Interpreter','latex')
        grid on
        title(' Reynolds Number vs Time','FontSize',20,'Interpreter','latex')
    end
    
    if sim_options.flight_path_angle_plot == true
        figure(12)
        plot(sim_data_out.t_out, sim_data_out.gamma)
        xlabel('time (s)','FontSize',20,'Interpreter','latex')
        ylabel(' $\gamma$ (deg) ','FontSize',20,'Interpreter','latex')
        grid on
        title('Flight Path Angle vs Time','FontSize',20,'Interpreter','latex')
    end

    if sim_options.aero_angle_plot == true
        figure(13)
        subplot(2,1,1)
        plot(sim_data_out.t_out, sim_data_out.alpha, sim_data_out.t_out, sim_data_out.beta)
        xlabel('time (s)','FontSize',20,'Interpreter','latex')
        ylabel(' $\alpha$, $\beta$ (deg)','FontSize',20, 'Interpreter','latex')
        grid on
        legend(' $\alpha$',' $\beta$','FontSize',20, 'Interpreter','latex','Location','southeast')
        
        subplot(2,1,2)
        plot(sim_data_out.t_out, sim_data_out.alpha_dot, sim_data_out.t_out, sim_data_out.beta_dot)
        xlabel('time (s)')
        ylabel(' $\dot{\alpha}$, $\dot{\beta}$ (deg/s)','FontSize',20,'Interpreter','latex')
        grid on
        legend(' $\alpha$',' $\beta$', 'Interpreter','latex','FontSize',20,'Location','southeast')
    end

    if sim_options.load_factor_plot == true
        figure(14)
        title('load factor')
        plot(sim_data_out.t_out, sim_data_out.nz_g)
        xlabel('time (s)','FontSize',20,'Interpreter','latex')
        ylabel('$n_z$ (g)','FontSize',20,'Interpreter','latex')
        grid on
    end

    if sim_options.controls_input_plot == true
        figure(15)
        title('Control Inputs')
        subplot(4,1,1)
        plot(sim_data_out.t_out, sim_data_out.delta_T)
        xlabel('time (s)','FontSize',20,'Interpreter','latex')
        ylabel('$\delta_T$','FontSize',20,'Interpreter','latex')
        grid on

        subplot(4,1,2)
        plot(sim_data_out.t_out, sim_data_out.delta_e)
        xlabel('time (s)','FontSize',20,'Interpreter','latex')
        ylabel('$\delta_e$ (deg)','FontSize',20,'Interpreter','latex')
        grid on

        subplot(4,1,3)
        plot(sim_data_out.t_out, sim_data_out.delta_a)
        xlabel('time (s)','FontSize',20,'Interpreter','latex')
        ylabel('$\delta_a$ (deg)','FontSize',20,'Interpreter','latex')
        grid on

        subplot(4,1,4)
        plot(sim_data_out.t_out, sim_data_out.delta_r)
        xlabel('time (s)','FontSize',20,'Interpreter','latex')
        ylabel('$\delta_r$ (deg)','FontSize',20,'Interpreter','latex')
        grid on
    end
    
    % visulaization
    if sim_options.attitude_movie_plot == true
        close all % matlab 2025 has issues when multiple figures are open
        controls_deflection_deg = [sim_data_out.delta_a, ...
            -sim_data_out.delta_a, ...
            zeros(length(sim_data_out.t_out),1), ...
            zeros(length(sim_data_out.t_out),1), ...
            sim_data_out.delta_r, ...
            sim_data_out.delta_e,...
            sim_data_out.delta_e]; % for f-16 model

        aircraft_3d_animation(sim_options.model_info_file,...
            sim_data_out.psi, ...            Heading angle [deg]
            sim_data_out.theta, ...              Pitch angle [deg]
            sim_data_out.phi, ...               Roll angle [deg]
            sim_data_out.delta_a/(aircraft.control.max_deflect(3,1)*180/pi), ... Roll  stick command [-1,+1] [-1 -> left,            +1 -> right]
            -sim_data_out.delta_e/(aircraft.control.max_deflect(2,2)*180/pi), ... Pitch stick command [-1,+1] [-1 -> full-back stick, +1 -> full-fwd stick]
            sim_data_out.alpha, ...    AoA [deg]
            sim_data_out.beta, ...  AoS [deg]
            sim_data_out.gamma, ...   Flight path angle [deg]
            sim_data_out.mach, ...                   Mach number
            sim_data_out.altitude/3.28084, ...            Altitude [ft]
            sim_data_out.nz_g,  ...                  Vertical load factor [g]
            controls_deflection_deg, ...Flight control deflection (each column is a control surface)
            sim_options.frame_sample_time, ...      Sample time [sec]
            sim_options.speedx, ...                 Reproduction speed
            sim_options.isave_movie, ...            Save the movie? 0-1
            sim_options.movie_file_name);           % Movie file name
    end
end