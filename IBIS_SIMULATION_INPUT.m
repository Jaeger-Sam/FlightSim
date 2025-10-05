%%%%%%%%%%%%%%%%%%%%%%% SIMULATION INPUT FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha_deg_0 = 10; % angle of attack
beta_deg_0 = 0;
gamma_deg_0 =  0; % flight path angle 
alt0 = 100; % ft
V0 = 50; % ft/s

%%%%%%%%%%%%%%% TIME %%%%%%%%%%%%%%%%%%%%%%%%%%
t_start = 0; % s
t_end = 80; % s
variable_dt = false; % not implemented in RK4
dt = 0.05; % s, will use this value if 'variable_dt' == false
n_tstep_disp = 1000; % number of timesteps to display trajectory time

%%%%%%%%%%%% INITAL CONDITONS %%%%%%%%%%%%%%%%%

% inital position with respect to flat earth ground (ft)
x0 = 0;
y0 = 0;
z0 = -alt0; % negative is positive altitude

% inital velocity (ft/s), approximate relation for small angles
u0 = cosd(alpha_deg_0)*cosd(beta_deg_0)*V0; 
v0 = sind(beta_deg_0)*V0;
w0 = sind(alpha_deg_0)*cosd(beta_deg_0)*V0;

% angular velocities start (deg/s)
p0 = 0;
q0 = 0;
r0 = 0;

% attitude start (deg)
theta0 = gamma_deg_0 + alpha_deg_0; % pitch attitude
phi0 = 0; % bank angle
psi0 = 0; % heading

% Latitude, Longitude, and Altitude (deg)
Lat_0 = 44.97; % Inital Latitude
Long_0 = -93.23; % Inital Longitude (negative west)

% corrections for lat, long outputs. If both false, no correction are made
ellipsoidal_earth = true;
spherical_earth = false;

%%%%%%%%%%%%%%% GRAVITY/ATMOSPHERE MODEL %%%%%%%%%%%%%%%%%
WGS84_gravity = false; 
const_grav_alt = 100000; % ft, fix gravity below this altitude
const_dens_alt = 3000; % ft, fix density below this altitude

%%%%%%%%%%%%%%% DISTURBANCES %%%%%%%%%%%%%%%%%%
wind_mag_fts = 0; % ft/s
wind_heading_deg = 0;
thermal_speed_fts = 0; %ft/s

gusts = true; % turn on/off turbulence

% standard deviation for Von Karman turbulence model
sigma_u = .1; 
sigma_v = .1;
sigma_w = .1;


%%%%%%%%%%%%%%%%% CONSTRAINTS %%%%%%%%%%%%%%%%%
% Can only constrain one at a time
constrain_roll = false;
constrain_pitch = false;
constrain_yaw = false;

pure_rolling = false;
pure_pitching = true;
pure_yawing = false;

%%%%%%%%%%%%% DEPARTURE %%%%%%%%%%%%%%%%%%%%%%%%
% hitting ground: z_f = 0 automatically
alpha_depart_deg = [-40,40]; % angle of attack
beta_depart_deg = [-20,20]; % true sideslip
theta_depart_deg = [-100,100]; % pitch attitude
phi_depart_deg = [-180,180]; % bank angle

%%%%%%%%%%%%% INTEGRATION OPTIONS %%%%%%%%%%%%%%%
time_integration = true;
plot_real_time = true; % display a plot of XYZ position
%int_options =
%odeset('RelTol',1e-4,'Stats','on','OutputFcn',@odeplot,'Events',@(t,x)ground_hit(t,x)); % for ode45
display_integration_time = true;
disp_alt_warn = false;

%%%%%%%%%%%%% POST PROCESSING %%%%%%%%%%%%%%%%%%%
% save data
%   default: t,x,control_vec,alpha,beta,alpha_dot,beta_dot,phi_theta_psi
save_atmos = true;
save_xdot = true;
save_forces = true;
save_Mach = false;
save_Re = false;
save_grav = true;
save_load = true;
save_coefs = true;
save_prop_forces = true;

% Plots
trajectory_3d_plot = true;
state_vec_plot = true;
mach_time_plot = false;
Reynolds_time_plot = false;
flight_path_angle_plot = true;
aero_angle_plot = true;
controls_input_plot = true;
load_factor_plot = false;
attitude_movie_plot = false;

% Attitude Movie Options
model_info_file = 'C:\Users\jaege\MATLAB Drive\tools\3d_animations_MATLAB\3d_models\f16_3d_model.mat';
frame_sample_time = 0.02;
speedx = 1; 
isave_movie = 0;
movie_file_name = 'IBIS_dive_turn.mp4';