% Input file for IBIS
%   Full Scale Ultra Stick

%%%%%%%%%%%%%%%%%%%%%%% AIRCRAFT INPUT FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m = 10; %kg

%%%%%%%% MASS PROPERTIES %%%%%%%%%%%%
W = m*2.205; % lb % weight of aircraft

% Inertia (slug*ft^2)
%   conversion from kgm^2 to slugs*ft^2 == 0.737562
%   Using inertias from orginal FASER paper
I_xx = 0.496;
I_yy = 0.656; % calc from swing test: 12.84 kg*m^2
I_zz = 1.164;
I_xy = 0;
I_xz = 0.560;
I_yz = 0;

% cg from orgin (ft)
%   unused?
xcg = 1.125;
ycg = 0;
zcg = 0.25;

%%%%%%% GEOMETRY PROPERTIES %%%%%%%%%%%
S_w = 7.92; % ft^2 % main wing area
b_w = 6.33; % ft % wing span
c_b_w = 1.375; % ft % mean aero chord
xnp = 1.33; % ft % neutral point, unused?

geom_from_stl_file = false;
stlfile = 'HGV_0.0.stl';

% specify rotation matrix if geomtry is not is standard flight dynamic
% coordinate system.
geom_rot_mat = [1,0,0;
                0,1,0;
                0,0,1];


%%%%%%% AERODYNAMIC COEFFICIENTS %%%%%

% load aero data... for lookup table
aero_data_mat = false; % if true then will load the mat file specified
%aero_data_matfile = 'NA_Aero_Coefs_delta_alpha_sweep.mat';

% interpolation method for matlabs built in interp2 function: 
%   'linear' or 'nearest' or 'cubic' or 'makima' or 'spline'
%method = 'linear'; 
                
% Newtoian Aerodynamics
Newtonian_Aerodynamics = false;
%plot_cp = false; % will plot the cp distribution of Newtoian Aerodynamics
%modNA = false;

unsteady = true;
nonlin = true;

% Coefficient matrix
S = [4.2015, 0, 0.4603, 0, 0, 1.4553;
     0, -0.1615, 0, -0.0545, 0.1123, 0;
     0, -0.0169, 0, -0.3275, 0.0022, 0;
     -1.45, 0, -1.3174, 0, 0, -1.2217;
     0, 0.1049, 0, 0.0276, -0.0779, 0]; % linear terms matrix (at firewall)

Sn = [0,0;
    -0.1285,-3.789;
    -0.4849, -0.4368;
     0,0;
    -0.1008, 0.5602]; % nonlinear terms matrix (at firewall)

Su = [0.0285,0,0,3.05,0;
     0,0,-0.0494,0,0.1274;
     0,0,-0.4437,0,0;
     0.0619,0,0,-7.96,0;
     0,0,0,0,-0.0645]; % unsteady terms matrix

CD_mat = diag([1.419,0.088,0.1363,0.3250,0.4583,0.2457]); % diagonal of quadratic drag terms

C_D_0 = 0.0353; % drag @ zero alpha
C_m_0 = 0.004; % pitching moment @ zero alpha
C_L_0 = 0.003;% lift @ zero alpha
alpha_stall  = 15*pi/180; % AoA at stall
alpha_b = alpha_stall + 5*pi/180;
M= 20; % blending function for poststall model


%%%%%%%  PROPULSION PROPERTIES %%%%%%%

% ==0 for glider, ==1 for prop, ==2 for jet, ==3 for rocket
glider_prop_jet = 2; 
z_prop = 0; 
x_prop = 1.25;
x_wing_CP = 1.55; % wing cp behind prop

Kmotor = 8000/60;
d = 1.5; % ft
p_CT = [-0.1714,-0.0105,0.1117]; % polynomial coef of thrust coefficient

prop_moments = false;
p_f = [-0.2008, 0.8512, 0.9744]; % polynomial coef of function for normal force

%max_thrust = 0; % lb


%%%%%%% CONTROL PROPERTIES %%%%%%%%%%%

full_state_info = false;
input_type = 3; % 0 for uncontrolled, 1 for step input, 2 for time history, 3 for feedback

% step input parameters
tstep = 0; % time for step to be applied
amp0 = []; % inital control amplitudes
amp = []; % step control control amplitudes

% input time history
thist = 0:0.01:100; % time history
uhist = zeros(5,length(thist));
uhist(1,:) = 1*ones(1,length(thist)); % input time history
uhist(2,:) = -25*pi/180*ones(1,length(thist));
ttol = 0.01; % time matching tolerance

% feedback, altitude holder + bank control
C = [0,0,0,  0,0,0,  0,0,-1, 0,0,0,0]; % sense negative zf
K  = [0,0;
    1/70,0;
    0, 1/10;
    0, 0;
    0, 0]; % Kp_e = 1/70 (longitudinal)
yref = [100;0*pi/180]; % [altitude to maintain (ft); bank to hold (rad) ]
ustat = [1;-10*pi/180;0;0;0];

% control size - currenly only for one control surface
%xcontrol = 0.25*c_b_w;
%ycontrol = 0.02*b_w;

max_deflect = [1, 0;
               25 , -30;
               19, -25;
               30, -30;
               30, 0]; % matrix of maximum deflections in degrees
%max_deflect(1,1) = 1; % for throttle

elevator = true;
aileron = true;
rudder = true;


%%%%%%%% GYRO PROPERTIES %%%%%%%%%%%%%

% gyroscopic moments (ft-lb)
h_x = 0;
h_y = 0;
h_z = 0;