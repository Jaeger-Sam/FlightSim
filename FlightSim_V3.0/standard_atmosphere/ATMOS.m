% ATMOS calculates atmospheric properties based on altitude. Based on
% McCormick Aerodynamics, Aeronautics, and Flight Mechanics.
%
% [rho, T, p, a, nu] = ATMOS(h,units)
%
% INPUTS:
%   h: altitude (ft or m)
%   units: 'SI' or 'US'
%
%
% OUTPUTS:
%   rho: density of air (slugs/ft^3 or kg/m^3)
%   T: temperature (R or K)
%   p: pressure (lb/ft^2 or Pa)
%   a: speed of sound (ft/s or m/s)
%   nu: kinematic viscosity (mu/rho) (ft^2/s or m^2/s)
%
% Sam Jaeger
% jaege246@umn.edu
% Revised: 5/18/2022
% Revised: 1/18/2024

function [rho, T, p, a, nu] = ATMOS(h,units)

if units == 'SI'
    
    p_0 = 101300; %(N/m^2)
    rho_0 = 1.225; %(kg/m^3)
    T_0 = 288.16; %(K)
    a_0 = 340.3; %(m/s)
    R = 286.97; %(m^2/s^2/K)
    delta_c = 0.225;
    theta_c = 0.752;
    g = 9.81; %(m/s^2)
    k = 1.4; % k=cp/cv (ratio of specific heats for air;
    
    if h < 11000 %below 11km or 36,000'
        %theta = T/T_0;
        %delta = p/p_0;
        %sigma = rho/rho_0;
        
        theta = 1 - 0.0226*(h/1000);
        
        delta = theta^5.256;
        sigma = theta^4.256;
        
        p = delta*p_0; % pressure at altitude
        rho = sigma*rho_0; % density at altitude
        T = theta*T_0; % temperature
        a = (k*R*T)^.5; %speed of sound
        
        A0 = 1.5723;
        A1 = 8.73065e-2;
        A2 = -1.18412e-2;
        A3 = 1.16978e-3;
        A4 = -5.27207e-5;
        A5 = 1.22466e-6;
        A6 = -1.369780e-8;
        A7 = 5.94238e-11;
        
        h=h*3.281; % convert to ft for nu calculation
        nu = A0 + A1*(h/1000) + A2*(h/1000)^2 + A3*(h/1000)^3 + A4*(h/1000)^4 + A5*(h/1000)^5 + A6*(h/1000)^6 + A7*(h/1000)^7;
        nu = (nu/(10^4))*0.3048^2;
        
    elseif h < 24.4e3
        h_c = 11000;
        T_c = theta_c*T_0;        
        
        delta = delta_c*exp( (-g/(R*T_c))*(h - h_c));
        sigma = delta / theta_c;
        
        p = delta*p_0; % pressure at altitude
        rho = sigma*rho_0; % density at altitude
        T = T_c; % temperature
        a = (k*R*T)^.5; %speed of sound
        
        A0 = 1.5723;
        A1 = 8.73065e-2;
        A2 = -1.18412e-2;
        A3 = 1.16978e-3;
        A4 = -5.27207e-5;
        A5 = 1.22466e-6;
        A6 = -1.369780e-8;
        A7 = 5.94238e-11;
        
        h=h*3.281; % convert to ft for nu calculation
        nu = A0 + A1*(h/1000) + A2*(h/1000)^2 + A3*(h/1000)^3 + A4*(h/1000)^4 + A5*(h/1000)^5 + A6*(h/1000)^6 + A7*(h/1000)^7;
        nu = (nu/(10^4))*0.3048^2;
    else
        error('altitude in excess of 80,000 ft / 24.2 km')
    end
        
elseif units == 'US'
    p_0 = 2116; %(psf)
    rho_0 = 0.002377; %(slugs/ft^3)
    T_0 = 518.7; %(R)
    a_0 = 1116; %(ft/s)
    R = 1716; %(ft^2/s^2/R)
    delta_c = 0.225;
    theta_c = 0.752;
    g = 32.2; %(ft/s^2)
    k = 1.4; % k=cp/cv (ratio of specific heats for air;
    
    if h < 36000 %below 11km or 36,000'
        %theta = T/T_0;
        %delta = p/p_0;
        %sigma = rho/rho_0;
        
        theta = 1 - 0.00688*(h/1000);
        
        delta = theta^5.256;
        sigma = theta^4.256;
        
        p = delta*p_0; % pressure at altitude
        rho = sigma*rho_0; % density at altitude
        T = theta*T_0; % temperature
        a = (k*R*T)^.5; %speed of sound
        
        A0 = 1.5723;
        A1 = 8.73065e-2;
        A2 = -1.18412e-2;
        A3 = 1.16978e-3;
        A4 = -5.27207e-5;
        A5 = 1.22466e-6;
        A6 = -1.369780e-8;
        A7 = 5.94238e-11;
        
        nu = A0 + A1*(h/1000) + A2*(h/1000)^2 + A3*(h/1000)^3 + A4*(h/1000)^4 + A5*(h/1000)^5 + A6*(h/1000)^6 + A7*(h/1000)^7;
        nu = (nu/(10^4));
        
    elseif h < 80000
        h_c = 36000;
        T_c = theta_c*T_0; 
        
        delta = delta_c*exp( (-g/(R*T_c))*(h - h_c));
        sigma = delta / theta_c;
        
        p = delta*p_0; % pressure at altitude
        rho = sigma*rho_0; % density at altitude
        T = T_c; % temperature
        a = (k*R*T)^.5; %speed of sound
        
        A0 = 1.5723;
        A1 = 8.73065e-2;
        A2 = -1.18412e-2;
        A3 = 1.16978e-3;
        A4 = -5.27207e-5;
        A5 = 1.22466e-6;
        A6 = -1.369780e-8;
        A7 = 5.94238e-11;
        
        nu = A0 + A1*(h/1000) + A2*(h/1000)^2 + A3*(h/1000)^3 + A4*(h/1000)^4 + A5*(h/1000)^5 + A6*(h/1000)^6 + A7*(h/1000)^7;
        nu = (nu/(10^4));
    else
        error('altitude in excess of 80,000 ft / 24.2 km')
    end
    
else
    error('units of SI or US must be specified')
end
    

end