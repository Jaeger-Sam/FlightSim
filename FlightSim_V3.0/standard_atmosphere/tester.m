clc;clear all;close all
%   Test Program for atmosphere with plots
%       Modifed to compare to ATMOS from McCormick
%
%       Several errors in the atmo code for US units
%           checked T, P, rho, c, nu agaist ATMOS function for both US and
%           SI units
%
%       code orginally from 1976 standard atmosphere on MATLAB file exchange
%   
%   Sam Jaeger
%   Last modified: 9/26/2024

% inputs
%altitude = 200000; % ft (final altitude to compute)
altitude = 80000; % m
km_to_ft = 3281;
psi_to_psf =  144;
units = 'SI';

tic
if units == 'US'
    [Z, Z_L, Z_U, T, P, rho, c, g, mu, nu, k, n, n_sum] = atmo(altitude/km_to_ft,0.5,2);
    Z = Z*1000; % convert to actual feet, not thousands of feet
    Z_L = Z_L*1000;
    for ii=1:length(Z)
        if Z(ii) <80000
            [rho_1(ii), T_1(ii), p_1(ii), a_1(ii), nu_1(ii)] = ATMOS(Z(ii),'US');
            Z_1(ii) = Z(ii);
        end
    end
else
    [Z, Z_L, Z_U, T, P, rho, c, g, mu, nu, k, n, n_sum] = atmo(altitude/1000,0.5,1);
    Z = Z*1000; % convert to actual feet, not thousands of feet
    Z_L = Z_L*1000;
    for ii=1:length(Z)
        if Z(ii) < 24200
            [rho_1(ii), T_1(ii), p_1(ii), a_1(ii), nu_1(ii)] = ATMOS(Z(ii),'SI');
            Z_1(ii) = Z(ii);
        end
    end

    beta = 0.1378/1000; %1/7163; %book % candler =  0.1378 /1000; % 1/km to 1/m
    rho_0 = 1.225;% 0.12492; % kg/m^3;
    rho_exp = rho_0*exp(-Z*beta);
end
toc
%% Plots

if units == 'US'
    figure(10)
    plot(T,Z,'.', T_1,Z_1,'o')
    ylabel('Z [ft]')
    xlabel('T [R]')
    legend('1976 atmo','McCormick ATMOS')
    grid on
    
    figure(11)
    semilogx(P*psi_to_psf, Z,'.', p_1, Z_1,'o')
    ylabel('Z [ft]')
    xlabel('P [lb/ft^2]')
    legend('1976 atmo','McCormick ATMOS')
    grid on
    
    figure(12)
    semilogx(rho/32.2,Z,'.', rho_1,Z_1,'o')
    ylabel('Z [ft]')
    xlabel('rho [slugs/ft^3]')
    legend('1976 atmo','McCormick ATMOS')
    grid on
    
    figure(13)
    plot(c,Z_L,'.', a_1, Z_1,'o')
    ylabel('Z [ft]')
    xlabel('c[ft/s]')
    legend('1976 atmo','McCormick ATMOS')
    grid on
    
    
    % figure
    % plot(g,Z)
    % ylabel('Z')
    % xlabel('g[m/s^2]')
    
    % figure
    % plot(mu,Z_L)
    % ylabel('Z')
    % xlabel('\mu')
    
    figure(14)
    semilogx(nu,Z_L,'.', nu_1, Z_1,'o')
    ylabel('Z [ft]')
    xlabel('\nu [ft^2/s]')
    legend('1976 atmo','McCormick ATMOS')
    grid on
    
    
    % figure
    % plot(k,Z_L)
    % ylabel('Z')
    % xlabel('k')
    % 
    % figure
    % semilogx(n,Z_U)
    % ylabel('Z')
    % xlabel('n')
    % 
    % figure
    % semilogx(n_sum,Z_U)
    % ylabel('Z')
    % xlabel('n_sum')
else
    figure(10)
    plot(T,Z,'.', T_1,Z_1,'o')
    ylabel('Z [m]')
    xlabel('T [K]')
    legend('1976 atmo','McCormick ATMOS')
    grid on
    
    figure(11)
    semilogx(P, Z,'.', p_1, Z_1,'o')
    ylabel('Z [m]')
    xlabel('P [Pa]')
    legend('1976 atmo','McCormick ATMOS')
    grid on
    
    figure(12)
    semilogx(rho,Z,'.', rho_1,Z_1,'o', rho_exp,Z,'*')
    ylabel('Z [m]')
    xlabel('rho [kg/m^3]')
    legend('1976 atmo','McCormick ATMOS','Exponential')
    grid on
    
    figure(13)
    plot(c,Z_L,'.', a_1, Z_1,'o')
    ylabel('Z [m]')
    xlabel('c[m/s]')
    legend('1976 atmo','McCormick ATMOS')
    grid on
    
    
    % figure
    % plot(g,Z)
    % ylabel('Z')
    % xlabel('g[m/s^2]')
    
    % figure
    % plot(mu,Z_L)
    % ylabel('Z')
    % xlabel('\mu')
    
    figure(14)
    semilogx(nu,Z_L,'.', nu_1, Z_1,'o')
    ylabel('Z [m]')
    xlabel('\nu [m^2/s]')
    legend('1976 atmo','McCormick ATMOS')
    grid on
    
    
    % figure
    % plot(k,Z_L)
    % ylabel('Z')
    % xlabel('k')
    % 
    % figure
    % semilogx(n,Z_U)
    % ylabel('Z')
    % xlabel('n')
    % 
    % figure
    % semilogx(n_sum,Z_U)
    % ylabel('Z')
    % xlabel('n_sum')
end