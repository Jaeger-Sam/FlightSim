% ATMOS calculates atmospheric properties based on altitude. Based on
% 1976 atmosphere. Calls revised function atmo.fcn and ATMOS.fcn. 
% SI units of m-s-kg, US units of ft-s-slug. 
% Same function call as ATMOS.fcn.
% The variable input argument can supress altitude warnings when the
% atmospheric model is likely incorrect (above 100 km and below 0 km).
%
% [rho, T, p, a, nu, g, mu] = ATMOS_1976(h,units,suppress_warnings)
%
% INPUTS:
%   h: altitude (ft or m)
%   units: 'SI' or 'US'
% OPTIONAL INPUT:
%   varargin: ==0, supress altitude warnings
%             ==1, show altitude warnings
%
% OUTPUTS:
%   rho: density of air (slugs/ft^3 or kg/m^3)
%   T: temperature (R or K)
%   p: pressure (lb/ft^2 or Pa)
%   a: speed of sound (ft/s or m/s)
%   nu: kinematic viscosity (mu/rho) (ft^2/s or m^2/s)
%   g: acceleration due to gravity (ft/s^2 or m/s^2)
%   mu: dynamic viscosity 
%
% Sam Jaeger
% jaege246@umn.edu
% 1/18/2024
%   Modified: 9/27/2024
%       logic to turn off warning

function [rho, T, p, a, nu, g, mu] = ATMOS_1976(h,units,varargin)
    if nargin<3
        varargin = {1};
    end
    if units == 'US'
        km_to_ft = 3281;
        alt_km = h/km_to_ft;
        if alt_km > 100 && cell2mat(varargin) == 1
            warning('h>100km Too high to use contiuum assumption.')
            alt_km = round(alt_km,1);
        else
            alt_km = round(alt_km,1);
        end
        psi_to_psf =  144;

        if alt_km == 0
            [rho, T, p, a, nu] = ATMOS(h,units);
            g = 32.2;
            mu = nu*rho;
        elseif alt_km < 0 && cell2mat(varargin) == 1
            warning('Altitude below sea level!')
            [rho, T, p, a, nu] = ATMOS(0,units);
            g = 32.2;
            mu = nu*rho;
        elseif alt_km < 0 
            [rho, T, p, a, nu] = ATMOS(0,units);
            g = 32.2;
            mu = nu*rho;
        else
            [~, ~, ~, T, p, rho_1, a, g, mu, nu, ~, ~, ~] = atmo(alt_km,alt_km,2);
            T = T(2);
            p = p(2)*psi_to_psf;
            g = g(2);
            rho = rho_1(2,1)/g;
            a = a(2);
            mu = mu(2);
            nu = nu(2);
        end
        

    elseif units == 'SI'
        km_to_m=1000;
        alt_km = h/km_to_m;

        if alt_km > 100 && cell2mat(varargin) == 1
            warning('h>100km Too high to use contiuum assumption.')
            alt_km = round(alt_km,2);
        else
            alt_km = round(alt_km,2);
        end

        if alt_km == 0
            [rho, T, p, a, nu] = ATMOS(h,units);
            g = 9.81;
            mu = nu*rho;
        elseif alt_km < 0 && cell2mat(varargin) == 1
            warning('Altitude below sea level!')
            [rho, T, p, a, nu] = ATMOS(0,units);
            g = 32.2;
            mu = nu*rho;
        elseif alt_km < 0 
            [rho, T, p, a, nu] = ATMOS(0,units);
            g = 32.2;
            mu = nu*rho;
        else
            [~, ~, ~, T, p, rho, a, g, mu, nu, ~, ~, ~] = atmo(alt_km,alt_km,1);
            T = T(2);
            p = p(2);
            rho = rho(2,1);
            a = a(2);
            g = g(2);
            mu = mu(2);
            nu = nu(2);
        end

    else
        error('units of SI or US must be specified')
    end

end