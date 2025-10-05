% grav_WGS84.fcn computes the acceleration due to gravity (and centrifugal)
% using the Department of Defense World Geodetic System 1984 (WGS-84). The 
% basic method is outlined in Section 1.6 in Stevens, Lewis Aicraft Control
% and Simulation. 
%
% g = grav_WGS84(h,psi, units)
%
% INPUTS:
%   h: altitude or geodetic height
%   psi: latitude in deg
%   units: 'SI' or 'US'
%
% OUTPUTS:
%   g: gravity vector in flat Earth North-East-Down frame
%
% Sam Jaeger
% jaege246@umn.edu
% 11/2/2024

function g = grav_WGS84(h,psi, units)
    %l = 0; % results should be invariant of longitude
    if units == 'US'
        h = h*3.28084;
    elseif units =='SI'
    else
        error('units must be SI or US!')
    end

    % constants
    J2 = 1.082626684e-3;
    a = 6378137; % Semimajor axis of Earth in m
    omega_E = 7.2921150e-5; %rad/s sidereal rate of rotation 
    GM = 3986004.418e8; % m^3/s^2, Earth's gravitational constant
    e = .0818191908426; % eccentricity
    b = 6356752; % semiminor axis of Earth in m

    % radius of oblate Earth
    r_c = sqrt((b^2)/(1 - (e*cosd(psi))^2 )); % 1.6-11
    
    r = r_c + h; %approximate radius
    phi = psi; % !!!! approximate relation, D=0
    
    % could fix this
%     % prime vertical radius of curvature
%     N = a/(1 - (e*sind(phi)^2))^.5; 
%     % normal
%     n = (e^2)*N/(N + h);
%     % radius to center of earth
%     r = (N+h)*(1 - n*(2-n)*sind(phi)^2)^.5;% 1.6-15

    % Earth-Centered-Earth-Fixed (ECEF) coord system
    %p = [r*cosd(psi)*cosd(l); r*cosd(psi)*sind(l); r*sind(psi)];
    p = [r*cosd(psi); 0; r*sind(psi)];

    % grav.
    G_ecef(1) = - GM/(r^2)*((1 + (3/2)*((a/r)^2)*J2*(1 - 5*sind(psi)^2)).*p(1)/r );
    G_ecef(2) = - GM/(r^2)*((1 + (3/2)*((a/r)^2)*J2*(1 - 5*sind(psi)^2)).*p(2)/r );
    G_ecef(3) = - GM/(r^2)*((1 + (3/2)*((a/r)^2)*J2*(3 - 5*sind(psi)^2)).*p(3)/r );

    % centrifugal
    omega_ei = omega_E*[0;0;1];
    g_cent = cross(omega_ei,cross(omega_ei,p));

    % gravity
    g = G_ecef' - g_cent;
    
    % convert to North-East-Down (NED) coordinate system
%     C_ned_ecef = [-sind(phi)*cosd(l), -sind(phi)*sind(l), cosd(phi);
%                   -sind(l),             cosd(l),            0;
%                   -cosd(phi)*cosd(l), -cosd(phi)*sind(l), -sind(phi)]; %   DCM
    C_ned_ecef = [-sind(phi), 0, cosd(phi);
                  0,          1,       0;
                  -cosd(phi), 0, -sind(phi)]; %   DCM
    g = C_ned_ecef*g;

    % convert to us units
    if units == 'US'
        g = g*3.28084;
    end
end