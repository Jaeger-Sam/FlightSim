% change_lift computes the change in lift for IBIS due to the
% slipstream effect of the propeller. Etkin B.7, Smelt and Davies (1937)
%
% DCL = change_lift(alpha,J, p_CT, CL_noprop,aircraft)
%
% INPUTS:
%   alpha: angle of attack (rad)
%   J: advance ratio (Vinf/n/d)
%   CL_noprop: lift coefficient with no propulsion
%   aircraft: data structure with aircraft properties
%
% OUTPUTS:
%   DCL: change in lift coefficient
%   
% Sam Jaeger
% jaege246@umn.edu
% 9/15/25

function DCL = change_lift(alpha,J, CL_noprop,aircraft)
    if abs(alpha)<0.01 % conditioning for numerics
        alpha = 0.01;
    end
    
    % x = 18.625; % inches, distance wing CP behind prop
    % S_w = 1140; % inches^2, wing area
    % c = 16.5; % inches, wing chord length
    % d = 18; %inches, prop diameter
    x = aircraft.propulsion.x_wing_CP;
    S_w = aircraft.geom.S_w;
    c = aircraft.geom.c_b_w;
    d = aircraft.propulsion.d;
    lambda = 1.0; % emperical factor, lookup in Etkin

    CT = polyval(aircraft.propulsion.p_CT,J); % evaluate thrust coefficient
    CL0 = CL_noprop;
    %CL0 = CL_alpha0*alpha + CL_0; % lift coefficent without prop

    Tc = CT/J^2;
    a = -0.5 + 0.5*(1 + 8*Tc/pi)^0.5;
    a0 = 2*pi; % 2d lift slope
    s = a + a*x/sqrt(0.25*d^2 + x^2);
    theta_0 = a*alpha/(1 + a);
    theta = ( 1/((0.16*x/d) + theta_0^(-0.8)))^1.25;

    D1 = d*sqrt( (1+a)/(1 + s) );

    DCL = D1*c/S_w*s*(lambda*CL0 - 0.6*a0*theta);
    if isnan(DCL)
        DCL=0;
    else 
        DCL = real(DCL);
    end
end