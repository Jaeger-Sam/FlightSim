% attitude_from_quats.fcn returns the bank angle, pitch, and yaw attitude
% angles from a quaternion vector. From Equation 11.7.11 in Phillips
% Mechanics of Flight. See also quat2eul.fcn.
%
% [phi_theta_psi] = attitude_from_quats(es)
%
% INPUTS:
%   es: vector or matrix of quaternions
%
% OUTPUTS:
%   phi: bank angle (rad)
%   theta: pitch angle (rad)
%   psi: heading (rad)
%
% Written by:
%   Sam Jaeger
%   jaege246@umn.edu
%   11/4/2023
%   Revised: 1/2/2024
%       Added comments and changed phi to NaN

function [phi_theta_psi] = attitude_from_quats(es)
    vect_matrix = size(es); % check if vector or matrix
    
    if vect_matrix(2) ~= 4 % logic 
        error('Quaternion vector must have 4 components!')
    end

    if vect_matrix(1)>1 % matrix of quaternions
        for ii=1:length(es(:,1)) % pointing straight up
            if es(ii,1)*es(ii,3) - es(ii,2)*es(ii,4) == 0.5
                phi_theta_psi(ii,1) = 2*asin(es(ii,2)/cos(pi/4)) + 0; % +phi but phi unknown?
                phi_theta_psi(ii,2) = pi/2;
                phi_theta_psi(ii,3) = NaN;
            elseif es(ii,1)*es(ii,3) - es(ii,2)*es(ii,4) == -0.5
                phi_theta_psi(ii,1) = 2*asin(es(2)/cos(pi/4)) - 0; % -phi but phi unknown?
                phi_theta_psi(ii,2) = -pi/2;
                phi_theta_psi(ii,3) = NaN;
            else 
                phi_theta_psi(ii,1) = atan2(2*(es(ii,1)*es(ii,2) + es(ii,3)*es(ii,4)), es(ii,1)^2 + es(ii,4)^2 - es(ii,2)^2 - es(ii,3)^2);
                phi_theta_psi(ii,2) = asin( 2*(es(ii,1)*es(ii,3) - es(ii,2)*es(ii,4) ) );
                phi_theta_psi(ii,3) = atan2(2*(es(ii,1)*es(ii,4) + es(ii,2)*es(ii,3) ), es(ii,1)^2 + es(ii,2)^2 - es(ii,3)^2 - es(ii,4)^2 );
            end
        end

    else % vector of quaternions
        if es(1)*es(3) - es(2)*es(4) == 0.5
            phi_theta_psi(1) = 2*asin(es(2)/cos(pi/4)) + 0; % +phi but phi unknown?
            phi_theta_psi(2) = pi/2;
            phi_theta_psi(3) = NaN;
        elseif es(1)*es(3) - es(2)*es(4) == -0.5
            phi_theta_psi(1) = 2*asin(es(2)/cos(pi/4)) - 0; % -phi but phi unknown?
            phi_theta_psi(2) = -pi/2;
            phi_theta_psi(3) = NaN;
        else 
            phi_theta_psi(1) = atan2(2*(es(1)*es(2) + es(3)*es(4)), es(1)^2 + es(4)^2 - es(2)^2 - es(3)^2);
            phi_theta_psi(2) = asin( 2*(es(1)*es(3) - es(2)*es(4) ) );
            phi_theta_psi(3) = atan2(2*(es(1)*es(4) + es(2)*es(3) ), es(1)^2 + es(2)^2 - es(3)^2 - es(4)^2 );
        end
    end
end