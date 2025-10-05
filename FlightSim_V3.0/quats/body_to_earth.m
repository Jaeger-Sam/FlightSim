% earth_to_body.fcn trasnforms vector v_b in body fixed coordinates to
% earth fixed coordinates using a quaternion, es, that describes the 
% rotation. This function uses the inverse of Equation 11.6.8 in Phillips 
% Mechanics of Flight.
%
%  [phi_theta_psi] = attitude_from_quats(es)
%
% INPUTS:
%   v_b 3x1 vector in body fixed frame
%   es: 4x1 quaterion describing orientation of the earth fixed frame
%
% OUTPUTS:
%   v_f: 3x1 vector in earth fixed fram
%
% Sam Jaeger
% jaege246@umn.edu
% 1/2/2024

function v_f = body_to_earth(es,v_b)
    v_f = quat_mult( es, quat_mult([0,v_b],quat_conj(es)) );
    v_f = v_f(2:4);
end