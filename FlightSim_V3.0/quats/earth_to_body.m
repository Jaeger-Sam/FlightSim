% earth_to_body.fcn trasnforms vector v_f in earth coordinates to body
% fixed coordinates using a quaternion, es, that describes the rotation.
% This function uses Equation 11.6.8 in Phillips Mechanics of Flight.
%
% v_b = earth_to_body(es,v_f)
%
% INPUTS:
%   v_f: 3x1 vector in earth fixed frame
%   es: 4x1 quaterion describing orientation of the body fixed frame
%
% OUTPUTS:
%   v_b: 3x1 vector in body fixed fram
%
% Sam Jaeger
% jaege246@umn.edu
% 1/2/2024

function v_b = earth_to_body(es,v_f)
    v_b = quat_mult( quat_conj(es), quat_mult([0,v_f],es));
    v_b = v_b(2:4);
end