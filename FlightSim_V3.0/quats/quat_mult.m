% quat_mult.fcn performs quaternion mutliplication defined by equation
% 11.6.3 in Phillips Mechanics of Flight.
%
% C = quat_mult(A,B)
%
% INPUTS:
%   A: quaternion with components A = [A_0, A_x, A_y, A_z]'
%   B: quaternion with components B = [B_0, B_x, B_y, B_z]'
%
% OUTPUTS:
%   C: quaternion with components C = [C_0, C_x, C_y, C_z]'
%
% Sam Jaeger
% jaege246@umn.edu
% 1/2/2024

function C = quat_mult(A,B)
    C(1) = A(1)*B(1) - A(2)*B(2) - A(3)*B(3) - A(4)*B(4);
    C(2) = A(1)*B(2) + A(2)*B(1) + A(3)*B(4) - A(4)*B(3);
    C(3) = A(1)*B(3) - A(2)*B(4) + A(3)*B(1) + A(4)*B(2);
    C(4) = A(1)*B(4) + A(2)*B(3) - A(3)*B(2) + A(4)*B(1);
end