% quat_conj.fcn computes the conjugate quaternion defined in Phillips
% Mechanics of Flight Equation 11.6.6.
%
% Qs = quat_conj(Q)
%
% INPUTS:
%   Q: Quaternion
%
% OUTPUTS:
%   Qs: Conjugate of quaternion Q
%
% Sam Jaeger
% jaege246@umn.edu
% 1/2/2024

function Qs = quat_conj(Q)
    Qs = [Q(1), -Q(2), -Q(3), -Q(4)];
end