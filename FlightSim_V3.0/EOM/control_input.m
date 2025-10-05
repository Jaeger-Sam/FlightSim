% This function retuns a vector of control inputs given the simulation time,
%  state vector, saturation limits, and input type. Currently, there are
%  four types of inputs: autonomous (zeros), step, predfined time history,
%  feedback control law in the form of a sensor mapping matrix C, and gain
%  matrix K.
%
% u = control_input(t,x,max_deflect,input_type)
%
% INPUTS:
%   t: simulation time
%   x: state vector
%   control: data structure with the following inputs
%       .max_deflect: matrix of saturation states (rad). First column is 
%           lower bound, second column is upper bound.
%       .input_type
%           integer between 1 and 3
%       .auto -> 0 (or anything else)
%           u == 0
%       .step -> 1
%           .tstep: time to apply step function
%           .amp0: u static values before step
%           .amp: u step values
%       .thist -> 2
%           .thist: timevector associated with input time history
%           .uhist: input history matirx (column times, row inputs)
%           .ttol: time tolerance for finding input number
%       .feedback -> 3
%           .C: matrix mapping state vector to sensor output
%           .K: gain matrix
%           .yref: yreference signal
%           .ustat: static control input
%
% OUTPUTS:
%   u: vector of control inputs
%       u(1): delta_T
%       u(2): delta_e
%       u(3): delta_a
%       u(4): delta_r
%       u(5): delta_f
%
% Sam Jaeger
% jaege246@umn.edu
% 1/3/2024
%   Revised: 1/23/2024
%   Revised: 3/15/2024
%   Revised: 9/21/2025

function u = control_input(t,x,control)
    if control.input_type == 1 % step
        if t < control.step.tstep
            u = control.step.amp0; % vector of 
        else % vector of amplitudes to apply step function to
            u = control.step.amp;
        end
        
    elseif control.input_type == 2 % predefined time history
        tindex = find(abs(t-control.thist.thist)<control.thist.ttol,1);
        %tindex = find(abs(t-control.thist.thist));
        u = control.thist.uhist(:,tindex);

    elseif control.input_type == 3 % feedback law
        C = control.feedback.C;
        K = control.feedback.K;
        yref = control.feedback.yref;
        ustat = control.feedback.ustat;

        %y = C*x;
        % proportional altitude & bank controller
        phi_theta_psi = attitude_from_quats(x(10:13)');
        y = [-x(9);phi_theta_psi(1)];
        u = K*(y-yref) + ustat;
        
        % proportional pitch attitude controller
        % [phi,theta,psi] = attitude_from_quats(x(10:13));
        % u(1) = 1;
        % u(2) = theta*K;
        % u(3:5) = zeros(3,1);
    else 
        u = zeros(5,1); % autonomous
    end

    % Saturation Limits
    for ii=1:length(u)
        if u(ii) > control.max_deflect(ii,1)
            u(ii) = control.max_deflect(ii,1);
        elseif u(ii) < control.max_deflect(ii,2)
            u(ii) = control.max_deflect(ii,2);
        end
    end
end