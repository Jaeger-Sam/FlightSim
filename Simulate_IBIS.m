% Simulates IBIS (Ultra Stick)
%
% Sam Jaeger
% jaege246@umn.edu
% 10/5/2025

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% INPUT
IBIS_input = 'IBIS_INPUT.m';
Sim_input = 'IBIS_SIMULATION_INPUT.m';

[aircraft,sim_options] = INITIALIZE_SIMULATION(IBIS_input, Sim_input);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Simulate IBIS
[sim_data_out] = RK4(aircraft,sim_options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot Integrated Quantities
SIMULATION_PLOTS(sim_data_out,aircraft,sim_options)