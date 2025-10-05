# FlightSim
6DOF flight simulator for flight control and estimation research.

Features:
  Runge-Kutta 4 numerical routine,
  quaternion attitude solution,
  Von Karman turbulence model,
  1976 standard atmosphere,
  WGS84 gravity model,
  propeller model including moments generated and lift enhancement,
  feedback control law module with saturation limits,
  flat Earth corrections,
  attitude animation plug in.

If attitude animations are desired, download and add this toolbox to your path.
https://github.com/Ro3code/aircraft_3d_animation

How to run:

  Assemble aircraft and simulation input data files, then call the following commands:
  
  [aircraft,sim_options] = INITIALIZE_SIMULATION(IBIS_input, Sim_input); % generates aircraft, simulation data structures
  
  [sim_data_out] = RK4(aircraft,sim_options); % integrates trajectory
  
  SIMULATION_PLOTS(sim_data_out,aircraft,sim_options) % plots data and generates animation
