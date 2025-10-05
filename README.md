# FlightSim
6DOF flight simulator for flight control and estimation research.

Features:
  Runge-Kutta 4 numerical routine,
  Quaterion attitude solution,
  Von Karman turbulence model,
  1976 standard atmosphere,
  WGS84 gravity model,
  Propeller model including moments generated and lift enhancement,
  Feedback control law module with saturation limits,
  Flat Earth corrections,
  Attitude animation plug in.

If attitude animations are desired, download and add this toolbox to your path.


How to run:

  Assemble aircraft and simulation input data files
  
  [aircraft,sim_options] = INITIALIZE_SIMULATION(IBIS_input, Sim_input); % generates aircraft, simulation data structures
  
  [sim_data_out] = RK4(aircraft,sim_options); % integrates trajectory
  
  SIMULATION_PLOTS(sim_data_out,aircraft,sim_options) % plots data
