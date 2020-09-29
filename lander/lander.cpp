// Mars lander simulator
// Version 1.11
// Mechanical simulation functions
// Gabor Csanyi and Andrew Gee, August 2019

// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation, to make use of it
// for non-commercial purposes, provided that (a) its original authorship
// is acknowledged and (b) no modified versions of the source code are
// published. Restriction (b) is designed to protect the integrity of the
// exercise for future generations of students. The authors would be happy
// to receive any suggested modifications by private correspondence to
// ahg@eng.cam.ac.uk and gc121@eng.cam.ac.uk.

#include "lander.h"

void autopilot_land (void)
  // Autopilot to adjust the engine throttle, parachute and attitude control
{
  // INSERT YOUR CODE HERE
  // Enable altitude stabilisation
  stabilized_attitude = true;

  // Orbital re-entry sequence: if periapsis of current orbit greater than re-entry alt, then thrust to decease speed
  if (r_p > (MARS_RADIUS+77500)) { // 99500
    stabilized_attitude_angle = -90;
    throttle = 1;
  }
  // Deploy parachute if lander under 70km and safe to deploy
  else if (parachute_status == NOT_DEPLOYED && safe_to_deploy_parachute() && altitude < 70000) {
    parachute_status = DEPLOYED;
    stabilized_attitude_angle = 0;
  }
  else {
    double throttle_vert, throttle_horiz;

    // Calculate error term
    double error = -(0.5 + K_h * altitude + (velocity * position.norm()));

    // Calculate P_out
    double P_out = K_p * error;

    // Calculate new downwards throttle level
    if (P_out <= (-delta)) {
      throttle_vert = 0;
    } else if (P_out < (1 - delta)) {
      throttle_vert = delta + P_out;
    } else {
      throttle_vert = 1;
    }

    if (altitude < 200 && ground_speed_wrt_mars > 1) { // wrt
      throttle_horiz = 0.2;
      stabilized_attitude_angle = -atan(throttle_horiz / throttle_vert) * 180 / M_PI;
      // cout << "vert = " << throttle_vert << "\thoriz = " << throttle_horiz <<  "\tangle = " << stabilized_attitude_angle <<  endl;
    } else {
      throttle_horiz = 0;
    }

    throttle = double (min(sqrt(pow(throttle_vert, 2) + pow(throttle_horiz, 2)), 1.0));
   }
}

void autopilot_inject(void) {
  
  stabilized_attitude = true;
  fuel_rate_at_max_thrust = 0;
  double target_semi_major = 0.5 * (periapsis + apoapsis);
  double target_eccentricity = 1 - (periapsis/target_semi_major);
//  double target_ground_speed = sqrt(GRAVITY * MARS_MASS / periapsis);
  
  cout  << "theta = " << (theta*180/M_PI) << "\te = " << eccentricity << "\ttarget_e = " << target_eccentricity << "\treached_circular = " << reached_circular_orbit << "\treached_final = " << reached_final_orbit << "\tperi = " << r_p << " " << to_string((r_p-periapsis)/periapsis*100) << "\tapo = " << r_a << " " << to_string((r_a-apoapsis)/apoapsis*100) << endl;

  if (!reached_circular_orbit && !reached_final_orbit) {
    if (position.abs() < (periapsis+2587894.3661202100)/1.6883873497) { 
      stabilized_attitude_angle = 0;
      throttle = 1;
    } else if ((theta*180/M_PI) < 179 && eccentricity > 0.001) {
      stabilized_attitude_angle = 0;
      throttle = 1;
    } else if ((theta*180/M_PI) >= 179 && (theta*180/M_PI) <= 181 && eccentricity > 0.001) {
      stabilized_attitude_angle = 90;
      throttle = 1;
    } else if ((theta*180/M_PI) > 181 && eccentricity > 0.001) {
      stabilized_attitude_angle = 180;
      throttle = 1;
    } else {
      stabilized_attitude_angle = 0;
      throttle = 0;
      reached_circular_orbit = true;
    }
//    if (eccentricity > 0.97 && position.abs() < periapsis) { // 2505579.999998770
//      stabilized_attitude_angle = 15;
//      throttle = 1;
//    } else if (eccentricity > 0.001 && climb_speed >  0) {
//      stabilized_attitude_angle = 90 - atan(control_function(climb_speed, 0, 1, 3, 0)/(eccentricity))*180/M_PI;
//      throttle = 1;
//    } else {
//      stabilized_attitude_angle = 0;
//      throttle = 0;
//      reached_circular_orbit = true;
//    }
  } else if (periapsis == apoapsis && reached_circular_orbit) {
    stabilized_attitude_angle = 0;
    throttle = 0;
    reached_final_orbit = true;
  } else if (reached_circular_orbit && !reached_final_orbit) {
    if (eccentricity < target_eccentricity) {
      stabilized_attitude_angle = 90;
      throttle = 1;
    } else {
      reached_final_orbit = true;
    }
  } else {
    stabilized_attitude_angle = 0;
    throttle = 0;
  }
}

void numerical_dynamics (void)
  // This is the function that performs the numerical integration to update the
  // lander's pose. The time step is delta_t (global variable).
{
  // INSERT YOUR CODE HERE
  // Define varibles for use in trajectory prediction
  double theta_dot;
  double omega = 2*M_PI / MARS_DAY;
  double sine_theta = (position^vector3d(0,0,1)).abs() / position.abs();
  vector3d planet_unit = (vector3d(0,0,1) ^ vector3d(position.x, position.y, 0)).norm();
  velocity_wrt_atm = (velocity - (omega * position.abs() * sine_theta * planet_unit));
      
  // Calculate current mass of lander with fuel
  tot_mass = UNLOADED_LANDER_MASS + (fuel * FUEL_CAPACITY * FUEL_DENSITY);
  
  // Calculate drag force on the lander and the parachute if it is deployed
  drag_force_lander = -0.5 * atmospheric_density(position) * DRAG_COEF_LANDER * (M_PI * pow(LANDER_SIZE, 2)) * pow(velocity_wrt_atm.abs(), 2) * velocity_wrt_atm.norm();
  drag_force_chute = vector3d(0,0,0);
  if (parachute_status == DEPLOYED) {
    drag_force_chute = -0.5 * atmospheric_density(position) * DRAG_COEF_CHUTE * (5 * pow((2*LANDER_SIZE) , 2)) * pow(velocity_wrt_atm.abs(), 2) * velocity_wrt_atm.norm();
  }
  
  // Calculate gravitational force
  grav_force = - (GRAVITY * MARS_MASS * tot_mass) / pow((position.abs()), 2) * position.norm();
  
  // Update acceleration
  acceleration = (thrust_wrt_world() + grav_force + drag_force_lander + drag_force_chute) / tot_mass;
  
  // COMMENT OUT VERLET OR EULER INTEGRATION CODE BELOW
  // VERLET INTEGRATION with Euler for first time step
  if (previous_position == vector3d(0 , 0 , 0)) {
    previous_position = position;
    // Euler integration for the first time step
    position = position + delta_t * velocity;
    velocity = velocity + delta_t * acceleration;
  } else {
    // Verlet integration for all the other time steps
    vector3d intermediary_position = position;
    position = 2 * intermediary_position - previous_position + pow(delta_t, 2) * acceleration;
    velocity = (1 / delta_t) * (position - intermediary_position);
    previous_position = intermediary_position;
  }
  
  // EULER INTEGRATION
//  position = position + delta_t * velocity;
//  velocity = velocity + delta_t * acceleration;

  // Calculation of periapsis for re-entry sequence
  ang_momentum = position ^ (tot_mass * velocity);
  theta_dot = ang_momentum.abs() / (tot_mass * position.abs2());
  orbit_energy = 0.5 * tot_mass * velocity.abs2() - (GRAVITY * MARS_MASS * tot_mass / position.abs());
  eccentricity = sqrt(1 + ((2 * orbit_energy * ang_momentum.abs2()) / (pow(tot_mass, 3) * pow((-GRAVITY * MARS_MASS), 2))));
  r_p = (ang_momentum.abs2() / (pow(tot_mass, 2) * (GRAVITY * MARS_MASS))) * (1 / (1 + eccentricity));
  semi_major = r_p / (1 - eccentricity); // a on the elipse diagram
  semi_minor = semi_major * sqrt(1 - pow(eccentricity, 2)); // b on the elipse diagram
  r_a = 2*semi_major - r_p;

  
  cos_theta = ((ang_momentum.abs2()/(pow(tot_mass, 2) * (GRAVITY * MARS_MASS) * position.abs())) - 1) / eccentricity;
  theta = acos(cos_theta); // theta in radians!
  if (climb_speed > 0) theta = 2*M_PI - theta;
    
//  cout << "theta = " << theta*180/M_PI << endl;
  
  // Here we can apply an autopilot to adjust the thrust, parachute and attitude
  if (autopilot_enabled && autopilot_mode == 0) autopilot_land();
  if (autopilot_enabled && autopilot_mode == 1) autopilot_inject();

  // Here we can apply 3-axis stabilization to ensure the base is always pointing downwards
  if (stabilized_attitude) attitude_stabilization();
}

void initialize_simulation (void)
  // Lander pose initialization - selects one of 10 possible scenarios
{
  // The parameters to set are:
  // position - in Cartesian planetary coordinate system (m)
  // velocity - in Cartesian planetary coordinate system (m/s)
  // orientation - in lander coordinate system (xyz Euler angles, degrees)
  // delta_t - the simulation time step
  // boolean state variables - parachute_status, stabilized_attitude, autopilot_enabled
  // scenario_description - a descriptive string for the help screen

  scenario_description[0] = "circular orbit";
  scenario_description[1] = "descent from 10km";
  scenario_description[2] = "elliptical orbit, thrust changes orbital plane";
  scenario_description[3] = "polar launch at escape velocity (but drag prevents escape)";
  scenario_description[4] = "elliptical orbit that clips the atmosphere and decays";
  scenario_description[5] = "descent from 200km";
  scenario_description[6] = "aerostationary orbit";
  scenario_description[7] = "";
  scenario_description[8] = "";
  scenario_description[9] = "";
  
  previous_position = vector3d(0,0,0);
  
  switch (scenario) {

  case 0:
    // a circular equatorial orbit
    position = vector3d(1.2*MARS_RADIUS, 0.0, 0.0);
    velocity = vector3d(0.0, -3247.087385863725, 0.0);
    orientation = vector3d(0.0, 90.0, 0.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    reached_circular_orbit = false;
    reached_final_orbit = false;
    break;

  case 1:
    // a descent from rest at 10km altitude
    position = vector3d(0.0, -(MARS_RADIUS + 10000.0), 0.0);
    velocity = vector3d(0.0, 0.0, 0.0);
    orientation = vector3d(0.0, 0.0, 90.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = true;
    autopilot_enabled = false;
    reached_circular_orbit = false;
    reached_final_orbit = false;
    break;

  case 2:
    // an elliptical polar orbit
    position = vector3d(0.0, 0.0, 1.2*MARS_RADIUS);
    velocity = vector3d(3500.0, 0.0, 0.0);
    orientation = vector3d(0.0, 0.0, 90.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    reached_circular_orbit = false;
    reached_final_orbit = false;
    break;

  case 3:
    // polar surface launch at escape velocity (but drag prevents escape)
    position = vector3d(0.0, 0.0, MARS_RADIUS + LANDER_SIZE/2.0);
    velocity = vector3d(0.0, 0.0, 5027.0);
    orientation = vector3d(0.0, 0.0, 0.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    reached_circular_orbit = false;
    reached_final_orbit = false;
    break;

  case 4:
    // an elliptical orbit that clips the atmosphere each time round, losing energy
    position = vector3d(0.0, 0.0, MARS_RADIUS + 100000.0);
    velocity = vector3d(4000.0, 0.0, 0.0);
    orientation = vector3d(0.0, 90.0, 0.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = false;
    autopilot_enabled = false;
    reached_circular_orbit = false;
    reached_final_orbit = false;
    break;

  case 5:
    // a descent from rest at the edge of the exosphere
    position = vector3d(0.0, -(MARS_RADIUS + EXOSPHERE), 0.0);
    velocity = vector3d(0.0, 0.0, 0.0);
    orientation = vector3d(0.0, 0.0, 90.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = true;
    autopilot_enabled = false;
    reached_circular_orbit = false;
    reached_final_orbit = false;
    break;

  case 6:
    // a circular equitorial orbit where the period of orbit is the same as the orbital period of Mars on its axis of rotation
    position = vector3d(20428000.0, 0.0, 0.0);
    velocity = vector3d(0.0, 1448.155346, 0.0);
    orientation = vector3d(0.0, 0.0, 90.0);
    delta_t = 0.1;
    parachute_status = NOT_DEPLOYED;
    stabilized_attitude = true;
    autopilot_enabled = false;
    reached_circular_orbit = false;
    reached_final_orbit = false;
    break;

  case 7:
    break;

  case 8:
    break;

  case 9:
    break;

  }
}
