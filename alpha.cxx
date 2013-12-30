/*
  alpha.cxx

  Copyright 2013 Mhoram Kerbin

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see
  <http://www.gnu.org/licenses/>.
*/

#include "psopt.h"
#include "alpha.hh"

static const double StageParameter[STAGES][SP_NUMBER] = STAGE_PARAM;
static double LaunchParameter[LP_NUMBER];

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase)
{
  if (iphase < STAGES) {
	return 0.0;
  } else {
	return tf;
	return 0.0;
  }
}

void dae(adouble* derivatives, adouble* path, adouble* states,
		 adouble* controls, adouble* parameters, adouble& time,
		 adouble* xad, int iphase)
{
  // Position derivatives
  derivatives[ST_POSX] = states[ST_VELX];
  derivatives[ST_POSY] = states[ST_VELY];
  derivatives[ST_POSZ] = states[ST_VELZ];

  // Velocity derivatives
  adouble Fx = controls[CO_THRX] * 1000;
  adouble Fy = controls[CO_THRY] * 1000;
  adouble Fz = controls[CO_THRZ] * 1000;

  derivatives[ST_VELX] = Fx / states[ST_MASS];
  derivatives[ST_VELY] = Fy / states[ST_MASS];
  derivatives[ST_VELZ] = Fz / states[ST_MASS];

  // Distance Path

  adouble pos[3]; pos[0] = states[ST_POSX]; pos[1] = states[ST_POSY]; pos[2] = states[ST_POSZ];
  adouble pos_norm2 = dot(pos, pos, 3);
  path[PA_DISTANCE] = pos_norm2;

  // Thrust Path
  adouble thrust_norm2 = norm2(controls[CO_THRX], controls[CO_THRY], controls[CO_THRZ]);
  path[PA_THRUST] = thrust_norm2;

  // Mass Derivative
  adouble altitude = sqrt(pos_norm2) - PLANET_RADIUS;
  adouble pressure = calc_pressure(altitude, PLANET_P_0, PLANET_SCALE_HEIGHT);
  adouble isp = get_isp(pressure, StageParameter[iphase-1][SP_ISP_0],
						StageParameter[iphase-1][SP_ISP_VAC]);
  derivatives[ST_MASS] = -sqrt(thrust_norm2) / (isp * G_0) * 1000;
}

adouble norm(adouble x, adouble y, adouble z)
{
  return sqrt(x*x + y*y + z*z);
}

adouble norm2(adouble x, adouble y, adouble z)
{
  return x*x + y*y + z*z;
}

adouble calc_pressure(adouble altitude, adouble p_0, adouble psh)
{
  return p_0 * exp(-altitude / psh);
}

adouble get_isp(adouble pressure, adouble isp_0, adouble isp_vac)
{
  adouble real_p = pressure;
  if (real_p > 1) {
	real_p = 1;
  }
  return isp_0 * real_p + isp_vac * (1 - real_p);
}

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase)
{
  if (iphase == 1) {
	e[E1_POSX] = initial_states[ST_POSX];
	e[E1_POSY] = initial_states[ST_POSY];
	e[E1_POSZ] = initial_states[ST_POSZ];
	e[E1_VELX] = initial_states[ST_VELX];
	e[E1_VELY] = initial_states[ST_VELY];
	e[E1_VELZ] = initial_states[ST_VELZ];
	//	cout << "E " << e[E1_POSX] << " " << e[E1_POSY] << " " << e[E1_POSZ] << endl;
	//	cout << "E " << e[E1_VELX] << " " << e[E1_VELY] << " " << e[E1_VELZ] << endl;
  }
}

void linkages( adouble* linkages, adouble* xad)
{
  adouble time_prev, time_next;
  adouble stat_prev[ST_NUMBER], stat_next[ST_NUMBER];

  int index = 0;
  int iphase;
  for(iphase=1;iphase<STAGES;iphase++) {
	// get states
	time_prev = get_final_time(xad, iphase);
	time_next = get_initial_time(xad, iphase+1);
	get_final_states(stat_prev, xad, iphase);
	get_initial_states(stat_next, xad, iphase+1);

	// time
	linkages[index++] = time_prev - time_next;

	// position
	linkages[index++] = stat_prev[ST_POSX] - stat_next[ST_POSX];
	linkages[index++] = stat_prev[ST_POSY] - stat_next[ST_POSY];
	linkages[index++] = stat_prev[ST_POSZ] - stat_next[ST_POSZ];

	// velocity
	linkages[index++] = stat_prev[ST_VELX] - stat_next[ST_VELX];
	linkages[index++] = stat_prev[ST_VELY] - stat_next[ST_VELY];
	linkages[index++] = stat_prev[ST_VELZ] - stat_next[ST_VELZ];

	// mass
	adouble mass_difference = StageParameter[iphase-1][SP_MASS]-StageParameter[iphase-1][SP_PROPELLANT] - StageParameter[iphase][SP_MASS];
	linkages[index++] = stat_prev[ST_MASS]-mass_difference-stat_next[ST_MASS];
  }

}


int main(void)
{

  Alg  algorithm;
  Sol  solution;
  Prob problem;

  init_launch_parameters();

  // Level 1 Setup
  problem.name = "KSP Launch Optimization Alpha";
  problem.outfilename = "alpha.txt";
  problem.nphases = STAGES;
  problem.nlinkages = (STAGES - 1) * (ST_NUMBER + 1);
  psopt_level1_setup(problem);

  // Level 2 Setup
  int iphase;
  for (iphase=1;iphase <= problem.nphases;iphase++) {
	int events = EN_NUMBER;
	if (iphase == 1) {
	  events += E1_NUMBER;
	}
	if (iphase == STAGES) {
	  events += EF_NUMBER;
	}

	problem.phases(iphase).nstates   = ST_NUMBER;
	problem.phases(iphase).ncontrols = CO_NUMBER;
	problem.phases(iphase).nevents   = events;
	problem.phases(iphase).npath     = PA_NUMBER;
	problem.phases(iphase).nodes     = NODES;
  }
  psopt_level2_setup(problem, algorithm);

  // Constraint Setup
  for (iphase=1;iphase <= problem.nphases;iphase++) {

	// State Constraints
	problem.phases(iphase).bounds.lower.states(BI(ST_POSX)) = -PLANET_SOI;
	problem.phases(iphase).bounds.upper.states(BI(ST_POSX)) =  PLANET_SOI;
	problem.phases(iphase).bounds.lower.states(BI(ST_POSY)) = -PLANET_SOI;
	problem.phases(iphase).bounds.upper.states(BI(ST_POSY)) =  PLANET_SOI;
	problem.phases(iphase).bounds.lower.states(BI(ST_POSZ)) = -PLANET_SOI;
	problem.phases(iphase).bounds.upper.states(BI(ST_POSZ)) =  PLANET_SOI;

	problem.phases(iphase).bounds.lower.states(BI(ST_VELX)) = -PLANET_MAX_V;
	problem.phases(iphase).bounds.upper.states(BI(ST_VELX)) =  PLANET_MAX_V;
	problem.phases(iphase).bounds.lower.states(BI(ST_VELY)) = -PLANET_MAX_V;
	problem.phases(iphase).bounds.upper.states(BI(ST_VELY)) =  PLANET_MAX_V;
	problem.phases(iphase).bounds.lower.states(BI(ST_VELZ)) = -PLANET_MAX_V;
	problem.phases(iphase).bounds.upper.states(BI(ST_VELZ)) =  PLANET_MAX_V;

	problem.phases(iphase).bounds.lower.states(BI(ST_MASS)) = StageParameter[iphase-1][SP_MASS] - StageParameter[iphase-1][SP_PROPELLANT];
	problem.phases(iphase).bounds.upper.states(BI(ST_MASS)) = StageParameter[iphase-1][SP_MASS];

	// Cotnrol Constraints
	problem.phases(iphase).bounds.lower.controls(BI(CO_THRX)) = -StageParameter[iphase-1][SP_THRUST];
	problem.phases(iphase).bounds.upper.controls(BI(CO_THRX)) =  StageParameter[iphase-1][SP_THRUST];
	problem.phases(iphase).bounds.lower.controls(BI(CO_THRY)) = -StageParameter[iphase-1][SP_THRUST];
	problem.phases(iphase).bounds.upper.controls(BI(CO_THRY)) =  StageParameter[iphase-1][SP_THRUST];
	problem.phases(iphase).bounds.lower.controls(BI(CO_THRZ)) = -StageParameter[iphase-1][SP_THRUST];
	problem.phases(iphase).bounds.upper.controls(BI(CO_THRZ)) =  StageParameter[iphase-1][SP_THRUST];

	// Path Constraints
	problem.phases(iphase).bounds.lower.path(BI(PA_DISTANCE)) = (double)PLANET_RADIUS * PLANET_RADIUS;
	problem.phases(iphase).bounds.upper.path(BI(PA_DISTANCE)) = (double)PLANET_SOI * PLANET_SOI;
	problem.phases(iphase).bounds.lower.path(BI(PA_THRUST)) = 1;
	problem.phases(iphase).bounds.upper.path(BI(PA_THRUST)) = StageParameter[iphase-1][SP_THRUST] * StageParameter[iphase-1][SP_THRUST];

  }

  // Event Constraints
  problem.phases(1).bounds.lower.events(BI(E1_POSX)) = LaunchParameter[LP_POSX];
  problem.phases(1).bounds.upper.events(BI(E1_POSX)) = LaunchParameter[LP_POSX];
  problem.phases(1).bounds.lower.events(BI(E1_POSY)) = LaunchParameter[LP_POSY];
  problem.phases(1).bounds.upper.events(BI(E1_POSY)) = LaunchParameter[LP_POSY];
  problem.phases(1).bounds.lower.events(BI(E1_POSZ)) = LaunchParameter[LP_POSZ];
  problem.phases(1).bounds.upper.events(BI(E1_POSZ)) = LaunchParameter[LP_POSZ];

  problem.phases(1).bounds.lower.events(BI(E1_VELX)) = LaunchParameter[LP_VELX];
  problem.phases(1).bounds.upper.events(BI(E1_VELX)) = LaunchParameter[LP_VELX];
  problem.phases(1).bounds.lower.events(BI(E1_VELY)) = LaunchParameter[LP_VELY];
  problem.phases(1).bounds.upper.events(BI(E1_VELY)) = LaunchParameter[LP_VELY];
  problem.phases(1).bounds.lower.events(BI(E1_VELZ)) = LaunchParameter[LP_VELZ];
  problem.phases(1).bounds.upper.events(BI(E1_VELZ)) = LaunchParameter[LP_VELZ];

  setup_time_constraints(problem);
  setup_linkage_constraints(problem);

  // Setup Guesses

  DMatrix x1(ST_NUMBER,SUBDIVISIONS);
  DMatrix x2(ST_NUMBER,SUBDIVISIONS);
  DMatrix u1(CO_NUMBER,SUBDIVISIONS);
  DMatrix u2(CO_NUMBER,SUBDIVISIONS);

  x1(BI(ST_POSX),colon()) = linspace(LaunchParameter[LP_POSX], PLANET_RADIUS + ATMOSPHERIC_HEIGHT, SUBDIVISIONS);
  x1(BI(ST_POSY),colon()) = linspace(LaunchParameter[LP_POSY], LaunchParameter[LP_POSY], SUBDIVISIONS);
  x1(BI(ST_POSZ),colon()) = linspace(LaunchParameter[LP_POSZ], LaunchParameter[LP_POSZ], SUBDIVISIONS);
  x1(BI(ST_VELX),colon()) = linspace(LaunchParameter[LP_VELX], LaunchParameter[LP_VELX], SUBDIVISIONS);
  x1(BI(ST_VELY),colon()) = linspace(LaunchParameter[LP_VELY], LaunchParameter[LP_VELY], SUBDIVISIONS);
  x1(BI(ST_VELZ),colon()) = linspace(LaunchParameter[LP_VELZ], LaunchParameter[LP_VELZ], SUBDIVISIONS);
  x1(BI(ST_MASS),colon()) = linspace(StageParameter[0][SP_MASS], StageParameter[0][SP_MASS] - StageParameter[0][SP_PROPELLANT], SUBDIVISIONS);

  problem.phases(1).guess.states = x1;

  x2(BI(ST_POSX),colon()) = linspace(PLANET_RADIUS + ATMOSPHERIC_HEIGHT, PLANET_RADIUS + ATMOSPHERIC_HEIGHT, SUBDIVISIONS);
  x2(BI(ST_POSY),colon()) = linspace(LaunchParameter[LP_POSY], LaunchParameter[LP_POSY] + ATMOSPHERIC_HEIGHT, SUBDIVISIONS);
  x2(BI(ST_POSZ),colon()) = linspace(LaunchParameter[LP_POSZ], LaunchParameter[LP_POSZ], SUBDIVISIONS);
  x2(BI(ST_VELY),colon()) = linspace(LaunchParameter[LP_VELY], LaunchParameter[LP_VELY], SUBDIVISIONS);
  x2(BI(ST_MASS),colon()) = linspace(StageParameter[1][SP_MASS], StageParameter[1][SP_MASS] - StageParameter[1][SP_PROPELLANT], SUBDIVISIONS);

  problem.phases(2).guess.states = x2;

  u1(BI(CO_THRX),colon()) = linspace(StageParameter[0][SP_THRUST], StageParameter[0][SP_THRUST], SUBDIVISIONS);
  u1(BI(CO_THRY),colon()) = linspace(0, 0, SUBDIVISIONS);
  u1(BI(CO_THRZ),colon()) = linspace(0, 0, SUBDIVISIONS);
  u2(BI(CO_THRX),colon()) = linspace(0, 0, SUBDIVISIONS);
  u2(BI(CO_THRY),colon()) = linspace(StageParameter[1][SP_THRUST], StageParameter[1][SP_THRUST], SUBDIVISIONS);
  u2(BI(CO_THRZ),colon()) = linspace(0, 0, SUBDIVISIONS);

  problem.phases(1).guess.controls = u1;
  problem.phases(2).guess.controls = u2;


  // Start
  problem.endpoint_cost = &endpoint_cost;
  problem.dae = &dae;
  problem.events = &events;
  problem.linkages = &linkages;

  algorithm.nlp_method                  	= "IPOPT";
  algorithm.scaling                     	= "automatic";
  algorithm.derivatives                 	= "automatic";
  algorithm.nlp_iter_max                	= 500;
  //algorithm.mesh_refinement                   = "automatic";
  //algorithm.collocation_method = "trapezoidal";
  algorithm.ode_tolerance			= 1.e-6;

  psopt(solution, problem, algorithm);

}

void init_launch_parameters()
{
  double dist = PLANET_RADIUS + LAUNCH_ELEVATION;

  double f = dist * 2 * M_PI * cos(LAUNCH_LATITUDE) / PLANET_ROT_PER;
  LaunchParameter[LP_POSX] = dist * cos(LAUNCH_LATITUDE) * cos(LAUNCH_LONGITUDE);
  LaunchParameter[LP_POSY] = dist * cos(LAUNCH_LATITUDE) * sin(LAUNCH_LONGITUDE);
  LaunchParameter[LP_POSZ] = dist * sin(LAUNCH_LATITUDE);
  LaunchParameter[LP_VELX] = -sin(LAUNCH_LONGITUDE) * f;
  LaunchParameter[LP_VELY] =  cos(LAUNCH_LONGITUDE) * f;
  LaunchParameter[LP_VELZ] = 0;

  cout << "LAUNCH " << LAUNCH_LATITUDE << " " << LAUNCH_LONGITUDE << " " << LaunchParameter[LP_POSX] << " " << LaunchParameter[LP_POSY] << " " << LaunchParameter[LP_POSZ] << " | " << LaunchParameter[LP_VELX] << " " << LaunchParameter[LP_VELY] << " " << LaunchParameter[LP_VELZ] << endl;

}

void setup_time_constraints(Prob problem)
{
  int iphase;
  double tim = 0;
  for (iphase = 1;iphase <= problem.nphases;iphase++) {
	problem.phases(iphase).bounds.lower.StartTime = tim;
	problem.phases(iphase).bounds.upper.StartTime = tim * STAGE_TIME_MAX_FACTOR;
	tim += StageParameter[iphase-1][SP_MIN_STAGE_TIME];
	problem.phases(iphase  ).bounds.lower.EndTime   = tim;
	problem.phases(iphase  ).bounds.upper.EndTime   = tim * STAGE_TIME_MAX_FACTOR;
  }
}

void setup_linkage_constraints(Prob problem)
{
  int i;
  for(i=1;i<= problem.nlinkages;i++) {
	problem.bounds.lower.linkage(i) = 0;
	problem.bounds.upper.linkage(i) = 0;
  }
}
