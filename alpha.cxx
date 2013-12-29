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
#include "setup.hh"

static const double StageParameter[STAGES][SP_NUMBER] = STAGE_PARAM;
static double LaunchParameter[LP_NUMBER];

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase)
{
  if (iphase < STAGES) {
	return 0.0;
  } else {
	return 0.0;
	return -tf;
  }
}

void dae(adouble* derivatives, adouble* path, adouble* states,
		 adouble* controls, adouble* parameters, adouble& time,
		 adouble* xad, int iphase)
{
  // derivatives of position
  derivatives[ST_POSX] = states[ST_VELX];
  derivatives[ST_POSY] = states[ST_VELY];
  derivatives[ST_POSZ] = states[ST_VELZ];

  // derivatives of velocity
  adouble Fx = controls[CO_THRX] * 1000;
  adouble Fy = controls[CO_THRY] * 1000;
  adouble Fz = controls[CO_THRZ] * 1000;

  derivatives[ST_VELX] = Fx / StageParameter[iphase-1][SP_MASS];
  derivatives[ST_VELY] = Fy / StageParameter[iphase-1][SP_MASS];
  derivatives[ST_VELZ] = Fz / StageParameter[iphase-1][SP_MASS];

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

	// Cotnrol Constraints
	problem.phases(iphase).bounds.lower.controls(BI(CO_THRX)) = -StageParameter[iphase-1][SP_THRUST];
	problem.phases(iphase).bounds.upper.controls(BI(CO_THRX)) =  StageParameter[iphase-1][SP_THRUST];
	problem.phases(iphase).bounds.lower.controls(BI(CO_THRY)) = -StageParameter[iphase-1][SP_THRUST];
	problem.phases(iphase).bounds.upper.controls(BI(CO_THRY)) =  StageParameter[iphase-1][SP_THRUST];
	problem.phases(iphase).bounds.lower.controls(BI(CO_THRZ)) = -StageParameter[iphase-1][SP_THRUST];
	problem.phases(iphase).bounds.upper.controls(BI(CO_THRZ)) =  StageParameter[iphase-1][SP_THRUST];

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

//  problem.phases(1).guess.states = zeros(ST_NUMBER,SUBDIVISIONS);
//  problem.phases(1).guess.states(BI(ST_POSX), colon()) = LaunchParameter[LP_POSX]* ones(1, SUBDIVISIONS);
//  //problem.phases(1).guess.states(2, colon()) = LaunchParameter[LP_POSY]* ones(1, SUBDIVISIONS);
//  problem.phases(1).guess.states(BI(ST_POSZ), colon()) = LaunchParameter[LP_POSZ]* ones(1, SUBDIVISIONS);
//
//  problem.phases(2).guess.states = zeros(ST_NUMBER,SUBDIVISIONS);
//  problem.phases(2).guess.states(BI(ST_POSX), colon()) = LaunchParameter[LP_POSX]* ones(1, SUBDIVISIONS);
//  //problem.phases(1).guess.states(2, colon()) = LaunchParameter[LP_POSY]* ones(1, SUBDIVISIONS);
//  problem.phases(2).guess.states(BI(ST_POSZ), colon()) = LaunchParameter[LP_POSZ]* ones(1, SUBDIVISIONS);



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
  algorithm.ode_tolerance			= 1.e-5;

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
