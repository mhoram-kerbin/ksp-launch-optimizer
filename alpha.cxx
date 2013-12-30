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

static const double StageParameters[STAGES][SP_NUMBER] = STAGE_PARAM;
static double LaunchParameters[LP_NUMBER];

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase)
{

  return 0.0;

}

void dae(adouble* derivatives, adouble* path, adouble* states,
		 adouble* controls, adouble* parameters, adouble& time,
		 adouble* xad, int iphase)
{

}

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase)
{

}

void linkages( adouble* linkages, adouble* xad)
{

}


int main(void)
{

  Alg  algorithm;
  Sol  solution;
  Prob problem;

  init_launch_parameters();

  level_1_setup(problem);
  level_2_setup();

  setup_phase_bounds_information(problem);


}

void init_launch_parameters()
{
  double dist = PLANET_RADIUS + LAUNCH_ELEVATION;

  double f = dist * 2 * M_PI * cos(LAUNCH_LATITUDE) / PLANET_ROT_PER;
  LaunchParameters[LP_POSX] = dist * cos(LAUNCH_LATITUDE) * cos(LAUNCH_LONGITUDE);
  LaunchParameters[LP_POSY] = dist * cos(LAUNCH_LATITUDE) * sin(LAUNCH_LONGITUDE);
  LaunchParameters[LP_POSZ] = dist * sin(LAUNCH_LATITUDE);
  LaunchParameters[LP_VELX] = -sin(LAUNCH_LONGITUDE) * f;
  LaunchParameters[LP_VELY] =  cos(LAUNCH_LONGITUDE) * f;
  LaunchParameters[LP_VELZ] = 0;

  cout << "LAUNCH " << LAUNCH_LATITUDE << " " << LAUNCH_LONGITUDE << " " << LaunchParameters[LP_POSX] << " " << LaunchParameters[LP_POSY] << " " << LaunchParameters[LP_POSZ] << " | " << LaunchParameters[LP_VELX] << " " << LaunchParameters[LP_VELY] << " " << LaunchParameters[LP_VELZ] << endl;

}

void level_1_setup(Prob problem)
{
  problem.name = "KSP Launch Optimization Alpha";
  problem.outfilename = "alpha.txt";
  problem.nphases = STAGES;
  problem.nlinkages = (STAGES - 1) * (ST_NUMBER + 1);
  psopt_level1_setup(problem);

}

void level_2_setup(Prob problem, Alg algorithm)
{
  int iphase;

  setup_phase_level_information(problem);
  psopt_level2_setup(problem, algorithm);
}

void setup_phase_level_information(Prob problem)
{
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

}

void setup_phase_bounds_information(Prob problem)
{
  int iphase;
  for (iphase=1;iphase <= problem.nphases;iphase++) {
	// States
	problem.phases(iphase).bounds.lower.states(ST_POSX) = 0;
	problem.phases(iphase).bounds.upper.states(ST_POSX) = PLANET_SOI;
	problem.phases(iphase).bounds.lower.states(ST_POSY) = 0;
	problem.phases(iphase).bounds.upper.states(ST_POSY) = PLANET_SOI;
	problem.phases(iphase).bounds.lower.states(ST_POSZ) = 0;
	problem.phases(iphase).bounds.upper.states(ST_POSZ) = PLANET_SOI;

	// Controls

	//Events

  }

}
