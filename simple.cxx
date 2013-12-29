/*
  launch.cxx

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

#include "simple.hh"
#include "setup.hh"

static const double StageParameters[STAGES][7] = STAGE_PARAM;
static double LaunchParameters[6];

adouble norm(adouble x, adouble y, adouble z)
{
  return sqrt(x*x + y*y + z*z);
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

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase)
{

    if (iphase == STAGES) {
	  return tf;
	  return -final_states[ST_MASS];
	} else {
	  return 0.0;
	}

}

void dae(adouble* derivatives, adouble* path, adouble* states,
		 adouble* controls, adouble* parameters, adouble& time,
		 adouble* xad, int iphase)
{

  //  adouble pos[3]; pos[0] = states[ST_POSX]; pos[1] = states[ST_POSY]; pos[2] = states[ST_POSZ];
  adouble pos_norm = PLANET_RADIUS + LAUNCH_ELEVATION; // sqrt(dot(pos, pos, 3));
  adouble altitude = pos_norm - PLANET_RADIUS;
  adouble pressure = calc_pressure(altitude, PLANET_P_0, PLANET_SCALE_HEIGHT);
  adouble isp = get_isp(pressure, StageParameters[iphase-1][SP_ISP_0],
						StageParameters[iphase-1][SP_ISP_VAC]);

  adouble thrust_norm = norm(controls[CO_THRX], controls[CO_THRY], controls[CO_THRZ]);

  derivatives[ST_MASS] = -thrust_norm / (isp * G_0) * 1000;

  //cout << "time " << time << "(" << iphase << ")" << " thrust " << thrust_norm << " dMass " << derivatives[ST_MASS] << " isp " << isp << endl;
  adouble time_start = get_initial_time(xad, iphase);
  // path restriction on thrust
  path[PA_THRUST] = thrust_norm;// * (1 + 1 / (1 + (time - time_start)/5) * sin(time/5));

}

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase)
{

  // Launch Parameters

  if (iphase == 1) {
	e[E1_L_MASS] = initial_states[ST_MASS];

	adouble controls[3];
	get_initial_controls(controls, xad, 1);
  }

  if (iphase == STAGES) {
	e[EF_MASS_F] = final_states[ST_MASS];
  }
}

void linkages( adouble* linkages, adouble* xad)
{
    int index=0;

	adouble time_prev, time_next;
	adouble stat_prev[7], stat_next[7];
	int i;
	for (i=1;i<STAGES;i++) {
	  // get states
	  time_prev = get_final_time(xad, i);
	  time_next = get_initial_time(xad, i+1);
	  get_final_states(stat_prev, xad, i);
	  get_initial_states(stat_next, xad, i+1);

	  // time
	  linkages[index++] = time_prev - time_next;
	  // mass
	  adouble mass_difference = StageParameters[i-1][SP_MASS]-StageParameters[i-1][SP_PROPELLANT] - StageParameters[i][SP_MASS];
	  linkages[index++] = stat_prev[ST_MASS]-mass_difference-stat_next[ST_MASS];
	}

}

int main(void)
{

   Alg  algorithm;
   Sol  solution;
   Prob problem;

   // Level 1 Setup

   problem.name = "KSP Launch Optimization";
   problem.outfilename = "launch.txt";
   problem.nphases = STAGES;
   problem.nlinkages = (STAGES - 1) * (ST_NUMBER + 1);
   psopt_level1_setup(problem);

   // Level 2 Setup

   assert(STAGES > 1);

   problem.phases(1).nstates   = ST_NUMBER;
   problem.phases(1).ncontrols = CO_NUMBER;
   problem.phases(1).nevents   = E1_NUMBER;
   problem.phases(1).npath     = PA_NUMBER;
   problem.phases(1).nodes     = "[5, 50]";

   int iphase;

   for (iphase=2;iphase<STAGES;iphase++) {
	 problem.phases(iphase).nstates   = ST_NUMBER;
	 problem.phases(iphase).ncontrols = CO_NUMBER;
	 problem.phases(iphase).nevents   = 0;
	 problem.phases(iphase).npath     = PA_NUMBER;
	 problem.phases(iphase).nodes     = "[5, 50]";
   }

   problem.phases(STAGES).nstates   = ST_NUMBER;
   problem.phases(STAGES).ncontrols = CO_NUMBER;
   problem.phases(STAGES).nevents   = EF_NUMBER;
   problem.phases(STAGES).npath     = PA_NUMBER;
   problem.phases(STAGES).nodes     = "[5, 50]";

   psopt_level2_setup(problem, algorithm);

   // Problem bounds


   // Times
   DMatrix ltM = DMatrix(STAGES + 1);
   DMatrix utM = DMatrix(STAGES + 1);
   double ltd = 0;
   ltM(1) = 0.0;
   utM(1) = 0.0;
   for (iphase=0;iphase<STAGES;iphase++) {
	 ltd += StageParameters[iphase][SP_MIN_STAGE_TIME];
	 cout << ltd << endl;
	 ltM(iphase+2) = ltd;
	 utM(iphase+2) = ltd * STAGE_TIME_MAX_FACTOR;
   }

   problem.bounds.lower.times = ltM;
   problem.bounds.upper.times = utM;


   for(iphase = 1;iphase <= STAGES;iphase++) {
	 problem.phases(iphase).bounds.lower.states(BI(ST_MASS)) = StageParameters[iphase-1][SP_MASS] - StageParameters[iphase-1][SP_PROPELLANT];
	 problem.phases(iphase).bounds.upper.states(BI(ST_MASS)) = StageParameters[iphase-1][SP_MASS];

	 problem.phases(iphase).bounds.lower.controls(BI(CO_THRX)) = 0;
	 problem.phases(iphase).bounds.upper.controls(BI(CO_THRX)) = StageParameters[iphase-1][SP_THRUST];
	 problem.phases(iphase).bounds.lower.controls(BI(CO_THRY)) = 0;
	 problem.phases(iphase).bounds.upper.controls(BI(CO_THRY)) = StageParameters[iphase-1][SP_THRUST];
	 problem.phases(iphase).bounds.lower.controls(BI(CO_THRZ)) = 0;
	 problem.phases(iphase).bounds.upper.controls(BI(CO_THRZ)) = StageParameters[iphase-1][SP_THRUST];

	 problem.phases(iphase).bounds.lower.path(BI(PA_THRUST)) = 0;
	 problem.phases(iphase).bounds.upper.path(BI(PA_THRUST)) = StageParameters[iphase-1][SP_THRUST];

	 if (iphase == STAGES) { // at the end of the final stage the mass can be larger than total mass - propellant mass
	 }

   }

   // Events

   problem.phases(1).bounds.lower.events(BI(E1_L_MASS)) = StageParameters[0][SP_MASS];
   problem.phases(1).bounds.upper.events(BI(E1_L_MASS)) = StageParameters[0][SP_MASS];

   problem.phases(STAGES).bounds.lower.events(BI(EF_MASS_F)) = StageParameters[STAGES-1][SP_MASS] - StageParameters[STAGES-1][SP_PROPELLANT];
   problem.phases(STAGES).bounds.upper.events(BI(EF_MASS_F)) = StageParameters[STAGES-1][SP_MASS] - StageParameters[STAGES-1][SP_PROPELLANT];

   // Guesses

   int i;
   for (i=0;i<STAGES;i++) {
	 problem.phases(i+1).guess.controls = zeros(3,5);
	 problem.phases(i+1).guess.controls(1,colon()) = StageParameters[i][SP_THRUST]*ones( 1, 5);
	 problem.phases(i+1).guess.controls(2,colon()) = zeros(1, 5);
	 problem.phases(i+1).guess.controls(3,colon()) = zeros(1, 5);
   }

   problem.phases(1).guess.states = zeros(ST_NUMBER,5);
   problem.phases(1).guess.states(BI(ST_MASS), colon()) = linspace(StageParameters[0][SP_MASS], StageParameters[0][SP_MASS] - StageParameters[0][SP_PROPELLANT] , 5);

   problem.phases(1).guess.time = linspace(0, 60, 5);
   problem.phases(2).guess.time = linspace(60, 180, 5);




   problem.endpoint_cost = &endpoint_cost;
   problem.dae = &dae;
   problem.events = &events;
   problem.linkages = &linkages;

   algorithm.nlp_method                  	= "IPOPT";
   algorithm.scaling                     	= "automatic";
   algorithm.derivatives                 	= "automatic";
   algorithm.nlp_iter_max                	= 500;
   //    algorithm.mesh_refinement                   = "automatic";
   //   algorithm.collocation_method = "trapezoidal";
   algorithm.ode_tolerance			= 2.e-4;

   psopt(solution, problem, algorithm);

   DMatrix x, u, t;

   x = (solution.get_states_in_phase(1) || solution.get_states_in_phase(2));
   u = (solution.get_controls_in_phase(1) || solution.get_controls_in_phase(2));
   t = (solution.get_time_in_phase(1) || solution.get_time_in_phase(2));
   x.Save("x.dat");
   u.Save("u.dat");
   t.Save("t.dat");

   DMatrix r = x(colon(1,3),colon());
   //   DMatrix mass = x(colon(1,1),colon());
   DMatrix v = x(colon(4,6),colon());
   DMatrix altitude = Sqrt(sum(elemProduct(r,r)))/1000.0;
   DMatrix speed = Sqrt(sum(elemProduct(v,v)));
   DMatrix um = u(colon(1,3),colon());
   DMatrix uma = Sqrt(sum(elemProduct(um,um)));

   plot(t,uma,problem.name, const_cast<char *>("time(s)"), const_cast<char *>("thrust (kN)"));
   plot(t,x,problem.name, const_cast<char *>("time(s)"), const_cast<char *>("Mass (kg)"));
//   plot(t,altitude,problem.name, "time (s)", "position (km)");
//   plot(t,speed,problem.name, "time (s)", "speed (m/s)");
//   plot(t,u,problem.name,"time (s)", "u");
//   plot(t,altitude,problem.name, "time (s)", "position (km)", "alt",
//		"pdf", "launch_position.pdf");
//   plot(t,speed,problem.name, "time (s)", "speed (m/s)", "speed",
//		"pdf", "launch_speed.pdf");
//   plot(t,u,problem.name,"time (s)", "u", "u1 u2 u3",
//		"pdf", "launch_control.pdf");

}
