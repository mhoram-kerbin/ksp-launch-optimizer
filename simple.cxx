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

static const double StageParameters[STAGES][SP_NUMBER] = STAGE_PARAM;
static double LaunchParameters[LP_NUMBER];

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
	  //	  return tf;
	  return -final_states[ST_POSX];
	  return -norm(final_states[ST_POSX], final_states[ST_POSY], final_states[ST_POSZ]);
	} else {
	  return 0.0;
	}

}

void dae(adouble* derivatives, adouble* path, adouble* states,
		 adouble* controls, adouble* parameters, adouble& time,
		 adouble* xad, int iphase)
{

  // Position Derivate
  derivatives[ST_POSX] = states[ST_VELX];
  derivatives[ST_POSY] = states[ST_VELY];
  derivatives[ST_POSZ] = states[ST_VELZ];

  // Velocity Derivative

  adouble Fx = controls[CO_THRX] * 1000;
  adouble Fy = controls[CO_THRY] * 1000;
  adouble Fz = controls[CO_THRZ] * 1000;

  derivatives[ST_VELX] = Fx / states[ST_MASS];
  derivatives[ST_VELY] = Fy / states[ST_MASS];
  derivatives[ST_VELZ] = Fz / states[ST_MASS];


  // Altitude Path

  adouble pos[3]; pos[0] = states[ST_POSX]; pos[1] = states[ST_POSY]; pos[2] = states[ST_POSZ];
  adouble pos_norm = sqrt(dot(pos, pos, 3));
  path[PA_ALTITUDE] = pos_norm;

  // Mass Derivative

  adouble altitude = pos_norm - PLANET_RADIUS;
  adouble pressure = calc_pressure(altitude, PLANET_P_0, PLANET_SCALE_HEIGHT);
  adouble isp = get_isp(pressure, StageParameters[iphase-1][SP_ISP_0],
						StageParameters[iphase-1][SP_ISP_VAC]);

  adouble thrust_norm = norm(controls[CO_THRX], controls[CO_THRY], controls[CO_THRZ]);

  derivatives[ST_MASS] = -thrust_norm / (isp * G_0) * 1000;

  // Thrust Path

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
	e[E1_L_POSX] = initial_states[ST_POSX];
	e[E1_L_POSY] = initial_states[ST_POSY];
	e[E1_L_POSZ] = initial_states[ST_POSZ];

	e[E1_L_VELX] = initial_states[ST_VELX];
	e[E1_L_VELY] = initial_states[ST_VELY];
	e[E1_L_VELZ] = initial_states[ST_VELZ];

	e[E1_L_MASS] = initial_states[ST_MASS];
  }

  if (iphase == STAGES) {
	e[EF_MASS_F] = final_states[ST_MASS];
  }
}

void linkages( adouble* linkages, adouble* xad)
{
    int index=0;

	adouble time_prev, time_next;
	adouble stat_prev[ST_NUMBER], stat_next[ST_NUMBER];
	int i;
	for (i=1;i<STAGES;i++) {
	  // get states
	  time_prev = get_final_time(xad, i);
	  time_next = get_initial_time(xad, i+1);
	  get_final_states(stat_prev, xad, i);
	  get_initial_states(stat_next, xad, i+1);

	  // time
	  linkages[index++] = time_prev - time_next;

	  // position
	  linkages[index++] = stat_prev[ST_POSX]-stat_next[ST_POSX];
	  linkages[index++] = stat_prev[ST_POSY]-stat_next[ST_POSY];
	  linkages[index++] = stat_prev[ST_POSZ]-stat_next[ST_POSZ];

	  // velocity
	  linkages[index++] = stat_prev[ST_VELX]-stat_next[ST_VELX];
	  linkages[index++] = stat_prev[ST_VELY]-stat_next[ST_VELY];
	  linkages[index++] = stat_prev[ST_VELZ]-stat_next[ST_VELZ];

	  // mass
	  adouble mass_difference = StageParameters[i-1][SP_MASS]-StageParameters[i-1][SP_PROPELLANT] - StageParameters[i][SP_MASS];
	  linkages[index++] = stat_prev[ST_MASS]-mass_difference-stat_next[ST_MASS];
	}

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

int main(void)
{

   Alg  algorithm;
   Sol  solution;
   Prob problem;

   init_launch_parameters();

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
   problem.phases(1).nodes     = "[5, 30]";

   int iphase;

   for (iphase=2;iphase<STAGES;iphase++) {
	 exit;
	 problem.phases(iphase).nstates   = ST_NUMBER;
	 problem.phases(iphase).ncontrols = CO_NUMBER;
	 problem.phases(iphase).nevents   = 0;
	 problem.phases(iphase).npath     = PA_NUMBER;
	 problem.phases(iphase).nodes     = "[5, 30]";
   }

   problem.phases(STAGES).nstates   = ST_NUMBER;
   problem.phases(STAGES).ncontrols = CO_NUMBER;
   problem.phases(STAGES).nevents   = EF_NUMBER;
   problem.phases(STAGES).npath     = PA_NUMBER;
   problem.phases(STAGES).nodes     = "[5, 30]";

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

   // Bound Constraints

   for(iphase = 1;iphase <= STAGES;iphase++) {
	 problem.phases(iphase).bounds.lower.states(BI(ST_MASS)) = StageParameters[iphase-1][SP_MASS] - StageParameters[iphase-1][SP_PROPELLANT];
	 problem.phases(iphase).bounds.upper.states(BI(ST_MASS)) = StageParameters[iphase-1][SP_MASS];

	 problem.phases(iphase).bounds.lower.states(BI(ST_POSX)) = 0;
	 problem.phases(iphase).bounds.upper.states(BI(ST_POSX)) = PLANET_SOI;
	 problem.phases(iphase).bounds.lower.states(BI(ST_POSY)) = 0;
	 problem.phases(iphase).bounds.upper.states(BI(ST_POSY)) = PLANET_SOI;
	 problem.phases(iphase).bounds.lower.states(BI(ST_POSZ)) = 0;
	 problem.phases(iphase).bounds.upper.states(BI(ST_POSZ)) = PLANET_SOI;

	 problem.phases(iphase).bounds.lower.states(BI(ST_VELX)) = -PLANET_MAX_V;
	 problem.phases(iphase).bounds.upper.states(BI(ST_VELX)) = PLANET_MAX_V;
	 problem.phases(iphase).bounds.lower.states(BI(ST_VELY)) = -PLANET_MAX_V;
	 problem.phases(iphase).bounds.upper.states(BI(ST_VELY)) = PLANET_MAX_V;
	 problem.phases(iphase).bounds.lower.states(BI(ST_VELZ)) = -PLANET_MAX_V;
	 problem.phases(iphase).bounds.upper.states(BI(ST_VELZ)) = PLANET_MAX_V;

	 problem.phases(iphase).bounds.lower.controls(BI(CO_THRX)) = 0;
	 problem.phases(iphase).bounds.upper.controls(BI(CO_THRX)) = StageParameters[iphase-1][SP_THRUST];
	 problem.phases(iphase).bounds.lower.controls(BI(CO_THRY)) = 0;
	 problem.phases(iphase).bounds.upper.controls(BI(CO_THRY)) = StageParameters[iphase-1][SP_THRUST];
	 problem.phases(iphase).bounds.lower.controls(BI(CO_THRZ)) = 0;
	 problem.phases(iphase).bounds.upper.controls(BI(CO_THRZ)) = StageParameters[iphase-1][SP_THRUST];

	 problem.phases(iphase).bounds.lower.path(BI(PA_THRUST)) = 0;
	 problem.phases(iphase).bounds.upper.path(BI(PA_THRUST)) = StageParameters[iphase-1][SP_THRUST];
	 problem.phases(iphase).bounds.lower.path(BI(PA_ALTITUDE)) = PLANET_RADIUS;
	 problem.phases(iphase).bounds.upper.path(BI(PA_ALTITUDE)) = PLANET_SOI;

   }

   // Events

   problem.phases(1).bounds.lower.events(BI(E1_L_POSX)) = LaunchParameters[LP_POSX];
   problem.phases(1).bounds.upper.events(BI(E1_L_POSX)) = LaunchParameters[LP_POSX];
   problem.phases(1).bounds.lower.events(BI(E1_L_POSY)) = LaunchParameters[LP_POSY];
   problem.phases(1).bounds.upper.events(BI(E1_L_POSY)) = LaunchParameters[LP_POSY];
   problem.phases(1).bounds.lower.events(BI(E1_L_POSZ)) = LaunchParameters[LP_POSZ];
   problem.phases(1).bounds.upper.events(BI(E1_L_POSZ)) = LaunchParameters[LP_POSZ];

   problem.phases(1).bounds.lower.events(BI(E1_L_VELX)) = LaunchParameters[LP_VELX];
   problem.phases(1).bounds.upper.events(BI(E1_L_VELX)) = LaunchParameters[LP_VELX];
   problem.phases(1).bounds.lower.events(BI(E1_L_VELY)) = LaunchParameters[LP_VELY];
   problem.phases(1).bounds.upper.events(BI(E1_L_VELY)) = LaunchParameters[LP_VELY];
   problem.phases(1).bounds.lower.events(BI(E1_L_VELZ)) = LaunchParameters[LP_VELZ];
   problem.phases(1).bounds.upper.events(BI(E1_L_VELZ)) = LaunchParameters[LP_VELZ];

   problem.phases(1).bounds.lower.events(BI(E1_L_MASS)) = StageParameters[0][SP_MASS];
   problem.phases(1).bounds.upper.events(BI(E1_L_MASS)) = StageParameters[0][SP_MASS];


   problem.phases(STAGES).bounds.lower.events(BI(EF_MASS_F)) = StageParameters[STAGES-1][SP_MASS] - StageParameters[STAGES-1][SP_PROPELLANT];
   problem.phases(STAGES).bounds.upper.events(BI(EF_MASS_F)) = StageParameters[STAGES-1][SP_MASS] - StageParameters[STAGES-1][SP_PROPELLANT];

   // Guesses

   int i;
   for (i=0;i<STAGES;i++) {
	 problem.phases(i+1).guess.controls = zeros(CO_NUMBER,SUBDIVISIONS);
	 problem.phases(i+1).guess.controls(BI(CO_THRX),colon()) = StageParameters[i][SP_THRUST]*ones( 1, SUBDIVISIONS);
	 problem.phases(i+1).guess.controls(BI(CO_THRY),colon()) = zeros(1, SUBDIVISIONS);
	 problem.phases(i+1).guess.controls(BI(CO_THRZ),colon()) = zeros(1, SUBDIVISIONS);
   }

   problem.phases(1).guess.states = zeros(ST_NUMBER,SUBDIVISIONS);
   problem.phases(1).guess.states(1, colon()) = LaunchParameters[LP_POSX]* ones(1, SUBDIVISIONS);
   problem.phases(1).guess.states(2, colon()) = LaunchParameters[LP_POSY]* ones(1, SUBDIVISIONS);
   problem.phases(1).guess.states(3, colon()) = LaunchParameters[LP_POSZ]* ones(1, SUBDIVISIONS);
   problem.phases(1).guess.states(4, colon()) = LaunchParameters[LP_VELX]* ones(1, SUBDIVISIONS);
   problem.phases(1).guess.states(5, colon()) = LaunchParameters[LP_VELY]* ones(1, SUBDIVISIONS);
   problem.phases(1).guess.states(6, colon()) = LaunchParameters[LP_VELZ]* ones(1, SUBDIVISIONS);
   problem.phases(1).guess.states(BI(ST_MASS), colon()) = linspace(StageParameters[0][SP_MASS], StageParameters[0][SP_MASS] - StageParameters[0][SP_PROPELLANT] , SUBDIVISIONS);

//   problem.phases(2).guess.states = zeros(ST_NUMBER,SUBDIVISIONS);
//   problem.phases(2).guess.states(1, colon()) = LaunchParameters[LP_POSX]* ones(1, SUBDIVISIONS);
//   problem.phases(2).guess.states(2, colon()) = LaunchParameters[LP_POSY]* ones(1, SUBDIVISIONS);
//   problem.phases(2).guess.states(3, colon()) = LaunchParameters[LP_POSZ]* ones(1, SUBDIVISIONS);
//   problem.phases(2).guess.states(4, colon()) = LaunchParameters[LP_VELX]* ones(1, SUBDIVISIONS);
//   problem.phases(2).guess.states(5, colon()) = LaunchParameters[LP_VELY]* ones(1, SUBDIVISIONS);
//   problem.phases(2).guess.states(6, colon()) = LaunchParameters[LP_VELZ]* ones(1, SUBDIVISIONS);
//   problem.phases(2).guess.states(BI(ST_MASS), colon()) = linspace(StageParameters[1][SP_MASS], StageParameters[1][SP_MASS] - StageParameters[1][SP_PROPELLANT] , SUBDIVISIONS);

   problem.phases(1).guess.time = linspace(0, 60, SUBDIVISIONS);
   problem.phases(2).guess.time = linspace(60, 180, SUBDIVISIONS);




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
   algorithm.ode_tolerance			= 2.e-5;

   psopt(solution, problem, algorithm);

   DMatrix x, u, t;

   x = (solution.get_states_in_phase(1) || solution.get_states_in_phase(2));
   u = (solution.get_controls_in_phase(1) || solution.get_controls_in_phase(2));
   t = (solution.get_time_in_phase(1) || solution.get_time_in_phase(2));
   x.Save("x.dat");
   u.Save("u.dat");
   t.Save("t.dat");

   DMatrix r = x(colon(1,3),colon());
   DMatrix mass = x(colon(7,7),colon());
   DMatrix v = x(colon(4,6),colon());
   DMatrix altitude = (Sqrt(sum(elemProduct(r,r)))-PLANET_RADIUS)/10000;
   DMatrix speed = Sqrt(sum(elemProduct(v,v)));
   DMatrix um = u(colon(1,3),colon());
   DMatrix uma = Sqrt(sum(elemProduct(um,um)));

   plot(t,u,problem.name, const_cast<char *>("time(s)"), const_cast<char *>("Thrust (kN)"));
   plot(t,uma,problem.name, const_cast<char *>("time(s)"), const_cast<char *>("Thrust Sum (kN)"));
   plot(t,mass,problem.name, const_cast<char *>("time(s)"), const_cast<char *>("Mass (kg)"));
   plot(t,v,problem.name, const_cast<char *>("time (s)"), const_cast<char *>("Velocity (m/s)"));
   plot(t,u,problem.name,"time (s)", "u", "u1 u2 u3", "pdf", "launch_control.pdf");
   plot(t,r,problem.name, "time (s)", "position (km)");
//   plot(t,speed,problem.name, "time (s)", "speed (m/s)");
//   plot(t,u,problem.name,"time (s)", "u");
//   plot(t,altitude,problem.name, "time (s)", "position (km)", "alt",
//		"pdf", "launch_position.pdf");
//   plot(t,speed,problem.name, "time (s)", "speed (m/s)", "speed",
//		"pdf", "launch_speed.pdf");
//   plot(t,u,problem.name,"time (s)", "u", "u1 u2 u3",
//		"pdf", "launch_control.pdf");

}
