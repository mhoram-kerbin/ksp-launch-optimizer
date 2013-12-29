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

#include "launch.hh"
#include "setup.hh"

static const double StageParameters[STAGES][SP_NUMBER] = STAGE_PARAM;
static double LaunchParameters[6];

adouble norm(adouble x, adouble y, adouble z)
{
  return sqrt(x*x + y*y + z*z);
}

adouble calc_pressure(adouble altitude, adouble p_0, adouble psh)
{
  return p_0 * exp(-altitude / psh);
}

adouble ground_velocity(adouble px, adouble py, adouble pz, adouble vx, adouble vy, adouble vz, adouble siderial_rotation_period)
{
  adouble dist = norm(px, py, pz);
  adouble longitude = atan2(py, px);
  adouble latitude = atan2(pz, sqrt(px*px + py*py));


  adouble f = dist * 2 * M_PI * cos(latitude) / siderial_rotation_period;
  adouble gvx = vx - (-sin(longitude) * f);
  adouble gvy = vy - cos(longitude) * f;

  //  cout << "dist " << dist << " px " << px << " py " << py << " pz " << pz << " long " << longitude << " lat " << latitude << " f " << f << " gvx " << gvx << " gvy " << gvy << " vz " << vz << " srp " << siderial_rotation_period << endl;

  return norm(gvx, gvy, vz);
}

adouble get_isp(adouble pressure, adouble isp_0, adouble isp_vac)
{
  adouble real_p = pressure;
  if (real_p > 1) {
	real_p = 1;
  }
  return isp_0 * real_p + isp_vac * (1 - real_p);
}

adouble get_periapsis(adouble* states)
{
	adouble p[3]; p[0] = states[ST_POSX]; p[1] = states[ST_POSY]; p[2] = states[ST_POSZ];
	adouble v[3]; v[0] = states[ST_VELX]; v[1] = states[ST_VELY]; v[2] = states[ST_VELZ];
	adouble h[3];
	cross(p, v, h);
	adouble ev[3];
	cross(v, h, ev);
	adouble pos_norm = norm(states[ST_POSX], states[ST_POSY], states[ST_POSZ]);
	adouble mu = GRAVITATIONAL_CONSTANT * PLANET_MASS;

	int i;
	for (i=0;i<3;i++) {
	  ev[i] = ev[i] / mu - states[ST_POSX + i] / pos_norm;
	}
	adouble e_norm = dot(ev, ev, 3);

	adouble a = 1 / (2 / pos_norm - 2 * pos_norm / mu);

	return a * (1 - e_norm);
}

adouble endpoint_cost(adouble* initial_states, adouble* final_states,
                      adouble* parameters,adouble& t0, adouble& tf,
                      adouble* xad, int iphase)
{

    if (iphase == STAGES) {
	  //	  if (final_states[ST_MASS] == parameters[SP_MASS] - parameters[SP_PROPELLANT]) {
	  //	adouble periapsis = get_periapsis(final_states);
	  //		return TARGET_PERIAPSIS - periapsis;
	  //} else {
		return -final_states[ST_MASS];
		//	  }
	} else {
	  return 0.0;
	}

}

void dae(adouble* derivatives, adouble* path, adouble* states,
		 adouble* controls, adouble* parameters, adouble& time,
		 adouble* xad, int iphase)
{

  // set change of position as velocity

  derivatives[ST_POSX] = states[ST_VELX]; // Px -> Vx
  derivatives[ST_POSY] = states[ST_VELY]; // Py -> Vy
  derivatives[ST_POSZ] = states[ST_VELZ]; // Pz -> Vz


  // calculate drag
  adouble pos[3]; pos[0] = states[ST_POSX]; pos[1] = states[ST_POSY]; pos[2] = states[ST_POSZ];
  adouble pos_norm = sqrt(dot(pos, pos, 3));
  adouble altitude = pos_norm - PLANET_RADIUS;
  adouble pressure = calc_pressure(altitude, PLANET_P_0, PLANET_SCALE_HEIGHT);

  adouble groundvelocity = ground_velocity
	(states[ST_POSX], states[ST_POSY], states[ST_POSZ],
	 states[ST_VELX], states[ST_VELY], states[ST_VELZ],
	 PLANET_ROT_PER);
  adouble vel[3]; vel[0] = states[ST_VELX]; vel[1] = states[ST_VELY]; vel[2] = states[ST_VELZ];
  adouble vel_mod = - 0.5 * pressure * groundvelocity *
	StageParameters[iphase-1][SP_DRAG_COEFFICIENT] * 0.008 * states[ST_MASS] /
	sqrt(dot(vel, vel, 3));
  adouble Dx = states[ST_VELX] * vel_mod;
  adouble Dy = states[ST_VELY] * vel_mod;
  adouble Dz = states[ST_VELZ] * vel_mod;

  // calculate gravity
  adouble grav_mod = - GRAVITATIONAL_CONSTANT * states[ST_MASS] * PLANET_MASS /pow(pos_norm, 3);

  adouble Gx = states[ST_POSX] * grav_mod;
  adouble Gy = states[ST_POSY] * grav_mod;
  adouble Gz = states[ST_POSZ] * grav_mod;

  // calculate change of velocity by application of thrust, drag and
  // gravity forces

  adouble Fx = controls[CO_THRX] + Dx + Gx;
  adouble Fy = controls[CO_THRY] + Dy + Gy;
  adouble Fz = controls[CO_THRZ] + Dz + Gz;


  derivatives[ST_VELX] = Fx / states[ST_MASS];
  derivatives[ST_VELY] = Fy / states[ST_MASS];
  derivatives[ST_VELZ] = Fz / states[ST_MASS];

  //cout << "Time(" << iphase << ") " << time << endl << "Pos     " << pos[0] << " " << pos[1] << " " << pos[2] << " norm " << pos_norm << endl << "Vel     " << states[3] << " " << states[4] << " " << states[5] << endl << "Mass    " << states[ST_MASS] << endl << "GV      " << groundvelocity << endl << "Thrust  " << controls[CO_THRX] << " " << controls[CO_THRY] << " " << controls[CO_THRZ] << endl << "Drag    " << Dx << " " << Dy << " " << Dz <<endl << "Gravity " << Gx << " " << Gy << " " << Gz << endl << "Res " << derivatives[ST_VELX] << " " << derivatives[ST_VELY] << " " << derivatives[ST_VELZ] << endl << endl;
  // calculate mass change

  adouble isp = get_isp(pressure, StageParameters[iphase-1][SP_ISP_0],
						StageParameters[iphase-1][SP_ISP_VAC]);

  adouble thrust_norm = norm(controls[CO_THRX], controls[CO_THRY], controls[CO_THRZ]);

  derivatives[ST_MASS] = -thrust_norm / (isp * G_0) * 1000;

  // path restriction on thrust
  path[PA_THRUST] = thrust_norm;

}

void init_launch_parameters()
{
  double dist = PLANET_RADIUS + LAUNCH_ELEVATION;

  double f = dist * 2 * M_PI * cos(LAUNCH_LATITUDE) / PLANET_ROT_PER;
  LaunchParameters[L_POSX] = dist * cos(LAUNCH_LATITUDE) * cos(LAUNCH_LONGITUDE);
  LaunchParameters[L_POSY] = dist * cos(LAUNCH_LATITUDE) * sin(LAUNCH_LONGITUDE);
  LaunchParameters[L_POSZ] = dist * sin(LAUNCH_LATITUDE);
  LaunchParameters[L_VELX] = -sin(LAUNCH_LONGITUDE) * f;
  LaunchParameters[L_VELY] =  cos(LAUNCH_LONGITUDE) * f;
  LaunchParameters[L_VELZ] = 0;

  //  cout << LAUNCH_LATITUDE << " " << LAUNCH_LONGITUDE << " " << LaunchParameters[0] << " " << LaunchParameters[1] << " " << LaunchParameters[2] << " " << LaunchParameters[3] << " " << LaunchParameters[4] << " " << LaunchParameters[5] << endl;

}

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase)
{
  // initial mass
  e[E_MASS_BEGIN] = initial_states[ST_MASS];

  // Launch Parameters
  if (iphase == 1) {
	e[E1_L_POSX] = initial_states[ST_POSX];
	e[E1_L_POSY] = initial_states[ST_POSY];
	e[E1_L_POSZ] = initial_states[ST_POSZ];
	e[E1_L_VELX] = initial_states[ST_VELX];
	e[E1_L_VELY] = initial_states[ST_VELY];
	e[E1_L_VELZ] = initial_states[ST_VELZ];

	// Scalarproduct of Pos and Vel to check if they are 90 deg apart
	e[E1_L_PVS] = initial_states[ST_POSX] * initial_states[ST_VELX] +
	  initial_states[ST_POSY] * initial_states[ST_VELY] +
	  initial_states[ST_POSZ] * initial_states[ST_VELZ];

	// Launchplace and thrust must point in the same direction
	adouble controls[3];
	get_initial_controls(controls, xad, 1);
	e[E1_L_THRUST_XY] = controls[CO_THRX] * initial_states[ST_POSY] - controls[CO_THRY] * initial_states[ST_POSX];
	e[E1_L_THRUST_XZ] = controls[CO_THRX] * initial_states[ST_POSZ] - controls[CO_THRZ] * initial_states[ST_POSX];
  }

  if (iphase == STAGES) {

	// Mass
	e[EF_MASS_F] = final_states[ST_MASS];

	// Periapsis
	e[EF_PERIAPSIS] = get_periapsis(final_states);
  }
}
void linkages( adouble* linkages, adouble* xad)
{
    int index=0;

	adouble time_prev, time_next;
	adouble stat_prev[7], stat_next[7];
	int i,j;
	for (i=1;i<STAGES;i++) {
	  // get states
	  get_final_states(stat_prev, xad, i);
	  get_initial_states(stat_next, xad, i+1);

	  // time
	  time_prev = get_final_time(xad, i);
	  time_next = get_initial_time(xad, i+1);
	  linkages[index++] = time_prev - time_next;

	  // position and velocity
	  for (j=0;j<6;j++) {
		linkages[index++] = stat_prev[j]-stat_next[j];
	  }
	  // mass
	  adouble mass_difference = StageParameters[i-1][SP_MASS]-StageParameters[i-1][SP_PROPELLANT] - StageParameters[i][SP_MASS];
	  linkages[index++] = stat_prev[ST_MASS]-mass_difference-stat_next[ST_MASS];
	  //cout << "md " << mass_difference << endl;
	}

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
   problem.nlinkages = (STAGES-1) * (ST_NUMBER + 1);
   psopt_level1_setup(problem);

   // Level 2 Setup

   problem.phases(1).nstates   = ST_NUMBER;
   problem.phases(1).ncontrols = CO_NUMBER;
   problem.phases(1).nevents   = E1_NUMBER;
   problem.phases(1).npath     = PA_NUMBER;
   problem.phases(1).nodes     = "[5, 15]";

   problem.phases(2).nstates   = ST_NUMBER;
   problem.phases(2).ncontrols = CO_NUMBER;
   problem.phases(2).nevents   = EF_NUMBER;
   problem.phases(2).npath     = PA_NUMBER;
   problem.phases(2).nodes     = "[5, 15]";

   psopt_level2_setup(problem, algorithm);

   // Problem bounds

   problem.bounds.lower.times = "[0,   60.0,  178.0]";
   problem.bounds.upper.times = "[0, 1000.0, 2000.0]";

   int iphase;

   for(iphase = 1;iphase <= STAGES;iphase++) {
	 problem.phases(iphase).bounds.lower.states(BI(ST_POSX)) = -PLANET_SOI;
	 problem.phases(iphase).bounds.upper.states(BI(ST_POSX)) = PLANET_SOI;
	 problem.phases(iphase).bounds.lower.states(BI(ST_POSY)) = -PLANET_SOI;
	 problem.phases(iphase).bounds.upper.states(BI(ST_POSY)) = PLANET_SOI;
	 problem.phases(iphase).bounds.lower.states(BI(ST_POSZ)) = -PLANET_SOI;
	 problem.phases(iphase).bounds.upper.states(BI(ST_POSZ)) = PLANET_SOI;
	 problem.phases(iphase).bounds.lower.states(BI(ST_VELX)) = -PLANET_MAX_V;
	 problem.phases(iphase).bounds.upper.states(BI(ST_VELX)) = PLANET_MAX_V;
	 problem.phases(iphase).bounds.lower.states(BI(ST_VELY)) = -PLANET_MAX_V;
	 problem.phases(iphase).bounds.upper.states(BI(ST_VELY)) = PLANET_MAX_V;
	 problem.phases(iphase).bounds.lower.states(BI(ST_VELZ)) = -PLANET_MAX_V;
	 problem.phases(iphase).bounds.upper.states(BI(ST_VELZ)) = PLANET_MAX_V;
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

   int i;
   for (i=1;i<=STAGES;i++) {
	 problem.phases(i).bounds.lower.events(BI(E_MASS_BEGIN)) = StageParameters[i-1][SP_MASS];
	 problem.phases(i).bounds.upper.events(BI(E_MASS_BEGIN)) = StageParameters[i-1][SP_MASS];
   }

   problem.phases(1).bounds.lower.events(BI(E1_L_POSX)) = LaunchParameters[L_POSX];
   problem.phases(1).bounds.upper.events(BI(E1_L_POSX)) = LaunchParameters[L_POSX];
   problem.phases(1).bounds.lower.events(BI(E1_L_POSY)) = LaunchParameters[L_POSY];
   problem.phases(1).bounds.upper.events(BI(E1_L_POSY)) = LaunchParameters[L_POSY];
   problem.phases(1).bounds.lower.events(BI(E1_L_POSZ)) = LaunchParameters[L_POSZ];
   problem.phases(1).bounds.upper.events(BI(E1_L_POSZ)) = LaunchParameters[L_POSZ];
   problem.phases(1).bounds.lower.events(BI(E1_L_VELX)) = LaunchParameters[L_VELX];
   problem.phases(1).bounds.upper.events(BI(E1_L_VELX)) = LaunchParameters[L_VELX];
   problem.phases(1).bounds.lower.events(BI(E1_L_VELY)) = LaunchParameters[L_VELY];
   problem.phases(1).bounds.upper.events(BI(E1_L_VELY)) = LaunchParameters[L_VELY];
   problem.phases(1).bounds.lower.events(BI(E1_L_VELZ)) = LaunchParameters[L_VELZ];
   problem.phases(1).bounds.upper.events(BI(E1_L_VELZ)) = LaunchParameters[L_VELZ];
   problem.phases(1).bounds.lower.events(BI(E1_L_PVS)) = 0;
   problem.phases(1).bounds.upper.events(BI(E1_L_PVS)) = 0;
   problem.phases(1).bounds.lower.events(BI(E1_L_THRUST_XY)) = 0;
   problem.phases(1).bounds.upper.events(BI(E1_L_THRUST_XY)) = 0;
   problem.phases(1).bounds.lower.events(BI(E1_L_THRUST_XZ)) = 0;
   problem.phases(1).bounds.upper.events(BI(E1_L_THRUST_XZ)) = 0;

   problem.phases(STAGES).bounds.lower.events(BI(EF_PERIAPSIS)) = TARGET_PERIAPSIS;
   problem.phases(STAGES).bounds.upper.events(BI(EF_PERIAPSIS)) = PLANET_SOI;
   problem.phases(STAGES).bounds.lower.events(BI(EF_MASS_F)) = StageParameters[STAGES-1][SP_MASS] - StageParameters[STAGES-1][SP_PROPELLANT];
   problem.phases(STAGES).bounds.upper.events(BI(EF_MASS_F)) = StageParameters[STAGES-1][SP_MASS];

   // Guesses

   double tim = 0;
   for (i=0;i<STAGES;i++) {
	 problem.phases(BI(i)).guess.controls = zeros(CO_NUMBER,SUBDIVISIONS);
	 problem.phases(BI(i)).guess.controls(1,colon()) = StageParameters[i][SP_THRUST]*ones( 1, SUBDIVISIONS);
	 problem.phases(BI(i)).guess.controls(2,colon()) = zeros(1, SUBDIVISIONS);
	 problem.phases(BI(i)).guess.controls(3,colon()) = zeros(1, SUBDIVISIONS);

	 problem.phases(BI(i)).guess.time = linspace(tim, tim + StageParameters[i][SP_MIN_STAGE_TIME], SUBDIVISIONS);
	 tim += StageParameters[i][SP_MIN_STAGE_TIME];
   }

   problem.phases(1).guess.states = zeros(ST_NUMBER,5);
   problem.phases(1).guess.states(1, colon()) = LaunchParameters[L_POSX]* ones(1, SUBDIVISIONS);
   problem.phases(1).guess.states(2, colon()) = LaunchParameters[L_POSY]* ones(1, SUBDIVISIONS);
   problem.phases(1).guess.states(3, colon()) = LaunchParameters[L_POSZ]* ones(1, SUBDIVISIONS);
   problem.phases(1).guess.states(4, colon()) = LaunchParameters[L_VELX]* ones(1, SUBDIVISIONS);
   problem.phases(1).guess.states(5, colon()) = LaunchParameters[L_VELY]* ones(1, SUBDIVISIONS);
   problem.phases(1).guess.states(6, colon()) = LaunchParameters[L_VELZ]* ones(1, SUBDIVISIONS);
   problem.phases(1).guess.states(7, colon()) = linspace(StageParameters[0][SP_MASS], StageParameters[0][SP_MASS] - StageParameters[0][SP_PROPELLANT] , SUBDIVISIONS);

   problem.phases(2).guess.states = zeros(ST_NUMBER,5);
   problem.phases(2).guess.states(1, colon()) = LaunchParameters[L_POSX]* ones(1, SUBDIVISIONS);
   problem.phases(2).guess.states(2, colon()) = LaunchParameters[L_POSY]* ones(1, SUBDIVISIONS);
   problem.phases(2).guess.states(3, colon()) = LaunchParameters[L_POSZ]* ones(1, SUBDIVISIONS);
   problem.phases(2).guess.states(4, colon()) = LaunchParameters[L_VELX]* ones(1, SUBDIVISIONS);
   problem.phases(2).guess.states(5, colon()) = LaunchParameters[L_VELY]* ones(1, SUBDIVISIONS);
   problem.phases(2).guess.states(6, colon()) = LaunchParameters[L_VELZ]* ones(1, SUBDIVISIONS);
   problem.phases(2).guess.states(7, colon()) = linspace(StageParameters[1][SP_MASS], StageParameters[1][SP_MASS] - StageParameters[1][SP_PROPELLANT], SUBDIVISIONS);




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
   DMatrix v = x(colon(4,6),colon());
   DMatrix mass = x(colon(7,7),colon());
   DMatrix altitude = Sqrt(sum(elemProduct(r,r)))/1000.0;
   DMatrix speed = Sqrt(sum(elemProduct(v,v)));
   DMatrix um = u(colon(1,3),colon());
   DMatrix uma = Sqrt(sum(elemProduct(um,um)));

   plot(t,uma,problem.name, const_cast<char *>("time(s)"), const_cast<char *>("thrust (kN)"));
   plot(t,mass,problem.name, const_cast<char *>("time(s)"), const_cast<char *>("mass (kg)"));
//   plot(t,altitude,problem.name, "time (s)", "position (km)");
//   plot(t,speed,problem.name, "time (s)", "speed (m/s)");
   plot(t,u,problem.name,const_cast<char *>("time (s)"), const_cast<char *>("u"));
//   plot(t,altitude,problem.name, "time (s)", "position (km)", "alt",
//		"pdf", "launch_position.pdf");
//   plot(t,speed,problem.name, "time (s)", "speed (m/s)", "speed",
//		"pdf", "launch_speed.pdf");
//   plot(t,u,problem.name,"time (s)", "u", "u1 u2 u3",
//		"pdf", "launch_control.pdf");

}
