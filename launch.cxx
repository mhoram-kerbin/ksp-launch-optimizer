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

static const double StageParameters[STAGES][5] = STAGE_PARAM;
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
	  if (final_states[ST_MASS] == parameters[S_MASS] - parameters[S_PROPELLANT]) {
		adouble periapsis = get_periapsis(final_states);
		return TARGET_PERIAPSIS - periapsis;
	  } else {
		return -final_states[ST_MASS];
	  }
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
  adouble pos_norm = norm(states[ST_POSX], states[ST_POSY], states[ST_POSZ]);
  adouble altitude = pos_norm - parameters[SP_PLANET_RADIUS];
  adouble pressure = calc_pressure(altitude, parameters[SP_PLANET_P_0],
	parameters[SP_PLANET_SCALE_HEIGHT]);

  adouble groundvelocity = ground_velocity
	(states[ST_POSX], states[ST_POSY], states[ST_POSZ],
	 states[ST_VELX], states[ST_VELY], states[ST_VELZ],
	 parameters[SP_PLANET_ROTATION_PERIOD]);
  adouble vel_mod = - 0.5 * pressure * groundvelocity *
	parameters[SP_ROCKET_DRAG] * 0.008 * states[ST_MASS] /
	norm(states[ST_VELX], states[ST_VELY], states[ST_VELZ]);
  adouble Dx = states[ST_VELX] * vel_mod;
  adouble Dy = states[ST_VELY] * vel_mod;
  adouble Dz = states[ST_VELZ] * vel_mod;

  // calculate gravity

  adouble grav_mod = - GRAVITATIONAL_CONSTANT * states[ST_MASS] * parameters[SP_PLANET_MASS] / pow(pos_norm, 3);

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

  // calculate mass change

  adouble isp = get_isp(pressure, parameters[SP_ROCKET_ISP_0],
						parameters[SP_ROCKET_ISP_VAC]);

  adouble thrust_norm = norm(controls[CO_THRX], controls[CO_THRY], controls[CO_THRZ]);

  derivatives[ST_MASS] = thrust_norm / (isp * G_0);

  // path restriction on thrust
  path[P_THRUST] = thrust_norm;

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

}

void events(adouble* e, adouble* initial_states, adouble* final_states,
            adouble* parameters,adouble& t0, adouble& tf, adouble* xad,
            int iphase)
{
  // Mass
  e[E_MASS_0] = initial_states[ST_MASS];
  e[E_MASS_F] = final_states[ST_MASS];

  // Launch Parameters

  if (iphase == 1) {
	e[E1_L_VELZ] = initial_states[ST_VELZ];
	e[E1_L_PVS] = LaunchParameters[L_POSX] * LaunchParameters[L_VELX] +
	  LaunchParameters[L_POSY] * LaunchParameters[L_VELY] +
	  LaunchParameters[L_POSZ] * LaunchParameters[L_VELZ];

	adouble controls[3];
	get_initial_controls(controls, xad, 1);
	e[E1_L_THRUST_XY] = controls[C_X] * initial_states[ST_POSY] - controls[C_Y] * initial_states[ST_POSX];
	e[E1_L_THRUST_XZ] = controls[C_X] * initial_states[ST_POSZ] - controls[C_Z] * initial_states[ST_POSX];
  }

  if (iphase == STAGES) {

	// calculate periapsis
	adouble p[3]; p[0] = final_states[ST_POSX]; p[1] = final_states[ST_POSY]; p[2] = final_states[ST_POSZ];
	adouble v[3]; v[0] = final_states[ST_VELX]; v[1] = final_states[ST_VELY]; v[2] = final_states[ST_VELZ];
	adouble h[3];
	cross(p, v, h);
	adouble ev[3];
	cross(v, h, ev);
	adouble pos_norm = norm(final_states[ST_POSX], final_states[ST_POSY], final_states[ST_POSZ]);
	adouble mu = GRAVITATIONAL_CONSTANT * PLANET_MASS;

	int i;
	for (i=0;i<3;i++) {
	  ev[i] = ev[i] / mu - final_states[ST_POSX + i] / pos_norm;
	}
	adouble e_norm = dot(ev, ev, 3);

	adouble a = 1 / (2 / pos_norm - 2 * pos_norm / mu);

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
	  time_prev = get_final_time(xad, i);
	  time_next = get_initial_time(xad, i+1);
	  get_final_states(stat_prev, xad, i);
	  get_initial_states(stat_next, xad, i+1);

	  // time
	  linkages[index++] = time_prev - time_next;
	  // position and velocity
	  for (j=0;j<6;j++) {
		linkages[index++] = stat_prev[j]-stat_next[j];
	  }
	  // mass
	  adouble mass_difference = StageParameters[i-1][S_MASS]-StageParameters[i-1][S_PROPELLANT] - StageParameters[i][S_MASS];
	  linkages[index++] = stat_prev[ST_MASS]-mass_difference-stat_next[ST_MASS];
	}

	/*
    double m_tot_first   = 104380.0;
    double m_prop_first  = 95550.0;
    double m_dry_first   = m_tot_first-m_prop_first;
    double m_tot_srb     = 19290.0;
    double m_prop_srb    = 17010.0;
    double m_dry_srb     = m_tot_srb-m_prop_srb;

    int index=0;

    auto_link(linkages, &index, xad, 1, 2 );
    linkages[index-2]-= 6*m_dry_srb;
    auto_link(linkages, &index, xad, 2, 3 );
    linkages[index-2]-= 3*m_dry_srb;
    auto_link(linkages, &index, xad, 3, 4 );
    linkages[index-2]-= m_dry_first;
	*/
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
   problem.nlinkages = 8;
   psopt_level1_setup(problem);

   // Level 2 Setup

   problem.phases(1).nstates   = 7;
   problem.phases(1).ncontrols = 3;
   problem.phases(1).nevents   = 6;
   problem.phases(1).npath     = 1;
   problem.phases(1).nodes     = "[5, 15]";

   problem.phases(2).nstates   = 7;
   problem.phases(2).ncontrols = 3;
   problem.phases(2).nevents   = 3;
   problem.phases(2).npath     = 1;
   problem.phases(2).nodes     = "[5, 15]";

   psopt_level2_setup(problem, algorithm);

   // Problem bounds

   problem.bounds.lower.times = "[0,   60.0,   178.0]";
   problem.bounds.upper.times = "[0, 1000.0, 2000.0]";

   int iphase;

   for(iphase = 1;iphase <= STAGES;iphase++) {
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
	 problem.phases(iphase).bounds.lower.states(BI(ST_MASS)) = StageParameters[iphase-1][S_MASS] - StageParameters[iphase-1][S_PROPELLANT];
	 problem.phases(iphase).bounds.upper.states(BI(ST_MASS)) = StageParameters[iphase-1][S_MASS];

	 problem.phases(iphase).bounds.lower.controls(BI(C_X)) = 0;
	 problem.phases(iphase).bounds.upper.controls(BI(C_X)) = StageParameters[iphase-1][S_THRUST];
	 problem.phases(iphase).bounds.lower.controls(BI(C_Y)) = 0;
	 problem.phases(iphase).bounds.upper.controls(BI(C_Y)) = StageParameters[iphase-1][S_THRUST];
	 problem.phases(iphase).bounds.lower.controls(BI(C_Z)) = 0;
	 problem.phases(iphase).bounds.upper.controls(BI(C_Z)) = StageParameters[iphase-1][S_THRUST];

	 problem.phases(iphase).bounds.lower.path(BI(P_THRUST)) = 0;
	 problem.phases(iphase).bounds.upper.path(BI(P_THRUST)) = StageParameters[iphase-1][S_THRUST];

	 problem.phases(iphase).bounds.lower.events(BI(E_MASS_0)) = StageParameters[iphase-1][S_MASS];
	 problem.phases(iphase).bounds.upper.events(BI(E_MASS_0)) = StageParameters[iphase-1][S_MASS];
	 problem.phases(iphase).bounds.lower.events(BI(E_MASS_F)) = StageParameters[iphase-1][S_MASS] - StageParameters[iphase-1][S_PROPELLANT];

	 if (iphase < STAGES) { // at the end of the final stage the mass can be larger than total mass - propellant mass
	   problem.phases(iphase).bounds.upper.events(BI(E_MASS_F)) = StageParameters[iphase-1][S_MASS] - StageParameters[iphase-1][S_PROPELLANT];
	 } else {
	   problem.phases(iphase).bounds.upper.events(BI(E_MASS_F)) = StageParameters[iphase-1][S_MASS];
	 }

   }

   problem.phases(1).bounds.lower.events(BI(E1_L_VELZ)) = 0;
   problem.phases(1).bounds.upper.events(BI(E1_L_VELZ)) = 0;
   problem.phases(1).bounds.lower.events(BI(E1_L_PVS)) = 0;
   problem.phases(1).bounds.upper.events(BI(E1_L_PVS)) = 0;
   problem.phases(1).bounds.lower.events(BI(E1_L_THRUST_XY)) = 0;
   problem.phases(1).bounds.upper.events(BI(E1_L_THRUST_XY)) = 0;
   problem.phases(1).bounds.lower.events(BI(E1_L_THRUST_XZ)) = 0;
   problem.phases(1).bounds.upper.events(BI(E1_L_THRUST_XZ)) = 0;

   problem.phases(STAGES).bounds.lower.events(BI(EF_PERIAPSIS)) = TARGET_PERIAPSIS;
   problem.phases(STAGES).bounds.upper.events(BI(EF_PERIAPSIS)) = PLANET_SOI;


   //   problem.integrand_cost = &integrand_cost;
   problem.endpoint_cost = &endpoint_cost;
   problem.dae = &dae;
   problem.events = &events;
   problem.linkages = &linkages;

   algorithm.nlp_method                  	= "IPOPT";
   algorithm.scaling                     	= "automatic";
   algorithm.derivatives                 	= "automatic";
   algorithm.nlp_iter_max                	= 500;
//    algorithm.mesh_refinement                   = "automatic";
//    algorithm.collocation_method = "trapezoidal";
   algorithm.ode_tolerance			= 1.e-2;

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
   DMatrix altitude = Sqrt(sum(elemProduct(r,r)))/1000.0;
   DMatrix speed = Sqrt(sum(elemProduct(v,v)));

   plot(t,altitude,problem.name, "time (s)", "position (km)");
   plot(t,speed,problem.name, "time (s)", "speed (m/s)");
   plot(t,u,problem.name,"time (s)", "u");
   plot(t,altitude,problem.name, "time (s)", "position (km)", "alt",
		"pdf", "launch_position.pdf");
   plot(t,speed,problem.name, "time (s)", "speed (m/s)", "speed",
		"pdf", "launch_speed.pdf");
   plot(t,u,problem.name,"time (s)", "u", "u1 u2 u3",
		"pdf", "launch_control.pdf");

}
