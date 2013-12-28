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

#include <math.h>
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
   problem.nlinkages = 0; // what does this mean?
   psopt_level1_setup(problem);

   // Level 2 Setup

   problem.phases(1).nstates   = 7;
   problem.phases(1).ncontrols = 3;
   problem.phases(1).nevents   = 0;
   problem.phases(1).npath     = 1;
   problem.phases(1).nodes     = "[5, 15]";

   problem.phases(2).nstates   = 7;
   problem.phases(2).ncontrols = 3;
   problem.phases(2).nevents   = 1;
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
	 problem.phases(iphase).bounds.lower.states(BI(ST_MASS)) = StageParameters[iphase-1][S_MASS] - StageParameters[iphase][S_PROPELLANT];
	 problem.phases(iphase).bounds.upper.states(BI(ST_MASS)) = StageParameters[iphase-1][S_MASS];

	 problem.phases(iphase).bounds.lower.controls(BI(C_X)) = 0;
	 problem.phases(iphase).bounds.upper.controls(BI(C_X)) = StageParameters[iphase-1][S_THRUST];
	 problem.phases(iphase).bounds.lower.controls(BI(C_Y)) = 0;
	 problem.phases(iphase).bounds.upper.controls(BI(C_Y)) = StageParameters[iphase-1][S_THRUST];
	 problem.phases(iphase).bounds.lower.controls(BI(C_Z)) = 0;
	 problem.phases(iphase).bounds.upper.controls(BI(C_Z)) = StageParameters[iphase-1][S_THRUST];

	 problem.phases(iphase).bounds.lower.path(BI(P_THRUST)) = 0;
	 problem.phases(iphase).bounds.upper.path(BI(P_THRUST)) = StageParameters[iphase-1][S_THRUST];

   }

   DMatrix x, u, t, H;


   return 0;
}
