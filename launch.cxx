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

  derivatives[ST_MASS] = norm(controls[CO_THRX], controls[CO_THRY], controls[CO_THRZ]) / (isp * G_0);
}

int main(void)
{

   Alg  algorithm;
   Sol  solution;
   Prob problem;

   // Configuration


   // Level 1 Setup

   problem.name = "KSP Launch Optimization";
   problem.outfilename = "launch.txt";
   problem.nphases = 2;
   problem.nlinkages = 0; // what does this mean?
   psopt_level1_setup(problem);

   // Level 2 Setup

    problem.phases(1).nstates   = 7;
    problem.phases(1).ncontrols = 3;
    problem.phases(1).nevents   = 0;
    problem.phases(1).npath     = 0;
    problem.phases(1).nodes     = "[5, 15]";

    problem.phases(2).nstates   = 7;
    problem.phases(2).ncontrols = 3;
    problem.phases(2).nevents   = 0;
    problem.phases(2).npath     = 0;
    problem.phases(2).nodes     = "[5, 15]";

	psopt_level2_setup(problem, algorithm);

	// Problem bounds

	problem.bounds.lower.times = "[0,   10.0,   20.0]";
    problem.bounds.upper.times = "[0, 1000.0, 2000.0]";

	int iphase;

	iphase = 1;

	problem.phases(iphase).bounds.lower.states(1) = 0;
	problem.phases(iphase).bounds.upper.states(1) = PLANET_SOI;
	problem.phases(iphase).bounds.lower.states(2) = 0;
	problem.phases(iphase).bounds.upper.states(2) = PLANET_SOI;
	problem.phases(iphase).bounds.lower.states(3) = 0;
	problem.phases(iphase).bounds.upper.states(3) = PLANET_SOI;
	problem.phases(iphase).bounds.lower.states(4) = 0;
	problem.phases(iphase).bounds.upper.states(4) = PLANET_MAX_V;
	problem.phases(iphase).bounds.lower.states(5) = 0;
	problem.phases(iphase).bounds.upper.states(5) = PLANET_MAX_V;
	problem.phases(iphase).bounds.lower.states(6) = 0;
	problem.phases(iphase).bounds.upper.states(6) = PLANET_MAX_V;
	problem.phases(iphase).bounds.lower.states(7) = 0;
	problem.phases(iphase).bounds.upper.states(7) = PLANET_MAX_V;




	DMatrix x, u, t, H;


   return 0;
}
