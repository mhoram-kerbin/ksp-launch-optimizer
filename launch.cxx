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

#include "Vector.hh"

#define GRAVITATIONAL_CONSTANT (6.674E-11)

#define ST_POSX 0
#define ST_POSY 1
#define ST_POSZ 2
#define ST_VELX 3
#define ST_VELY 4
#define ST_VELZ 5
#define ST_MASS 6

#define SP_PLANET_MASS 0
#define SP_PLANET_RADIUS 1
#define SP_PLANET_SCALE_HEIGHT 2
#define SP_PLANET_P_0 3
#define SP_PLANET_ROTATION_PERIOD 4
#define SP_PLANET_SOI 5
#define SP_ROCKET_DRAG 6

#define CO_THRX 0
#define CO_THRY 1
#define CO_THRZ 2

adouble norm(adouble x, adouble y, adouble z)
{
  return sqrt(x*x + y*y + z*z);
}

adouble pressure(adouble altitude, adouble p_0, adouble psh)
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
  adouble pressure = pressure(altitude, parameters[SP_PLANET_P_0],
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

  adouble grav_mod = - GRAVITATIONAL_CONSTANT * states[ST_MASS] * paremeters[SP_PLANET_MASS] / pos_norm ** 3;

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

  derivatives[ST_MASS] = 0;
}

int main(void)
{

  return 0;
}
