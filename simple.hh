/*
  launch.hh

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

#pragma once

#include <math.h>
#include <string>

#define BI(i) ((1 + i))

#define ST_MASS 0
#define ST_NUMBER 1

#define CO_THRX 0
#define CO_THRY 1
#define CO_THRZ 2
#define CO_NUMBER 3

#define SP_PROPELLANT 0
#define SP_MASS 1
#define SP_THRUST 2
#define SP_ISP_0 3
#define SP_ISP_VAC 4
#define SP_DRAG_COEFFICIENT 5
#define SP_MIN_STAGE_TIME 6
#define SP_NUMBER 7

#define E1_L_MASS 0
#define E1_NUMBER 1

#define EF_MASS_F 0
#define EF_NUMBER 1

#define PA_THRUST 0
#define PA_NUMBER 1


// as global as constants can get

#define GRAVITATIONAL_CONSTANT (6.674E-11)

// conversion constant for fuelconsumption of engines
#define G_0 (9.82)
