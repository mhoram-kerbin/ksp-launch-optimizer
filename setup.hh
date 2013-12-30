/*
  setup.hh

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

// PSOPT Configuration

#define SUBDIVISIONS 5
#define NODES "[5, 30]"


#define PLANET_MASS (5.2915793E22)
#define PLANET_RADIUS (600000)
#define PLANET_SCALE_HEIGHT (5000)
#define PLANET_P_0 (1)
#define PLANET_ROT_PER (21600)
#define PLANET_SOI (84159286)
#define PLANET_MAX_V (10000)
#define LAUNCH_LATITUDE (-(0 + 5.0/60 + 49.0/3600)/180 * M_PI)
#define LAUNCH_LONGITUDE 0 // +(-(74.0 + 33.0/60 + 27.0/3600)/180 * M_PI)
#define LAUNCH_ELEVATION (77.1)
#define TARGET_PERIAPSIS (75000)

#define STAGES (2)

#define STAGE_PARAM {{4000, 9640, 200.0, 320, 370, 0.2, 58}, {2000, 3590,  50.0, 300, 390, 0.2, 117}}


#define STAGE_TIME_MAX_FACTOR 10
