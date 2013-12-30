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

#define SUBDIVISIONS 10
#define NODES "[10, 40]"

// PROBLEM Configuration

#define PLANET_MASS (5.2915793E22) // kg
#define PLANET_RADIUS (600000) // m
#define PLANET_SCALE_HEIGHT (5000) // m
#define PLANET_P_0 (1) // atm
#define PLANET_ROT_PER (21600) // s
#define PLANET_SOI (84159286) // m
#define PLANET_MAX_V (10000) // m/s
#define LAUNCH_LATITUDE (-(0 + 5.0/60 + 49.0/3600)/180 * M_PI) // dimensionless
#define LAUNCH_LONGITUDE 0 // +(-(74.0 + 33.0/60 + 27.0/3600)/180 * M_PI) // dimensionless
#define LAUNCH_ELEVATION (77.1) // m
#define TARGET_PERIAPSIS (75000) // m as altitude above sealevel

#define STAGES (2) // dimensionless number of stages

// STAGE_PARAM contains an array of array of doubles
// the first array contains information about the first stage
//
// The elements of the inner arrays containf the following information
// 1. mass of fuel in kg
// 2. total mass of ship at begin of the stage in kg
// 3. thrust of the engines in kN
// 4. ISP at 1atm in s
// 5. ISP in vacuun in s
// 6. drag coefficient (dimensionless)
// 7. minimal duration of stage in s (will probably be removed)
#define STAGE_PARAM {{4000, 9640, 200.0, 320, 370, 0.2, 58}, {2000, 3590,  50.0, 300, 390, 0.2, 117}}

#define STAGE_TIME_MAX_FACTOR 10 // factor for maximum duration of stages
