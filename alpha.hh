/*
  alpha.hh

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

#include "setup.hh"

#define BI(x) ((x)+1)

#define SP_PROPELLANT 0
#define SP_MASS 1
#define SP_THRUST 2
#define SP_ISP_0 3
#define SP_ISP_VAC 4
#define SP_DRAG_COEFFICIENT 5
#define SP_MIN_STAGE_TIME 6
#define SP_NUMBER 7

#define LP_POSX 0
#define LP_POSY 1
#define LP_POSZ 2
#define LP_VELX 3
#define LP_VELY 4
#define LP_VELZ 5
#define LP_NUMBER 6

#define ST_POSX   0
#define ST_POSY   1
#define ST_POSZ   2
#define ST_VELX   3
#define ST_VELY   4
#define ST_VELZ   5
#define ST_MASS   6
#define ST_NUMBER 7

#define PA_DISTANCE2     0
#define PA_THRUST2       1
#define PA_ECCENTRICITY2 2
#define PA_NUMBER        3

#define CO_THRX   0
#define CO_THRY   1
#define CO_THRZ   2
#define CO_NUMBER 3

#define EN_NUMBER       (0)

#define E1_POSX         (0 + EN_NUMBER)
#define E1_POSY         (1 + EN_NUMBER)
#define E1_POSZ         (2 + EN_NUMBER)
#define E1_VELX         (3 + EN_NUMBER)
#define E1_VELY         (4 + EN_NUMBER)
#define E1_VELZ         (5 + EN_NUMBER)
#define E1_NUMBER       (6 + EN_NUMBER)

#define EF_DISTANCE2     (0 + EN_NUMBER)
#define EF_PERIAPSIS     (1 + EN_NUMBER)
#define EF_ECCENTRICITY2 (2 + EN_NUMBER)
#define EF_NUMBER        (3 + EN_NUMBER)

void init_launch_parameters();
void setup_state_constraints(Prob problem, int iphase);
void setup_control_constraints(Prob problem, int iphase);
void setup_event_constraints(Prob problem, int iphase);
void setup_time_constraints(Prob problem);
void setup_linkage_constraints(Prob problem);
adouble norm(adouble x, adouble y, adouble z);
adouble norm2(adouble x, adouble y, adouble z);
adouble calc_pressure(adouble altitude);
adouble get_isp(adouble pressure, adouble isp_0, adouble isp_vac);
adouble get_periapsis(adouble* states);
void calc_eccentricity_vector(adouble* states, adouble* ev);
adouble get_latitude(adouble* pos);
adouble get_longitude(adouble* pos);
DMatrix extend_dmatrix_row(DMatrix matrix, DMatrix row);
adouble get_orientation_pitch(adouble* pos, adouble* dir);
void calc_ground_velocity_vector(adouble* states, adouble* gvv);
//void myplot(DMatrix x, DMatrix y, string name, string xt, string yt, string coord=NULL, string filename);

// as global as constants can get

#define GRAVITATIONAL_CONSTANT (6.674E-11)

// conversion factor for pressure and density
#define CONVERSION_FACTOR (1.2230948554874)

// conversion constant for fuelconsumption of engines
#define G_0 (9.82)

// Constant calculations

#define ATMOSPHERIC_HEIGHT (-log(1 / 1000000.0 *) * PLANET_SCALE_HEIGHT)
#define PLANET_MU (GRAVITATIONAL_CONSTANT * PLANET_MASS)

#define SEMI_MAJOR(POS_NORM, VEL_NORM) (1 / (2 / (POS_NORM) - (VEL_NORM) * (VEL_NORM) / PLANET_MU))

#define CALC_DENSITY(altitude) (calc_pressure(altitude) * CONVERSION_FACTOR)

// Code stuff
#define CC(x) const_cast<char *>(x)
