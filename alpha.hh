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
#define ST_NUMBER 6

#define PA_NUMBER 0

#define CO_THRX   0
#define CO_THRY   1
#define CO_THRZ   2
#define CO_NUMBER 3

#define EN_NUMBER (0)

#define E1_POSX   (0 + EN_NUMBER)
#define E1_POSY   (1 + EN_NUMBER)
#define E1_POSZ   (2 + EN_NUMBER)
#define E1_VELX   (3 + EN_NUMBER)
#define E1_VELY   (4 + EN_NUMBER)
#define E1_VELZ   (5 + EN_NUMBER)
#define E1_NUMBER (6 + EN_NUMBER)

#define EF_NUMBER (0 + EN_NUMBER)

void init_launch_parameters();
void setup_state_constraints(Prob problem, int iphase);
void setup_control_constraints(Prob problem, int iphase);
void setup_event_constraints(Prob problem, int iphase);
void setup_time_constraints(Prob problem);
void setup_linkage_constraints(Prob problem);
