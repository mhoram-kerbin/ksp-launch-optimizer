#pragma once

#include <math.h>

#define BI(i) ((1 + i))

#define ST_POSX 0
#define ST_POSY 1
#define ST_POSZ 2
#define ST_VELX 3
#define ST_VELY 4
#define ST_VELZ 5
#define ST_MASS 6

//#define SP_PLANET_MASS 0
//#define SP_PLANET_RADIUS 1
//#define SP_PLANET_SCALE_HEIGHT 2
//#define SP_PLANET_P_0 3
//#define SP_PLANET_ROTATION_PERIOD 4
//#define SP_PLANET_SOI 5
//#define SP_ROCKET_DRAG 6
//#define SP_ROCKET_ISP_0 7
//#define SP_ROCKET_ISP_VAC 8
//#define SP_ROCKET_THRUST 9

#define CO_THRX 0
#define CO_THRY 1
#define CO_THRZ 2

#define S_PROPELLANT 0
#define S_MASS 1
#define S_THRUST 2
#define S_ISP_0 3
#define S_ISP_VAC 4
#define S_DRAG_COEFFICIENT 5

#define C_X 0
#define C_Y 1
#define C_Z 2

#define P_THRUST 0

#define L_POSX 0
#define L_POSY 1
#define L_POSZ 2
#define L_VELX 3
#define L_VELY 4
#define L_VELZ 5

#define E1_L_POSX 0
#define E1_L_POSY 1
#define E1_L_POSZ 2
#define E1_L_VELX 3
#define E1_L_VELY 4
#define E1_L_VELZ 5
#define E1_L_MASS 6
#define E1_L_PVS 7
#define E1_L_THRUST_XY 8
#define E1_L_THRUST_XZ 9
#define EF_MASS_F 2
#define EF_PERIAPSIS 3
// as global as constants can get

#define GRAVITATIONAL_CONSTANT (6.674E-11)

// conversion constant for fuelconsumption of engines
#define G_0 (9.82)
