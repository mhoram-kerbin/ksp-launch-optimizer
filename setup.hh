#pragma once

#include <math.h>

#define PLANET_MASS (5.2915793E22)
#define PLANET_RADIUS (600000)
#define PLANET_SCALE_HEIGHT (5000)
#define PLANET_P_0 (1)
#define PLANET_ROT_PER (21600)
#define PLANET_SOI (84159286)
#define PLANET_MAX_V (10000)
#define LAUNCH_LATITUDE (-(0 + 5/60 + 49/3600)/180 * M_PI)
#define LAUNCH_LONGITUDE (-(74 + 33/60 + 27/3600)/180 * M_PI)
#define LAUNCH_ELEVATION (77.1)

#define STAGES (2)

#define STAGE_PARAM {{4000, 9640, 200.0, 320, 370}, {2000, 3590,  50.0, 300, 390}}
