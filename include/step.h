#ifndef STEP
#define STEP

#include "constants.h"
#include "basic_linalg.h"
#include "general_ellipse_properties.h"

struct ellipse_parameters planar_next_step(struct ellipse_parameters* prev_ep);

struct ellipse_parameters_3D space_next_step(struct ellipse_parameters_3D* prev_ep);

#endif