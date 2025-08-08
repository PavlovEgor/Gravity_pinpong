#ifndef INIT
#define INIT

#include "constants.h"
#include "general_ellipse_properties.h"

struct ellipse_parameters planar_initialize_first_ep(
    double vx, 
    double vy, 
    double x,
    double y);

    
struct ellipse_parameters_3D space_initialize_first_ep(
    double vx, 
    double vy,
    double vz, 
    double x,
    double y, 
    double z);


#endif