#include "initialize.h"


struct ellipse_parameters planar_initialize_first_ep(double vx, double vy, double x, double y){
    struct ellipse_parameters first_ep;

    first_ep.x_finish = x;
    first_ep.y_finish = y;
    first_ep.vx_finish = vx;
    first_ep.vy_finish = -vy;

    first_ep.phi_finish = phi_from_xy(x, y);
    first_ep.r_finish = hypot(x, y);
    first_ep.v_finish = hypot(vx, vy);

    first_ep.Energy = M * first_ep.v_finish * first_ep.v_finish / 2 - G / first_ep.r_finish;

    return first_ep;
}


struct ellipse_parameters_3D space_initialize_first_ep(
    double vx, 
    double vy,
    double vz, 
    double x,
    double y, 
    double z){
    
    struct ellipse_parameters_3D first_ep;

    first_ep.x_finish     =   x;
    first_ep.y_finish     =   y;
    first_ep.z_finish     =   z;
    first_ep.vx_finish    =   vx;
    first_ep.vy_finish    =   vy;
    first_ep.vz_finish    =   -vz;
    first_ep.v_finish     =   sqrt(vx*vx + vy*vy + vz*vz);

    first_ep.ep.Energy    =   M * (vx*vx + vy*vy + vz*vz) / 2 - G / sqrt(x*x + y*y + z*z);
    
    return first_ep;
}