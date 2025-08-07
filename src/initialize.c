#include "initialize.h"


struct ellipse_parameters planar_initialize_first_ep(double vx, double vy, double phi){
    struct ellipse_parameters first_ep;

    first_ep.phi_finish = phi;
    first_ep.vx_finish = vx;
    first_ep.vy_finish = -vy;
    first_ep.y_finish = 1;

    first_ep.r_finish = 1 / sin(phi);
    first_ep.x_finish = 1 / tan(phi);
    first_ep.v_finish = pow(vx * vx + vy * vy, 0.5);


    first_ep.Energy = M * first_ep.v_finish * first_ep.v_finish / 2 - G / first_ep.r_finish;

    return first_ep;
}