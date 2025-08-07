#include "initialize.h"


struct ellipse_parameters planar_initialize_first_ep(double vx, double vy, double x, double y){
    struct ellipse_parameters first_ep;

    first_ep.x_finish = x;
    first_ep.y_finish = y;
    first_ep.vx_finish = vx;
    first_ep.vy_finish = -vy;

    first_ep.phi_finish = phi_from_xy(x, y);
    first_ep.r_finish = pow(x * x + y * y, 0.5);
    first_ep.v_finish = pow(vx * vx + vy * vy, 0.5);

    first_ep.Energy = M * first_ep.v_finish * first_ep.v_finish / 2 - G / first_ep.r_finish;

    return first_ep;
}