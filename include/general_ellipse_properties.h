#ifndef EP
#define EP

#include "constants.h"

struct ellipse_parameters
{
    // r = p / (1 + e cos (phi - theta))
    double p;
    double e;
    double theta;

    // Ax^2 + 2Bxy + Cy^2 + Dx + Ey + F = 0
    double A;
    double B;
    double C;
    double D;
    double E;
    double F;

    double Energy;
    double Moment;

    double phi_start;
    double r_start;
    double x_start;
    double y_start;
    double vx_start;
    double vy_start;
    double v_start;

    double phi_finish;
    double r_finish;
    double x_finish;
    double y_finish;
    double vx_finish;
    double vy_finish;
    double v_finish;
};


void def_general_param(struct ellipse_parameters* ep);

double phi_from_xy(double x, double y);

double tangent(struct ellipse_parameters* ep, double x, double y);



#endif