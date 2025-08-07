#include "step.h"


struct ellipse_parameters planar_next_step(struct ellipse_parameters* prev_ep){
    struct ellipse_parameters next_ep;
    double phi, r, vx, vy, v;

    phi = prev_ep -> phi_finish;
    r = prev_ep -> r_finish;
    vx = prev_ep -> vx_finish;
    vy = -(prev_ep -> vy_finish);
    v = prev_ep -> v_finish;

    next_ep.phi_start = phi;
    next_ep.r_start = r;
    next_ep.x_start = (prev_ep -> x_finish);
    next_ep.y_start = (prev_ep -> y_finish);
    next_ep.vx_start = vx;
    next_ep.vy_start =  vy;
    next_ep.v_start = v;

    next_ep.y_finish = next_ep.y_start;

    next_ep.Energy = (prev_ep -> Energy);
    next_ep.Moment = M * r * v * sin(atan(vy / vx) - phi);

    next_ep.p = (next_ep.Moment) * (next_ep.Moment) / G;
    next_ep.e = pow(1 + (2 * (next_ep.Energy) * (next_ep.p) / G), 0.5);

    next_ep.theta = phi - acos((((next_ep.p) / r) - 1) / (next_ep.e));
    
    def_general_param(&next_ep);
    // printf("%3.5f %3.5f %3.5f\n",tangent(&next_ep, next_ep.x_start, next_ep.y_start),  atan(vy/vx), fabs(tangent(&next_ep, next_ep.x_start, next_ep.x_start) - atan(vy/vx)));
    if (fabs(tangent(&next_ep, next_ep.x_start, next_ep.y_start) - atan(vy/vx)) > TOL){
        next_ep.theta = phi + acos((((next_ep.p) / r) - 1) / (next_ep.e));
        def_general_param(&next_ep);
    }

    // y = 1: Ax^2 + (2B + D)x + C + E + F = 0 => a_1x^2 + b_1x + c_1 = 0

    double a_1, b_1, c_1, d, x1, x2, alpha;

    a_1 = next_ep.A;
    b_1 = 2 * next_ep.B + next_ep.D;
    c_1 = next_ep.C + next_ep.E + next_ep.F;

    d = b_1 * b_1 - 4 * a_1 * c_1;

    x1 = (-b_1 + pow(d, 0.5)) / (2 * a_1);
    x2 = (-b_1 - pow(d, 0.5)) / (2 * a_1);

    if (next_ep.x_start * next_ep.vy_start - next_ep.vx_start * next_ep.y_start < 0){
        next_ep.x_finish = x1;
    } else {
        next_ep.x_finish = x2;
    }

    next_ep.phi_finish = phi_from_xy(next_ep.x_finish, next_ep.y_finish);

    next_ep.r_finish = pow(next_ep.y_finish * next_ep.y_finish + next_ep.x_finish * next_ep.x_finish, 0.5);
    next_ep.v_finish = pow(next_ep.v_start * next_ep.v_start - 2 * G * ((1 / next_ep.r_start) - (1 / next_ep.r_finish)), 0.5);

    alpha = tangent(&next_ep, next_ep.x_finish, next_ep.y_finish);
    next_ep.vx_finish = -next_ep.v_finish * cos(alpha) * alpha / fabs(alpha);
    next_ep.vy_finish = -fabs(next_ep.v_finish * sin(alpha));

    return next_ep;
}