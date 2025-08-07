#include "stdio.h"
#include "math.h"
#include "stdlib.h"

#define G 1.0
#define M 1.0
#define TOL 1e-12
#define PI 3.141592653589793
#define NUM_OF_ELLIPSE 10000000
#define VX0 0.2
#define VY0 0.5
#define PHI0 2.475859696
#define NUM_OF_POINT 10

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


void def_general_param(struct ellipse_parameters* ep){
    double a, b, c, theta;
    
    a = (ep -> p) / (1 - (ep -> e) * (ep -> e));
    b = pow(a * (ep -> p), 0.5); 
    c = (ep -> e) * a;
    theta = ep -> theta;

    ep -> A = (cos(theta) * cos(theta) / (a * a)) + (sin(theta) * sin(theta) / (b * b));
    ep -> C = (sin(theta) * sin(theta) / (a * a)) + (cos(theta) * cos(theta) / (b * b));
    ep -> B = -cos(theta) * sin(theta) * (c * c / ((a * a) * (b * b)));
    ep -> D = 2 * cos(theta) * c / (a * a);
    ep -> E = 2 * sin(theta) * c / (a * a);
    ep -> F = (c * c / (a * a)) - 1;
}


double phi_from_xy(double x, double y){
    if (x > 0 && y > 0){
        return atan(y / x);
    }else if ((x < 0 && y > 0) || (x > 0 && y < 0))
    {
        return PI + atan(y / x);
    } else if (x < 0 && y < 0)
    {
        return 2 * PI + atan(y / x);
    }
    
    
}

double tangent(struct ellipse_parameters* ep, double x, double y){
    return atan(- ((ep -> D) + 2 * (ep -> A) * x + 2 * (ep -> B) * y) / ((ep -> E) + 2 * (ep -> C) * y + 2 * (ep -> B) * x));
}


struct ellipse_parameters next_step(struct ellipse_parameters* prev_ep){
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

struct ellipse_parameters init_first_ep(double vx, double vy, double phi){
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

void save_m_points_of_ellipse(struct ellipse_parameters* ep, int m, FILE *file){
    double r, phi;

    for (int i = 0; i < m; i++)
    {
        phi = (ep -> phi_start) + i * ((ep -> phi_finish) - (ep -> phi_start)) / (m-1);
        r = (ep -> p) / (1 + (ep -> e) * cos(phi - (ep -> theta)));
        fprintf(file, "%.6f %.6f \n", r * cos(phi), r * sin(phi));
    }
}

void save_max_of_ellipse(struct ellipse_parameters* ep, FILE *file1, FILE *file2){

    double a, b, c, d;
    double x1, x2, y1, y2;
    double A, B, C, D, E, F;
    A = ep -> A; B = ep -> B; C = ep -> C; D = ep -> D;
    E = ep -> E; F = ep -> F;

    // y' = 0 => y = - (2Ax + D) / 2B => x^2((CA^2/B^2) - A) + x((A/B)((CD/B) - E)) + F + (D/2B)((CD/2B) - E) = 0
    // => x^2 a + x b + c = 0

    a = (C * A * A / (B * B)) - A;
    b = (A / B) * ((C * D / B) - E);
    c = F + (D / (2 * B)) * ((C * D / (2 * B)) - E);
    d = b * b - 4 * a * c;

    x1 = (-b + pow(d, 0.5)) / (2 * a);
    x2 = (-b - pow(d, 0.5)) / (2 * a);

    y1 = - (2 * A * x1 + D) / (2 * B);
    y2 = - (2 * A * x2 + D) / (2 * B);

    if (y1 > 0){
        fprintf(file1, "%.6f %.6f \n", x1, y1);
        fprintf(file2, "%.6f %.6f \n", pow(x1 * x1 + y1 * y1, 0.5), phi_from_xy(x1, y1));
    } else{
        fprintf(file1, "%.6f %.6f \n", x2, y2);
        fprintf(file2, "%.6f %.6f \n", pow(x2 * x2 + y2 * y2, 0.5), phi_from_xy(x2, y2));
    }
}


int main(){

    FILE *file1, *file2, *file3, *file4;
    file1 = fopen("ellipses.txt", "w");
    file2 = fopen("x_finish.txt", "w");
    file3 = fopen("max.txt", "w");
    file4 = fopen("max_polar.txt", "w");

    double E;

    struct ellipse_parameters tmp;
    tmp = init_first_ep(VX0, VY0, PHI0);

    printf("%.10f \n", tmp.Energy);

    for (int i = 0; i < NUM_OF_ELLIPSE; i++)
    {

        fprintf(file2, "%.6f \n", tmp.x_finish);
        tmp = next_step(&tmp);
        save_max_of_ellipse(&tmp, file3, file4);
        save_m_points_of_ellipse(&tmp, NUM_OF_POINT, file1);
    }
    
    return 0;
}