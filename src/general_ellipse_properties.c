#include "general_ellipse_properties.h"
#include "constants.h"


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

void def_polar_param(struct ellipse_parameters* ep){
    ep -> p = (ep -> Moment) * (ep->Moment) / G;
    ep ->e = pow(1 + (2 * (ep ->Energy) * (ep ->p) / G), 0.5);

    ep ->theta = ep -> phi_start - acos((((ep ->p) / ep ->r_start) - 1) / (ep ->e));

    def_general_param(ep);
    if (fabs(tangent(ep, ep ->x_start, ep ->y_start) - atan( ep -> vy_start / ep -> vx_start)) > TOL){
        ep ->theta = ep -> phi_start  + acos((((ep ->p) / ep ->r_start) - 1) / (ep ->e));
        def_general_param(ep);
    }
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

void cross_prod(double cross[SPACE_DIM], double x1[SPACE_DIM], double x2[SPACE_DIM]){
    cross[0] = x1[1] * x2[2] - x1[2] * x2[1];
    cross[1] = x1[2] * x2[0] - x1[0] * x2[2];
    cross[2] = x1[0] * x2[1] - x1[1] * x2[0];

    double norm = pow(cross[0] * cross[0] + cross[1] * cross[1] + cross[2] * cross[2], 0.5);

    for (int i = 0; i < SPACE_DIM; i++)
    {
        cross[i] /= norm;
    }
    
}

double det(double A[SPACE_DIM][SPACE_DIM]){
    return  A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) - 
            A[0][1] * (A[1][0] * A[2][2] - A[1][2] * A[2][0]) + 
            A[0][2] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]);
}

void inv(double invA[SPACE_DIM][SPACE_DIM], double A[SPACE_DIM][SPACE_DIM]){
    double detM = det(A);

    invA[0][0] = (A[1][1] * A[2][2] - A[1][2] * A[2][1]) / detM;
    invA[0][1] = -(A[1][0] * A[2][2] - A[1][2] * A[2][0]) / detM;
    invA[0][2] = (A[1][0] * A[2][1] - A[1][1] * A[2][0]) / detM;

    invA[1][0] = -(A[0][1] * A[2][2] - A[0][2] * A[2][1]) / detM;
    invA[1][1] = (A[0][0] * A[2][2] - A[0][2] * A[2][0]) / detM;
    invA[1][2] = -(A[1][0] * A[2][1] - A[1][1] * A[2][0]) / detM;

    invA[2][0] = (A[0][1] * A[1][2] - A[0][2] * A[2][1]) / detM;
    invA[2][1] = -(A[0][0] * A[1][2] - A[0][2] * A[1][0]) / detM;
    invA[2][2] = (A[0][0] * A[1][1] - A[0][1] * A[1][0]) / detM;
    
}

void mat_vec_prod(double y[SPACE_DIM], double A[SPACE_DIM][SPACE_DIM], double x[SPACE_DIM]){
    for (int i = 0; i < SPACE_DIM; i++)
    {
        for (int j = 0; j < SPACE_DIM; j++)
        {
            y[i] += A[i][j] * x[j];
        }
    }
}

void transposition(double A[SPACE_DIM][SPACE_DIM]){
    double tmp;

    for (int i = 0; i < SPACE_DIM; i++)
    {
        for (int j = 0; j < i; j++)
        {
            tmp = A[i][j];
            A[i][j] = A[j][i];
            A[j][i] = tmp;
        }
    }
    
}

double find_intersection(struct ellipse_parameters* ep){
    // y = y0: Ax^2 + (2By0 + D)x + Cy0^2 + Ey0 + F = 0 => a_1x^2 + b_1x + c_1 = 0

    double a_1, b_1, c_1, d, x1, x2, alpha;

    a_1 = ep -> A;
    b_1 = 2 * ep -> B * ep -> y_finish  + ep -> D;
    c_1 = ep -> C * ep -> y_finish * ep -> y_finish + ep -> E * ep -> y_finish + ep -> F;

    d = b_1 * b_1 - 4 * a_1 * c_1;

    x1 = (-b_1 + pow(d, 0.5)) / (2 * a_1);
    x2 = (-b_1 - pow(d, 0.5)) / (2 * a_1);

    if (ep -> x_start * ep -> vy_start - ep -> vx_start * ep -> y_start < 0){
        ep -> x_finish = x1;
    } else {
        ep -> x_finish = x2;
    }

    ep -> phi_finish = phi_from_xy(ep -> x_finish, ep -> y_finish);

    ep -> r_finish = pow(ep -> y_finish * ep -> y_finish + ep -> x_finish * ep -> x_finish, 0.5);
    ep -> v_finish = pow(ep -> v_start * ep -> v_start - 2 * G * ((1 / ep -> r_start) - (1 / ep -> r_finish)), 0.5);

    alpha = tangent(ep, ep -> x_finish, ep -> y_finish);
    ep -> vx_finish = -ep -> v_finish * cos(alpha) * alpha / fabs(alpha);
    ep -> vy_finish = -fabs(ep -> v_finish * sin(alpha));

}