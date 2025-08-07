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