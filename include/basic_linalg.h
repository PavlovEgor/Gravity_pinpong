#ifndef LA
#define LA

#include "math.h"

double phi_from_xy(double x, double y);

void cross_prod(double cross[SPACE_DIM], double x1[SPACE_DIM], double x2[SPACE_DIM]);

void inv(double invA[SPACE_DIM][SPACE_DIM], double A[SPACE_DIM][SPACE_DIM]);

void mat_vec_prod(double y[SPACE_DIM], double A[SPACE_DIM][SPACE_DIM], double x[SPACE_DIM]);

void transposition(double A[SPACE_DIM][SPACE_DIM]);

double mag_vec(double vec[SPACE_DIM]);

#endif