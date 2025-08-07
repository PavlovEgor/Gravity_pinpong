#include "basic_linalg.h"


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


double mag_vec(double vec[SPACE_DIM]){
    return pow(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2], 0,5);
}