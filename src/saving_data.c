#include "saving_data.h"

void save_m_points_of_ellipse_2D(struct ellipse_parameters* ep, int m, FILE *ellipse_point_file){
    double r, phi;

    for (int i = 0; i < m; i++)
    {
        phi = (ep -> phi_start) + i * ((ep -> phi_finish) - (ep -> phi_start)) / (m-1);
        r = (ep -> p) / (1 + (ep -> e) * cos(phi - (ep -> theta)));
        fprintf(ellipse_point_file, "%.6f %.6f \n", r * cos(phi), r * sin(phi));
    }
}

void save_max_of_ellipse_2D(
    struct ellipse_parameters* ep, 
    FILE *cartesian_point_file, 
    FILE *polar_point_file,
    char output_type){

    double a, b, c, d;
    double x1, x2, y1, y2, x, y;
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

    x = x1 ? y1>0 : x2;
    y = y1 ? y1>0 : y2;

    if (output_type == 0){
        fprintf(cartesian_point_file, "%.6f %.6f \n", x, y);
    } else if (output_type == 1)
    {
        fprintf(polar_point_file, "%.6f %.6f \n", pow(x1 * x1 + y1 * y1, 0.5), phi_from_xy(x, y));
    } else {
        fprintf(cartesian_point_file, "%.6f %.6f \n", x, y);
        fprintf(polar_point_file, "%.6f %.6f \n", pow(x1 * x1 + y1 * y1, 0.5), phi_from_xy(x, y));
    }
    
}


void save_m_points_of_ellipse_3D(
    struct ellipse_parameters_3D* ep, 
    int m, 
    FILE *ellipse_point_file){

    double S[SPACE_DIM][SPACE_DIM];
    double x_vec_pc[SPACE_DIM];
    double x_vec_sc[SPACE_DIM];

    x_vec_pc[2] = 0;

    for (int i = 0; i < SPACE_DIM; i++)
    {
        S[i][0] = ep -> ex[i];
        S[i][1] = ep -> ey[i];
        S[i][2] = ep -> ez[i];
    }

    double r, phi;


    for (int i = 0; i < m; i++)
    {
        phi = (ep->ep.phi_start) + i * ((ep->ep.phi_finish) - (ep->ep.phi_start)) / (m-1);
        r = (ep->ep.p) / (1 + (ep->ep.e) * cos(phi - (ep->ep.theta)));
        
        x_vec_pc[0] = r * cos(phi);
        x_vec_pc[1] = r * sin(phi);
        
        mat_vec_prod(x_vec_sc, S, x_vec_pc);
        printf("%.6f %.6f %.6f \n", x_vec_pc[0], x_vec_pc[1], x_vec_pc[2]);


        fprintf(ellipse_point_file, "%.6f %.6f %.6f \n", x_vec_sc[0], x_vec_sc[1], x_vec_sc[2]);
    }

}