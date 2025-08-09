#include "stdio.h"
#include "math.h"
#include "stdlib.h"
#include "constants.h"
#include "general_ellipse_properties.h"
#include "initialize.h"
#include "saving_data.h"
#include "step.h"





int main(int argc, char *argv[]){

    double VX0, VY0, VZ0, X0, Y0, Z0;
    int NUM_OF_ELLIPSE, NUM_OF_POINT;

    if (argc == 1){
        VX0 = 0.2;
        VY0 = 0.5;
        VZ0 = 0.7;
        X0  = 0.5;
        Y0  = 0.2;
        Z0  = 1;

        NUM_OF_ELLIPSE = 10000;
        NUM_OF_POINT = 100;

    } else
    {

        VX0 = atof(argv[1]);
        VY0 = atof(argv[2]);
        VZ0 = atof(argv[3]);

        X0 = atof(argv[4]);
        Y0 = atof(argv[5]);
        Z0 = atof(argv[6]);

        NUM_OF_ELLIPSE  = atoi(argv[7]);
        NUM_OF_POINT    = atoi(argv[8]);
    }
    

    FILE *ellipse_point_file, *x_finish_point_file;
    ellipse_point_file = fopen("../data/ellipses.txt", "w");
    x_finish_point_file = fopen("../data/x_finish.txt", "w");

    struct ellipse_parameters_3D tmp;
    tmp = space_initialize_first_ep(VX0, VY0, VZ0, X0, Y0, Z0);

    printf("Energy = %.10f, must be < 0 \n", tmp.ep.Energy);

    for (int i = 0; i < NUM_OF_ELLIPSE; i++)
    {

        tmp = space_next_step(&tmp);

        fprintf(x_finish_point_file, "%.6f %.6f %.6f %.6f\n", tmp.x_finish, tmp.y_finish, hypot(tmp.x_finish, tmp.y_finish), phi_from_xy(tmp.x_finish, tmp.y_finish));
        // fprintf(x_finish_point_file, "%.6f %.6f\n", tmp.x_finish, tmp.y_finish);
        save_m_points_of_ellipse_3D(&tmp, NUM_OF_POINT, ellipse_point_file);
    }
    
    return 0;
}