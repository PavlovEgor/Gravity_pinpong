#include "stdio.h"
#include "math.h"
#include "stdlib.h"
#include "constants.h"
#include "general_ellipse_properties.h"
#include "initialize.h"
#include "saving_data.h"
#include "step.h"


#define VX0 0.2
#define VY0 0.5
#define VZ0 0.7
#define X0 0.5
#define Y0 0.2
#define Z0 1

#define NUM_OF_ELLIPSE 100000
#define NUM_OF_POINT 100


int main(){

    FILE *ellipse_point_file, *x_finish_point_file;
    ellipse_point_file = fopen("../data/ellipses.txt", "w");
    x_finish_point_file = fopen("../data/x_finish.txt", "w");

    struct ellipse_parameters_3D tmp;
    tmp = space_initialize_first_ep(VX0, VY0, VZ0, X0, Y0, Z0);

    printf("Energy = %.10f, must be < 0 \n", tmp.ep.Energy);

    for (int i = 0; i < NUM_OF_ELLIPSE; i++)
    {

        tmp = space_next_step(&tmp);

        fprintf(x_finish_point_file, "%.6f %.6f\n", tmp.x_finish, tmp.y_finish);
        save_m_points_of_ellipse_3D(&tmp, NUM_OF_POINT, ellipse_point_file);
    }
    
    return 0;
}