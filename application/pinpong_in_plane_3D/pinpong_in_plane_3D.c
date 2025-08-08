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
#define VZ0 0.2
#define X0 0.5
#define Y0 0.2
#define Z0 1

#define NUM_OF_ELLIPSE 10000
#define NUM_OF_POINT 10


int main(){

    FILE *ellipse_point_file, *x_finish_point_file, *cartesian_point_file, *polar_point_file;
    ellipse_point_file = fopen("../data/ellipses.txt", "w");
    x_finish_point_file = fopen("../data/x_finish.txt", "w");

    char output_type = 0; // 0 - XY, 1 - RPhi, != 0 or != 1 both

    struct ellipse_parameters_3D tmp;
    tmp = space_initialize_first_ep(VX0, VY0, VZ0, X0, Y0, Z0);

    printf("%.10f \n", tmp.ep.Energy);

    for (int i = 0; i < NUM_OF_ELLIPSE; i++)
    {

        tmp = space_next_step(&tmp);

        fprintf(x_finish_point_file, "%.6f %.6f\n", tmp.x_finish, tmp.y_finish);
        save_m_points_of_ellipse_2D(&tmp, NUM_OF_POINT, ellipse_point_file);
    }
    
    return 0;
}