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
#define X0 -1
#define Y0 2
#define PHI0 2.475859696
#define NUM_OF_ELLIPSE 100
#define NUM_OF_POINT 100


int main(){

    FILE *ellipse_point_file, *x_finish_point_file, *cartesian_point_file, *polar_point_file;
    ellipse_point_file = fopen("../data/ellipses.txt", "w");
    x_finish_point_file = fopen("../data/x_finish.txt", "w");
    cartesian_point_file = fopen("../data/max.txt", "w");
    polar_point_file = fopen("../data/max_polar.txt", "w");

    char output_type = 2; // 0 - XY, 1 - RPhi, != 0 or != 1 both

    struct ellipse_parameters tmp;
    tmp = planar_initialize_first_ep(VX0, VY0, X0, Y0);

    printf("%.10f \n", tmp.Energy);

    for (int i = 0; i < NUM_OF_ELLIPSE; i++)
    {
        fprintf(x_finish_point_file, "%.6f \n", tmp.x_finish);
        tmp = planar_next_step(&tmp);
        save_max_of_ellipse_2D(&tmp, cartesian_point_file, polar_point_file, output_type);
        save_m_points_of_ellipse_2D(&tmp, NUM_OF_POINT, ellipse_point_file);
    }
    
    return 0;
}