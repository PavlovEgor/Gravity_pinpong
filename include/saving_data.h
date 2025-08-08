#ifndef SAVE
#define SAVE

#include "constants.h"
#include "general_ellipse_properties.h"


void save_m_points_of_ellipse_2D(
    struct ellipse_parameters* ep, 
    int m, 
    FILE *ellipse_point_file);

void save_m_points_of_ellipse_3D(
    struct ellipse_parameters_3D* ep, 
    int m, 
    FILE *ellipse_point_file);

void save_max_of_ellipse_2D(
    struct ellipse_parameters* ep, 
    FILE *cartesian_point_file, 
    FILE *polar_point_file,
    char output_type);

#endif