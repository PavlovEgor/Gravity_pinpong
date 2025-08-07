#include "step.h"


struct ellipse_parameters planar_next_step(struct ellipse_parameters* prev_ep){
    struct ellipse_parameters next_ep;

    next_ep.x_start     = (prev_ep -> x_finish);
    next_ep.y_start     = (prev_ep -> y_finish);
    next_ep.vx_start    = prev_ep -> vx_finish;
    next_ep.vy_start    = -(prev_ep -> vy_finish);
    next_ep.v_start     = prev_ep -> v_finish;

    next_ep.phi_start   = prev_ep -> phi_finish;
    next_ep.r_start     = prev_ep -> r_finish;

    next_ep.y_finish    = next_ep.y_start;

    next_ep.Energy      = (prev_ep -> Energy);
    next_ep.Moment      = M * (next_ep.r_start) * (next_ep.v_start) * sin(atan(next_ep.vy_start / next_ep.vx_start) - next_ep.phi_start);

    def_polar_param(&next_ep);

    find_intersection(&next_ep);

    return next_ep;
}

struct ellipse_parameters_3D space_next_step(struct ellipse_parameters_3D* prev_ep){
    struct ellipse_parameters_3D next_ep;

    next_ep.x_start     =   (prev_ep -> x_finish);
    next_ep.y_start     =   (prev_ep -> y_finish);
    next_ep.z_start     =   (prev_ep -> z_finish);
    next_ep.vx_start    =   (prev_ep -> vx_finish);
    next_ep.vy_start    =   (prev_ep -> vy_finish);
    next_ep.vz_start    =  -(prev_ep -> vz_finish);
    next_ep.v_start     =   (prev_ep -> v_finish);

    double vec_x[SPACE_DIM] = {next_ep.x_start, next_ep.y_start, next_ep.z_start};
    double vec_v[SPACE_DIM] = {next_ep.vx_start, next_ep.vy_start, next_ep.vz_start};
    double n[SPACE_DIM] = {0, 0, 1};

    cross_prod(next_ep.ez, vec_x, vec_v);
    if (next_ep.ez[2] < 0){
        next_ep.ez[0] = -next_ep.ez[0];
        next_ep.ez[1] = -next_ep.ez[1];
        next_ep.ez[2] = -next_ep.ez[2];
    }

    cross_prod(next_ep.ex, next_ep.ez, n);
    
    cross_prod(next_ep.ey, next_ep.ez, next_ep.ex);
    if (next_ep.ey[2] < 0){
        next_ep.ex[0] = -next_ep.ex[0];
        next_ep.ex[1] = -next_ep.ex[1];
        next_ep.ex[2] = -next_ep.ex[2];

        next_ep.ey[0] = -next_ep.ey[0];
        next_ep.ey[1] = -next_ep.ey[1];
        next_ep.ey[2] = -next_ep.ey[2];
    }
    double S[SPACE_DIM][SPACE_DIM];
    double invS[SPACE_DIM][SPACE_DIM];

    

    for (int i = 0; i < SPACE_DIM; i++)
    {
        S[i][0] = next_ep.ex[i];
        S[i][1] = next_ep.ey[i];
        S[i][2] = next_ep.ez[i];
    }
    inv(invS, S);

    double vec_x_in_plane_coord[SPACE_DIM], vec_v_in_plane_coord[SPACE_DIM];
    mat_vec_prod(vec_x_in_plane_coord, invS, vec_x);
    mat_vec_prod(vec_v_in_plane_coord, invS, vec_v);

    next_ep.ep.x_start = vec_x_in_plane_coord[0];
    next_ep.ep.y_start = vec_x_in_plane_coord[1];
    next_ep.ep.vx_start = vec_v_in_plane_coord[0];
    next_ep.ep.vy_start = vec_v_in_plane_coord[1];
    next_ep.ep.v_start = next_ep.v_start;

    next_ep.ep.phi_start = phi_from_xy(next_ep.ep.x_start, next_ep.ep.y_start);
    next_ep.ep.r_start = pow(next_ep.ep.x_start * next_ep.ep.x_start + next_ep.ep.y_start * next_ep.ep.y_start, 0.5);

    next_ep.ep.Energy = prev_ep -> ep.Energy;
    next_ep.ep.Moment = M * next_ep.ep.r_start * next_ep.ep.v_start * sin(atan(next_ep.ep.vy_start / next_ep.ep.vx_start) - next_ep.ep.phi_start);

    def_polar_param(&next_ep.ep);

    find_intersection(&next_ep.ep);

    double  new_vec_x[SPACE_DIM], new_vec_v[SPACE_DIM];
    double  new_vec_x_ps[SPACE_DIM] = {next_ep.ep.x_finish, next_ep.ep.y_finish, 0}, 
            new_vec_v_ps[SPACE_DIM] = {next_ep.ep.vx_finish, next_ep.ep.vy_finish, 0};

    mat_vec_prod(new_vec_x, S, new_vec_x_ps);
    mat_vec_prod(new_vec_v, S, new_vec_v_ps);

    next_ep.x_finish = new_vec_x[0];
    next_ep.y_finish = new_vec_x[1];
    next_ep.z_finish = new_vec_x[2];

    next_ep.vx_finish = new_vec_v[0];
    next_ep.vy_finish = new_vec_v[1];
    next_ep.vz_finish = new_vec_v[2];

    next_ep.v_finish = mag_vec(new_vec_v);
}