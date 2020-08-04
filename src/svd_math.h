#ifndef svd_math_h
#define svd_math_h
#include <inttypes.h>

#define FIXED_POINT_SCALE_FACTOR 17
#define DOUBLE_PRECISION_FIXED_POINT_SCALE_FACTOR FIXED_POINT_SCALE_FACTOR * 2 + 1
#define N_ROWS 4
#define N_COLS 4

typedef float floating_point_t;
typedef int32_t fixed_point_t;
typedef uint32_t unsigned_fixed_point_t;
typedef int64_t fixed_point_double_t;

fixed_point_t convert_to_fixed_point(floating_point_t floating);
floating_point_t convert_to_floating_point(floating_point_t fixed);
void convert_mat_to_fixed_point(floating_point_t M_in[N_ROWS][N_COLS], fixed_point_t M_out[N_ROWS][N_COLS]);
void convert_mat_to_floating_point(fixed_point_t M_in[N_ROWS][N_COLS], floating_point_t M_out[N_ROWS][N_COLS]);
floating_point_t fixed_point_mult(fixed_point_t LHS, fixed_point_t RHS);
#endif