#ifndef svd_math_h
#define svd_math_h
#include <inttypes.h>

#define FIXED_POINT_SCALE_FACTOR 17
#define DOUBLE_PRECISION_FIXED_POINT_SCALE_FACTOR FIXED_POINT_SCALE_FACTOR * 2 + 1

typedef float floating_point_t;
typedef int32_t fixed_point_t;
typedef uint32_t unsigned_fixed_point_t;
typedef int64_t fixed_point_double_t;

fixed_point_t convert_to_fixed_point(floating_point_t floating);
floating_point_t fixed_point_mult(floating_point_t LHS, floating_point_t RHS);
#endif