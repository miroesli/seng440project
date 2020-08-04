#include "svd_math.h"
#include "stdio.h"

fixed_point_t convert_to_fixed_point(floating_point_t floating)
{
    // printf("Floating Point input: %f\n", floating);
    fixed_point_t fixed = floating * (floating_point_t)(1 << FIXED_POINT_SCALE_FACTOR);
    // printf("Fixed Point output: %d\n", fixed);
    return fixed;
}

floating_point_t fixed_point_mult(floating_point_t x, floating_point_t y)
{
    fixed_point_t LHS = convert_to_fixed_point(x);
    fixed_point_t RHS = convert_to_fixed_point(y);

    int num_bits = sizeof(fixed_point_t) * 8;
    int num_bits_double = sizeof(fixed_point_double_t) * 8;

    // Calculate the sign of the product
    char sign_x = 1 - 2 * ((unsigned)LHS >> (num_bits - 1));
    char sign_y = 1 - 2 * ((unsigned)RHS >> (num_bits - 1));
    char sign_pdp = sign_x * sign_y;
    // Calculate the magnitude of the products
    fixed_point_double_t mag_x = LHS * sign_x;
    fixed_point_double_t mag_y = RHS * sign_y;
    fixed_point_double_t mag_xy = mag_x * mag_y * 2;

    fixed_point_double_t pdp = mag_xy >> (num_bits_double - num_bits);

    floating_point_t pdp_floating = pdp / (floating_point_t)(1 << (DOUBLE_PRECISION_FIXED_POINT_SCALE_FACTOR - (num_bits_double - num_bits))); // Apply scale factor of 2^3
    return pdp_floating * sign_pdp;
}