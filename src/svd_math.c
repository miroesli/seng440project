#include "svd_math.h"
#include "stdio.h"
#include "tables.h"

/**
 * @brief Convert a floating_point_t to a fixed_point_t with a given scale factor
 *
 * You have to manually track the scale factors after converting to fixed point.
 *
 * @param floating - The floating point number x
 * @param scale_factor - The scale factor
 * @return fixed_point_t - X = x * SF
 */
fixed_point_t convert_to_fixed(floating_point_t floating, size_t scale_factor)
{
    return floating * (floating_point_t)(1 << scale_factor);
}

/**
 * @brief Multiplies two fixed_point_t to produce a fixed_point_double_t
 *
 * You have to manually track the scale factors when doing this multiplication.
 *
 * @param LHS
 * @param RHS
 * @return fixed_point_u_dp_t
 */
fixed_point_u_dp_t fixed_point_mul(fixed_point_u_t LHS, fixed_point_u_t RHS)
{
    return (fixed_point_double_t)LHS * (fixed_point_double_t)RHS * 2;
}

/**
 * @brief Truncates a fixed_point_double_t to a fixed_point_t
 *
 * @param x - The value to be truncated
 * @return fixed_point_m_tmp_t - The truncated value
 */
fixed_point_t truncate(fixed_point_double_t x)
{
    return x >> 32;
}

/**
 * @brief Converts a fixed_point_double_t into a floating_point_t
 *
 * You have to manually track the scale factor to make sure that this conversion is correct
 *
 * @param f - The fixed point number
 * @param scale_factor - The scale factor
 * @return floating_point_t
 */
floating_point_t convert_to_floating(fixed_point_double_t f, size_t scale_factor)
{
    return f / (floating_point_t)((fixed_point_double_t)1 << scale_factor);
}

/**
 * @brief Obtains the arctan value of a floating_point_t from a lookup table
 *
 * //TODO consider not handling out of bounds?
 *
 * @param frac - The floating point number to calcualte artan with
 * @return floating_point_t
 */
floating_point_t arctan_lookup(floating_point_t frac)
{
    floating_point_t theta;
    printf("FRAC: %f\n", frac);
    if (frac < 0) {
        theta = arctan_lookup_table_old[(uint32_t)(-frac)];
        theta = -theta;
    }
    else {
        theta = arctan_lookup_table_old[(uint32_t)frac];
    }
    return theta;
    // return convert_to_floating(theta, SCALE_FACTOR_ARCTAN);
}

/**
 * @brief Obtains the sin value of a floating_point_t from a lookup table
 *
 * //TODO need to use fixedpoint modulo?
 *
 * @param theta - The floating point number to calcualte artan with
 * @return floating_point_t
 */
floating_point_t sin_lookup(floating_point_t theta)
{
    // theta = theta % (2*M_PI);
    // fixed_point_t frac = sin_lookup_table[(uint32_t)(theta)];
    // return convert_to_floating(frac, SCALE_FACTOR_ARCTAN);
}