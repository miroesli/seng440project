#include "svd_math.h"

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
fixed_point_t arctan_lookup(floating_point_t frac)
{
    floating_point_t neg = 1;
    frac = frac * VALUES_IN_RANGE / ARCTAN_RANGE;
    // printf("FRAC: %f\n", frac);
    if (frac < 0) {
        frac = -frac; //or abs(frac)?
        neg = -1;
    }
    // if out of bounds, return signed pi/2
    if (frac >= VALUES_IN_RANGE) {
        return (M_PI / 2) * neg;
    }
    // lookup value in lookup table
    // floating_point_t theta = arctan_lookup_table_old[(uint32_t)(frac)];
    fixed_point_t theta_fixed = arctan_lookup_table[(uint32_t)(frac)];

    // printf("\ntheta float: %f, theta fixed: %d\n", theta, theta_fixed);
    // fixed_point_t theta2 = convert_to_fixed(theta, SCALE_FACTOR_ARCTAN);
    // printf("theta fixed: %d\n", theta2);
    // floating_point_t theta3 = convert_to_floating(theta2, SCALE_FACTOR_ARCTAN);
    // printf("theta float: %f\n\n", theta3);
    return theta_fixed*neg;
    // return convert_to_floating(theta, SCALE_FACTOR_ARCTAN);
}

/**
 * @brief Obtains the sin value of a floating_point_t from a lookup table
 *
 * //TODO need to use fixedpoint modulo?
 *
 * @param frac_fixed - The floating point number to calcualte artan with
 * @return fixed_point_t
 */
fixed_point_t sin_lookup(fixed_point_t theta)
{
    // mod to 0-2pi
    fixed_point_t theta_mod = theta % convert_to_fixed(M_PI*2, SCALE_FACTOR_MOD);
    // convert theta from fixed to floating
    floating_point_t theta_float = convert_to_floating(theta_mod, SCALE_FACTOR_ARCTAN);
    // floating_point_t frac_float = sin_lookup_table_old[(uint32_t)(theta_float)];
    // fixed_point_t frac_converted = convert_to_floating(frac_float, SCALE_FACTOR_U);
    fixed_point_t frac_fixed = sin_lookup_table[(uint32_t)(theta_float)];
    printf("sin: %d, %d\n", theta_mod, frac_fixed);
    return frac_fixed;
}

/**
 * @brief Obtains the cos value of a floating_point_t from a lookup table
 *
 * //TODO need to use fixedpoint modulo?
 *
 * @param frac_fixed - The floating point number to calcualte artan with
 * @return fixed_point_t
 */
fixed_point_t cos_lookup(fixed_point_t theta)
{
    // mod to 0-2pi
    // fixed_point_t theta_mod = theta % convert_to_fixed(M_PI*2);
    // // TODO convert theta from fixed to floating
    // floating_point_t theta_float = convert_to_floating(theta_mod, SCALE_FACTOR_ARCTAN);
    // floating_point_t frac_float = sin_lookup_table_old[(uint32_t)(theta_float)];
    // fixed_point_t frac_converted = convert_to_floating(frac_float, SCALE_FACTOR_U);
    // fixed_point_t frac_fixed = sin_lookup_table[(uint32_t)(theta_float)];
    // printf("sin: %f, %d, %d", frac_float, frac_converted, frac_fixed);
    return theta;
}