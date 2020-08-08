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
        frac = -frac; //TODO or abs(frac)?
        neg = -1;
    }
    // if out of bounds, return signed pi/2
    if (frac >= VALUES_IN_RANGE) {
        return (M_PI / 2) * neg;
    }
    // lookup value in lookup table
    fixed_point_t theta_fixed = arctan_lookup_table[(uint32_t)(frac)];
    return theta_fixed*neg;
}

/**
 * @brief Obtains the sin value of a floating_point_t from a lookup table
 *
 * @param frac_fixed - The floating point number to calcualte artan with
 * @return fixed_point_t
 */
fixed_point_t sin_lookup(fixed_point_t theta)
{
    floating_point_t neg = 1;
    if (theta < 0) {
        theta = -theta; //TODO or abs(frac)?
        neg = -1;
    }
    // convert theta from fixed to floating
    floating_point_t theta_float = convert_to_floating(theta, SCALE_FACTOR_SINCOS);
    theta_float = theta_float * VALUES_IN_RANGE / SINCOS_RANGE;
    printf("sin theta float: %f\n", theta_float);
    // fixed_point_t frac_fixed = sin_lookup_table[(uint32_t)(theta_float)]>>1;
    fixed_point_t frac_fixed = sin_lookup_table[(uint32_t)(theta_float)];
    return frac_fixed*neg;
}

/**
 * @brief Obtains the cos value of a floating_point_t from a lookup table
 *
 * @param frac_fixed - The floating point number to calcualte artan with
 * @return fixed_point_t
 */
fixed_point_t cos_lookup(fixed_point_t theta)
{
    floating_point_t neg = 1;
    if (theta < 0) {
        theta = -theta; //TODO or abs(frac)?
        neg = -1;
    }
    floating_point_t theta_float = convert_to_floating(theta, SCALE_FACTOR_ARCTAN);
    theta_float = theta_float * VALUES_IN_RANGE / SINCOS_RANGE;
    printf("cos theta float: %f\n", theta_float);
    fixed_point_t frac_fixed = cos_lookup_table[(uint32_t)(theta_float)];
    return frac_fixed*neg;
}