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
    return (fixed_point_double_t)LHS * (fixed_point_double_t)RHS * 8;
}

/**
 * @brief Divides two fixed point_t after scaling the LHS to a fixed_point_double_t
 * 
 * @param LHS 
 * @param RHS 
 * @return fixed_point_double_t 
 */
fixed_point_double_t fixed_point_div(fixed_point_t LHS, fixed_point_t RHS)
{
    if (RHS == 0)
    {
        int sign_LHS = LHS < 0 ? -1 : 1,
            sigh_RHS = RHS < 0 ? -1 : 1;

        return sign_LHS * sigh_RHS == 1 ? 2147483647 : -2147483648;
    }
    return ((fixed_point_double_t)LHS << 32) / RHS;
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
fixed_point_t arctan_lookup(fixed_point_double_t x)
{
    int sign_x = x < 0 ? -1 : 1;
    fixed_point_double_t abs_x = x * sign_x;
    fixed_point_double_t limit = (fixed_point_double_t)ARCTAN_RANGE << 32;

    if (abs_x > limit)
        return convert_to_fixed(M_PI / 2.0, SCALE_FACTOR_ARCTAN) * sign_x;

    size_t idx = (abs_x * VALUES_IN_RANGE / ARCTAN_RANGE) >> 32;
    return arctan_lookup_table[idx] * sign_x;
}

/**
 * @brief Obtains the sin value of a floating_point_t from a lookup table
 *
 * @param frac_fixed - The floating point number to calcualte artan with
 * @return fixed_point_t
 */
fixed_point_t sin_lookup(fixed_point_double_t x)
{
    // printf("x (flt): %f\n", convert_to_floating(x, SCALE_FACTOR_ARCTAN));
    int sign_x = x < 0 ? -1 : 1;
    // abs_X = abs_x * SF
    fixed_point_double_t abs_x = x * sign_x;

    // abx_X / SF = abs_x

    fixed_point_double_t PI_FIXED = M_PI * (floating_point_t)(1 << SCALE_FACTOR_ARCTAN);
    int idx = abs_x * VALUES_IN_RANGE / PI_FIXED;
    return sin_lookup_table[idx] * sign_x;
}

/**
 * @brief Obtains the cos value of a floating_point_t from a lookup table
 *
 * @param frac_fixed - The floating point number to calcualte artan with
 * @return fixed_point_t
 */
fixed_point_t cos_lookup(fixed_point_double_t x)
{
    // printf("x (flt): %f\n", convert_to_floating(x, SCALE_FACTOR_ARCTAN));
    int sign_x = x < 0 ? -1 : 1;
    // abs_X = abs_x * SF
    fixed_point_double_t abs_x = x * sign_x;

    // abx_X / SF = abs_x

    fixed_point_double_t PI_FIXED = M_PI * (floating_point_t)(1 << SCALE_FACTOR_ARCTAN);

    // printf("abs x: %ld\n", abs_x);
    // printf("abs x * N: %ld\n", abs_x * VALUES_IN_RANGE);
    // printf("PI_FIXED: %ld\n", PI_FIXED);

    int idx = abs_x * VALUES_IN_RANGE / PI_FIXED;
    // printf("idx: %d\n", idx);
    return cos_lookup_table[idx];
}