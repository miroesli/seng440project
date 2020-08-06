#include "svd_math.h"
#include "stdio.h"

fixed_point_t convert_to_fixed(floating_point_t floating, size_t scale_factor)
{
    fixed_point_t fixed = floating * (floating_point_t)(1 << scale_factor);
    return fixed;
}

fixed_point_u_dp_t fixed_point_mul(fixed_point_u_t LHS, fixed_point_u_t RHS)
{
    return (fixed_point_double_t)LHS * (fixed_point_double_t)RHS * 2;
}
fixed_point_u_dp_t fixed_point_mul_u_x_u(fixed_point_u_t LHS, fixed_point_u_t RHS) __attribute__((alias("fixed_point_mul")));
fixed_point_m_tmp_dp_t fixed_point_mul_u_x_m(fixed_point_u_t LHS, fixed_point_m_t RHS) __attribute__((alias("fixed_point_mul")));
fixed_point_m_dp_t fixed_point_mul_m_x_v(fixed_point_m_tmp_t LHS, fixed_point_v_t RHS) __attribute__((alias("fixed_point_mul")));
fixed_point_v_dp_t fixed_point_mul_v_x_v(fixed_point_v_t LHS, fixed_point_v_t RHS) __attribute__((alias("fixed_point_mul")));
fixed_point_m_tmp_t truncate_m_tmp(fixed_point_m_tmp_dp_t m)
{
    return m >> 32;
}

floating_point_t convert_to_floating(fixed_point_double_t f, size_t scale_factor)
{
    return f / (floating_point_t)((fixed_point_double_t)1 << scale_factor);
}