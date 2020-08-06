#include "config.h"
#include "svd.h"
#include "svd_math.h"

int main()
{
    floating_point_t m = 1234.56;
    floating_point_t u = 0.654321;
    floating_point_t v = -0.123456;

    fixed_point_m_t m_sp = convert_to_fixed(m, SCALE_FACTOR_M);
    fixed_point_u_t u_sp = convert_to_fixed(u, SCALE_FACTOR_U);
    fixed_point_v_t v_sp = convert_to_fixed(v, SCALE_FACTOR_V);
    printf("m_sp: %d\n", m_sp);
    printf("u_sp: %d\n", u_sp);
    printf("v_sp: %d\n\n", v_sp);

    fixed_point_u_dp_t u_x_u_dp = fixed_point_mul_u_x_u(u_sp, u_sp);
    fixed_point_m_tmp_dp_t u_x_m_dp = fixed_point_mul_u_x_m(u_sp, m_sp);
    fixed_point_m_tmp_t u_x_m_sp = truncate_m_tmp(u_x_m_dp);
    fixed_point_m_dp_t u_x_m_x_v_dp = fixed_point_mul_m_x_v(u_x_m_sp, v_sp);
    fixed_point_v_dp_t v_x_v_dp = fixed_point_mul_v_x_v(v_sp, v_sp);

    printf("u_x_u_dp: %ld\n", u_x_u_dp);
    printf("u_x_m_dp: %ld\n", u_x_m_dp);
    printf("u_x_m_sp: %d\n", u_x_m_sp);
    printf("u_x_m_x_v_dp: %ld\n", u_x_m_x_v_dp);
    printf("v_x_v_dp: %ld\n\n", v_x_v_dp);

    floating_point_t u_x_u_f = convert_to_floating(u_x_u_dp, SCALE_FACTOR_U_DP);
    floating_point_t u_x_m_f = convert_to_floating(u_x_m_dp, SCALE_FACTOR_M_tmp_DP);
    floating_point_t u_x_m_x_v_f = convert_to_floating(u_x_m_x_v_dp, SCALE_FACTOR_M_DP);
    floating_point_t v_x_v_f = convert_to_floating(v_x_v_dp, SCALE_FACTOR_V_DP);

    printf("u_x_u_f: %f u*u: %f\n", u_x_u_f, u * u);
    printf("u_x_m_f: %f u*m: %f\n", u_x_m_f, u * m);
    printf("m_x_v_f: %f u*m*v: %f\n", u_x_m_x_v_f, u * m * v);
    printf("v_x_v_f: %f v*v: %f\n", v_x_v_f, v * v);
}