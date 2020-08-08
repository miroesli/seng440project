#include "config.h"
#include "svd.h"
#include "svd_math.h"

int main()
{
    // floating_point_t m = 1234.56;
    // floating_point_t u = 0.654321;
    // floating_point_t v = -0.123456;

    // fixed_point_m_t m_sp = convert_to_fixed(m, SCALE_FACTOR_M);
    // fixed_point_u_t u_sp = convert_to_fixed(u, SCALE_FACTOR_U);
    // fixed_point_v_t v_sp = convert_to_fixed(v, SCALE_FACTOR_V);
    // printf("m_sp: %d\n", m_sp);
    // printf("u_sp: %d\n", u_sp);
    // printf("v_sp: %d\n\n", v_sp);

    // fixed_point_u_dp_t u_x_u_dp = fixed_point_mul(u_sp, u_sp);
    // fixed_point_m_tmp_dp_t u_x_m_dp = fixed_point_mul(u_sp, m_sp);
    // fixed_point_m_tmp_t u_x_m_sp = truncate(u_x_m_dp);
    // fixed_point_m_dp_t u_x_m_x_v_dp = fixed_point_mul(u_x_m_sp, v_sp);
    // fixed_point_v_dp_t v_x_v_dp = fixed_point_mul(v_sp, v_sp);

    // printf("u_x_u_dp: %ld\n", u_x_u_dp);
    // printf("u_x_m_dp: %ld\n", u_x_m_dp);
    // printf("u_x_m_sp: %d\n", u_x_m_sp);
    // printf("u_x_m_x_v_dp: %ld\n", u_x_m_x_v_dp);
    // printf("v_x_v_dp: %ld\n\n", v_x_v_dp);

    // floating_point_t u_x_u_f = convert_to_floating(u_x_u_dp, SCALE_FACTOR_U_DP);
    // floating_point_t u_x_m_f = convert_to_floating(u_x_m_dp, SCALE_FACTOR_M_tmp_DP);
    // floating_point_t u_x_m_x_v_f = convert_to_floating(u_x_m_x_v_dp, SCALE_FACTOR_M_DP);
    // floating_point_t v_x_v_f = convert_to_floating(v_x_v_dp, SCALE_FACTOR_V_DP);

    // printf("u_x_u_f: %f u*u: %f\n", u_x_u_f, u * u);
    // printf("u_x_m_f: %f u*m: %f\n", u_x_m_f, u * m);
    // printf("m_x_v_f: %f u*m*v: %f\n", u_x_m_x_v_f, u * m * v);
    // printf("v_x_v_f: %f v*v: %f\n", v_x_v_f, v * v);

    // floating_point_t theta = arctan_lookup_table_old[(uint32_t)(frac)];
    // fixed_point_t theta_fixed = arctan_lookup_table[(uint32_t)(theta_sum)];
    // printf("\ntheta float: %f, theta fixed: %d\n", theta, theta_fixed);
    // fixed_point_t theta2 = convert_to_fixed(theta, SCALE_FACTOR_ARCTAN);
    // printf("theta fixed: %d\n", theta2);
    // floating_point_t theta3 = convert_to_floating(theta2, SCALE_FACTOR_ARCTAN);
    // printf("theta float: %f\n\n", theta3);


    fixed_point_t theta_sum = arctan_lookup((floating_point_t)0.3);
    fixed_point_t theta_diff = arctan_lookup((floating_point_t)0.4);
    fixed_point_t theta_l = (theta_sum - theta_diff) >> 1;
    fixed_point_t theta_r = theta_sum - theta_l;
    fixed_point_t sin_theta_l_fixed = sin_lookup(theta_l);
    fixed_point_t cos_theta_l_fixed = cos_lookup(theta_l);
    fixed_point_t sin_theta_r_fixed = sin_lookup(theta_r);
    fixed_point_t cos_theta_r_fixed = cos_lookup(theta_r);
}
