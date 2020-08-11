
#include <assert.h>
#include <math.h>
#include "../config.h"
#include "../svd.h"
#include "../svd_math.h"
#include "../tables/arctan_lookup_table.h"
#include "../tables/sin_lookup_table.h"
#include "../tables/cos_lookup_table.h"

floating_point_t m[4][4] ={
    { 31, 77, -11, 26 },
    { -42, 14, 79, -53 },
    { -68, -10, 45, 90 },
    { 34, 16, 38, -19 },
};

fixed_point_t m_fixed[4][4];

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

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            m_fixed[i][j] = convert_to_fixed(m[i][j], SCALE_FACTOR_M);
        }
    }

    int size = 4;
    for (int sweep = 0; sweep < 4; sweep++)
    {
        for (int i = 0; i < size - 1; i++)
        {
            for (int j = i + 1; j < size; j++)
            {
                // i = 0; j = 1;
                // Replicating the first iteration of the sweep
                fixed_point_t m_ij = m_fixed[i][j], m_ji = m_fixed[j][i], m_ii = m_fixed[i][i], m_jj = m_fixed[j][j];

                floating_point_t theta_sum_fixed, theta_diff_fixed;

                // Calculate theta_sum_fixed
                fixed_point_double_t x = fixed_point_div(m_ji + m_ij, m_jj - m_ii);
                theta_sum_fixed = arctan_lookup(x);

                floating_point_t theta_sum_exp = atan((m[j][i] + m[i][j]) / (m[j][j] - m[i][i]));

                // printf("x (flt): %f\n", convert_to_floating(x, 32));
                // printf("x (exp): %f\n", (m[j][i] + m[i][j]) / (m[j][j] - m[i][i]));
                // printf("theta_sum: %f\n", convert_to_floating(theta_sum_fixed, SCALE_FACTOR_ARCTAN));
                // printf("theta_sum (exp): %f", theta_sum_exp);

                // Calculate theta_diff_fixed
                x = fixed_point_div(m_ji - m_ij, m_jj + m_ii);
                theta_diff_fixed = arctan_lookup(x);

                floating_point_t theta_diff_exp = atan((m[j][i] - m[i][j]) / (m[j][j] + m[i][i]));

                // printf("x (flt): %f\n", convert_to_floating(x, 32));
                // printf("x (exp): %f\n", (m[j][i] - m[i][j]) / (m[j][j] + m[i][i]));
                // printf("theta_diff: %f\n", convert_to_floating(theta_diff_fixed, SCALE_FACTOR_ARCTAN));
                // printf("theta_diff (exp): %f\n\n", theta_diff_exp);

                // Ok up to here
                fixed_point_double_t theta_l_fixed, theta_r_fixed;
                theta_l_fixed = ((fixed_point_double_t)theta_sum_fixed - (fixed_point_double_t)theta_diff_fixed) >> 1;
                theta_r_fixed = theta_sum_fixed - theta_l_fixed;

                floating_point_t theta_l_exp = (theta_sum_exp - theta_diff_exp) / 2;
                floating_point_t theta_r_exp = (theta_sum_exp - theta_l_exp);
                // printf("theta_l      : %f\n", convert_to_floating(theta_l_fixed, SCALE_FACTOR_ARCTAN));
                // printf("theta_l (exp): %f\n", theta_l_exp);
                // printf("theta_r      : %f\n", convert_to_floating(theta_r_fixed, SCALE_FACTOR_ARCTAN));
                // printf("theta_r (exp): %f\n", theta_r_exp);

                fixed_point_t sin_theta_l, cos_theta_l, sin_theta_r, cos_theta_r;
                sin_theta_l = sin_lookup(theta_l_fixed);
                sin_theta_r = sin_lookup(theta_r_fixed);
                cos_theta_l = cos_lookup(theta_l_fixed);
                cos_theta_r = cos_lookup(theta_r_fixed);

                floating_point_t sin_theta_l_exp = sin(theta_l_exp);
                floating_point_t sin_theta_r_exp = sin(theta_r_exp);
                floating_point_t cos_theta_l_exp = cos(theta_l_exp);
                floating_point_t cos_theta_r_exp = cos(theta_r_exp);

                printf("sin_theta_l (flt): %f\n", convert_to_floating(sin_theta_l, SCALE_FACTOR_SINCOS));
                printf("sin_theta_l (exp): %f\n", sin_theta_l_exp);
                printf("sin_theta_r (flt): %f\n", convert_to_floating(sin_theta_r, SCALE_FACTOR_SINCOS));
                printf("sin_theta_r (exp): %f\n", sin_theta_r_exp);
                printf("cos_theta_l (flt): %f\n", convert_to_floating(cos_theta_l, SCALE_FACTOR_SINCOS));
                printf("cos_theta_l (exp): %f\n", cos_theta_l_exp);
                printf("cos_theta_r (flt): %f\n", convert_to_floating(cos_theta_r, SCALE_FACTOR_SINCOS));
                printf("cos_theta_r (exp): %f\n", cos_theta_r_exp);
                printf("\n");
            }
        }
    }
}
