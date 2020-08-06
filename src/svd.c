/*
 *
 * svd.c
 *
 */
#include "svd.h"
#include "svd_math.h"

/**
 * @brief Perfoms a fixed point matrix multiplication
 * 
 * Result is placed in out[][]
 * 
 * @param size 
 * @param LHS 
 * @param RHS 
 * @param out 
 */
void mat_mul(int size, fixed_point_t LHS[size][size], fixed_point_t RHS[size][size], fixed_point_double_t out[size][size])
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            out[i][j] = 0;
            for (int k = 0; k < size; k++)
            {
                out[i][j] += fixed_point_mul(LHS[i][k], RHS[k][j]);
            }
        }
    }
}

void mat_mul_u_x_u(
    int size,
    fixed_point_u_t LHS[size][size],
    fixed_point_u_t RHS[size][size],
    fixed_point_u_dp_t out[size][size]) __attribute__((alias("mat_mul")));

void mat_mul_u_x_m(
    int size,
    fixed_point_u_t LHS[size][size],
    fixed_point_m_t RHS[size][size],
    fixed_point_m_tmp_dp_t out[size][size]) __attribute__((alias("mat_mul")));

void mat_mul_m_x_v(
    int size,
    fixed_point_m_tmp_t LHS[size][size],
    fixed_point_v_t RHS[size][size],
    fixed_point_m_dp_t out[size][size]) __attribute__((alias("mat_mul")));

void mat_mul_v_x_v(
    int size,
    fixed_point_v_t LHS[size][size],
    fixed_point_v_t RHS[size][size],
    fixed_point_v_dp_t out[size][size]) __attribute__((alias("mat_mul")));

/**
 * @brief Performes a single sweep of the svd algorithm
 * 
 */
void sweep(floating_point_t m[4][4], floating_point_t u[4][4], floating_point_t v_trans[4][4])
{

    for (int i = 0; i < 3; i++)
    {
        for (int j = i + 1; j < 4; j++)
        {
            fixed_point_m_t m_fixed[4][4];
            fixed_point_u_t u_fixed[4][4];
            fixed_point_v_t v_trans_fixed[4][4];

            // Convert the input matricies to fixed point.
            for (int row = 0; row < 4; row++)
            {
                for (int col = 0; col < 4; col++)
                {
                    m_fixed[row][col] = convert_to_fixed(m[row][col], SCALE_FACTOR_M);
                    u_fixed[row][col] = convert_to_fixed(u[row][col], SCALE_FACTOR_U);
                    v_trans_fixed[row][col] = convert_to_fixed(v_trans[row][col], SCALE_FACTOR_V);
                }
            }

            /**
             * Do all of the angle calculations
             * 
             * TODO: Implement all of these functions.
             */
            floating_point_t theta_sum = atan((m[j][i] + m[i][j]) / (m[j][j] - m[i][i]));
            floating_point_t theta_diff = atan((m[j][i] - m[i][j]) / (m[j][j] + m[i][i]));
            floating_point_t theta_l = (theta_sum - theta_diff) / 2;
            floating_point_t theta_r = theta_sum - theta_l;
            floating_point_t sin_theta_l = sin(theta_l);
            floating_point_t cos_theta_l = cos(theta_l);
            floating_point_t sin_theta_r = sin(theta_r);
            floating_point_t cos_theta_r = cos(theta_r);

            /**
             * @brief Create temporary matricies for u_ij, u_ij_trans and v_ij_trans
             * 
             */
            fixed_point_u_t u_ij[4][4] = {
                {one_u, 0, 0, 0},
                {0, one_u, 0, 0},
                {0, 0, one_u, 0},
                {0, 0, 0, one_u},
            };

            fixed_point_u_t u_ij_trans[4][4] = {
                {one_u, 0, 0, 0},
                {0, one_u, 0, 0},
                {0, 0, one_u, 0},
                {0, 0, 0, one_u},
            };

            fixed_point_v_t v_ij_trans[4][4] = {
                {one_v, 0, 0, 0},
                {0, one_v, 0, 0},
                {0, 0, one_v, 0},
                {0, 0, 0, one_v},
            };

            /**
             * U_ij = [ cos(θl) -sin(θl) ]
             *        [ sin(θl)  cos(θl) ]
             */
            u_ij[i][i] = convert_to_fixed(cos_theta_l, SCALE_FACTOR_U);
            u_ij[j][j] = convert_to_fixed(cos_theta_l, SCALE_FACTOR_U);
            u_ij[i][j] = convert_to_fixed(-sin_theta_l, SCALE_FACTOR_U);
            u_ij[j][i] = convert_to_fixed(sin_theta_l, SCALE_FACTOR_U);

            /**
             * U_ij_Trans = [  cos(θl) sin(θl) ]
             *              [ -sin(θl) cos(θl) ]
             */
            u_ij_trans[i][i] = convert_to_fixed(cos_theta_l, SCALE_FACTOR_U);
            u_ij_trans[j][j] = convert_to_fixed(cos_theta_l, SCALE_FACTOR_U);
            u_ij_trans[i][j] = convert_to_fixed(sin_theta_l, SCALE_FACTOR_U);
            u_ij_trans[j][i] = convert_to_fixed(-sin_theta_l, SCALE_FACTOR_U);

            /**
             * V_ij_Trans = [  cos(θr) sin(θr) ]
             *              [ -sin(θr) cos(θr) ]
             */
            v_ij_trans[i][i] = convert_to_fixed(cos_theta_r, SCALE_FACTOR_V);
            v_ij_trans[j][j] = convert_to_fixed(cos_theta_r, SCALE_FACTOR_V);
            v_ij_trans[i][j] = convert_to_fixed(sin_theta_r, SCALE_FACTOR_V);
            v_ij_trans[j][i] = convert_to_fixed(-sin_theta_r, SCALE_FACTOR_V);

            /**
             * Create temporary matricies for claculations.
             */
            fixed_point_u_dp_t u_prime[4][4];
            fixed_point_v_dp_t v_trans_prime[4][4];
            fixed_point_m_tmp_dp_t m_prime_tmp[4][4];
            fixed_point_m_tmp_t m_prime_tmp_trunc[4][4];
            fixed_point_m_dp_t m_prime[4][4];

            // Do the calculations
            mat_mul_u_x_u(4, u_fixed, u_ij_trans, u_prime); // [U][U_ij_T] = [U']
            mat_mul_u_x_m(4, u_ij, m_fixed, m_prime_tmp);   // [U_ij][M] = [M'_tmp]
            for (int row = 0; row < 4; row++)
            {
                for (int col = 0; col < 4; col++)
                {
                    m_prime_tmp_trunc[row][col] = truncate_m_tmp(m_prime_tmp[row][col]);
                }
            }
            mat_mul_m_x_v(4, m_prime_tmp_trunc, v_ij_trans, m_prime);   // [M_tmp][V_ij_T] = [M']
            mat_mul_v_x_v(4, v_ij_trans, v_trans_fixed, v_trans_prime); // [V_ij][V_T] = [V'_T] <- I need to do this wrong to get it to work?????

            /**
             * Copy the values into U, V, and M and convert them back to floating point.
             * I think there we can avoid doing this for every ij pair, but for now this works.
             */
            for (int row = 0; row < 4; row++)
            {
                for (int col = 0; col < 4; col++)
                {
                    u[row][col] = convert_to_floating(u_prime[row][col], SCALE_FACTOR_U_DP);
                    v_trans[row][col] = convert_to_floating(v_trans_prime[row][col], SCALE_FACTOR_V_DP);
                    m[row][col] = convert_to_floating(m_prime[row][col], SCALE_FACTOR_M_DP);
                }
            }
        }
    }
}