/*
 *
 * svd.c
 *
 */
#include "svd.h"
#include "svd_math.h"
#include "config.h"

fixed_point_t *access(fixed_point_t *arr, size_t size, size_t row, size_t col)
{
    return arr + size * row + col;
}

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
void mat_mul(int size, fixed_point_t *LHS, fixed_point_t *RHS, fixed_point_t *out)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            *access(out, size, i, j) = 0;
            for (int k = 0; k < size; k++)
            {
                *access(out, size, i, j) += truncate(
                    fixed_point_mul(
                        *access(LHS, size, i, k),
                        *access(RHS, size, k, j)));
            }
        }
    }
}

void mat_mul_u_x_u(
    int size,
    fixed_point_u_t *LHS,
    fixed_point_u_t *RHS,
    fixed_point_u_t *out) __attribute__((alias("mat_mul")));

void mat_mul_u_x_m(
    int size,
    fixed_point_u_t *LHS,
    fixed_point_m_t *RHS,
    fixed_point_m_tmp_t *out) __attribute__((alias("mat_mul")));

void mat_mul_m_x_v(
    int size,
    fixed_point_m_tmp_t *LHS,
    fixed_point_v_t *RHS,
    fixed_point_m_t *out) __attribute__((alias("mat_mul")));

void mat_mul_v_x_v(
    int size,
    fixed_point_v_t *LHS,
    fixed_point_v_t *RHS,
    fixed_point_v_t *out) __attribute__((alias("mat_mul")));

/**
 * @brief Performes a single sweep of the svd algorithm
 * 
 */
void sweep(floating_point_t m[4][4], floating_point_t u[4][4], floating_point_t v_trans[4][4])
{
    /**
     * Create temporary matricies for claculations.
     */
    fixed_point_u_t u_prime_1[4][4], u_prime_2[4][4];
    fixed_point_v_t v_trans_prime_1[4][4], v_trans_prime_2[4][4];
    fixed_point_m_t m_prime_1[4][4], m_prime_2[4][4];

    // Convert the input matricies to fixed point.
    for (int row = 0; row < 4; row++)
    {
        for (int col = 0; col < 4; col++)
        {
            m_prime_1[row][col] = convert_to_fixed(m[row][col], SCALE_FACTOR_M);
            u_prime_1[row][col] = convert_to_fixed(u[row][col], SCALE_FACTOR_U);
            v_trans_prime_1[row][col] = convert_to_fixed(v_trans[row][col], SCALE_FACTOR_V);
        }
    }

    int count = 0;

    for (int i = 0; i < 3; i++)
    {
        for (int j = i + 1; j < 4; j++)
        {

            /**
             * Do all of the angle calculations
             * 
             * TODO: Implement all of these functions.
             */
            floating_point_t theta_sum, theta_diff;
            if (count % 2 == 0)
            {
                theta_sum = atan((m_prime_1[j][i] + m_prime_1[i][j]) / (floating_point_t)(m_prime_1[j][j] - m_prime_1[i][i]));
                theta_diff = atan((m_prime_1[j][i] - m_prime_1[i][j]) / (floating_point_t)(m_prime_1[j][j] + m_prime_1[i][i]));
            }
            else
            {
                theta_sum = atan((m_prime_2[j][i] + m_prime_2[i][j]) / (floating_point_t)(m_prime_2[j][j] - m_prime_2[i][i]));
                theta_diff = atan((m_prime_2[j][i] - m_prime_2[i][j]) / (floating_point_t)(m_prime_2[j][j] + m_prime_2[i][i]));
            }

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

            fixed_point_m_tmp_t m_prime_tmp[4][4];

            if (count++ % 2 == 0)
            {                                                                   // Do the calculations
                mat_mul_u_x_u(4, u_prime_1, u_ij_trans, u_prime_2);             // [U][U_ij_T] = [U']
                mat_mul_u_x_m(4, u_ij, m_prime_1, m_prime_tmp);                 // [U_ij][M] = [M'_tmp]
                mat_mul_m_x_v(4, m_prime_tmp, v_ij_trans, m_prime_2);           // [M_tmp][V_ij_T] = [M']
                mat_mul_v_x_v(4, v_ij_trans, v_trans_prime_1, v_trans_prime_2); // [V_ij][V_T] = [V'_T] <- I need to do this wrong to get it to work?????
            }
            else
            {
                mat_mul_u_x_u(4, u_prime_2, u_ij_trans, u_prime_1);             // [U][U_ij_T] = [U']
                mat_mul_u_x_m(4, u_ij, m_prime_2, m_prime_tmp);                 // [U_ij][M] = [M'_tmp]
                mat_mul_m_x_v(4, m_prime_tmp, v_ij_trans, m_prime_1);           // [M_tmp][V_ij_T] = [M']
                mat_mul_v_x_v(4, v_ij_trans, v_trans_prime_2, v_trans_prime_1); // [V_ij][V_T] = [V'_T] <- I need to do this wrong to get it to work?????
            }
        }
    }

    for (int row = 0; row < 4; row++)
    {
        for (int col = 0; col < 4; col++)
        {
            u[row][col] = convert_to_floating(u_prime_1[row][col], SCALE_FACTOR_U);
            v_trans[row][col] = convert_to_floating(v_trans_prime_1[row][col], SCALE_FACTOR_V);
            m[row][col] = convert_to_floating(m_prime_1[row][col], SCALE_FACTOR_M);
        }
    }
}