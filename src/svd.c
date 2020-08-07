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
     * 
     * Two matricies are used to avoid copying. One matrix is used for input, and the 
     * other is used for output. Every iteration, the input and output matricies switch.
     * 
     * To avoid copying during the switch, just the pointers to the matricies are flipped.
     */
    fixed_point_u_t u_prime_1[4][4], u_prime_2[4][4];
    fixed_point_v_t v_trans_prime_1[4][4], v_trans_prime_2[4][4];
    fixed_point_m_t m_prime_1[4][4], m_prime_2[4][4];

    /**
     * Create a table of pointers to the matricies for calculations.
     */
    fixed_point_u_t *u_mats[] = {&u_prime_1[0][0], &u_prime_2[0][0]};
    fixed_point_u_t *v_mats[] = {&v_trans_prime_1[0][0], &v_trans_prime_2[0][0]};
    fixed_point_m_t *m_mats[] = {&m_prime_1[0][0], &m_prime_2[0][0]};

    // Variables to track which matrix to use for input vs. output
    int input = 0, output = 1;

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

    /**
     * Start of iterations
     */
    for (int i = 0; i < 3; i++)
    {
        for (int j = i + 1; j < 4; j++)
        {

            /**
             * Do all of the angle calculations
             * 
             * TODO: Implement all of these functions.
             */
            floating_point_t theta_sum = atan((*access(m_mats[input], 4, j, i) + *access(m_mats[input], 4, i, j)) / (floating_point_t)(*access(m_mats[input], 4, j, j) - *access(m_mats[input], 4, i, i)));
            floating_point_t theta_diff = atan((*access(m_mats[input], 4, j, i) - *access(m_mats[input], 4, i, j)) / (floating_point_t)(*access(m_mats[input], 4, j, j) + *access(m_mats[input], 4, i, i)));
            floating_point_t theta_l = (theta_sum - theta_diff) / 2;
            floating_point_t theta_r = theta_sum - theta_l;
            floating_point_t sin_theta_l = sin(theta_l);
            floating_point_t cos_theta_l = cos(theta_l);
            floating_point_t sin_theta_r = sin(theta_r);
            floating_point_t cos_theta_r = cos(theta_r);

            /**
             * Create temporary matricies for u_ij, u_ij_trans and v_ij_trans
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

            // Do the calculations
            mat_mul_u_x_u(4, u_mats[input], &u_ij_trans[0][0], u_mats[output]);      // [U][U_ij_T] = [U']
            mat_mul_u_x_m(4, &u_ij[0][0], m_mats[input], &m_prime_tmp[0][0]);        // [U_ij][M] = [M'_tmp]
            mat_mul_m_x_v(4, &m_prime_tmp[0][0], &v_ij_trans[0][0], m_mats[output]); // [M_tmp][V_ij_T] = [M']
            mat_mul_v_x_v(4, &v_ij_trans[0][0], v_mats[input], v_mats[output]);      // [V_ij][V_T] = [V'_T] <- I need to do this wrong to get it to work?????

            // swap input and output matricies.
            int tmp = input;
            input = output;
            output = tmp;
        }
    }

    // Convert the fixed point matricies back into floating point.
    for (int row = 0; row < 4; row++)
    {
        for (int col = 0; col < 4; col++)
        {
            u[row][col] = convert_to_floating(*access(u_mats[input], 4, row, col), SCALE_FACTOR_U);
            v_trans[row][col] = convert_to_floating(*access(v_mats[input], 4, row, col), SCALE_FACTOR_V);
            m[row][col] = convert_to_floating(*access(m_mats[input], 4, row, col), SCALE_FACTOR_M);
        }
    }
}