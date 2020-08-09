/*
 *
 * svd.c
 *
 */
#include "svd.h"
#include "arm_neon.h"

/**
 * Create matricies for claculations.
 *
 * Two matricies are used to avoid copying. One matrix is used for input, and the
 * other is used for output. Every iteration, the input and output matricies switch.
 *
 * To avoid copying during the switch, just the pointers to the matricies are flipped.
 */
static volatile fixed_point_double_t u_prime_1[SIZE][SIZE], u_prime_2[SIZE][SIZE];
static volatile fixed_point_double_t v_trans_prime_1[SIZE][SIZE], v_trans_prime_2[SIZE][SIZE];
static volatile fixed_point_double_t m_prime_1[SIZE][SIZE], m_prime_2[SIZE][SIZE];

/**
 * Create matricies for internal calculations
 */
static volatile fixed_point_double_t u_ij[SIZE][SIZE];
static volatile fixed_point_double_t u_ij_trans[SIZE][SIZE];
static volatile fixed_point_double_t v_ij_trans[SIZE][SIZE];
static volatile fixed_point_double_t m_prime_tmp[SIZE][SIZE];

// Variable to track which matrix to use for input vs. output
static int input;

/**
 * This is the algorithm for SIMD matrix multiplication.
 * Because we are multiplying 2 fixed_point_t numbers, the result is a 
 * fixed_point_double_t. We to truncate the least significant bits to turn
 * the result into a fixed_point_t for the next round of computation.
 * 
 * Let X deonote the matrix on the LHS, and Y deonte the matrix on the RHS.
 * Let M_in denote the the nth element in the ith row of a matrix M.
 * Let M_i denote the vector [M_i0, M_i1, M_i2, M_i3] for a matrix M
 * let B denote the number of bits to shift
 *  
 * Each row of the square matrix is caluclated as such
 * 
 * X_i = (X_i0 * Y_0 + X_i1 * Y_1 + X_i2 * Y_2 + X_i3 * Y_3) >> B;
 * 
 * In this function X = U, Y = U_ij_transform, and the result is placed in
 * U` (where U and U` are either U_prime_1 or U_prime_2 depending on input)
 */
static void mat_mul_u_x_u_ij_trans_NEON()
{
    int32x4_t row_0, row_1, row_2, row_3, out_neon;

    // Read the rows of u_ij_trans
    row_0 = vld1q_s32((const fixed_point_double_t *)&u_ij_trans[0][0]); // Y_0
    row_1 = vld1q_s32((const fixed_point_double_t *)&u_ij_trans[1][0]); // Y_1
    row_2 = vld1q_s32((const fixed_point_double_t *)&u_ij_trans[2][0]); // Y_2
    row_3 = vld1q_s32((const fixed_point_double_t *)&u_ij_trans[3][0]); // Y_3

    for (int i = 0; i < SIZE; i++)
    {
        // Prevent copies by switching between matricies.
        if (input == 0)
        {
            out_neon = vmulq_n_s32(row_0, u_prime_1[i][0]);                      // X_i0 * Y_0
            out_neon = vaddq_s32(vmulq_n_s32(row_1, u_prime_1[i][1]), out_neon); // X_i1 * Y_1
            out_neon = vaddq_s32(vmulq_n_s32(row_2, u_prime_1[i][2]), out_neon); // X_i2 * Y_2
            out_neon = vaddq_s32(vmulq_n_s32(row_3, u_prime_1[i][3]), out_neon); // X_i3 * Y_3
            out_neon = vshrq_n_s32(out_neon, 13);                                // X_i >> B
            vst1q_s32((fixed_point_double_t *)&u_prime_2[i][0], out_neon);
        }
        else // input == 1
        {
            out_neon = vmulq_n_s32(row_0, u_prime_2[i][0]);                      // X_i0 * Y_0
            out_neon = vaddq_s32(vmulq_n_s32(row_1, u_prime_2[i][1]), out_neon); // X_i1 * Y_1
            out_neon = vaddq_s32(vmulq_n_s32(row_2, u_prime_2[i][2]), out_neon); // X_i2 * Y_2
            out_neon = vaddq_s32(vmulq_n_s32(row_3, u_prime_2[i][3]), out_neon); // X_i3 * Y_3
            out_neon = vshrq_n_s32(out_neon, 13);                                // X_i >> B
            vst1q_s32((fixed_point_double_t *)&u_prime_1[i][0], out_neon);
        }
    }
}

/**
 * This is the algorithm for SIMD matrix multiplication.
 * Because we are multiplying 2 fixed_point_t numbers, the result is a 
 * fixed_point_double_t. We to truncate the least significant bits to turn
 * the result into a fixed_point_t for the next round of computation.
 * 
 * Let X deonote the matrix on the LHS, and Y deonte the matrix on the RHS.
 * Let M_in denote the the nth element in the ith row of a matrix M.
 * Let M_i denote the vector [M_i0, M_i1, M_i2, M_i3] for a matrix M
 * let B denote the number of bits to shift
 *  
 * Each row of the square matrix is caluclated as such
 * 
 * X_i = (X_i0 * Y_0 + X_i1 * Y_1 + X_i2 * Y_2 + X_i3 * Y_3) >> B;
 * 
 * In this function X = U_ij, Y = M, and the result is placed in m_tmp (where
 * M is either M_prime_1 or M_prime_2 depending on input)
 */
static void mat_mul_u_ij_x_m_NEON()
{
    int32x4_t row_0, row_1, row_2, row_3, out_neon;

    if (input == 0)
    {
        row_0 = vld1q_s32((const fixed_point_double_t *)&m_prime_1[0][0]); // Y_0
        row_1 = vld1q_s32((const fixed_point_double_t *)&m_prime_1[1][0]); // Y_1
        row_2 = vld1q_s32((const fixed_point_double_t *)&m_prime_1[2][0]); // Y_2
        row_3 = vld1q_s32((const fixed_point_double_t *)&m_prime_1[3][0]); // Y_3
    }
    else // Input == 1
    {
        row_0 = vld1q_s32((const fixed_point_double_t *)&m_prime_2[0][0]); // Y_0
        row_1 = vld1q_s32((const fixed_point_double_t *)&m_prime_2[1][0]); // Y_1
        row_2 = vld1q_s32((const fixed_point_double_t *)&m_prime_2[2][0]); // Y_2
        row_3 = vld1q_s32((const fixed_point_double_t *)&m_prime_2[3][0]); // Y_3
    }

    for (int i = 0; i < SIZE; i++)
    {

        out_neon = vmulq_n_s32(row_0, u_ij[i][0]);                      // X_i0 * Y_0
        out_neon = vaddq_s32(vmulq_n_s32(row_1, u_ij[i][1]), out_neon); // X_i1 * Y_1
        out_neon = vaddq_s32(vmulq_n_s32(row_2, u_ij[i][2]), out_neon); // X_i2 * Y_2
        out_neon = vaddq_s32(vmulq_n_s32(row_3, u_ij[i][3]), out_neon); // X_i3 * Y_3
        out_neon = vshrq_n_s32(out_neon, 13);                           // X_i >> B
        vst1q_s32((fixed_point_double_t *)&m_prime_tmp[i][0], out_neon);
    }
}

/**
 * This is the algorithm for SIMD matrix multiplication.
 * Because we are multiplying 2 fixed_point_t numbers, the result is a 
 * fixed_point_double_t. We to truncate the least significant bits to turn
 * the result into a fixed_point_t for the next round of computation.
 * 
 * Let X deonote the matrix on the LHS, and Y deonte the matrix on the RHS.
 * Let M_in denote the the nth element in the ith row of a matrix M.
 * Let M_i denote the vector [M_i0, M_i1, M_i2, M_i3] for a matrix M
 * let B denote the number of bits to shift
 *  
 * Each row of the square matrix is caluclated as such
 * 
 * X_i = (X_i0 * Y_0 + X_i1 * Y_1 + X_i2 * Y_2 + X_i3 * Y_3) >> B;
 * 
 * In this function X = m_tmp, Y = v_ij_trans, and the result is placed in 
 * M` (where M` is either M_prime_1 or M_prime_2 depending on input)
 */
static void mat_mul_m_x_v_ij_trans_NEON()
{
    int32x4_t row_0, row_1, row_2, row_3, out_neon;

    row_0 = vld1q_s32((const fixed_point_double_t *)&v_ij_trans[0][0]); // Y_0
    row_1 = vld1q_s32((const fixed_point_double_t *)&v_ij_trans[1][0]); // Y_1
    row_2 = vld1q_s32((const fixed_point_double_t *)&v_ij_trans[2][0]); // Y_2
    row_3 = vld1q_s32((const fixed_point_double_t *)&v_ij_trans[3][0]); // Y_3

    for (int i = 0; i < SIZE; i++)
    {
        out_neon = vmulq_n_s32(row_0, m_prime_tmp[i][0]);                      // X_i0 * Y_0
        out_neon = vaddq_s32(vmulq_n_s32(row_1, m_prime_tmp[i][1]), out_neon); // X_i1 * Y_1
        out_neon = vaddq_s32(vmulq_n_s32(row_2, m_prime_tmp[i][2]), out_neon); // X_i2 * Y_2
        out_neon = vaddq_s32(vmulq_n_s32(row_3, m_prime_tmp[i][3]), out_neon); // X_i3 * Y_3
        out_neon = vshrq_n_s32(out_neon, 13);                                  // X_i >> B
        if (input == 0)
        {
            vst1q_s32((fixed_point_double_t *)&m_prime_2[i][0], out_neon);
        }
        else // input = 1;
        {
            vst1q_s32((fixed_point_double_t *)&m_prime_1[i][0], out_neon);
        }
    }
}

/**
 * This is the algorithm for SIMD matrix multiplication.
 * Because we are multiplying 2 fixed_point_t numbers, the result is a 
 * fixed_point_double_t. We to truncate the least significant bits to turn
 * the result into a fixed_point_t for the next round of computation.
 * 
 * Let X deonote the matrix on the LHS, and Y deonte the matrix on the RHS.
 * Let M_in denote the the nth element in the ith row of a matrix M.
 * Let M_i denote the vector [M_i0, M_i1, M_i2, M_i3] for a matrix M
 * let B denote the number of bits to shift
 *  
 * Each row of the square matrix is caluclated as such
 * 
 * X_i = (X_i0 * Y_0 + X_i1 * Y_1 + X_i2 * Y_2 + X_i3 * Y_3) >> B;
 * 
 * In this function X = v_ij_trans, Y = V_trans, and the result is placed in 
 * v_trans` (where v_trans and v_trans` are either v_trans_prime_1 or 
 * v_trans_prime_2 depending on input)
 */
static void mat_mul_v_ij_trans_x_v_trans_NEON()
{
    int32x4_t row_0, row_1, row_2, row_3, out_neon;

    if (input == 0)
    {
        row_0 = vld1q_s32((const fixed_point_double_t *)&v_trans_prime_1[0][0]); // Y_0
        row_1 = vld1q_s32((const fixed_point_double_t *)&v_trans_prime_1[1][0]); // Y_1
        row_2 = vld1q_s32((const fixed_point_double_t *)&v_trans_prime_1[2][0]); // Y_2
        row_3 = vld1q_s32((const fixed_point_double_t *)&v_trans_prime_1[3][0]); // Y_3
    }
    else
    {
        row_0 = vld1q_s32((const fixed_point_double_t *)&v_trans_prime_2[0][0]); // Y_0
        row_1 = vld1q_s32((const fixed_point_double_t *)&v_trans_prime_2[1][0]); // Y_1
        row_2 = vld1q_s32((const fixed_point_double_t *)&v_trans_prime_2[2][0]); // Y_2
        row_3 = vld1q_s32((const fixed_point_double_t *)&v_trans_prime_2[3][0]); // Y_3
    }

    for (int i = 0; i < SIZE; i++)
    {
        out_neon = vmulq_n_s32(row_0, v_ij_trans[i][0]);                      // X_i0 * Y_0
        out_neon = vaddq_s32(vmulq_n_s32(row_1, v_ij_trans[i][1]), out_neon); // X_i1 * Y_1
        out_neon = vaddq_s32(vmulq_n_s32(row_2, v_ij_trans[i][2]), out_neon); // X_i2 * Y_2
        out_neon = vaddq_s32(vmulq_n_s32(row_3, v_ij_trans[i][3]), out_neon); // X_i3 * Y_3
        out_neon = vshrq_n_s32(out_neon, 13);                                 // X_i >> B
        if (input == 0)
        {
            vst1q_s32((fixed_point_double_t *)&v_trans_prime_2[i][0], out_neon);
        }
        else // input = 1;
        {
            vst1q_s32((fixed_point_double_t *)&v_trans_prime_1[i][0], out_neon);
        }
    }
}

/**
 * This is a function to set the values of u_ij, u_ij_trans and v_ij_trans to 
 * zero (or one), because memset doesn't work with volatile pointers (as far 
 * as I know).
 */
static void reset_mats()
{
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            u_ij[i][j] = i == j ? one_u : 0;
            u_ij_trans[i][j] = i == j ? one_u : 0;
            v_ij_trans[i][j] = i == j ? one_v : 0;
        }
    }
}

/**
 * @brief Performes a single sweep of the svd algorithm
 *
 */
void sweep(floating_point_t m[SIZE][SIZE], floating_point_t u[SIZE][SIZE], floating_point_t v_trans[SIZE][SIZE])
{
    input = 0;

    // Convert the input matricies to fixed point.
    for (int row = 0; row < SIZE; row++)
    {
        for (int col = 0; col < SIZE; col++)
        {
            m_prime_1[row][col] = convert_to_fixed(m[row][col], SCALE_FACTOR_M);
            u_prime_1[row][col] = convert_to_fixed(u[row][col], SCALE_FACTOR_U);
            v_trans_prime_1[row][col] = convert_to_fixed(v_trans[row][col], SCALE_FACTOR_V);
        }
    }

    /**
     * Start of iterations
     */
    for (int i = 0; i < SIZE - 1; i++)
    {
        for (int j = i + 1; j < SIZE; j++)
        {
            fixed_point_t m_ij, m_ji, m_ii, m_jj;

            // Read m_ij, m_ji, m_ii, m_jj from M.
            if (input == 0)
            {
                m_ij = m_prime_1[i][j];
                m_ji = m_prime_1[j][i];
                m_ii = m_prime_1[i][i];
                m_jj = m_prime_1[j][j];
            }
            else
            {
                m_ij = m_prime_2[i][j];
                m_ji = m_prime_2[j][i];
                m_ii = m_prime_2[i][i];
                m_jj = m_prime_2[j][j];
            }

            /**
             * Do all of the angle calculations
             * 
             * θ_sum  = θr + θl = atan((m_ji + m_ij) / (m_jj - m_ii))
             * θ_diff = θr - θl= atan((m_ji - m_ij) / (m_jj + m_ii))
             * 
             * θr = θ_diff + θl
             * θl = θ_sum - θr = θ_sum - (θ_diff + θ_l) = (θ_sum - θ_diff) /2
             * 
             */
            fixed_point_double_t x = fixed_point_div(m_ji + m_ij, m_jj - m_ii);
            fixed_point_t theta_sum_fixed = arctan_lookup(x);

            x = fixed_point_div(m_ji - m_ij, m_jj + m_ii);
            fixed_point_t theta_diff_fixed = arctan_lookup(x);

            fixed_point_double_t theta_l_fixed, theta_r_fixed;
            theta_l_fixed = ((fixed_point_double_t)theta_sum_fixed - (fixed_point_double_t)theta_diff_fixed) >> 1;
            theta_r_fixed = theta_sum_fixed - theta_l_fixed;

            fixed_point_double_t sin_theta_l_fixed = sin_lookup(theta_l_fixed);
            fixed_point_double_t cos_theta_l_fixed = cos_lookup(theta_l_fixed);
            fixed_point_double_t sin_theta_r_fixed = sin_lookup(theta_r_fixed);
            fixed_point_double_t cos_theta_r_fixed = cos_lookup(theta_r_fixed);

            // Reset the matricies to unit matricies.
            reset_mats();

            /**
             * U_ij = [ cos(θl) -sin(θl) ]
             *        [ sin(θl)  cos(θl) ]
             */
            u_ij[i][i] = cos_theta_l_fixed;
            u_ij[j][j] = cos_theta_l_fixed;
            u_ij[i][j] = -sin_theta_l_fixed;
            u_ij[j][i] = sin_theta_l_fixed;

            /**
             * U_ij_Trans = [  cos(θl) sin(θl) ]
             *              [ -sin(θl) cos(θl) ]
             */
            u_ij_trans[i][i] = cos_theta_l_fixed;
            u_ij_trans[j][j] = cos_theta_l_fixed;
            u_ij_trans[i][j] = sin_theta_l_fixed;
            u_ij_trans[j][i] = -sin_theta_l_fixed;

            /**
             * V_ij_Trans = [  cos(θr) sin(θr) ]
             *              [ -sin(θr) cos(θr) ]
             */
            v_ij_trans[i][i] = cos_theta_r_fixed;
            v_ij_trans[j][j] = cos_theta_r_fixed;
            v_ij_trans[i][j] = sin_theta_r_fixed;
            v_ij_trans[j][i] = -sin_theta_r_fixed;

            // Do the calculations
            mat_mul_u_x_u_ij_trans_NEON();       // [U][U_ij_T] = [U']
            mat_mul_u_ij_x_m_NEON();             // [U_ij][M] = [M'_tmp]
            mat_mul_m_x_v_ij_trans_NEON();       // [M_tmp][V_ij_T] = [M']
            mat_mul_v_ij_trans_x_v_trans_NEON(); // [V_ij][V_T] = [V'_T] <- I need to do this wrong to get it to work?????

            // swap input and output matricies.
            input = input == 0 ? 1 : 0;
        }
    }

    // Convert the fixed point matricies back into floating point.
    for (int row = 0; row < SIZE; row++)
    {
        for (int col = 0; col < SIZE; col++)
        {
            if (input == 0)
            {
                u[row][col] = convert_to_floating(u_prime_1[row][col], SCALE_FACTOR_U);
                v_trans[row][col] = convert_to_floating(v_trans_prime_1[row][col], SCALE_FACTOR_V);
                m[row][col] = convert_to_floating(m_prime_1[row][col], SCALE_FACTOR_M);
            }
            else
            {
                u[row][col] = convert_to_floating(u_prime_2[row][col], SCALE_FACTOR_U);
                v_trans[row][col] = convert_to_floating(v_trans_prime_2[row][col], SCALE_FACTOR_V);
                m[row][col] = convert_to_floating(m_prime_2[row][col], SCALE_FACTOR_M);
            }
        }
    }
}