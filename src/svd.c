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
 * Create a table of pointers to the matricies for calculations.
 */
static volatile fixed_point_double_t *u_mats[] = {&u_prime_1[0][0], &u_prime_2[0][0]};
static volatile fixed_point_double_t *v_mats[] = {&v_trans_prime_1[0][0], &v_trans_prime_2[0][0]};
static volatile fixed_point_double_t *m_mats[] = {&m_prime_1[0][0], &m_prime_2[0][0]};

/**
 * Create matricies for internal calculations
 */
static volatile fixed_point_double_t u_ij[SIZE][SIZE];
static volatile fixed_point_double_t u_ij_trans[SIZE][SIZE];
static volatile fixed_point_double_t v_ij_trans[SIZE][SIZE];
static volatile fixed_point_double_t m_prime_tmp[SIZE][SIZE];

// Variables to track which matrix to use for input vs. output
static int input = 0, output = 1;

void print_matrix(const fixed_point_double_t *m)
{
    for (int i = 0; i < 4; i++)
    {
        printf("[");
        for (int j = 0; j < 4 - 1; j++)
        {
            printf(" %d", *(m + i * 4 + j));
        }
        printf(" %d ]\n", *(m + i * 4 + 4 - 1));
    }
    printf("\n");
}

// U * U_ij_trans -> u_prime
static void mat_mul_u_x_u_ij_trans_NEON()
{
    int32x4_t row_0, row_1, row_2, row_3, out_neon;

    print_matrix((const fixed_point_double_t *)&u_ij_trans[0][0]);
    if (input == 0)
        print_matrix((const fixed_point_double_t *)&u_prime_1[0][0]);
    else
        print_matrix((const fixed_point_double_t *)&u_prime_2[0][0]);

    row_0 = vld1q_s32((const fixed_point_double_t *)&u_ij_trans[0][0]);
    row_1 = vld1q_s32((const fixed_point_double_t *)&u_ij_trans[1][0]);
    row_2 = vld1q_s32((const fixed_point_double_t *)&u_ij_trans[2][0]);
    row_3 = vld1q_s32((const fixed_point_double_t *)&u_ij_trans[3][0]);

    for (int i = 0; i < SIZE; i++)
    {
        if (input == 0)
        {
            out_neon = vmulq_n_s32(row_0, u_prime_1[i][0]);
            out_neon = vaddq_s32(vmulq_n_s32(row_1, u_prime_1[i][1]), out_neon);
            out_neon = vaddq_s32(vmulq_n_s32(row_2, u_prime_1[i][2]), out_neon);
            out_neon = vaddq_s32(vmulq_n_s32(row_3, u_prime_1[i][3]), out_neon);
            vst1q_s32((fixed_point_double_t *)&u_prime_2[i][0], out_neon);
        }
        else // input == 1
        {
            out_neon = vmulq_n_s32(row_0, u_prime_2[i][0]);
            out_neon = vaddq_s32(vmulq_n_s32(row_1, u_prime_2[i][1]), out_neon);
            out_neon = vaddq_s32(vmulq_n_s32(row_2, u_prime_2[i][2]), out_neon);
            out_neon = vaddq_s32(vmulq_n_s32(row_3, u_prime_2[i][3]), out_neon);
            vst1q_s32((fixed_point_double_t *)&u_prime_1[i][0], out_neon);
        }
    }

    if (input == 0)
        print_matrix((const fixed_point_double_t *)&u_prime_2[0][0]);
    else
        print_matrix((const fixed_point_double_t *)&u_prime_1[0][0]);
}

// U_ij x M -> m_prime_tmp
static void mat_mul_u_ij_x_m_NEON()
{
    int32x4_t row_0, row_1, row_2, row_3, out_neon;

    // print_matrix((const fixed_point_double_t *)&u_ij_trans[0][0]);
    // if (input == 0)
    //     print_matrix((const fixed_point_double_t *)&u_prime_1[0][0]);
    // else
    //     print_matrix((const fixed_point_double_t *)&u_prime_2[0][0]);
    if (input == 0)
    {
        row_0 = vld1q_s32((const fixed_point_double_t *)&m_prime_1[0][0]);
        row_1 = vld1q_s32((const fixed_point_double_t *)&m_prime_1[1][0]);
        row_2 = vld1q_s32((const fixed_point_double_t *)&m_prime_1[2][0]);
        row_3 = vld1q_s32((const fixed_point_double_t *)&m_prime_1[3][0]);
    }
    else // Input == 1
    {
        row_0 = vld1q_s32((const fixed_point_double_t *)&m_prime_2[0][0]);
        row_1 = vld1q_s32((const fixed_point_double_t *)&m_prime_2[1][0]);
        row_2 = vld1q_s32((const fixed_point_double_t *)&m_prime_2[2][0]);
        row_3 = vld1q_s32((const fixed_point_double_t *)&m_prime_2[3][0]);
    }

    for (int i = 0; i < SIZE; i++)
    {

        out_neon = vmulq_n_s32(row_0, u_ij[i][0]);
        out_neon = vaddq_s32(vmulq_n_s32(row_1, u_ij[i][1]), out_neon);
        out_neon = vaddq_s32(vmulq_n_s32(row_2, u_ij[i][2]), out_neon);
        out_neon = vaddq_s32(vmulq_n_s32(row_3, u_ij[i][3]), out_neon);
        vst1q_s32((fixed_point_double_t *)&m_prime_tmp[i][0], out_neon);
    }

    // if (input == 0)
    //     print_matrix((const fixed_point_double_t *)&u_prime_2[0][0]);
    // else
    //     print_matrix((const fixed_point_double_t *)&u_prime_1[0][0]);
}

// m_prime_tmp * v_ij_trans -> m_prime
static void mat_mul_m_x_v_ij_trans_NEON()
{
    int32x4_t row_0, row_1, row_2, row_3, out_neon;

    // print_matrix((const fixed_point_double_t *)&u_ij_trans[0][0]);
    // if (input == 0)
    //     print_matrix((const fixed_point_double_t *)&u_prime_1[0][0]);
    // else
    //     print_matrix((const fixed_point_double_t *)&u_prime_2[0][0]);
    row_0 = vld1q_s32((const fixed_point_double_t *)&v_ij_trans[0][0]);
    row_1 = vld1q_s32((const fixed_point_double_t *)&v_ij_trans[1][0]);
    row_2 = vld1q_s32((const fixed_point_double_t *)&v_ij_trans[2][0]);
    row_3 = vld1q_s32((const fixed_point_double_t *)&v_ij_trans[3][0]);

    for (int i = 0; i < SIZE; i++)
    {
        out_neon = vmulq_n_s32(row_0, m_prime_tmp[i][0]);
        out_neon = vaddq_s32(vmulq_n_s32(row_1, m_prime_tmp[i][1]), out_neon);
        out_neon = vaddq_s32(vmulq_n_s32(row_2, m_prime_tmp[i][2]), out_neon);
        out_neon = vaddq_s32(vmulq_n_s32(row_3, m_prime_tmp[i][3]), out_neon);
        if (input == 0)
        {
            vst1q_s32((fixed_point_double_t *)&m_prime_1[i][0], out_neon);
        }
        else // input = 1;
        {
            vst1q_s32((fixed_point_double_t *)&m_prime_2[i][0], out_neon);
        }
    }

    // if (input == 0)
    //     print_matrix((const fixed_point_double_t *)&u_prime_2[0][0]);
    // else
    //     print_matrix((const fixed_point_double_t *)&u_prime_1[0][0]);
}

// v_ij_trans * v_trans -> v_prime
static void mat_mul_v_ij_trans_x_v_trans_NEON()
{
    int32x4_t row_0, row_1, row_2, row_3, out_neon;

    // print_matrix((const fixed_point_double_t *)&u_ij_trans[0][0]);
    // if (input == 0)
    //     print_matrix((const fixed_point_double_t *)&u_prime_1[0][0]);
    // else
    //     print_matrix((const fixed_point_double_t *)&u_prime_2[0][0]);

    if (input == 0)
    {
        row_0 = vld1q_s32((const fixed_point_double_t *)&v_trans_prime_1[0][0]);
        row_1 = vld1q_s32((const fixed_point_double_t *)&v_trans_prime_1[1][0]);
        row_2 = vld1q_s32((const fixed_point_double_t *)&v_trans_prime_1[2][0]);
        row_3 = vld1q_s32((const fixed_point_double_t *)&v_trans_prime_1[3][0]);
    }
    else
    {
        row_0 = vld1q_s32((const fixed_point_double_t *)&v_trans_prime_2[0][0]);
        row_1 = vld1q_s32((const fixed_point_double_t *)&v_trans_prime_2[1][0]);
        row_2 = vld1q_s32((const fixed_point_double_t *)&v_trans_prime_2[2][0]);
        row_3 = vld1q_s32((const fixed_point_double_t *)&v_trans_prime_2[3][0]);
    }

    for (int i = 0; i < SIZE; i++)
    {
        out_neon = vmulq_n_s32(row_0, v_ij_trans[i][0]);
        out_neon = vaddq_s32(vmulq_n_s32(row_1, v_ij_trans[i][1]), out_neon);
        out_neon = vaddq_s32(vmulq_n_s32(row_2, v_ij_trans[i][2]), out_neon);
        out_neon = vaddq_s32(vmulq_n_s32(row_3, v_ij_trans[i][3]), out_neon);
        if (input == 0)
        {
            vst1q_s32((fixed_point_double_t *)&v_trans_prime_2[i][0], out_neon);
        }
        else // input = 1;
        {
            vst1q_s32((fixed_point_double_t *)&v_trans_prime_1[i][0], out_neon);
        }
    }

    // if (input == 0)
    //     print_matrix((const fixed_point_double_t *)&u_prime_2[0][0]);
    // else
    //     print_matrix((const fixed_point_double_t *)&u_prime_1[0][0]);
}

static void zero_mats()
{
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            u_ij[i][j] = 0;
            u_ij_trans[i][j] = 0;
            v_ij_trans[i][j] = 0;
        }
    }
}

static inline volatile fixed_point_double_t *access(volatile fixed_point_double_t *arr, size_t row, size_t col)
{
    return arr + SIZE * row + col;
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
void mat_mul(volatile fixed_point_double_t *LHS, volatile fixed_point_double_t *RHS, volatile fixed_point_double_t *out)
{

    // printf("NEON result: \n");
    // print_matrix(out);

    for (int i = 0; i < SIZE; i++)
    {
        for (int j = 0; j < SIZE; j++)
        {
            *access(out, i, j) = 0;
            for (int k = 0; k < SIZE; k++)
            {
                *access(out, i, j) += truncate(
                    fixed_point_mul(
                        *access(LHS, i, k),
                        *access(RHS, k, j)));
            }
        }
    }
    // printf("OUR result: \n");
    // print_matrix(out);
}

/**
 * @brief Performes a single sweep of the svd algorithm
 *
 */
void sweep(floating_point_t m[SIZE][SIZE], floating_point_t u[SIZE][SIZE], floating_point_t v_trans[SIZE][SIZE])
{

    // Convert the input matricies to fixed point.
    for (int row = 0; row < SIZE; row++)
    {
        for (int col = 0; col < SIZE; col++)
        {
            *access(m_mats[input], row, col) = convert_to_fixed(m[row][col], SCALE_FACTOR_M);
            *access(u_mats[input], row, col) = convert_to_fixed(u[row][col], SCALE_FACTOR_U);
            *access(v_mats[input], row, col) = convert_to_fixed(v_trans[row][col], SCALE_FACTOR_V);
        }
    }

    /**
     * Start of iterations
     */
    for (int i = 0; i < SIZE - 1; i++)
    {
        for (int j = i + 1; j < SIZE; j++)
        {
            fixed_point_t m_ij = *access(m_mats[input], i, j),
                          m_ji = *access(m_mats[input], j, i),
                          m_ii = *access(m_mats[input], i, i),
                          m_jj = *access(m_mats[input], j, j);

            /**
             * Do all of the angle calculations
             *
             * TODO: Implement all of these functions.
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

            // Reset the matricies to zero.
            zero_mats();

            for (int k = 0; k < SIZE; k++)
            {
                u_ij[k][k] = one_u;
                u_ij_trans[k][k] = one_u;
                v_ij_trans[k][k] = one_v;
            }

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
            int tmp = input;
            input = output;
            output = tmp;
        }
    }

    // Convert the fixed point matricies back into floating point.
    for (int row = 0; row < SIZE; row++)
    {
        for (int col = 0; col < SIZE; col++)
        {
            u[row][col] = convert_to_floating(*access(u_mats[input], row, col), SCALE_FACTOR_U);
            v_trans[row][col] = convert_to_floating(*access(v_mats[input], row, col), SCALE_FACTOR_V);
            m[row][col] = convert_to_floating(*access(m_mats[input], row, col), SCALE_FACTOR_M);
        }
    }
}