/*
 *
 * svd.c
 *
 */
#include "svd.h"

/**
 * @brief Helper function to multiply a (size x size) matrix
 * 
 * Result is placed in out[][]
 * 
 * @param size 
 * @param LHS 
 * @param RHS 
 * @param out 
 */
void mat_mul(int size, double LHS[size][size], double RHS[size][size], double out[size][size])
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            out[i][j] = 0;
            for (int k = 0; k < size; k++)
            {
                out[i][j] += LHS[i][k] * RHS[k][j];
            }
        }
    }
}

/**
 * @brief Function to obtain arctan using integer / fixed-point arithmetic
 * and taylor series expansion
 * 
 * Result is placed in out[][]
 * 
 * @param approximation
 * @param out
 * @return out
 */
double arctan(int approximation, double out)
{
    for (int i = 0; i < approximation; i++)
    {
        out += pow(-1, (i + 1)) * (1 / (i * 2 + 1)) * pow(out, (i * 2 + 1));
    }
    return out;
}

/**
 * @brief Performes a single sweep of the svd algorithm
 * 
 */
void sweep(double m[4][4], double u[4][4], double v_trans[4][4])
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = i + 1; j < 4; j++)
        {
            /**
             * Do all of the angle calculations
             * 
             * TODO: Implement all of these functions.
             * 
             */
            double theta_sum = arctan(1, (m[j][i] + m[i][j]) / (m[j][j] - m[i][i]));
            double theta_diff = arctan(1, (m[j][i] - m[i][j]) / (m[j][j] + m[i][i]));
            double theta_l = (theta_sum - theta_diff) / 2;
            double theta_r = theta_sum - theta_l;
            double sin_theta_l = sin(theta_l);
            double cos_theta_l = cos(theta_l);
            double sin_theta_r = sin(theta_r);
            double cos_theta_r = cos(theta_r);

            /**
             * @brief Create temporary matricies for u_ij, u_ij_trans and v_ij_trans
             * 
             */
            double u_ij[4][4] = {
                {1, 0, 0, 0},
                {0, 1, 0, 0},
                {0, 0, 1, 0},
                {0, 0, 0, 1},
            };

            double u_ij_trans[4][4] = {
                {1, 0, 0, 0},
                {0, 1, 0, 0},
                {0, 0, 1, 0},
                {0, 0, 0, 1},
            };

            double v_ij_trans[4][4] = {
                {1, 0, 0, 0},
                {0, 1, 0, 0},
                {0, 0, 1, 0},
                {0, 0, 0, 1},
            };

            /**
             * U_ij = [ cos(θl) -sin(θl) ]
             *        [ sin(θl)  cos(θl) ]
             */
            u_ij[i][i] = cos_theta_l;
            u_ij[j][j] = cos_theta_l;
            u_ij[i][j] = -sin_theta_l;
            u_ij[j][i] = sin_theta_l;

            /**
             * U_ij_Trans = [  cos(θl) sin(θl) ]
             *              [ -sin(θl) cos(θl) ]
             */
            u_ij_trans[i][i] = cos_theta_l;
            u_ij_trans[j][j] = cos_theta_l;
            u_ij_trans[i][j] = sin_theta_l;
            u_ij_trans[j][i] = -sin_theta_l;

            /**
             * V_ij_Trans = [  cos(θr) sin(θr) ]
             *              [ -sin(θr) cos(θr) ]
             */
            v_ij_trans[i][i] = cos_theta_r;
            v_ij_trans[j][j] = cos_theta_r;
            v_ij_trans[i][j] = sin_theta_r;
            v_ij_trans[j][i] = -sin_theta_r;

            /**
             * Create temporary matricies for claculations.
             */
            double u_prime[4][4];
            double v_trans_prime[4][4];
            double m_prime_tmp[4][4];
            double m_prime[4][4];

            // Do the calculations
            mat_mul(4, u, u_ij_trans, u_prime);             // [U][U_ij_T] = [U']
            mat_mul(4, u_ij, m, m_prime_tmp);               // [U_ij][M] = [M'_tmp]
            mat_mul(4, m_prime_tmp, v_ij_trans, m_prime);   // [M_tmp][V_ij_T] = [M']
            mat_mul(4, v_ij_trans, v_trans, v_trans_prime); // [V_ij][V_T] = [V'_T] <- I need to do this wrong to get it to work?????

            /**
             * Copy the values into U, V, and M.
             * I think there we can avoid doing this for every ij pair, but for now this works.
             */
            for (int row = 0; row < 4; row++)
            {
                for (int col = 0; col < 4; col++)
                {
                    u[row][col] = u_prime[row][col];
                    v_trans[row][col] = v_trans_prime[row][col];
                    m[row][col] = m_prime[row][col];
                }
            }
        }
    }
}