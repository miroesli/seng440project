/*
 *
 * main.c
 *
 */

#include "config.h"
#include "svd.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double m[4][4] = {
    {31, 77, -11, 26},
    {-42, 14, 79, -53},
    {-68, -10, 45, 90},
    {34, 16, 38, -19},
};

double u[4][4] = {
    {1, 0, 0, 0},
    {0, 1, 0, 0},
    {0, 0, 1, 0},
    {0, 0, 0, 1},
};

double v_trans[4][4] = {
    {1, 0, 0, 0},
    {0, 1, 0, 0},
    {0, 0, 1, 0},
    {0, 0, 0, 1},
};

void mat_print(double mat[4][4])
{
    for (int row = 0; row < 4; row++)
    {
        printf("[ ");
        for (int col = 0; col < 3; col++)
        {
            printf("%f, ", mat[row][col]);
        }
        printf("%f ]\n", mat[row][3]);
    }
}

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

void sweep()
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = i + 1; j < 4; j++)
        {

            printf("************************%d & %d************************\n", i, j);
            double m_tmp[2][2] = {
                {m[i][i], m[i][j]},
                {m[j][i], m[j][j]},
            };

            printf("m[%d][%d]: %f\n", i, i, m[i][i]);
            printf("m[%d][%d]: %f\n", i, j, m[i][j]);
            printf("m[%d][%d]: %f\n", j, i, m[j][i]);
            printf("m[%d][%d]: %f\n", j, j, m[j][j]);

            double theta_sum = atan((m[j][i] + m[i][j]) / (m[j][j] - m[i][i]));
            double theta_diff = atan((m[j][i] - m[i][j]) / (m[j][j] + m[i][i]));

            printf("Theta SUM: %f\n", theta_sum);
            printf("Theta DIFF: %f\n", theta_diff);

            double theta_l = (theta_sum - theta_diff) / 2;
            double theta_r = theta_sum - theta_l;

            printf("Theta L: %f\n", theta_l);
            printf("Theta R: %f\n", theta_r);

            double sin_theta_l = sin(theta_l);
            double cos_theta_l = cos(theta_l);

            printf("cos(Theta L): %f\n", cos_theta_l);
            printf("sin(Theta L): %f\n", sin_theta_l);

            double sin_theta_r = sin(theta_r);
            double cos_theta_r = cos(theta_r);

            printf("cos(Theta R): %f\n", cos_theta_r);
            printf("sin(Theta R): %f\n", sin_theta_r);

            double u_tmp[4][4] = {
                {1, 0, 0, 0},
                {0, 1, 0, 0},
                {0, 0, 1, 0},
                {0, 0, 0, 1},
            };

            double u_tmp_trans[4][4] = {
                {1, 0, 0, 0},
                {0, 1, 0, 0},
                {0, 0, 1, 0},
                {0, 0, 0, 1},
            };

            u_tmp[i][i] = cos_theta_l;
            u_tmp[j][j] = cos_theta_l;
            u_tmp[i][j] = -sin_theta_l;
            u_tmp[j][i] = sin_theta_l;

            u_tmp_trans[i][i] = cos_theta_l;
            u_tmp_trans[j][j] = cos_theta_l;
            u_tmp_trans[i][j] = sin_theta_l;
            u_tmp_trans[j][i] = -sin_theta_l;

            printf("--------------u_tmp----------\n");
            mat_print(u_tmp);
            printf("--------------u_tmp_T--------\n");
            mat_print(u_tmp_trans);

            double v_tmp[4][4] = {
                {1, 0, 0, 0},
                {0, 1, 0, 0},
                {0, 0, 1, 0},
                {0, 0, 0, 1},
            };

            double v_tmp_trans[4][4] = {
                {1, 0, 0, 0},
                {0, 1, 0, 0},
                {0, 0, 1, 0},
                {0, 0, 0, 1},
            };

            v_tmp[i][i] = cos_theta_r;
            v_tmp[j][j] = cos_theta_r;
            v_tmp[i][j] = -sin_theta_r;
            v_tmp[j][i] = sin_theta_r;

            v_tmp_trans[i][i] = cos_theta_r;
            v_tmp_trans[j][j] = cos_theta_r;
            v_tmp_trans[i][j] = sin_theta_r;
            v_tmp_trans[j][i] = -sin_theta_r;

            printf("--------------v_tmp----------\n");
            mat_print(v_tmp);
            printf("--------------v_tmp_T--------\n");
            mat_print(v_tmp_trans);

            double u_prime[4][4];
            double v_trans_prime[4][4];
            double m_prime_tmp[4][4];
            double m_prime[4][4];

            mat_mul(4, u, u_tmp_trans, u_prime);             // [U][U_ij_T] = [U']
            mat_mul(4, u_tmp, m, m_prime_tmp);               // [U_ij][M] = [M'_tmp]
            mat_mul(4, m_prime_tmp, v_tmp_trans, m_prime);   // [M_tmp][V_ij_T] = [M']
            mat_mul(4, v_tmp_trans, v_trans, v_trans_prime); // [V_ij][V_T] = [V'_T] <- I need to do this wrong to get it to work?????

            printf("--------------u'-------------\n");
            mat_print(u_prime);
            printf("-------------v_T'------------\n");
            mat_print(v_trans_prime);
            printf("--------------m'-------------\n");
            mat_print(m_prime);

            for (int row = 0; row < 4; row++)
            {
                for (int col = 0; col < 4; col++)
                {
                    u[row][col] = u_prime[row][col];
                    v_trans[row][col] = v_trans_prime[row][col];
                    m[row][col] = m_prime[row][col];
                }
            }

            printf("--------------u--------------\n");
            mat_print(u);
            printf("-------------v_T-------------\n");
            mat_print(v_trans);
            printf("--------------m--------------\n");
            mat_print(m);
            printf("**********************************************************************\n");
        }
    }
}

int main(void)
{
    for (int i = 0; i < 4; i++)
    {
        sweep();
        mat_print(m);
    }
}