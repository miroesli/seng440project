#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "arm_neon.h"
#include "svd_math.h"

#define N 4
#define M 4

void print_matrix(const fixed_point_double_t *m)
{
    for (int i = 0; i < N; i++)
    {
        printf("[");
        for (int j = 0; j < M - 1; j++)
        {
            printf(" %d", *(m + i * N + j));
        }
        printf(" %d ]\n", *(m + i * N + M - 1));
    }
    printf("\n");
}

static volatile fixed_point_double_t X[N][M] = {
    {0, 1, 2, 3},
    {4, 5, 6, 7},
    {8, 9, 10, 11},
    {12, 13, 14, 15},
};

static volatile fixed_point_double_t Y[N][M] = {
    {0, 1, 2, 3},
    {4, 5, 6, 7},
    {8, 9, 10, 11},
    {12, 13, 14, 15},
};

static volatile fixed_point_double_t OUT[N][M];
static volatile fixed_point_double_t OUT_NEON[N][M];

void matrix_multiply(
    volatile fixed_point_double_t LHS[N][M],
    volatile fixed_point_double_t RHS[N][M],
    volatile fixed_point_double_t RESULT[N][M])
{
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            RESULT[i][j] = 0;
            for (int k = 0; k < 4; k++)
            {
                RESULT[i][j] += LHS[i][k] * RHS[k][j];
            }
        }
    }
}

void matrix_multiply_NEON()
{
    int32x4_t Y_row_0, Y_row_1, Y_row_2, Y_row_3, out_neon;

    Y_row_0 = vld1q_s32((const fixed_point_double_t *)&Y[0][0]);
    Y_row_1 = vld1q_s32((const fixed_point_double_t *)&Y[1][0]);
    Y_row_2 = vld1q_s32((const fixed_point_double_t *)&Y[2][0]);
    Y_row_3 = vld1q_s32((const fixed_point_double_t *)&Y[3][0]);

    for (int i = 0; i < N; i++)
    {
        out_neon = vmulq_n_s32(Y_row_0, X[i][0]);
        out_neon = vaddq_s32(vmulq_n_s32(Y_row_1, X[i][1]), out_neon);
        out_neon = vaddq_s32(vmulq_n_s32(Y_row_2, X[i][2]), out_neon);
        out_neon = vaddq_s32(vmulq_n_s32(Y_row_3, X[i][3]), out_neon);
        // out_neon = vshrq_n_s32(out_neon, SHIFT_AMOUNT);
        vst1q_s32((fixed_point_double_t *)&OUT_NEON[i][0], out_neon);
    }
}

int main(void)
{

    matrix_multiply(X, Y, OUT);
    matrix_multiply_NEON();
    print_matrix((const fixed_point_double_t *)&OUT[0][0]);
    print_matrix((const fixed_point_double_t *)&OUT_NEON[0][0]);

    exit(0);
}