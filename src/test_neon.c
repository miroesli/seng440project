#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "arm_neon.h"
#include "svd_math.h"

#define N 4
#define M 4

void print_matrix(const fixed_point_t *m)
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

volatile fixed_point_t X[N][M] = {
    {0x1, 0x2, 0x3, 0x4},
    {0x5, 0x6, 0x7, 0x8},
    {0x9, 0xA, 0xB, 0xC},
    {0xD, 0xE, 0xF, 0x0},
};

volatile fixed_point_t Y[N][M] = {
    {0x1, 0x2, 0x3, 0x4},
    {0x5, 0x6, 0x7, 0x8},
    {0x9, 0xA, 0xB, 0xC},
    {0xD, 0xE, 0xF, 0x0},
};

volatile fixed_point_t OUT[N][M];

int main(void)
{
    int i, j;
    int32x4_t Y_row_0, Y_row_1, Y_row_2, Y_row_3, out_neon;

    print_matrix((const fixed_point_t *)&X[0][0]);
    print_matrix((const fixed_point_t *)&Y[0][0]);

    Y_row_0 = vld1q_s32((const fixed_point_t *)&Y[0][0]);
    Y_row_1 = vld1q_s32((const fixed_point_t *)&Y[1][0]);
    Y_row_2 = vld1q_s32((const fixed_point_t *)&Y[2][0]);
    Y_row_3 = vld1q_s32((const fixed_point_t *)&Y[3][0]);

    for (int i = 0; i < N; i++)
    {
        out_neon = vmulq_n_s32(Y_row_0, X[i][0]);
        out_neon = vaddq_s32(vmulq_n_s32(Y_row_1, X[i][1]), out_neon);
        out_neon = vaddq_s32(vmulq_n_s32(Y_row_2, X[i][2]), out_neon);
        out_neon = vaddq_s32(vmulq_n_s32(Y_row_3, X[i][3]), out_neon);
        vst1q_s32((fixed_point_t *)&OUT[i][0], out_neon);
    }

    print_matrix((const fixed_point_t *)&OUT[0][0]);

    exit(0);
}