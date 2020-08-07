#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "arm_neon.h"
#include "svd_math.h"

#define N 4
#define M 4

volatile fixed_point_t A[N][M], B[N][M], SUM[N][M];

void print_matrix(fixed_point_t m[N][M])
{
    for (int i = 0; i < N; i++)
    {
        printf("[");
        for (int j = 0; j < M - 1; j++)
        {
            printf(" %d", m[i][j]);
        }
        printf(" ]\n");
    }
}

int main(void)
{
    int i, j;
    fixed_point_t A_neon, B_neon, SUM_neon;

    for (i = 0; i < N; i++)
        for (j = 0; j < M; j++)
        {
            A[i][j] = 7 * i;
            B[i][j] = 15 * j;
        }

    print_matrix(A);
    print_matrix(B);

    for (i = 0; i < N; i += 4)
        for (j = 0; j < M; j++)
        {
            A_neon = vld1_s32((const fixed_point_t *)&A[i][j]);
            B_neon = vld1_s32((const fixed_point_t *)&B[i][j]);
            SUM_neon = vadd_s32(A_neon, B_neon);
            vst1_s32((fixed_point_t *)&SUM[i][j], SUM_neon);
        }

    print_matrix(SUM);

    exit(0);
}