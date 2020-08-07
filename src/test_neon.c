#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "arm_neon.h"

#define N 16
#define M 32

volatile int16_t A[N][M], B[N][M], SUM[N][M];

int main(void)
{
    int i, j;
    int16x4_t A_neon, B_neon, SUM_neon;

    for (i = 0; i < N; i++)
        for (j = 0; j < M; j++)
        {
            A[i][j] = 7 * i;
            B[i][j] = 15 * j;
        }

    for (i = 0; i < 4; i++)
        for (j = 0; j < 2; j++)
            printf("A[%i][%i] = %i\n", i, j, A[i][j]);
    printf("\n");

    for (i = 0; i < 4; i++)
        for (j = 0; j < 2; j++)
            printf("B[%i][%i] = %i\n", i, j, B[i][j]);
    printf("\n");

    for (i = 0; i < N; i += 4)
        for (j = 0; j < M; j++)
        {
            A_neon = vld1_s16((const int16_t *)&A[i][j]);
            B_neon = vld1_s16((const int16_t *)&B[i][j]);
            SUM_neon = vadd_s16(A_neon, B_neon);
            vst1_s16((int16_t *)&SUM[i][j], SUM_neon);
        }

    for (i = 0; i < 4; i++)
        for (j = 0; j < 2; j++)
            printf("SUM[%i][%i] = %i\n", i, j, SUM[i][j]);
    printf("\n");

    exit(0);
}