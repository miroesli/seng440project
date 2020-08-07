/*
 *
 * main.c
 *
 */

#include "config.h"
#include "svd.h"
#include "svd_math.h"

floating_point_t m[4][4] = {
    {31, 77, -11, 26},
    {-42, 14, 79, -53},
    {-68, -10, 45, 90},
    {34, 16, 38, -19},
};

floating_point_t u[4][4] = {
    {1, 0, 0, 0},
    {0, 1, 0, 0},
    {0, 0, 1, 0},
    {0, 0, 0, 1},
};

floating_point_t v_trans[4][4] = {
    {1, 0, 0, 0},
    {0, 1, 0, 0},
    {0, 0, 1, 0},
    {0, 0, 0, 1},
};

/**
 * @brief Helper function to print a 4x4 matrix
 * 
 * @param mat 
 */
void mat_print(size_t size, floating_point_t mat[size][size])
{
    for (int row = 0; row < size; row++)
    {
        printf("[ ");
        for (int col = 0; col < size - 1; col++)
        {
            printf("%f, ", mat[row][col]);
        }
        printf("%f ]\n", mat[row][size - 1]);
    }
}

int main(void)
{

    for (int i = 0; i < 4; i++)
    {
        sweep(4, m, u, v_trans);
        mat_print(4, m);
        printf("\n");
    }
}