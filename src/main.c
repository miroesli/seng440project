/*
 *
 * main.c
 *
 */

#include "config.h"
#include "svd.h"

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
void mat_print(floating_point_t mat[4][4])
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

int main(void)
{

    for (int i = 0; i < 4; i++)
    {
        sweep(m, u, v_trans);
        mat_print(m);
        printf("\n");
    }
}