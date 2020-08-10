/*
 *
 * svd.h
 *
 */
#ifndef svd_h
#define svd_h

#include <math.h>
#include "svd_math.h"
#include "memory.h"
#include "config.h"

#define SIZE 4

/**
  * @brief Performes a single sweep of the svd algorithm
  *
  */
void sweep(floating_point_t m[SIZE][SIZE], floating_point_t u[SIZE][SIZE], floating_point_t v_trans[SIZE][SIZE]);

#endif