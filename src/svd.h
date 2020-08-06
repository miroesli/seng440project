/*
 *
 * svd.h
 *
 */
#ifndef svd_h
#define svd_h

#include <math.h>
#include "svd_math.h"

/**
 * @brief Performes a single sweep of the svd algorithm
 * 
 */
void sweep(floating_point_t m[4][4], floating_point_t u[4][4], floating_point_t v_trans[4][4]);

#endif