/*
 *
 * svd.h
 *
 */
#ifndef svd_h
#define svd_h

#include <math.h>

/**
 * @brief Performes a single sweep of the svd algorithm
 * 
 */
void sweep(double m[4][4], double u[4][4], double v_trans[4][4]);

#endif