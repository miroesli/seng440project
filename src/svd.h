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

 /**
  * @brief Performes a single sweep of the svd algorithm
  *
  */
void sweep(const size_t size, floating_point_t m[size][size], floating_point_t u[size][size], floating_point_t v_trans[size][size]);

#endif