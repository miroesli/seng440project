#ifndef svd_math_h
#define svd_math_h
#include <inttypes.h>
#include <stddef.h>
#include <math.h>

typedef float floating_point_t;
typedef int32_t fixed_point_t;
typedef int64_t fixed_point_double_t;

typedef fixed_point_t fixed_point_u_t;
typedef fixed_point_t fixed_point_v_t;
typedef fixed_point_t fixed_point_m_t;
typedef fixed_point_t fixed_point_m_tmp_t;

typedef fixed_point_double_t fixed_point_u_dp_t;
typedef fixed_point_double_t fixed_point_v_dp_t;
typedef fixed_point_double_t fixed_point_m_dp_t;
typedef fixed_point_double_t fixed_point_m_tmp_dp_t;

#include "config.h"
// #include "tables.h"
#include "arctan_lookup_table.h"
#include "sin_lookup_table.h"
#include "cos_lookup_table.h"


// Table definitions
// exclusive end
#define ARCTAN_RANGE 10 
#define VALUES_IN_RANGE 10000
#define SINCOS_RANGE M_PI * 2
// Scale Factors
// TODO check if factor u can be used instead?
#define SCALE_FACTOR_ARCTAN 30
#define SCALE_FACTOR_MOD 28
#define SCALE_FACTOR_U 31
#define SCALE_FACTOR_U_DP 2 * SCALE_FACTOR_U + 1
#define SCALE_FACTOR_V 31
#define SCALE_FACTOR_V_DP 2 * SCALE_FACTOR_V + 1
#define SCALE_FACTOR_M 17
#define SCALE_FACTOR_M_tmp_DP SCALE_FACTOR_M + SCALE_FACTOR_V + 1
#define SCALE_FACTOR_M_tmp SCALE_FACTOR_M_tmp_DP - (sizeof(fixed_point_double_t) - sizeof(fixed_point_t)) * 8
#define SCALE_FACTOR_M_DP SCALE_FACTOR_M_tmp + SCALE_FACTOR_V + 1

fixed_point_t convert_to_fixed(floating_point_t f, size_t scale_factor);
fixed_point_double_t fixed_point_mul(fixed_point_t, fixed_point_t);
fixed_point_t truncate(fixed_point_m_tmp_dp_t);
floating_point_t convert_to_floating(fixed_point_double_t f, size_t scale_factor);
fixed_point_t arctan_lookup(floating_point_t frac);
fixed_point_t sin_lookup(fixed_point_t theta);
fixed_point_t cos_lookup(fixed_point_t theta);

static const fixed_point_u_t one_u = (1 << SCALE_FACTOR_U);
static const fixed_point_v_t one_v = (1 << SCALE_FACTOR_V);

#endif