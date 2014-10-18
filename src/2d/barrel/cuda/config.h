#pragma once

// no time-of-flight MC
#define NO_TOF 0

// fast aprox intrinsics
// dont use fpow and sqrtf intrinsics -> bad output, problem with precision
#define SINGLE_PRECISION_INTRINSIC 1

// test warp divergence path in "up", "down" intersection methods
#define WARP_DIVERGENCE_TEST 0

// calculate number of cycles per PRNG and number of execution of intersection
// methods
#define CLOCK_TEST 0

// FIXME: this must be dynamic!
#define NUMBER_OF_DETECTORS 204
#define LORS (NUMBER_OF_DETECTORS * (NUMBER_OF_DETECTORS + 1)) / 2

#define TP 6.28318530f
#define ITP 0.15915494f
