/**
 * Parallel Sobol' sequence generator. To use, simply
 *
 * $include "sobolp.h"
 *
 * struct sobolp sp;
 * double rands[dimensions];
 * sobolp_init(&sp,dimensions,seed);
 * sobolp_generateSamples(&sp,rands);
 *
 * and the array 'rands' contains the random double values in the interval
 * [0,1].
 *
 * Massaged from the original into parallel version by Austen McDonald,
 * <austen@cc.gatech.edu>
 */
#ifndef SOBOLP_H
#define SOBOLP_H

#define MAXDIM 5
#define VMAX 30

struct sobolp {
	double sequence[MAXDIM];
	int x[MAXDIM];
	int v[MAXDIM][VMAX];
	double RECIPD;
	int _dim;
	int _skip;
	unsigned long _nextn;
	unsigned long cur_seed;
};

void sobolp_init(struct sobolp*,int,unsigned long);
void sobolp_generateSamples(struct sobolp*,double*);
#endif
