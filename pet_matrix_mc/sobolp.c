/**
 * Parallel sobol sequence generator. Altered from original by Austen McDonald
 * <austen@cc.gatech.edu>
 */

/*==========================================================================
//
// Filename: QMCSobol.cc
//
// Creator: Simon Premoze
// Date:    February 12, 1999
// Revision: $Id$
//
// Description:
//
// Report problems and suggestions to premoze@cs.utah.edu
*/


#include "sobolp.h"
#include <math.h>

/*
// Return GRAY code
*/

#define GRAY(n) (n ^ ( n >> 1 ))

static void
nextSobol(struct sobolp* config, int seed )
{
	int c,gray,i;
  seed += config->_skip + 1;
  for( i = 0; i < config->_dim; i++ ) {
    config->x[i] = 0;
    c = 1;
    gray = GRAY( seed );
    while( gray ) {
      if( gray & 1 ) {
        config->x[i] = config->x[i] ^ ( config->v[i][c-1] << ( VMAX - c ) );
      }
      c++;
      gray >>= 1;
    }
    
    config->sequence[i] = config->x[i] * config->RECIPD;
  }
	
}

/*
// [in] numSamples    Number of samples desired
// [in] numDimensions The dimension of the sample spaces
// [out] realNum      Number actually created (<= numSamples)
*/

void
sobolp_generateSamples( struct sobolp* config, double* samples)
{
	int i;
	nextSobol(config, config->cur_seed);
	config->cur_seed++;
	for(i = 0; i < config->_dim; i++ )
		samples[i] = config->sequence[i];
}

static void nextSobolNoSeed(struct sobolp* config) {
	int c = 1;
	int i;
	int save = config->_nextn;
	while((save %2) == 1) {
			c += 1;
			save  = save /2;
	}
	for(i=0;i<config->_dim;i++) {
		config->x[i] = config->x[i]^(config->v[i][c-1]<< (VMAX-c));
		config->sequence[i] = config->x[i]*config->RECIPD;
	}
	config->_nextn += 1;
}

void
sobolp_init(struct sobolp* config, int dim, unsigned long seed)
{
  int d[MAXDIM], POLY[MAXDIM];
  int save;
  int m,i,j,k;
 
  config->_dim = dim;
  config->_nextn = 0;
  config->RECIPD = 1.0 / pow( 2.0, VMAX );
  config->cur_seed = seed;
  
  POLY[0] = 3;  d[0] = 1;  /* x + 1         */
  POLY[1] = 7;  d[1] = 2;  /* x^2 + x + 1   */
  POLY[2] = 11; d[2] = 3;  /* x^3 + x + 1   */
  POLY[3] = 19; d[3] = 4;  /* x^4 + x  + 1  */
  POLY[4] = 37; d[4] = 5;  /* x^5 + x^2 + 1 */
  
  for(i=0; i < config->_dim; i++ )
    for(j = 0; j < d[i]; j++ )
      config->v[i][j] = 1;
  
  for( i = 0; i < config->_dim; i++ )
    for( j = d[i]; j < VMAX; j++ ) {
      config->v[i][j] = config->v[i][j-d[i]];
      save = POLY[i];
      m = pow( 2, d[i] );
      for( k = d[i]; k > 0; k-- ) {
        config->v[i][j] = config->v[i][j] ^ m*(save%2)*config->v[i][j-k];
        save = save/2;
	m = m/2;
      }
    }
  
  for( i = 0; i < config->_dim; i++ )
    config->x[i]=0;
  config->_skip = pow( 2, 6 );
  for( i = 1; i <= config->_skip; i++ )
    nextSobolNoSeed(config);
}
