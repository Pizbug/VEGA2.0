/***************************************************/
/*												   */
/*		MESHGEN HEADER							   */
/*												   */
/*		JG LEE									   */
/*		28/DEC/2015								   */
/***************************************************/

#ifndef	_MESHGEN_H
#define	_MESHGEN_H

#include "vega.h"

void MESHPAR(SYSALL& sys);

TESTPOINT GETPOINT(SYSALL& sys, double r1, double z1,
	double r2, double z2, double psi, double ratio);

TESTPOINT NORVECTRC(SYSALL& sys, double r, double z,
	double r2, double z2, int LR, double psi, double w);

void STRETCHPAR(dvector dist, double L, int n,
	double E, double D, double S);

void STRETCHING(SYSALL& sys);

void MESHGEN(SYSALL& sys);

TESTPOINT BVEC(SYSALL &sys, int k, int LR, int YCo);
void PARTMESHGEN(SYSALL &sys, int &k, int YCo, int nDi, int nL, int nR, double *psiL, double *psiR);
void MESHGEN2(SYSALL& sys);


#endif

