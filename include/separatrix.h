/***************************************************/
/*												   */
/*		SEPARATRIX HEADER						   */
/*												   */
/*		JG LEE									   */
/*		28/DEC/2015								   */
/***************************************************/

#ifndef	_SEPARATRIX_H
#define	_SEPARATRIX_H

#include "vega.h"

void FIELD(SYSALL& sys);
void BILINEAR(SYSALL& sys);
TESTPOINT NULLSEARCH(SYSALL& sys, int i, int j);

TESTPOINT CROSSLIMT(SYSALL& sys, double r1, double z1);
TESTPOINT CROSSDIVT(SYSALL& sys, double r1, double z1, double r2, double z2);

void GRAD0(SYSALL& sys);

void MULTIX(SYSALL& sys);

void FIRSTSTEP(SYSALL& sys, int iXp, dvector r, dvector z, double psi);
void EXEMPTION(SYSALL& sys, int iXp, double r, double z,
	dvector newr, dvector newz, double dmin);
void GETFIELD(SYSALL& sys, double r, double z, dvector dr, dvector dz, double d);
void VECTORFOLLOW(SYSALL& sys, double psi, double& r, double& z, 
	double r0, double z0, double d, int drct);

void FRAME(SYSALL& sys);

#endif