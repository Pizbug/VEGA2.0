/***************************************************/
/*												   */
/*		VEGA FUNCTION							   */
/*												   */
/*		JG LEE									   */
/*		14/JAN/2016								   */
/***************************************************/

#define	_CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

#include "parameter.h"
#include "svc.h"
#include "vega.h"
#include "memalloc.h"
#include "separatrix.h"


//***********************************************************NEED TO BE CHANGE!!!*********************//
//	Derivative function
//	Central difference schemes are used
double DERIV_R(int rN, dvector r, dvector z, dmatrix V, int i, int j) {
	double	dVr;

	if (i == 0)
		dVr = (V[i + 1][j] - V[i][j]) / (r[i + 1] - r[i]);
	else if (i == rN - 1)
		dVr = (V[i][j] - V[i - 1][j]) / (r[i] - r[i - 1]);
	else if (i == 1 || i == rN - 2)
		dVr = (V[i + 1][j] - V[i - 1][j]) / (r[i + 1] - r[i - 1]);
	else
		dVr = (-V[i + 2][j] + 8 * V[i + 1][j] - 8 * V[i - 1][j] + V[i - 2][j]) / 6.0 / (r[i + 1] - r[i - 1]);

	return	dVr;
}

double DERIV_Z(int zN, dvector r, dvector z, dmatrix V, int i, int j) {
	double	dVz;

	if (j == 0)
		dVz = (V[i][j + 1] - V[i][j]) / (z[j + 1] - z[j]);
	else if (j == zN - 1)
		dVz = (V[i][j] - V[i][j - 1]) / (z[j] - z[j - 1]);
	else if (j == 1 || j == zN - 2)
		dVz = (V[i][j + 1] - V[i][j - 1]) / (z[j + 1] - z[j - 1]);
	else
		dVz = (-V[i][j + 2] + 8 * V[i][j + 1] - 8 * V[i][j - 1] + V[i][j - 2]) / 6.0 / (z[j + 1] - z[j - 1]);

	return	dVz;
}

//	Derivative Function for whole regime
void DERIV(int rN, int zN, dvector r, dvector z, dmatrix V, dmatrix dVr, dmatrix dVz) {
	int	i, j;
	DoAllDomain(i, rN) {
		DoAllDomain(j, zN) {
			if (i == 0)
				dVr[i][j] = (V[i + 1][j] - V[i][j]) / (r[i + 1] - r[i]);
			else if (i == rN - 1)
				dVr[i][j] = (V[i][j] - V[i - 1][j]) / (r[i] - r[i - 1]);
			else if (i == 1 || i == rN - 2)
				dVr[i][j] = (V[i + 1][j] - V[i - 1][j]) / (r[i + 1] - r[i - 1]);
			else
				dVr[i][j] = (-V[i + 2][j] + 8 * V[i + 1][j] - 8 * V[i - 1][j] + V[i - 2][j]) / 6.0 / (r[i + 1] - r[i - 1]);

			if (j == 0)
				dVz[i][j] = (V[i][j + 1] - V[i][j]) / (z[j + 1] - z[j]);
			else if (j == rN - 1)
				dVz[i][j] = (V[i][j] - V[i][j - 1]) / (z[j] - z[j - 1]);
			else if (j == 1 || j == rN - 2)
				dVz[i][j] = (V[i][j + 1] - V[i][j - 1]) / (z[j + 1] - z[j - 1]);
			else
				dVz[i][j] = (-V[i][j + 2] + 8 * V[i][j + 1] - 8 * V[i][j - 1] + V[i][j - 2]) / 6.0 / (z[j + 1] - z[j - 1]);
		}
	}
}
//***********************************************************NEED TO BE CHANGE!!!*********************//



/*	CROSSFUNC
*/
TESTPOINT CROSSFUNC(double r1, double z1, double r2, double z2,
					double br1, double bz1, double br2, double bz2) {
	TESTPOINT	V;
	double		det, d1, d2, mult1, mult2;

	det = (-(r2 - r1)) * (bz2 - bz1) + (z2 - z1) * (br2 - br1);
	if (det != 0) {
		d1 = (-(br1 - r1) * (bz2 - bz1) + (bz1 - z1) * (br2 - br1));
		mult1 = d1 / det;
		if ((mult1 >= 0) && (mult1 <= 1)) {
			d2 = ((r2 - r1) * (bz1 - z1) - (z2 - z1) * (br1 - r1));
			mult2 = d2 / det;
			if ((mult2 >= 0) && (mult2 <= 1)) {
				V.id = 1;							//	Test Point is out of Limiter
				V.r = r1 + (r2 - r1) * mult1;		//	Cross point
				V.z = z1 + (z2 - z1) * mult1;

				return V;
			}
		}
	}
	V.id = 0;
	V.r = r2;
	V.z = z2;
	V.i = 0;

	return V;
}


/*	GETPOINT
Find the desired psi point in the given line
*/
TESTPOINT GETPOINT(SYSALL & sys, double r1, double z1,
	double r2, double z2, double psi, double ratio) {
	usealias_geo(g, sys);
	usealias_fld(f, sys);
	usealias_mes(m, sys);
	TESTPOINT	V;
	int i, n = 10;
	int i1, j1, i2, j2;
	double p1, p2;
	double r, z, dr, dz;

	//	ratio: enlength the line
	dr = ((r2 - (r1 - r2) * ratio) - (r1 - (r2 - r1) * ratio)) / (double)(n);
	dz = ((z2 - (z1 - z2) * ratio) - (z1 - (z2 - z1) * ratio)) / (double)(n);
	r = r1 - (r2 - r1) * ratio;
	z = z1 - (z2 - z1) * ratio;

	i1 = (int)((r - g.r[0]) / (g.r[1] - g.r[0]));	//	Fine index of r, z
	j1 = (int)((z - g.z[0]) / (g.z[1] - g.z[0]));	//	Assume Equi-distribution
	p1 = f.Psi00[i1][j1] + f.Psi10[i1][j1] * r + f.Psi01[i1][j1] * z + f.Psi11[i1][j1] * r * z;
	DoAllDomain(i, n) {
		r += dr;
		z += dz;

		i2 = (int)((r - g.r[0]) / (g.r[1] - g.r[0]));	//	Fine index of r, z
		j2 = (int)((z - g.z[0]) / (g.z[1] - g.z[0]));	//	Assume Equi-distribution
		if (i2 < 0)		i2 = 0;
		if (j2 < 0)		j2 = 0;
		if (i2 > g.rN - 2)	i2 = g.rN - 2;
		if (j2 > g.zN - 2)	j2 = g.zN - 2;
		p2 = f.Psi00[i2][j2] + f.Psi10[i2][j2] * r + f.Psi01[i2][j2] * z + f.Psi11[i2][j2] * r * z;

		if (BTWN(psi, p1, p2)) {
			V.id = 1;
			V.r = LINEAR1D(r - dr, r, p1, p2, psi);
			V.z = LINEAR1D(z - dz, z, p1, p2, psi);
			V.v = psi;
			V.i = i1;
			V.j = j1;
			return V;
		}
		p1 = p2;
	}
	V.id = 0;
	V.r = 0;
	V.z = 0;
	V.v = p2;
	V.i = 0;
	V.j = 0;

	return V;
}

double BiliInterp(double **f, double *x, double *y, 
	double xx, double yy, int i, int j) {
	return (f[i][j] * (x[i + 1] - xx) * (y[j + 1] - yy)
		+ f[i + 1][j] * (xx - x[i]) * (y[j + 1] - yy)
		+ f[i][j + 1] * (x[i + 1] - xx) * (yy - y[j])
		+ f[i + 1][j + 1] * (xx - x[i]) * (yy - y[j]))
		/ ((x[i + 1] - x[i]) * (y[j + 1] - y[j]));
}

