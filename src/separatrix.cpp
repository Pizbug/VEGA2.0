/***************************************************/
/*												   */
/*		SEPARATRIX FUNCTION							   */
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


//	Get field data
void FIELD(SYSALL & sys)
{
	usealias_geo(g, sys);
	usealias_fld(f, sys);
	usealias_mes(m, sys);
	int 	i, j;
	dmatrix dVr, dVz;
	double	Sq;
	//FILE *out;

	dVr = alloc_matrix<double>(0, rN - 1, 0, zN - 1);
	dVz = alloc_matrix<double>(0, rN - 1, 0, zN - 1);

	DERIV(rN, zN, g.r, g.z, g.Psi, dVr, dVz);

	DoAllDomain(i, rN) {
		DoAllDomain(j, zN) {
			f.Br[i][j] = -1.0 / g.r[i] * dVz[i][j];// *(2.0 * PI);
			f.Bz[i][j] = 1.0 / g.r[i] * dVr[i][j];// *(2.0 * PI);

			Sq = sqrt((f.Br[i][j])*(f.Br[i][j]) + (f.Bz[i][j])*(f.Bz[i][j]));
			f.Bp[i][j] = Sq;
		}
	}
	free_matrix<double>(dVr, 0, rN - 1, 0, zN - 1);
	free_matrix<double>(dVz, 0, rN - 1, 0, zN - 1);

//	printf("(FIELD) -> ");
	/*
	fclose(out);
	out = fopen("Br.txt", "w");
	DoAllDomain(j, zN) {
		DoAllDomain(i, rN) {
			fprintf(out, "%lf\t", f.Br[i][j]);
		}
		fprintf(out, "\n");
	}
	fclose(out);
	out = fopen("Bz.txt", "w");
	DoAllDomain(j, zN) {
		DoAllDomain(i, rN) {
			fprintf(out, "%lf\t", f.Bz[i][j]);
		}
		fprintf(out, "\n");
	}
	fclose(out);
	out = fopen("Bp.txt", "w");
	DoAllDomain(j, zN) {
		DoAllDomain(i, rN) {
			fprintf(out, "%lf\t", f.Bp[i][j]);
		}
		fprintf(out, "\n");
	}
	fclose(out);*/
}


//	Get Bilinear coefficients of all regime
/*	Associated with
(i,j+1)		(i+1,j+1)

(i,j)		 (i+1,j)
*/
void BILINEAR(SYSALL& sys) {
	usealias_geo(g, sys);
	usealias_fld(f, sys);
	usealias_mes(m, sys);
	int		i, j;
	double	dxy;

	DoAllDomain(i, rN - 1) {
		DoAllDomain(j, zN - 1) {
			dxy = (g.r[i + 1] - g.r[i]) * (g.z[j + 1] - g.z[j]);
			f.Psi00[i][j] = (g.Psi[i][j] * g.r[i + 1] * g.z[j + 1] - g.Psi[i + 1][j] * g.r[i] * g.z[j + 1]
				- g.Psi[i][j + 1] * g.r[i + 1] * g.z[j] + g.Psi[i + 1][j + 1] * g.r[i] * g.z[j]) / dxy;
			f.Br00[i][j] = (f.Br[i][j] * g.r[i + 1] * g.z[j + 1] - f.Br[i + 1][j] * g.r[i] * g.z[j + 1]
				- f.Br[i][j + 1] * g.r[i + 1] * g.z[j] + f.Br[i + 1][j + 1] * g.r[i] * g.z[j]) / dxy;
			f.Bz00[i][j] = (f.Bz[i][j] * g.r[i + 1] * g.z[j + 1] - f.Bz[i + 1][j] * g.r[i] * g.z[j + 1]
				- f.Bz[i][j + 1] * g.r[i + 1] * g.z[j] + f.Bz[i + 1][j + 1] * g.r[i] * g.z[j]) / dxy;

			f.Psi10[i][j] = ((g.Psi[i + 1][j] - g.Psi[i][j])*g.z[j + 1] + (g.Psi[i][j + 1] - g.Psi[i + 1][j + 1])*g.z[j]) / dxy;
			f.Br10[i][j] = ((f.Br[i + 1][j] - f.Br[i][j])*g.z[j + 1] + (f.Br[i][j + 1] - f.Br[i + 1][j + 1])*g.z[j]) / dxy;
			f.Bz10[i][j] = ((f.Bz[i + 1][j] - f.Bz[i][j])*g.z[j + 1] + (f.Bz[i][j + 1] - f.Bz[i + 1][j + 1])*g.z[j]) / dxy;

			f.Psi01[i][j] = ((g.Psi[i + 1][j] - g.Psi[i + 1][j + 1])*g.r[i] + (g.Psi[i][j + 1] - g.Psi[i][j])*g.r[i + 1]) / dxy;
			f.Br01[i][j] = ((f.Br[i + 1][j] - f.Br[i + 1][j + 1])*g.r[i] + (f.Br[i][j + 1] - f.Br[i][j])*g.r[i + 1]) / dxy;
			f.Bz01[i][j] = ((f.Bz[i + 1][j] - f.Bz[i + 1][j + 1])*g.r[i] + (f.Bz[i][j + 1] - f.Bz[i][j])*g.r[i + 1]) / dxy;

			f.Psi11[i][j] = (g.Psi[i][j] - g.Psi[i + 1][j] - g.Psi[i][j + 1] + g.Psi[i + 1][j + 1]) / dxy;
			f.Br11[i][j] = (f.Br[i][j] - f.Br[i + 1][j] - f.Br[i][j + 1] + f.Br[i + 1][j + 1]) / dxy;
			f.Bz11[i][j] = (f.Bz[i][j] - f.Bz[i + 1][j] - f.Bz[i][j + 1] + f.Bz[i + 1][j + 1]) / dxy;
		}
	}
//	printf("(BILINEAR) -> ");
}



//	Find the Gradient 0 point in given regime
/*	The Given Regime is
(i,j+1)		(i+1,j+1)

(i,j)		 (i+1,j)
*/
TESTPOINT NULLSEARCH(SYSALL& sys, int i, int j) {
	usealias_geo(g, sys);
	usealias_fld(f, sys);
	usealias_mes(m, sys);
	TESTPOINT	N;
	double	a, b, c;
	double	r1, r2, z1, z2;

	N.id = 0;
	N.r = 0;
	N.z = 0;
	N.i = 0;
	N.j = 0;
	N.v = 0;

	a = f.Br11[i][j] * f.Bz01[i][j] - f.Br01[i][j] * f.Bz11[i][j];
	b = f.Br11[i][j] * f.Bz00[i][j] - f.Br01[i][j] * f.Bz10[i][j]
		+ f.Br10[i][j] * f.Bz01[i][j] - f.Br00[i][j] * f.Bz11[i][j];
	c = f.Br10[i][j] * f.Bz00[i][j] - f.Br00[i][j] * f.Bz10[i][j];

	if ((4 * a*c) < (b*b)) {
		if ((b > 0) && (a != 0)) {
			z1 = -2 * c / (sqrt(b*b - 4 * a*c) + b);
			z2 = (-b - sqrt(b*b - 4 * a*c)) / (2 * a);
		}
		else if (a != 0) {
			z1 = 2 * c / (sqrt(b*b - 4 * a*c) - b);
			z2 = (-b + sqrt(b*b - 4 * a*c)) / (2 * a);
		}
		else if (b != 0) {
			z1 = -c / b;
		}
		else {
			z1 = 2 * g.z[0] - g.z[1];
		}
		if ((z1 >= g.z[j]) && (z1 <= g.z[j + 1])) {
			if ((f.Br10[i][j] + f.Br11[i][j] * z1) != 0)
				r1 = (-f.Br00[i][j] - f.Br01[i][j] * z1) / (f.Br10[i][j] + f.Br11[i][j] * z1);
			else
				r1 = 2 * g.r[0] - g.r[1];
			if ((r1 >= g.r[i]) && (r1 <= g.r[i + 1])) {
				N.id = 1;
				N.r = r1;
				N.z = z1;
				N.i = i;
				N.j = j;
				N.v = f.Psi00[i][j] + f.Psi10[i][j] * g.r[i] + f.Psi01[i][j] * g.z[j] + f.Psi11[i][j] * g.r[i] * g.z[j];
			}
		}

		if ((z2 >= g.z[j]) && (z2 <= g.z[j + 1]) && (a != 0)) {
			if ((f.Br10[i][j] + f.Br11[i][j] * z2) != 0)
				r2 = (-f.Br00[i][j] - f.Br01[i][j] * z2) / (f.Br10[i][j] + f.Br11[i][j] * z2);
			else
				r2 = 2 * g.r[0] - g.r[1];
			if ((r2 >= g.r[i]) && (r2 <= g.r[i + 1])) {
				N.id = 1;
				N.r = r2;
				N.z = z2;
				N.i = i;
				N.j = j;
				N.v = f.Psi00[i][j] + f.Psi10[i][j] * g.r[i] + f.Psi01[i][j] * g.z[j] + f.Psi11[i][j] * g.r[i] * g.z[j];
			}
		}
	}

	return N;
}



//	Determine whether the test point(x1,y1: Null point) located out of limiter or not
//	And give Cross point if it is
TESTPOINT CROSSLIMT(SYSALL& sys, double r1, double z1) {
	usealias_geo(g, sys);	//	r1, z1 -> test point	r2,z2 -> Center Point
	usealias_fld(f, sys);
	usealias_mes(m, sys);
	double 	r2 = g.r[((int)(rN / 2.0))], z2 = g.z[((int)(zN / 2.0))];
	TESTPOINT	V;
	int		i;

	DoAllDomain(i, g.lnum - 1) {
		V = CROSSFUNC(r1, z1, r2, z2, g.Lmt_r[i], g.Lmt_z[i], g.Lmt_r[i + 1], g.Lmt_z[i + 1]);
		if (V.id == 1) {
			V.i = i;
			return V;
		}
	}

	return V;
}


//	Determine whether the test line located out of Divertor or not
//	And give Cross point if it exists
TESTPOINT CROSSDIVT(SYSALL& sys, double r1, double z1, double r2, double z2) {
	usealias_geo(g, sys);
	usealias_fld(f, sys);
	usealias_mes(m, sys);
	TESTPOINT	V;
	int			i;

	DoAllDomain(i, g.ndiv) {
		V = CROSSFUNC(r1, z1, r2, z2,
			(2 * g.Divt_r[i][0] - g.Divt_r[i][1]), (2 * g.Divt_z[i][0] - g.Divt_z[i][1]),
			(2 * g.Divt_r[i][1] - g.Divt_r[i][0]), (2 * g.Divt_z[i][1] - g.Divt_z[i][0]));
		if (V.id == 1) 	return V;
	}
	return V;
}


//	Magnetic Field Configurator
//	Find Gradient = 0 points and group thoses whether O / X point
void GRAD0(SYSALL& sys) {
	usealias_geo(g, sys);
	usealias_fld(f, sys);
	usealias_mes(m, sys);
	int 	i, j, ii, jj, iii;
	int 	ntotal = 0;
	int		xpn = 0;
	int		opn = 0;
	int		Condi;
	int		ovlap = 0;
	dmatrix	dVr, dVz;
	double	ddVr, ddVz, ddVrz, S;
	TESTPOINT	N;	//	Null Test point
	TESTPOINT*	V;	//	Storage for Null points
	TESTPOINT*	T;	//	Storage for Cross points

	V = alloc_vector<TESTPOINT>(0, 11);	//	(ASSUMPTION)Maximum Null point number -> 12
	T = alloc_vector<TESTPOINT>(0, 11);

	// Find every Null point
	DoDomain(j, 2, zN - 1) {
		DoDomain(i, 2, rN - 1) {
			Condi = ((f.Bp[i][j] < f.Bp[i - 1][j]) && (f.Bp[i][j] < f.Bp[i][j - 1]) &&
				(f.Bp[i][j] < f.Bp[i + 1][j]) && (f.Bp[i][j] < f.Bp[i][j + 1]))
				||	/* If Center of Bp value is MIN of MAX -> Condi = 1 */
				((f.Bp[i][j] > f.Bp[i - 1][j]) && (f.Bp[i][j] > f.Bp[i][j - 1]) &&
					(f.Bp[i][j] > f.Bp[i + 1][j]) && (f.Bp[i][j] > f.Bp[i][j + 1]));
			if (Condi) {
				DoDomain(ii, i - 1, i + 1) {
					DoDomain(jj, j - 1, j + 1) {
						ovlap = 0;
						N = NULLSEARCH(sys, ii, jj);	//	Find Null point
						if (N.id == 1) {
							T[ntotal] = CROSSLIMT(sys, N.r, N.z);	//	Check the
							if (T[ntotal].id == 1)
								continue;

							if (ntotal > 0) {
								DoAllDomain(iii, ntotal) {
									if (N.i == V[iii].i && N.j == V[iii].j) {
										ovlap = 1;
										break;
									}
								}
							}
							if (ovlap == 1)
								continue;

							V[ntotal] = N;	// Store test point's data
							ntotal++;
						}
					}
				}
			}
		}
	}

	dVr = alloc_matrix<double>(0, rN - 1, 0, zN - 1);
	dVz = alloc_matrix<double>(0, rN - 1, 0, zN - 1);
	f.Xp = alloc_vector<TESTPOINT>(0, ntotal - 1);
	f.Op = alloc_vector<TESTPOINT>(0, ntotal - 1);

	//	Distinguish the Null Points whether O or X
	DERIV(rN, zN, g.r, g.z, g.Psi, dVr, dVz);
	DoDomain(i, 0, ntotal) {
		ddVr = DERIV_R(rN, g.r, g.z, dVr, V[i].i, V[i].j);
		ddVz = DERIV_Z(zN, g.r, g.z, dVz, V[i].i, V[i].j);
		ddVrz = DERIV_R(rN, g.r, g.z, dVz, V[i].i, V[i].j);

		S = ddVr*ddVz - ddVrz*ddVrz;

		if (S < 0) {
			f.Xp[xpn] = V[i];
			f.Xp[xpn].id = xpn + 1;
			xpn++;
		}
		else {// if(S > 0){
			f.Op[opn] = V[i];
			f.Op[opn].id = opn + 1;
			opn++;
			//	Plasma Current -> Counter-Clockwise
			if (ddVr > 0)	f.ro = 1;	//	O-point is Local Minimum value
			else 			f.ro = -1;	//	O-point is Local Maximum value
		}
	}
	f.xpn = xpn;
	f.opn = opn;

	free_matrix<double>(dVr, 0, rN - 1, 0, zN - 1);
	free_matrix<double>(dVz, 0, rN - 1, 0, zN - 1);
	free_vector<TESTPOINT>(T, 0, 11);
	free_vector<TESTPOINT>(V, 0, 11);

	printf("  # of X points: %d\t# of O points: %d\n", f.xpn, f.opn);
	DoAllDomain(i, f.xpn){
		printf("  %d: (%.4f,%.4f)\t", i + 1, f.Xp[i].r, f.Xp[i].z);
		if(i < f.opn)
			printf("%d: (%.4f,%.4f)\n", i + 1, f.Op[i].r, f.Op[i].z);
		else
			printf("\n");
	}
}

/*	MULTIX
	Especially in Snowflake divertor configuration,
	X points could be are located very close unintentionally.
	This function murge those X points to one, So it
	fits(looks like) to the SN or DN cases.

	Assume 2 points are located each positive/negative
	region
*/
void MULTIX(SYSALL& sys) {
	usealias_geo(g, sys);
	usealias_fld(f, sys);
	usealias_mes(m, sys);
	int			i, j;
	int			NCNT, PCNT;
	double		Nval, Pval;
	double		d, Nd, Pd;
	TESTPOINT	NXP[2], PXP[2];	//	Dummy for many xpoint.
	TESTPOINT	N, P;			//	Assume only 2 points exist.

	NCNT = PCNT = 0;
	Nval = Pval = 0;

	d = distance(g.r[0], g.z[0], g.r[1], g.z[1]);

	if (f.xpn > 2) {
		DoAllDomain(i, f.xpn) {
			if (f.Xp[i].z < 0) {
				NXP[NCNT] = f.Xp[i];
				NCNT++;
			}
			else {
				PXP[PCNT] = f.Xp[i];
				PCNT++;
			}
		}

		DoAllDomain(i, NCNT)
			Nval += NXP[i].v / ((double)(NCNT));

		DoAllDomain(i, PCNT)
			Pval += PXP[i].v / ((double)(PCNT));

		N = GETPOINT(sys, NXP[0].r, NXP[0].z, NXP[1].r, NXP[1].z, Nval, 0.);
		P = GETPOINT(sys, PXP[0].r, PXP[0].z, PXP[1].r, PXP[1].z, Pval, 0.);
		Nd = distance(NXP[0].r, NXP[0].z, NXP[1].r, NXP[1].z);
		Pd = distance(PXP[0].r, PXP[0].z, PXP[1].r, PXP[1].z);

		printf("There are more than TWO X points found : %d\n", f.xpn);
		printf("Assuming :::SNOWFLAKE DIVERTOR:::\n");
		printf("%d negative, %d positive X points in z axis\n", NCNT, PCNT);
		printf("Old Negative\t(r : %f, z : %f)\t(r : %f, z : %f)\n",
			NXP[0].r, NXP[0].z, NXP[1].r, NXP[1].z);
		printf("Distance : %f\n", Nd);
		printf("Old Positive\t(r : %f, z : %f)\t(r : %f, z : %f)\n",
			PXP[0].r, PXP[0].z, PXP[1].r, PXP[1].z);
		printf("Distance : %f\n", Pd);

		if (Nd > 5 * d || Pd > 5 * d) {
			Errmsg("The magnetic configuration FAILED", -99);
		}

		printf("New Negative\t(r : %f, z : %f, Psi : %f)\n", N.r, N.z, N.v);
		printf("New Positive\t(r : %f, z : %f, Psi : %f)\n", P.r, P.z, P.v);

		f.Xp[0] = N;
		f.Xp[0].id = 1;
		f.Xp[1] = P;
		f.Xp[1].id = 2;

		f.xpn = 2;

		f.dXp[0].id = 1;
		f.dXp[0].v = Nd;
		f.dXp[0].i = (int)(fabs(NXP[0].r - NXP[1].r) / (g.r[1] - g.r[0]) * 2.5);	//	Grid distance
		f.dXp[0].j = (int)(fabs(NXP[0].z - NXP[1].z) / (g.z[1] - g.z[0]) * 2.);	//	in FIRSTSTEP

		f.dXp[1].id = 2;
		f.dXp[1].v = Pd;
		f.dXp[1].i = (int)(fabs(PXP[0].r - PXP[1].r) / (g.r[1] - g.r[0]) * 2.5);
		f.dXp[1].j = (int)(fabs(PXP[0].z - PXP[1].z) / (g.z[1] - g.z[0]) * 2.);
	}
}

/*	FIRST STEP
Get all segments of separatrix from the certain X-point
Searching Area
(i-2)(j+2) - 11 - (i)(j+2) - 10 - (i+1)(j+2) - 9 - (i+3)(j+2)
	|													|
	0													8
	|													|
(i-2)(j+1)		  (i)(j+1)--------(i+1)(j+1)	   (i+3)(j+1)
	|				|				|					|
	1				|		X		|					7
	|				|				|					|
(i-2)(j)		  (i)(j)----------(i+1)(j)		   (i+3)(j)
	|													|
	2													6
	|													|
(i-2)(j-1) -  3 - (i)(j-1) -  4 - (i+1)(j-1) - 5 - (i+3)(j-1)
*/
void FIRSTSTEP(SYSALL& sys, int iXp, dvector r, dvector z, double psi) {
	usealias_geo(g, sys);
	usealias_fld(f, sys);
	usealias_mes(m, sys);
	TESTPOINT	&Xp = f.Xp[iXp];
	int		i, j, ii, jj;
	int		di[4];
	int		dj[4];
	double	BpMin, CP;	//	Cross Product

	DoAllDomain(j, 12) {
		r[j] = 0;
		z[j] = 0;
	}

	/* Find Max Bp value */
	i = (int)((f.Xp[iXp].r - g.r[0]) / (g.r[1] - g.r[0]));	//	r from X point
	j = (int)((f.Op[0].z - g.z[0]) / (g.z[1] - g.z[0]));	//	z from O point

	if (i < 0)		i = 0;
	if (j < 0)		j = 0;
	if (i > g.rN - 1)	i = g.rN - 1;
	if (j > g.zN - 1)	j = g.zN - 1;

	BpMin = f.Bp[i][j] * 0.03;	//	5% of Max Bp value

	di[0] = -2;	di[1] = 0;	di[2] = 1;	di[3] = 3;
	dj[0] = -1;	dj[1] = 0;	dj[2] = 1;	dj[3] = 2;

	for (i = -1; i <= 1;i += 2) {
		for (j = -1; j <= 1; j += 2) {
			ii = f.Xp[iXp].i + di[(int)(1.5*i + 1.5)];
			jj = f.Xp[iXp].j + dj[(int)(1.5*j + 1.5)];
			while (f.Bp[ii][jj] < BpMin	&& ((ii - f.Xp[iXp].i) < rN * 0.2)
				&& ((jj - f.Xp[iXp].j) < zN * 0.2)) {
				CP = (g.r[ii] - f.Xp[iXp].r) * i - (g.z[jj] - f.Xp[iXp].z) * j;
				if (CP > 0)	jj += j;
				else		ii += i;
			}
			di[(int)(1.5*i + 1.5)] = ii - f.Xp[iXp].i;	//	0 or 3
			di[(int)(0.5*i + 1.5)] = (int)(0.5 * (ii - f.Xp[iXp].i + 0.5));	//	1 or 2
			dj[(int)(1.5*j + 1.5)] = jj - f.Xp[iXp].j;
			dj[(int)(0.5*j + 1.5)] = (int)(0.5 * (jj - f.Xp[iXp].j));
		}
	}

	/*
	di[0] =-2 - (int)(f.dXp[iXp].i / 1.);
	di[1] = 0 - (int)(f.dXp[iXp].i / 1.);
	di[2] = 1 + (int)(f.dXp[iXp].i / 1.);
	di[3] = 3 + (int)(f.dXp[iXp].i / 1.);

	dj[0] =-1 - (int)(f.dXp[iXp].j / 1.);
	dj[1] = 0 - (int)(f.dXp[iXp].j / 2.);
	dj[2] = 1 + (int)(f.dXp[iXp].j / 2.);
	dj[3] = 2 + (int)(f.dXp[iXp].j / 1.);
	*/

	DoAllDomain(j, 3) {
		DoDomain(jj, dj[j], dj[j + 1]) {
			//	0, 1, 2
			if (BTWN(psi, g.Psi[Xp.i + di[0]][Xp.j + jj], g.Psi[Xp.i + di[0]][Xp.j + jj + 1])) {
				r[2 - j] = g.r[Xp.i + di[0]];
				z[2 - j] = LINEAR1D(g.z[Xp.j + jj], g.z[Xp.j + jj + 1],
					g.Psi[Xp.i + di[0]][Xp.j + jj], g.Psi[Xp.i + di[0]][Xp.j + jj + 1], psi);
			}
			//	6, 7, 8
			if (BTWN(psi, g.Psi[Xp.i + di[3]][Xp.j + jj], g.Psi[Xp.i + di[3]][Xp.j + jj + 1])) {
				r[j + 6] = g.r[Xp.i + di[3]];
				z[j + 6] = LINEAR1D(g.z[Xp.j + jj], g.z[Xp.j + jj + 1],
					g.Psi[Xp.i + di[3]][Xp.j + jj], g.Psi[Xp.i + di[3]][Xp.j + jj + 1], psi);
			}
		}
	}
	DoAllDomain(i, 3) {
		DoDomain(ii, di[i], di[i + 1]) {
			//	3, 4, 5
			if (BTWN(psi, g.Psi[Xp.i + ii][Xp.j + dj[0]], g.Psi[Xp.i + ii + 1][Xp.j + dj[0]])) {
				r[i + 3] = LINEAR1D(g.r[Xp.i + ii], g.r[Xp.i + ii + 1],
					g.Psi[Xp.i + ii][Xp.j + dj[0]], g.Psi[Xp.i + ii + 1][Xp.j + dj[0]], psi);
				z[i + 3] = g.z[Xp.j + dj[0]];
			}
			//	9, 10, 11
			if (BTWN(psi, g.Psi[Xp.i + ii][Xp.j + dj[3]], g.Psi[Xp.i + ii + 1][Xp.j + dj[3]])) {
				r[11 - i] = LINEAR1D(g.r[Xp.i + ii], g.r[Xp.i + ii + 1],
					g.Psi[Xp.i + ii][Xp.j + dj[3]], g.Psi[Xp.i + ii + 1][Xp.j + dj[3]], psi);
				z[11 - i] = g.z[Xp.j + dj[3]];
			}
		}
	}
}


/*	EXEMPTION
Exemination overlaped flux of X-point by determining their distances
Least distance
!!!!!At least one point have to be exempted
*/
void EXEMPTION(SYSALL& sys, int iXp, double r, double z,
	dvector newr, dvector newz, double dmin) {
	usealias_geo(g, sys);
	usealias_fld(f, sys);
	usealias_mes(m, sys);
	int		i, j = 0, ii;
	int		i_min = 0;
	ivector	di;
	dvector	d;
	double	A, Ar, Az, B, Br, Bz, C;
	double	a, b, c;

	di = alloc_vector<int>(0, 11);
	d = alloc_vector<double>(0, 11);

	DoAllDomain(i, 12)
		di[i] = d[i] = 100;

	a = (f.Xp[iXp].z - z) / (f.Xp[iXp].r - r);
	b = -1.;
	c = z - a*r;

	Ar = r - f.Xp[iXp].r;
	Az = z - f.Xp[iXp].z;
	A = sqrt(Ar*Ar + Az*Az);
	DoAllDomain(i, 12) {
		if (newr[i] != 0) {
			Br = newr[i] - f.Xp[iXp].r;
			Bz = newz[i] - f.Xp[iXp].z;
			B = sqrt(Br*Br + Bz*Bz);
			C = (Ar * Br + Az * Bz) / A / B;
			if (C > 0) {	//	Need to be more rigid ( One more Condition )
				//d[j] = distance(r, z, newr[i], newz[i]);	//	Should Check : pointer <-> array
				d[j] = fabs(a*newr[i] + b*newz[i] + c) / sqrt(a*a + b*b);
				di[j] = i;
				j++;
			}
		}
	}
	DoAllDomain(i, j) {
		if (d[i_min] > d[i]) {
			i_min = i;
		}
	}
	if (d[i_min] < dmin) {
		newr[di[i_min]] = 0;
		newz[di[i_min]] = 0;
	}

	free_vector<int>(di, 0, 11);
	free_vector<double>(d, 0, 11);
}


/*	GETFIELD
*/
void GETFIELD(SYSALL& sys, double r, double z, dvector dr, dvector dz, double d) {
	usealias_geo(g, sys);
	usealias_fld(f, sys);
	usealias_mes(m, sys);
	int 	i, j;
	double 	Br, Bz, Bp;
	double	ratio = 1.0;

	i = (int)((r - g.r[0]) / (g.r[1] - g.r[0]));	//	Fine index of r, z
	j = (int)((z - g.z[0]) / (g.z[1] - g.z[0]));	//	Assume Equi-distribution

	if (i < 0)		i = 0;
	if (j < 0)		j = 0;
	if (i > g.rN - 1)	i = g.rN - 1;
	if (j > g.zN - 1)	j = g.zN - 1;

	Br = f.Br00[i][j] + f.Br10[i][j] * r + f.Br01[i][j] * z + f.Br11[i][j] * r*z;
	Bz = f.Bz00[i][j] + f.Bz10[i][j] * r + f.Bz01[i][j] * z + f.Bz11[i][j] * r*z;
	Bp = sqrt(Br*Br + Bz*Bz);

	*dr = Br / Bp*d * ratio;
	*dz = Bz / Bp*d * ratio;
}


/*	VECTORFOLLOW
Drow line by vector following method -> RUNGE-KUTTA 4th
*/
void VECTORFOLLOW(SYSALL& sys, double psi, double& r, double& z,
	double r0, double z0, double d, int drct) {
	usealias_geo(g, sys);
	usealias_fld(f, sys);
	usealias_mes(m, sys);
	int 	i, j, k = 0;
	double 	dr, dz;
	double 	k1r, k1z, k2r, k2z, k3r, k3z, k4r, k4z;
	double	p1, p2, p3, p4;
	double	w = 0.5;
	double	dr0, dz0;
	double	drctr, drctz;

	if (r0 == r)
		drctr = drct;
	else if (z0 == z)
		drctz = drct;
	else {
		drctr = (r - r0) / fabs(r - r0);
		drctz = (z - z0) / fabs(z - z0);
	}

	while (1) {
		//	Runge-Kutta 4th method  ----->>>> WIERD!
		GETFIELD(sys, r, z, &k1r, &k1z, d);
		GETFIELD(sys, r + k1r / 2, z + k1z / 2, &k2r, &k2z, d);
		GETFIELD(sys, r + k2r / 2, z + k2z / 2, &k3r, &k3z, d);
		GETFIELD(sys, r + k3r, z + k3z, &k4r, &k4z, d);

		dr = fabs(1.0 / 6.0 * (k1r + 2 * k2r + 2 * k3r + k4r)) * drctr;		//	One step further
		dz = fabs(1.0 / 6.0 * (k1z + 2 * k2z + 2 * k3z + k4z)) * drctz;

		i = (int)((r + dr - g.r[0]) / (g.r[1] - g.r[0]));	//	Fine index of r, z
		j = (int)((z + dz - g.z[0]) / (g.z[1] - g.z[0]));	//	Assume Equi-distribution

		if (i < 0)		i = 0;
		if (j < 0)		j = 0;
		if (i > g.rN - 1)	i = g.rN - 1;
		if (j > g.zN - 1)	j = g.zN - 1;

		p1 = f.Psi00[i][j] + f.Psi10[i][j] * g.r[i] + f.Psi01[i][j] * (z + dz) + f.Psi11[i][j] * g.r[i] * (z + dz);
		p2 = f.Psi00[i][j] + f.Psi10[i][j] * g.r[i + 1] + f.Psi01[i][j] * (z + dz) + f.Psi11[i][j] * g.r[i + 1] * (z + dz);
		p3 = f.Psi00[i][j] + f.Psi10[i][j] * (r + dr) + f.Psi01[i][j] * g.z[j] + f.Psi11[i][j] * (r + dr) * g.z[j];
		p4 = f.Psi00[i][j] + f.Psi10[i][j] * (r + dr) + f.Psi01[i][j] * g.z[j + 1] + f.Psi11[i][j] * (r + dr) * g.z[j + 1];

		if (BTWN(psi, p1, p2)) {
			r = LINEAR1D(g.r[i], g.r[i + 1], p1, p2, psi);
			z = z + dz;
			break;
		}
		else if (BTWN(psi, p3, p4)) {
			r = r + dr;
			z = LINEAR1D(g.z[j], g.z[j + 1], p3, p4, psi);
			break;
		}
		else if (k > 3) {
			printf("Rough estimation at  i = %d, j = %d (%.4f, %.4f)\n", i, j, r, z);

			dr0 = (r - r0) / (distance(r, z, r0, z0) ) * d;
			dz0 = (z - z0) / (distance(r, z, r0, z0)) * d;
			r += (w * dr + (1 - w) * (-dr0)) * 1.0;
			z += (w * dz + (1 - w) * dz0) * 1.0;
			break;
		}
		else {
			d *= 0.7;
			k++;
		}
	}
}

/*	FRAME
Set frame of entire mesh through lines from X-points
*/
void FRAME(SYSALL& sys) {
	usealias_geo(g, sys);
	usealias_fld(f, sys);
	usealias_mes(m, sys);
	const int limt = KMAX;
	const int prestep = 10;
	int		i, j, k = 0, l, nl = 0;	//	k: index of overlap point, nl: # of whole flux line
	int		iXp, jXp;	//	iXp: index of Xpoint
	int 	brk = 0;			// 	break condition
	int		drct = 1;			//	direction indicator
	int		LGi;
	double	newr[12] = { 0 };		//	first segment of separatrix
	double	newz[12] = { 0 };
	double	a, b, dummy1 = 0, dummy2;
	dmatrix Dumr, Dumz;		//	Dummy r, z matrix
	double	r, z, r0, z0, testd;	//	testd: Distance from first segment to stretched point
	double 	d, dmin;		//	dmin: Distance of cell diagonal
	TESTPOINT	DumV;	//	Dummy Point
	TESTPOINT	*PreV;	//	Save Point to prevent repetition

	d = distance(g.r[0], g.z[0], g.r[1], g.z[1]) * 1.0;	//	d = radial differece
	Dumr = alloc_matrix<double>(0, 12 * f.xpn - 1, 0, limt);	//	Allocate 12x300 matrix
	Dumz = alloc_matrix<double>(0, 12 * f.xpn - 1, 0, limt);	//	to each X-point
	PreV = alloc_vector<TESTPOINT>(0, 12 * f.xpn - 1);
	m.MaxGr = alloc_matrix<double>(0, 11, 0, f.xpn - 1);
	m.MaxGz = alloc_matrix<double>(0, 11, 0, f.xpn - 1);

	DoAllDomain(i, 12 * f.xpn) {
		DoAllDomain(j, limt) {
			Dumr[i][j] = Dumz[i][j] = 0;
		}
	}
	DoAllDomain(i, 12) {
		DoAllDomain(j, f.xpn) {
			m.MaxGr[i][j] = m.MaxGz[i][j] = 0;
		}
	}
	DoAllDomain(iXp, f.xpn) {
		FIRSTSTEP(sys, iXp, newr, newz, f.Xp[iXp].v);		//	Get first-vector points

		dmin = distance(g.r[0], g.z[0], g.r[2 + (int)(f.dXp[iXp].i / 2.)],
			g.z[1 + (int)(f.dXp[iXp].j / 2.)]) * 1.0;	//	Assume equi-distribution

		if (k > 0) {						//	To Escape overlaped point
			DoAllDomain(l, k) {
				if (iXp == PreV[l].id) {
					EXEMPTION(sys, iXp, PreV[l].r, PreV[l].z, newr, newz, dmin);
				}					//	Overlap Flux location -> 0
			}
		}
		DoAllDomain(i, 12) {
			if (newr[i] != 0) {
				//	Considering plasma rotation direction
				if ((f.Xp[iXp].r > newr[i] && f.Xp[iXp].z > newz[i])
					|| (f.Xp[iXp].r < newr[i] && f.Xp[iXp].z < newz[i])) {
					drct = f.ro*(1);
				}
				else {
					drct = f.ro*(-1);
				}
				brk = 0;
				j = 0;
				nl++;					//	# of Separatrix line
				Dumr[12 * iXp + i][j] = r = newr[i];
				Dumz[12 * iXp + i][j] = z = newz[i];

				DoDomain(j, 1, 3) {
					if (j == 1) {
						r0 = f.Xp[iXp].r;
						z0 = f.Xp[iXp].z;
					}
					else {
						r0 = Dumr[12 * iXp + i][j - 2];
						z0 = Dumz[12 * iXp + i][j - 2];
					}
					VECTORFOLLOW(sys, f.Xp[iXp].v, r, z, r0, z0, d, drct);
					Dumr[12 * iXp + i][j] = r;
					Dumz[12 * iXp + i][j] = z;
				}
				while (1) {
					r0 = Dumr[12 * iXp + i][j - 2];
					z0 = Dumz[12 * iXp + i][j - 2];
					VECTORFOLLOW(sys, f.Xp[iXp].v, r, z, r0, z0, d, drct);
					DoDomain(jXp, iXp, f.xpn) {	//	Checking loop
						testd = distance(r, z, f.Xp[jXp].r, f.Xp[jXp].z);
						if (testd < dmin) {		//	Could not be Correct
							Dumr[12 * iXp + i][j] = f.Xp[jXp].r;
							Dumz[12 * iXp + i][j] = f.Xp[jXp].z;
							PreV[k].id = jXp;			//	Which X-point
							PreV[k].r = Dumr[12 * iXp + i][j - 1];		//	And its location
							PreV[k].z = Dumz[12 * iXp + i][j - 1];
							if (jXp == iXp) {
								EXEMPTION(sys, iXp, PreV[k].r, PreV[k].z, newr, newz, dmin);
								nl++;
								brk = 1;
								break;
							}
							k++;
							brk = 1;			//	to escape while loop
							break;
						}
					}
					if (brk == 1)
						break;		//	End of flux

					DumV = CROSSDIVT(sys, r, z, r0, z0);
					if (DumV.id == 1) {
						Dumr[12 * iXp + i][j] = DumV.r;	//	point at limiter
						Dumz[12 * iXp + i][j] = DumV.z;
						DumV.id = 0;
						break;
					}

					if ((r < g.r[0]) || (r > g.r[g.rN - 1]) || (z < g.z[0]) || (z > g.z[g.zN - 1]))
						break;

					Dumr[12 * iXp + i][j] = r;	//	Save
					Dumz[12 * iXp + i][j] = z;
					j++;
					if (j > limt - 1) {
						//Errmsg("Not Enough Matrix in FRAME");
						break;
					}
				}
			}
		}
	}
	m.Sr = alloc_matrix<double>(0, nl - 1, 0, limt);
	m.Sz = alloc_matrix<double>(0, nl - 1, 0, limt);
	DoAllDomain(i, nl)
		DoAllDomain(j, limt)
			m.Sr[i][j] = m.Sz[i][j] = 0;

	nl = 0;

	//	Allocate Core Separatrix First
	DoAllDomain(i, 12 * f.xpn) {
		if (Dumr[i][0] != 0) {
			DoDomain(j, 1, limt) {
				if (Dumr[i][j] == 0)
					break;
				DoAllDomain(k, f.xpn) {
					if ((Dumr[i][j] == f.Xp[k].r) && (Dumz[i][j] == f.Xp[k].z)) {
						if ((Dumr[i][j] == f.Xp[((int)(i / 12))].r) && (Dumz[i][j] == f.Xp[((int)(i / 12))].z)) {
							if (f.xpn < 2) {	//	SN
								r0 = Dumr[i][j];	//	Dummy X point for SN case
								z0 = -Dumz[i][j];
							}
							else {				//	DDN	(Only the case of Max(xpn) == 2)
								r0 = f.Xp[1 - k].r;
								z0 = f.Xp[1 - k].z;
							}
							m.Sr[nl][0] = m.Sr[nl + 1][0] = f.Xp[((int)(i / 12))].r;
							m.Sz[nl][0] = m.Sz[nl + 1][0] = f.Xp[((int)(i / 12))].z;
							DoDomain(l, 1, j + 1) {
								DumV = CROSSFUNC(Dumr[i][l - 1], Dumz[i][l - 1], Dumr[i][l], Dumz[i][l],
									f.Op[0].r, f.Op[0].z, r0, z0);
								if (DumV.id == 1) {
									DumV.j = l;
									break;
								}
							}
							DoDomain(l, 1, DumV.j) {
								m.Sr[nl][l] = Dumr[i][l - 1];
								m.Sz[nl][l] = Dumz[i][l - 1];
							}
							m.Sr[nl][l] = DumV.r;
							m.Sz[nl][l] = DumV.z;
							DoDomain(l, 1, j + 1 - DumV.j) {
								m.Sr[nl + 1][l] = Dumr[i][j - l];
								m.Sz[nl + 1][l] = Dumz[i][j - l];
							}
							m.Sr[nl + 1][l] = DumV.r;
							m.Sz[nl + 1][l] = DumV.z;
							Dumr[i][0] = 0;
							Dumz[i][0] = 0;
							nl = nl + 2;
						}
						else {
							m.Sr[nl][0] = f.Xp[((int)(i / 12))].r;
							m.Sz[nl][0] = f.Xp[((int)(i / 12))].z;
							DoDomain(l, 1, j + 2) {
								m.Sr[nl][l] = Dumr[i][l - 1];	//	0 ~ j
								m.Sz[nl][l] = Dumz[i][l - 1];
							}
							Dumr[i][0] = 0;
							Dumz[i][0] = 0;
							nl++;
						}
					}
				}
			}
		}
	}

	//	Allocate DDN separatrix
	DoAllDomain(i, 12 * f.xpn) {
		if (Dumr[i][0] != 0) {
			DoDomain(j, 1, limt) {
				if (Dumr[i][j] == 0)
					break;
			}
			if ((Dumz[i][j - 1] - f.Op[0].z) * (Dumz[i][0] - f.Op[0].z) < 0) {
				m.Sr[nl][0] = f.Xp[((int)(i / 12.0))].r;
				m.Sz[nl][0] = f.Xp[((int)(i / 12.0))].z;
				DoDomain(j, 1, limt) {
					m.Sr[nl][j] = Dumr[i][j - 1];
					m.Sz[nl][j] = Dumz[i][j - 1];
					if (m.Sr[nl][j] != m.Sr[nl][j])	m.Sr[nl][j] = 0;	//	erase NAN
					if (m.Sz[nl][j] != m.Sz[nl][j])	m.Sz[nl][j] = 0;
				}
				Dumr[i][0] = 0;
				Dumz[i][0] = 0;
				nl++;
			}
		}
	}

	//	Other separatrix
	DoAllDomain(i, 12 * f.xpn) {
		if (Dumr[i][0] != 0) {
			m.Sr[nl][0] = f.Xp[((int)(i / 12.0))].r;
			m.Sz[nl][0] = f.Xp[((int)(i / 12.0))].z;
			DoDomain(j, 1, limt) {
				m.Sr[nl][j] = Dumr[i][j - 1];
				m.Sz[nl][j] = Dumz[i][j - 1];
				if (m.Sr[nl][j] != m.Sr[nl][j])	m.Sr[nl][j] = 0;	//	erase NAN
				if (m.Sz[nl][j] != m.Sz[nl][j])	m.Sz[nl][j] = 0;
			}
			nl++;
		}
	}
	m.Sn = nl;

	free_matrix<double>(Dumr, 0, 12 * f.xpn - 1, 0, limt);
	free_matrix<double>(Dumz, 0, 12 * f.xpn - 1, 0, limt);
	free_vector<TESTPOINT>(PreV, 0, 12 * f.xpn - 1);

	printf("# of segments: %d\n", m.Sn);
}
