/***************************************************/
/*												   */
/*		MESHGEN FUNCTION						   */
/*												   */
/*		JG LEE									   */
/*		14/JAN/2016								   */
/***************************************************/

#define	_CRT_SECURE_NO_WARNINGS
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "math.h"

#include "parameter.h"
#include "svc.h"
#include "vega.h"
#include "memalloc.h"
#include "separatrix.h"
#include "meshgenerator.h"

using namespace std;

/*	MESHPAR
Calculate parameters of mesh
*/
void MESHPAR(SYSALL& sys) {
	usealias_geo(g, sys);
	usealias_fld(f, sys);
	usealias_mes(m, sys);
	int i;

	/*	Number allocation			*/
	/*m.nD[Co]	= 61;
	m.nD[IBdy]	= 41;
	m.nD[OBdy]	= 41;
	m.nD[IDt]	= 10;	//	Even number
	m.nD[ODt]	= 20;
	m.nD[AS]	= 0;
	m.nD[SOL]	= 10;
	m.nD[Pr]	= 10;
	m.nD[APr]	= 0;*/	//	Move into read_par
	
	printf("RE-processing Mesh parameters\n");
	/*	Type Define					*/
	switch (f.xpn)
	{
	case 0:
		m.Type = f.Op[0];
		m.Type.id = LMT;				//	3 -> LMT
		break;
	case 1:
		m.Type = f.Xp[0];
		m.Type.id = SN;					//	0 -> SN
		break;
	case 2:
		//	Disconnected Double Null Case
		if (fabs(f.Xp[0].v - f.Xp[1].v) > (fabs(f.Op[0].v*0.9 - f.Xp[0].v) / m.nD[Co])) {
			if (max(fabs(f.Xp[0].v - f.Op[0].v), fabs(f.Xp[1].v - f.Op[0].v)) == fabs(f.Xp[1].v - f.Op[0].v)) {
				m.Type = f.Xp[1];	//	save a distant X point
				m.Type.id = DDN;
			}
			else {
				m.Type = f.Xp[0];
				m.Type.id = DDN;		//	2 -> DDN
			}
		}
		else {
			m.Type = f.Xp[0];
			m.Type.id = CDN;			//	1 -> CDN
		}
		break;
	default:
		Errmsg("There are more than Two X-points", -3);
		break;
	}

	/*	Determination of Mesh Parameter	*/
	if (m.Type.id == DDN) {

		m.Psi[APr] = (f.Xp[0].v + f.Xp[1].v - m.Type.v)	+ 
			(f.Op[0].v - (f.Xp[0].v + f.Xp[1].v - m.Type.v)) * PRIP;
		m.nD[AS] = (int)fabs((m.nD[Pr] * (m.Type.v - (f.Xp[0].v + f.Xp[1].v - m.Type.v)) 
			/ (m.Type.v - m.Psi[APr]))) + 2;
		m.nD[Co] -= (m.nD[AS] - 2);
		m.nD[APr] = m.nD[Pr] - m.nD[AS] + 2;
		m.dP[Co] = (f.Op[0].v - (f.Xp[0].v + f.Xp[1].v - m.Type.v)) * 
			PENET / ((double)(m.nD[Co]));
		m.dP[AS] = -((f.Xp[0].v + f.Xp[1].v - m.Type.v) - m.Type.v) / ((double)(m.nD[AS] - 1));
		m.dP[APr] = (m.Psi[APr] - (f.Xp[0].v + f.Xp[1].v - m.Type.v)) /
			((double)(m.nD[APr])) ;
	}
	else {
		m.dP[Co] = (f.Op[0].v - m.Type.v) * PENET / ((double)(m.nD[Co]));	// Exclude separatrix
	}
	m.Psi[SOL] = m.Type.v - (f.Op[0].v - m.Type.v) * SOLP;
	m.Psi[Pr] = m.Type.v + (f.Op[0].v - m.Type.v) * PRIP;	
	m.dP[SOL] = (m.Psi[SOL] - m.Type.v) / ((double)(m.nD[SOL] - 1));	//	Include X point
	m.dP[Pr] = (m.Psi[Pr] - m.Type.v) / ((double)(m.nD[Pr]));		//	Exclude X point

	printf("  Configuration Type: %d\t", m.Type.id);
	printf("  (0: SN, 1: CDN, 2: DDN, 3: Limiter Plasma)\n");
	printf("  # distribution of grid NOW:\t");
	printf("(Co, IB, OB, ID, OD,SOL, AS, PR,APR)\n\t\t\t\t");
	DoAllDomain(i, MGD)	printf(" %d ", m.nD[i]);
	printf("\n");
}

/*	NORVECTRC
	Normal Vector Tracing method
*/
TESTPOINT NORVECTRC(SYSALL& sys, double r, double z, double fr, double fz,
	int LR, double psi, double w) {
	usealias_geo(g, sys);
	usealias_fld(f, sys);
	usealias_mes(m, sys);
	int 	i;
	double 	dr0, dz0, dr, dz;
	double 	k1r, k1z, k2r, k2z, k3r, k3z, k4r, k4z;
	double	d;
	TESTPOINT V;

	V.id = 0;	V.r = 0;	V.z = 0;
	V.v = 0;	V.i = 0;	V.j = 0;

	//	Get Tangential Vector
	d = distance(g.r[0], g.z[0], g.r[1], g.z[1]);

	i = 0;
	while (i < 10) {
		GETFIELD(sys, r, z, &k1r, &k1z, d);
		GETFIELD(sys, r + k1r / 2, z + k1z / 2, &k2r, &k2z, d);
		GETFIELD(sys, r + k2r / 2, z + k2z / 2, &k3r, &k3z, d);
		GETFIELD(sys, r + k3r, z + k3z, &k4r, &k4z, d);

		dr0 = 1.0 / 6.0 * (k1r + 2 * k2r + 2 * k3r + k4r);
		dz0 = 1.0 / 6.0 * (k1z + 2 * k2z + 2 * k3z + k4z);

		dr = -(1 - w) * dz0 * LR * f.ro + w * fr;		//	LR -> -1: SOL(left?), 1: CORE, Private(right?)
		dz =  (1 - w) * dr0 * LR * f.ro + w * fz;
	
	//DoAllDomain(i, f.xpn) {
	//	if (r == f.Xp[i].r && z == f.Xp[i].z) {
	//		GETFIELD(sys, r2, z2, &k1r, &k1z, d);
	//		GETFIELD(sys, r2 + k1r / 2, z2 + k1z / 2, &k2r, &k2z, d);
	//		GETFIELD(sys, r2 + k2r / 2, z2 + k2z / 2, &k3r, &k3z, d);
	//		GETFIELD(sys, r2 + k3r, z2 + k3z, &k4r, &k4z, d);

	//		dr0 = 1.0 / 6.0 * (k1r + 2 * k2r + 2 * k3r + k4r);
	//		dz0 = 1.0 / 6.0 * (k1z + 2 * k2z + 2 * k3z + k4z);

	//		if (LR == -1) {		//	for SOL
	//			drct = (f.Xp[i].r - f.Op[0].r) * dr0 + (f.Xp[i].z - f.Op[0].z) * dz0;	//	Inner product
	//			drct = drct / fabs(drct);
	//			dr = -(f.Xp[i].z - f.Op[0].z) / dxo * d * drct;
	//			dz =  (f.Xp[i].r - f.Op[0].r) / dxo * d * drct;
	//		}
	//		else {		//	for CORE, Private
	//			drct = (f.Xp[i].r - f.Op[0].r) * (dz0)+(f.Xp[i].z - f.Op[0].z) * (-dr0);	//	Inner product
	//			drct = drct / fabs(drct);
	//			dr = (f.Xp[i].r - f.Op[0].r) / dxo * d * drct;
	//			dz = (f.Xp[i].z - f.Op[0].z) / dxo * d * drct;
	//		}
	//		break;
	//	}
	//}	

	
		/*V = CROSSDIVT(sys, r, z, r + dr, z + dz);
		if (V.id == 1) {
			return V;
		}*/
		V = GETPOINT(sys, r, z, r + dr, z + dz, psi, 0.2);
		if (V.id == 1) {
			return V;
		}
		r += dr;
		z += dz;
		i++;
	}
	return V;
}



void STRETCHPAR(dvector dist, double L, int n, double E, double D, double S){
	int		i;
	double	phi;
	dvector s, ds;
	s = alloc_vector<double>(0, n - 1);
	ds = alloc_vector<double>(0, n - 1);

	DoAllDomain(i, n)
		s[i] = 0;

	DoAllDomain(i, n) {
		phi = (((double)(i)-1.0) / ((double)(n)) - D) * S;
		ds[i] = pow((2.0 / (exp(phi) + exp(-phi))), E);
		if (i == 0)
			s[i] = ds[i];
		else
			s[i] = s[i - 1] + ds[i];
	}
	DoAllDomain(i, n)
		dist[i] = s[i] * L / s[n - 1];

	free_vector<double>(s, 0, n - 1);
	free_vector<double>(ds, 0, n - 1);
}


/*	STRETCHING
	Stretch out the separatrix
*/
void STRETCHING(SYSALL & sys){
	usealias_geo(g, sys);
	usealias_fld(f, sys);
	usealias_mes(m, sys);
	int i, j, k;
	dvector l, dl, dist;
	dmatrix Dumr, Dumz;
	int n, nl, dst = 0;
	double d, z;
	FILE *out;

	//	Print out
	out = fopen("Sr2.txt", "w");
	DoAllDomain(j, NMAX) {
		DoAllDomain(i, m.Sn) {
			fprintf(out, "%lf\t", m.Sr[i][j]);
		}
		fprintf(out, "\n");
	}
	fclose(out);
	out = fopen("Sz2.txt", "w");
	DoAllDomain(j, NMAX) {
		DoAllDomain(i, m.Sn) {
			fprintf(out, "%lf\t", m.Sz[i][j]);
		}
		fprintf(out, "\n");
	}
	fclose(out);

	printf("Distribution Separatrix points\t");

	Dumr = alloc_matrix<double>(0, m.Sn - 1, 0, NMAX - 1);
	Dumz = alloc_matrix<double>(0, m.Sn - 1, 0, NMAX - 1);
	DoAllDomain(i, m.Sn) {
		DoAllDomain(j, NMAX) {
			Dumr[i][j] = Dumz[i][j] = 0;
		}
	}

	d = distance(m.Type.r, m.Type.z, f.Op[0].r, f.Op[0].z);
	
	l = alloc_vector<double>(0, KMAX - 1);
	dl = alloc_vector<double>(0, KMAX - 1);

	DoAllDomain(k, m.Sn){
		DoAllDomain(i, KMAX)
			l[i] = dl[i] = 0;

		l[0] = dl[0] = distance(m.Sr[k][0], m.Sz[k][0], m.Sr[k][1], m.Sz[k][1]);
		DoDomain(i, 1, KMAX - 1) {
			if (m.Sr[k][i + 1] == 0)
				break;
			dl[i] = distance(m.Sr[k][i], m.Sz[k][i], m.Sr[k][i + 1], m.Sz[k][i + 1]);
			l[i] = l[i - 1] + dl[i];
		}
		nl = i;
		
		z = inorout(f.Op[0].r, f.Op[0].z, m.Sr[k][0], m.Sz[k][0],
			m.Sr[k][2], m.Sz[k][2]);

		if (k < 2) {	//	CORE
			if (z > 0.0){
				n = m.nD[IBdy] - 1;
				dst = IBdy;
			}
			else{
				n = m.nD[OBdy] - 1;
				dst = OBdy;
			}
		}
		else if (m.Type.id == DDN && k < 4) {	//	
			if (z > 0){
				n = m.nD[IBdy] + m.nD[IDt] - 1;
				dst = IBdy;
			}
			else{
				n = m.nD[OBdy] + m.nD[ODt] - 1;
				dst = OBdy;
			}
		}
		else {							//	Dt
			if (z > 0){
				n = m.nD[IDt] - 1;
				dst = IDt;
			}
			else{
				n = m.nD[ODt] - 1;
				dst = ODt;
			}
		}
		dist = alloc_vector<double>(0, n - 1);
		STRETCHPAR(dist, l[nl - 1], n, m.EDS[dst][0], m.EDS[dst][1], m.EDS[dst][2]);

		DoAllDomain(i, n - 1) {
			DoAllDomain(j, nl) {
				if (l[j] > dist[i]) {
					Dumr[k][i + 1] = m.Sr[k][j + 1] - (l[j] - dist[i]) / dl[j] * (m.Sr[k][j + 1] - m.Sr[k][j]);
					Dumz[k][i + 1] = m.Sz[k][j + 1] - (l[j] - dist[i]) / dl[j] * (m.Sz[k][j + 1] - m.Sz[k][j]);
					break;
				}
			}
		}
		Dumr[k][0] = m.Sr[k][0];
		Dumz[k][0] = m.Sz[k][0];
		Dumr[k][n] = m.Sr[k][nl];
		Dumz[k][n] = m.Sz[k][nl];
		free_vector<double>(dist, 0, n - 1);
	}

	m.Sr = alloc_matrix<double>(0, m.Sn - 1, 0, NMAX - 1);
	m.Sz = alloc_matrix<double>(0, m.Sn - 1, 0, NMAX - 1);
	DoAllDomain(i, m.Sn) {
		DoAllDomain(j, NMAX) {
			m.Sr[i][j] = m.Sz[i][j] = 0;
		}
	}
	DoAllDomain(i, m.Sn) {
		DoAllDomain(j, NMAX) {
			if (Dumr[i][j] == 0)
				break;
			m.Sr[i][j] = Dumr[i][j];
			m.Sz[i][j] = Dumz[i][j];
		}
	}
	free_matrix<double>(Dumr, 0, m.Sn - 1, 0, NMAX - 1);
	free_matrix<double>(Dumz, 0, m.Sn - 1, 0, NMAX - 1);
	free_vector<double>(l, 0, KMAX - 1);
	free_vector<double>(dl, 0, KMAX - 1);

	printf("Done!\n");

	//	Print out
	out = fopen("Sr.txt", "w");
	DoAllDomain(j, NMAX) {
		DoAllDomain(i, m.Sn) {
			fprintf(out, "%lf\t", m.Sr[i][j]);
		}
		fprintf(out, "\n");
	}
	fclose(out);
	out = fopen("Sz.txt", "w");
	DoAllDomain(j, NMAX) {
		DoAllDomain(i, m.Sn) {
			fprintf(out, "%lf\t", m.Sz[i][j]);
		}
		fprintf(out, "\n");
	}
	fclose(out);
}


/*	MESHGEN
	NOT USED --> Backed up
*/
void MESHGEN(SYSALL& sys) {
	
}



/*  Boundary Vector  */
TESTPOINT BVEC(SYSALL &sys, int k, int LR, int YCo) {
	usealias_geo(g, sys);
	usealias_fld(f, sys);
	usealias_mes(m, sys);
	TESTPOINT	B;
	TESTPOINT	L[3];
	double		k1r, k1z, k2r, k2z, k3r, k3z, k4r, k4z;
	double		r, z, dr0, dz0, dr1, dz1, d;
	double		Nr, Nz;
	double		MDr = 0, MDz = 0;
	double		dummy, Bsize;
	double		A, Ar, Az, C, Cr, Cz;
	int			i, j, iXp;

	B.id = 0;

	r = m.Mr[HALF - 1][k];
	z = m.Mz[HALF - 1][k];

	if(YCo == 1){
		DoAllDomain(i, f.xpn) {
			if ((m.Mr[HALF - 1][k] == f.Xp[i].r) && (m.Mz[HALF - 1][k] == f.Xp[i].z)) {
				iXp = i;
				r = m.Mr[HALF - 1][k + 2];
				z = m.Mz[HALF - 1][k + 2];
				DoAllDomain(j, f.xpn) {
					if ((m.Mr[HALF - 1][k + 1] == f.Xp[j].r) && (m.Mz[HALF - 1][k + 1] == f.Xp[j].z)) {
						//	"==" <- Dangerous
						r = m.Mr[HALF - 1][k - 2];
						z = m.Mz[HALF - 1][k - 2];
						break;
					}
				}
			}
		}
		dr1 = r - f.Xp[iXp].r;
		dz1 = z - f.Xp[iXp].z;
	}
	d = distance(g.r[0], g.z[0], g.r[1], g.z[1]);

	GETFIELD(sys, r, z, &k1r, &k1z, d);
	GETFIELD(sys, r + k1r / 2, z + k1z / 2, &k2r, &k2z, d);
	GETFIELD(sys, r + k2r / 2, z + k2z / 2, &k3r, &k3z, d);
	GETFIELD(sys, r + k3r, z + k3z, &k4r, &k4z, d);

	dr0 = 1.0 / 6.0 * (k1r + 2 * k2r + 2 * k3r + k4r);
	dz0 = 1.0 / 6.0 * (k1z + 2 * k2z + 2 * k3z + k4z);

	if (YCo == 1) {		//	For Near X point
		if (dr0*dr1 + dz0*dz1 > 0) {
			Nr = -dz1 * LR * f.ro;
			Nz =  dr1 * LR * f.ro;
		}
		else {
			Nr =  dz1 * LR * f.ro;
			Nz = -dr1 * LR * f.ro;
		}
		Ar = r - f.Xp[iXp].r;
		Az = z - f.Xp[iXp].z;
		A = sqrt(Ar*Ar + Az*Az);
		if (m.MaxGr[0][iXp] != 0) {
			Cr = m.MaxGr[0][iXp] - f.Xp[iXp].r;
			Cz = m.MaxGz[0][iXp] - f.Xp[iXp].z;
			C = sqrt(Cr*Cr + Cz*Cz);
			L[0].v = (-1) * (Ar*Cr + Az*Cz) / A / C;
		}
		L[0].id = 0;
		L[0].r = m.MaxGr[0][iXp];
		L[0].z = m.MaxGz[0][iXp];

		DoDomain(i, 1, 12) {		//	Two least distance boundary point
			if (L[0].v < L[1].v) {
				L[2] = L[1];	//	L[2] : dummy
				L[1] = L[0];
				L[0] = L[2];
			}
			if (m.MaxGr[i][iXp] != 0) {
				Cr = m.MaxGr[i][iXp] - f.Xp[iXp].r;
				Cz = m.MaxGz[i][iXp] - f.Xp[iXp].z;
				C = sqrt(Cr*Cr + Cz*Cz);
				dummy = (-1) * (Ar*Cr + Az*Cz) / A / C;
				if (dummy <= L[0].v) {
					L[0].v = dummy;
					L[0].id = i;
					L[0].r = m.MaxGr[i][iXp];
					L[0].z = m.MaxGz[i][iXp];
				}
			}
		}
		DoAllDomain(i, 2) {			//	Choose One
			dummy = Nr * (L[i].r - f.Xp[iXp].r) + Nz * (L[i].z - f.Xp[iXp].z);
			if (dummy > 0) {
				B.r		= L[i].r - f.Xp[iXp].r;
				B.z		= L[i].z - f.Xp[iXp].z;
				Bsize	= sqrt(B.r * B.r + B.z * B.z);
				B.r		= B.r / Bsize * d;		//	Move it to outward
				B.z		= B.z / Bsize * d;
				B.id	= L[i].id + 1;
				return B;
			}
		}
	}
	else if (YCo == 2) {		// For SN, DDN core
		Nr = -dz0 * LR * f.ro;
		Nz = dr0 * LR * f.ro;
		
		if ((f.Op[0].r - r)*Nr + (f.Op[0].z - z)*Nz > 0) {
			B.r = f.Op[0].r - r;
			B.z = f.Op[0].z - z;
		}
		else {
			B.r = r - f.Op[0].r;
			B.z = z - f.Op[0].z;
		}
		Bsize = sqrt(B.r * B.r + B.z * B.z);
		B.r = B.r / Bsize * d;
		B.z = B.z / Bsize * d;
		B.id = 22;
		return B;
	}
	else {				//	For Divertor
		Nr = -dz0 * LR * f.ro;
		Nz =  dr0 * LR * f.ro;
		DoAllDomain(i, g.ndiv) {
			Ar = g.Divt_r[i][0] - r;
			Az = g.Divt_z[i][0] - z;
			A = sqrt(Ar*Ar + Az*Az);
			Cr = g.Divt_r[i][1] - r;
			Cz = g.Divt_z[i][1] - z;
			C = sqrt(Cr*Cr + Cz*Cz);

			if ((Ar*Cr + Az*Cz) / A / C < -0.8) {
				if ((Nr * Ar + Nz * Az) > 0) {
					B.r = g.Divt_r[i][0] - r;
					B.z = g.Divt_z[i][0] - z;
				}
				else if ((Nr * Cr + Nz * Cz) > 0) {
					B.r = g.Divt_r[i][1] - r;
					B.z = g.Divt_z[i][1] - z;
				}
				else
					Errmsg("Something is WRONG in BVEC", -25);

				Bsize = sqrt(B.r * B.r + B.z * B.z);
				B.r = B.r / Bsize * d;
				B.z = B.z / Bsize * d;
				B.id = 13;
				return B;
			}
		}
	}
	return B;
}


/*  Part Mesh Generator	*/	
void PARTMESHGEN(SYSALL &sys, int &k, int YCo, int nDi, int nL, int nR, double *psiL, double *psiR) {
	usealias_geo(g, sys);
	usealias_fld(f, sys);
	usealias_mes(m, sys);
	TESTPOINT	BL, BR, V;
	double		w, wY = 1.0, wX = 2.5;
	int			i, j;
	
	//	Half upper
	BL = BVEC(sys, k, -1, 1);	//	Start from X point
	if (BL.id == 0)
		Errmsg("NORMAL VECTOR TRACE ERROR IN PARTMESHGEN", -k * 10 - 5);
	BR = BVEC(sys, k,  1, 1);
	if (BR.id == 0)
		Errmsg("NORMAL VECTOR TRACE ERROR IN PARTMESHGEN", -k * 10 - 6);

	DoAllDomain(i, (nDi + 1) / 2) {
		DoAllDomain(j, nL - 1) {	//	Left side of Matrix (SOL, AS)
			w = (erf(-(i + j - f.dXp[0].j) / wX) + 1) * wY;
			V = NORVECTRC(sys, m.Mr[HALF - 1 - j][k], m.Mz[HALF - 1 - j][k],
				BL.r, BL.z, -1, psiL[j], w);
			if (V.id == 0)
				Errmsg("NORMAL VECTOR TRACE ERROR IN PARTMESHGEN", -k * 10 - 1);
			m.Mr[HALF - j - 2][k] = V.r;
			m.Mz[HALF - j - 2][k] = V.z;
		}
		DoAllDomain(j, nR - 1) {	//	Right side of Matrix (Pr, Core)
			w = (erf(-(i + j - f.dXp[0].j) / wX) + 1) * wY;
			V = NORVECTRC(sys, m.Mr[HALF - 1 + j][k], m.Mz[HALF - 1 + j][k],
				BR.r, BR.z,  1, psiR[j], w);
			if (V.id == 0)
				Errmsg("NORMAL VECTOR TRACE ERROR IN PARTMESHGEN", -k * 10 - 2);
			m.Mr[HALF + j][k] = V.r;
			m.Mz[HALF + j][k] = V.z;
		}
		k++;
	}

	//	Half under
	BL = BVEC(sys, k + nDi / 2 - 1, -1, YCo);
	if (BL.id == 0)
		Errmsg("NORMAL VECTOR TRACE ERROR IN PARTMESHGEN", -k * 10 - 7);
	BR = BVEC(sys, k + nDi / 2 - 1,  1, YCo);
	if (BR.id == 0)
		Errmsg("NORMAL VECTOR TRACE ERROR IN PARTMESHGEN", -k * 10 - 8);
	DoAllDomain(i, nDi / 2) {
		DoAllDomain(j, nL - 1) {	//	Left side of Matrix (SOL, AS)
			if (YCo == 2 && i == nDi / 2 - 1){
				w = 1;
				if (m.Type.id == DDN && j == m.nD[AS] - 2) {
					m.Mr[HALF - j - 2][k] = m.Type.r;
					m.Mz[HALF - j - 2][k] = m.Type.z;
					break;
				}
			}
			else if (YCo != 0) 
				w = (erf((i - j + f.dXp[0].j - (nDi / 2 - 1)) / wX) + 1) * wY;
			else
				w = (erf((i + f.dXp[0].j - (nDi / 2 - 1)) / wX / 2.5) + 1) * wY;

			V = NORVECTRC(sys, m.Mr[HALF - 1 - j][k], m.Mz[HALF - 1 - j][k],
				BL.r, BL.z, -1, psiL[j], w);
			if (V.id == 0)
				Errmsg("NORMAL VECTOR TRACE ERROR IN PARTMESHGEN", -k * 10 - 3);
			m.Mr[HALF - j - 2][k] = V.r;
			m.Mz[HALF - j - 2][k] = V.z;
		}
		DoAllDomain(j, nR - 1) {	//	Right side of Matrix (Pr, Core)
			if (YCo != 0)
				w = (erf((i - j + f.dXp[0].j - (nDi / 2 - 1)) / wX) + 1) * wY;
			else
				w = (erf((i + f.dXp[0].j - (nDi / 2 - 1)) / wX / 2.5) + 1) * wY;

			V = NORVECTRC(sys, m.Mr[HALF - 1 + j][k], m.Mz[HALF - 1 + j][k], 
				BR.r, BR.z, 1, psiR[j], w);
			if (V.id == 0)
				Errmsg("NORMAL VECTOR TRACE ERROR IN PARTMESHGEN", -k * 10 - 4);
			m.Mr[HALF + j][k] = V.r;
			m.Mz[HALF + j][k] = V.z;
		}
		k++;
	}
}


/*  MESH GENERATOR2  */
void MESHGEN2(SYSALL &sys) {
	usealias_geo(g, sys);
	usealias_fld(f, sys);
	usealias_mes(m, sys);
	int i, j, k, l, ii, jj;
	double r, z, psiL[HALF], psiR[HALF], psiDUM[HALF];
	double d, zV;
	int nBdy, nM, nS;
	double dumVr, dumVz, dumV;
	TESTPOINT V;
	FILE *out;

	MESHPAR(sys);

	STRETCHING(sys);

	nBdy = 2;
	
	DoAllDomain(i, 2 * HALF) {
		DoAllDomain(j, KMAX) {
			m.Mr[i][j] = 0;
			m.Mz[i][j] = 0;
		}
	}

	k = 0;
	DoAllDomain(i, m.Sn) {
		j = 0;
		while (1) {
			if (m.Sr[i][j] == 0)
				break;
			m.Mr[HALF - 1][k] = m.Sr[i][j];
			m.Mz[HALF - 1][k] = m.Sz[i][j];
			j++;
			k++;
		}
	}
	nM = k;

	d = distance(m.Type.r, m.Type.z, f.Op[0].r, f.Op[0].z);

	/////////////////////////////////////////////
	k = 0;
	DoAllDomain(i, HALF)
		psiL[i] = psiR[i] = 0;
	/*  Core  */
	if (m.Type.id == DDN) {
		STRETCHPAR(psiL, m.dP[AS] * (m.nD[AS] - 1), m.nD[AS] - 1,
			m.EDS[AS][0], m.EDS[AS][1], m.EDS[AS][2]);
		DoAllDomain(i, m.nD[AS] - 1)	psiL[i] += f.Xp[0].v + f.Xp[1].v - m.Type.v;

		STRETCHPAR(psiDUM, m.dP[SOL] * m.nD[SOL], m.nD[SOL] - 1,
			m.EDS[SOL][0], m.EDS[SOL][1], m.EDS[SOL][2]);
		DoAllDomain(i, m.nD[SOL] - 1)	psiL[i + m.nD[AS] - 1] = m.Type.v + psiDUM[i];

		STRETCHPAR(psiR, m.dP[Co] * m.nD[Co], m.nD[Co] - 1,
			m.EDS[Co][0], m.EDS[Co][1], m.EDS[Co][2]);
		DoAllDomain(i, m.nD[Co] - 1)	psiR[i] += f.Xp[0].v + f.Xp[1].v - m.Type.v;

		DoAllDomain(l, nBdy) {
			zV = inorout(f.Op[0].r, f.Op[0].z, m.Mr[HALF - 1][k], m.Mz[HALF - 1][k],
				m.Mr[HALF - 1][k + 2], m.Mz[HALF - 1][k + 2]);
			if (zV > 0)
				nS = m.nD[IBdy];
			else
				nS = m.nD[OBdy];
		
			PARTMESHGEN(sys, k, 2, nS, m.nD[AS] + m.nD[SOL] - 1, m.nD[Co], psiL, psiR);
		}
	}
	else {
		STRETCHPAR(psiL, m.dP[SOL] * m.nD[SOL], m.nD[SOL] - 1,
			m.EDS[SOL][0], m.EDS[SOL][1], m.EDS[SOL][2]);
		DoAllDomain(i, m.nD[SOL] - 1)	psiL[i] += m.Type.v;

		STRETCHPAR(psiR, m.dP[Co] * m.nD[Co], m.nD[Co] - 1,
			m.EDS[Co][0], m.EDS[Co][1], m.EDS[Co][2]);	//	1 0.7 5.0
		DoAllDomain(i, m.nD[Co] - 1)	psiR[i] += m.Type.v;

		DoAllDomain(l, nBdy) {
			zV = inorout(f.Op[0].r, f.Op[0].z, m.Mr[HALF - 1][k], m.Mz[HALF - 1][k],
				m.Mr[HALF - 1][k + 2], m.Mz[HALF - 1][k + 2]);
			if (zV > 0)
				nS = m.nD[IBdy];
			else
				nS = m.nD[OBdy];
			if (m.Type.id == SN)
				PARTMESHGEN(sys, k, 2, nS, m.nD[SOL], m.nD[Co], psiL, psiR);
			else
				PARTMESHGEN(sys, k, 1, nS, m.nD[SOL], m.nD[Co], psiL, psiR);
		}
	}
	printf("(Body) ");

	/*  Other Separatrix  */
	DoAllDomain(j, m.Sn - nBdy) {
		DoAllDomain(i, HALF)
			psiL[i] = psiR[i] = 0;

		dumVr = m.Mr[HALF - 1][k + 5] - m.Mr[HALF - 1][k];
		dumVz = m.Mz[HALF - 1][k + 5] - m.Mz[HALF - 1][k];
		dumV = dumVr * (f.Op[0].r - m.Type.r) + dumVz * (f.Op[0].z - m.Type.z); 
		if ((m.Type.id == DDN) && (dumV > 0) && 
			(m.Mr[HALF - 1][k] == m.Type.r && m.Mz[HALF - 1][k] == m.Type.z)) {
			do
				k++;
			while ((m.Mr[HALF - 1][k] != f.Xp[0].r && m.Mz[HALF - 1][k] != f.Xp[0].z)
				&& (m.Mr[HALF - 1][k] != f.Xp[1].r && m.Mz[HALF - 1][k] != f.Xp[1].z));
		}
		else {
			//	Separatirx line inside of others
			zV = inorout(f.Op[0].r, f.Op[0].z, m.Mr[HALF - 1][k], m.Mz[HALF - 1][k],
				m.Mr[HALF - 1][k + 2], m.Mz[HALF - 1][k + 2]);
			if (zV > 0)
				nS = m.nD[IDt];
			else
				nS = m.nD[ODt];

			if ((m.Type.id == DDN) &&
				fabs(m.Mr[HALF - 1][k] - (f.Xp[0].r + f.Xp[1].r - m.Type.r)) < EPS &&
				fabs(m.Mz[HALF - 1][k] - (f.Xp[0].z + f.Xp[1].z - m.Type.z)) < EPS) {

				STRETCHPAR(psiL, m.dP[AS] * (m.nD[AS] - 1), m.nD[AS] - 1,
					m.EDS[AS][0], m.EDS[AS][1], m.EDS[AS][2]);
				DoAllDomain(i, m.nD[AS] - 1)	psiL[i] += f.Xp[0].v + f.Xp[1].v - m.Type.v;

				STRETCHPAR(psiDUM, m.dP[SOL] * m.nD[SOL], m.nD[SOL] - 1,
					m.EDS[SOL][0], m.EDS[SOL][1], m.EDS[SOL][2]);
				DoAllDomain(i, m.nD[SOL] - 1)	psiL[i + m.nD[AS] - 1] = m.Type.v + psiDUM[i];

				STRETCHPAR(psiR, m.dP[APr] * m.nD[APr], m.nD[APr] - 1,
					m.EDS[APr][0], m.EDS[APr][1], m.EDS[APr][2]);
				DoAllDomain(i, m.nD[APr] - 1)	psiR[i] += f.Xp[0].v + f.Xp[1].v - m.Type.v;

				PARTMESHGEN(sys, k, 0, nS, m.nD[AS] + m.nD[SOL] - 1, m.nD[APr], psiL, psiR);
			}
			else {
				STRETCHPAR(psiL, m.dP[SOL] * m.nD[SOL], m.nD[SOL] - 1,
					m.EDS[SOL][0], m.EDS[SOL][1], m.EDS[SOL][2]);
				DoAllDomain(i, m.nD[SOL] - 1)	psiL[i] += m.Type.v;

				STRETCHPAR(psiR, m.dP[Pr] * m.nD[Pr], m.nD[Pr] - 1,
					m.EDS[Pr][0], m.EDS[Pr][1], m.EDS[Pr][2]);
				DoAllDomain(i, m.nD[Pr] - 1)	psiR[i] += m.Type.v;

				PARTMESHGEN(sys, k, 0, nS, m.nD[SOL], m.nD[Pr], psiL, psiR);
			}
			printf("-> (Divt%d)", j);
		}
		
		//	Print out
		out = fopen("Mr.txt", "w");
		DoAllDomain(jj, KMAX) {
			DoAllDomain(ii, 2 * HALF) {
				fprintf(out, "%lf\t", m.Mr[ii][jj]);
			}
			fprintf(out, "\n");
		}
		fclose(out);
		out = fopen("Mz.txt", "w");
		DoAllDomain(jj, KMAX) {
			DoAllDomain(ii, 2 * HALF) {
				fprintf(out, "%lf\t", m.Mz[ii][jj]);
			}
			fprintf(out, "\n");
		}
		fclose(out);
	}
	if (k != nM)
		Errmsg("SOMETHING is MISSED in MESH MATRIX", -14);

	//	Print out
	out = fopen("Mr.txt", "w");
	DoAllDomain(j, KMAX) {
		DoAllDomain(i, 2 * HALF) {
			fprintf(out, "%lf\t", m.Mr[i][j]);
		}
		fprintf(out, "\n");
	}
	fclose(out);
	out = fopen("Mz.txt", "w");
	DoAllDomain(j, KMAX) {
		DoAllDomain(i, 2 * HALF) {
			fprintf(out, "%lf\t", m.Mz[i][j]);
		}
		fprintf(out, "\n");
	}
	fclose(out);

	printf("\tDone!\n");
}


/*  MESH REFINEMENT  
void MESHREFINE(SYSALL &sys) {
	usealias_geo(g, sys);
	usealias_fld(f, sys);
	usealias_mes(m, sys);
	int i, j, ii, jj;
	
	FILE *out;

	





}*/
