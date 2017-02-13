/***************************************************/
/*												   */
/*		POST PROCESSING FUNCTION				   */
/*					for C2 Code					   */
/*		JG LEE									   */
/*		20/JUL/2016								   */
/***************************************************/

#define	_CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <string>

#include "parameter.h"
#include "svc.h"
#include "vega.h"
#include "memalloc.h"
#include "separatrix.h"
#include "postprocess.h"

void MESHPART(SYSALL &sys) {
	usealias_geo(g, sys);
	usealias_fld(f, sys);
	usealias_mes(m, sys);
	int i, j, ii, jj, k, nUO, nUI, nLO, nLI;
	int iXp;
	int nS, nBdy;
	double	OX, OP;
	double	d, z;
	double	dummy;
	char	title[255];
	FILE*	out;

	/*	0 : Outer Body, 1 : Inner Body
		2 ~ 5 : Upper Separatrix (Outer to Inner)
		6 ~ 9 : Lower Separatirx (Outer to Inner)	*/
	m.SMr = alloc_tensor<double>(0, 10 - 1, 0, 2 * HALF - 1, 0, 100 - 1);
	m.SMz = alloc_tensor<double>(0, 10 - 1, 0, 2 * HALF - 1, 0, 100 - 1);

	printf("Separation of mesh matrix\t");

	DoAllDomain(i, 10) {
		DoAllDomain(j, 2 * HALF) {
			DoAllDomain(k, 100)
				m.SMr[i][j][k] = m.SMz[i][j][k] = 0;
		}
	}

	d = distance(m.Type.r, m.Type.z, f.Op[0].r, f.Op[0].z);

	k = 0;
	DoAllDomain(ii, 2) {
		z = inorout(f.Op[0].r, f.Op[0].z, m.Mr[HALF - 1][k], m.Mz[HALF - 1][k],
			m.Mr[HALF - 1][k + 2], m.Mz[HALF - 1][k + 2]);
		if (z < 0) {			//	Outer
			DoAllDomain(jj, m.nD[OBdy]) {
				DoAllDomain(i, 2 * HALF) {
					m.SMr[0][i][jj] = m.Mr[i][k];
					m.SMz[0][i][jj] = m.Mz[i][k];
				}
				k++;
			}
		}
		else {					//	Inner
			DoAllDomain(jj, m.nD[IBdy]) {
				DoAllDomain(i, 2 * HALF) {
					m.SMr[1][i][jj] = m.Mr[i][k];
					m.SMz[1][i][jj] = m.Mz[i][k];
				}
				k++;
			}
		}
	}	

	if (m.Type.id == DDN) {
		k += m.nD[IBdy] + m.nD[OBdy] + m.nD[ODt] + m.nD[IDt];
	}

	nUO = 2;	nUI = 4;
	nLO = 6;	nLI = 8;
	while (m.Mr[HALF - 1][k] != 0) {
		/*z = ((m.Mz[HALF - 1][k + 2] > 0) - (m.Mz[HALF - 1][k + 2] < 0)) * 
			(f.Op[0].r - m.Mr[HALF - 1][k + 2]) * (f.Op[0].z - m.Type.z) / d
			+ (f.Op[0].z - m.Mz[HALF - 1][k + 2]) * (f.Op[0].r - m.Type.r) / d;*/
		z = inorout(f.Op[0].r, f.Op[0].z, m.Mr[HALF - 1][k], m.Mz[HALF - 1][k],
			m.Mr[HALF - 1][k + 2], m.Mz[HALF - 1][k + 2]);
		if (m.Mz[HALF - 1][k] > f.Op[0].z) {
			if (z < 0) {
				DoAllDomain(jj, m.nD[ODt]) {
					DoAllDomain(i, 2 * HALF) {
						m.SMr[nUO][i][jj] = m.Mr[i][k];
						m.SMz[nUO][i][jj] = m.Mz[i][k];
					}
					k++;
				}
				nUO++;
			}
			else {
				DoAllDomain(jj, m.nD[IDt]) {
					DoAllDomain(i, 2 * HALF) {
						m.SMr[nUI][i][jj] = m.Mr[i][k];
						m.SMz[nUI][i][jj] = m.Mz[i][k];
					}
					k++;
				}
				nUI++;
			}
		}
		else {
			if (z < 0) {
				DoAllDomain(jj, m.nD[ODt]) {
					DoAllDomain(i, 2 * HALF) {
						m.SMr[nLO][i][jj] = m.Mr[i][k];
						m.SMz[nLO][i][jj] = m.Mz[i][k];
					}
					k++;
				}
				nLO++;
			}
			else {
				DoAllDomain(jj, m.nD[IDt]) {
					DoAllDomain(i, 2 * HALF) {
						m.SMr[nLI][i][jj] = m.Mr[i][k];
						m.SMz[nLI][i][jj] = m.Mz[i][k];
					}
					k++;
				}
				nLI++;
			}
		}
	}

	/*   SORT	*/
	if (nUO == 4) {
		if (m.SMz[nUO - 2][HALF - 1][2] > m.SMz[nUO - 1][HALF - 1][2]) {
			DoAllDomain(jj, m.nD[ODt]) {
				DoAllDomain(i, 2 * HALF) {
					dummy = m.SMr[nUO - 2][i][jj];
					m.SMr[nUO - 2][i][jj] = m.SMr[nUO - 1][i][jj];
					m.SMr[nUO - 1][i][jj] = dummy;
					dummy = m.SMz[nUO - 2][i][jj];
					m.SMz[nUO - 2][i][jj] = m.SMz[nUO - 1][i][jj];
					m.SMz[nUO - 1][i][jj] = dummy;
				}
			}
		}
	}
	if (nUI == 6) {
		if (m.SMz[nUI - 2][HALF - 1][2] > m.SMz[nUI - 1][HALF - 1][2]) {
			DoAllDomain(jj, m.nD[IDt]) {
				DoAllDomain(i, 2 * HALF) {
					dummy = m.SMr[nUI - 2][i][jj];
					m.SMr[nUI - 2][i][jj] = m.SMr[nUI - 1][i][jj];
					m.SMr[nUI - 1][i][jj] = dummy;
					dummy = m.SMz[nUI - 2][i][jj];
					m.SMz[nUI - 2][i][jj] = m.SMz[nUI - 1][i][jj];
					m.SMz[nUI - 1][i][jj] = dummy;
				}
			}
		}
	}
	if (nLO == 8) {
		if (m.SMz[nLO - 2][HALF - 1][2] < m.SMz[nLO - 1][HALF - 1][2]) {
			DoAllDomain(jj, m.nD[ODt]) {
				DoAllDomain(i, 2 * HALF) {
					dummy = m.SMr[nLO - 2][i][jj];
					m.SMr[nLO - 2][i][jj] = m.SMr[nLO - 1][i][jj];
					m.SMr[nLO - 1][i][jj] = dummy;
					dummy = m.SMz[nLO - 2][i][jj];
					m.SMz[nLO - 2][i][jj] = m.SMz[nLO - 1][i][jj];
					m.SMz[nLO - 1][i][jj] = dummy;
				}
			}
		}
	}
	if (nLI == 10) {
		if (m.SMz[nLI - 2][HALF - 1][2] < m.SMz[nLI - 1][HALF - 1][2]) {
			DoAllDomain(jj, m.nD[IDt]) {
				DoAllDomain(i, 2 * HALF) {
					dummy = m.SMr[nLI - 2][i][jj];
					m.SMr[nLI - 2][i][jj] = m.SMr[nLI - 1][i][jj];
					m.SMr[nLI - 1][i][jj] = dummy;
					dummy = m.SMz[nLI - 2][i][jj];
					m.SMz[nLI - 2][i][jj] = m.SMz[nLI - 1][i][jj];
					m.SMz[nLI - 1][i][jj] = dummy;
				}
			}
		}
	}
	DoAllDomain(i, 10) {
		sprintf(title, "%d.dat", i);
		out = fopen(title,"w");
		//fprintf(out, "ZONE T=\"ZONE%d\", I=\t%d, J=\t%d, F=POINT\n", i, 120, 100);
		DoAllDomain(jj, 100) {
			DoAllDomain(ii, 2 * HALF) {
				if(m.SMr[i][ii][jj] != 0)
					fprintf(out, "%e %e\n", m.SMr[i][ii][jj], m.SMz[i][ii][jj]);
			}
		}
		fclose(out);
	}
	printf("Done!\n");
}

void ASSEMBLER(SYSALL &sys, int id, int nS, int nW, int nH, int &k, int dW, int dH, int jW, int jH) {
	usealias_geo(g, sys);
	usealias_fld(f, sys);
	usealias_mes(m, sys);
	int i, j, ii, jj;
	int iii, jjj;

	DoAllDomain(i, nH) {
		DoAllDomain(j, nW) {
			iii = jH + dH * i;
			jjj = HALF - 1 + jW + dW  * j;
			m.SM[id][k][0] = m.SMr[nS][jjj][iii];
			m.SM[id][k][1] = m.SMz[nS][jjj][iii];

			ii = (int)((m.SM[id][k][0] - g.r[0]) / (g.r[1] - g.r[0]));
			jj = (int)((m.SM[id][k][1] - g.z[0]) / (g.z[1] - g.z[0]));
			if (ii < 0)		ii = 0;
			if (jj < 0)		jj = 0;
			if (ii > g.rN - 1)	ii = g.rN - 1;
			if (jj > g.zN - 1)	jj = g.zN - 1;

			m.SM[id][k][2] = BiliInterp(f.Br, g.r, g.z, m.SM[id][k][0], m.SM[id][k][1], ii, jj);
			m.SM[id][k][3] = BiliInterp(f.Bz, g.r, g.z, m.SM[id][k][0], m.SM[id][k][1], ii, jj);
			m.SM[id][k][4] = LINEAR1D(g.Bt[ii], g.Bt[ii + 1], g.r[ii], g.r[ii + 1], m.SM[id][k][0]);

			k++;
		}
	}
}

void PRINTMESH(SYSALL& sys, int id, int dom, int k, int nH, int nW, int *b, int *t) {
	//usealias_mes(m, sys);
	MESHDATA&	m = sys.mes;
	int		i, j, ii;
	char 	title[256];
	FILE	*fp;

	sprintf(title, "geo.%d.inp", id);
	//std::cfp << std::to_string(id) << '\n';
	//sprintf(title, "D:\\JAEGON\\Projects\\VEGA\\VEGA\\geo.%d.inp", id);
	fp = fopen(title, "w");
	fprintf(fp, "@DomainID\n%d\n", dom);
	fprintf(fp, "@NodeNum\n%d %d\n", nH, nW);
	fprintf(fp, "@Est\n%d\n", 0);
	fprintf(fp, "@Wst\n%d\n", 0);

	fprintf(fp, "@Top\n");
	if (t[0] == 0) fprintf(fp, "%d\n", 0);
	else {
		i = 0;
		while (t[i] != 0) {
			fprintf(fp, "%d ", t[i]);
			i++;
		}
		fprintf(fp, "\n");
	}

	fprintf(fp, "@Bot\n");
	if (b[0] == 0) fprintf(fp, "%d\n", 0);
	else {
		i = 0;
		while (b[i] != 0) {
			fprintf(fp, "%d ", b[i]);
			i++;
		}
		fprintf(fp, "\n");
	}

	fprintf(fp, "@XY\n");
	DoAllDomain (i, k) {
		DoAllDomain(j, 5) {
			fprintf(fp, "%e ", m.SM[id][i][j]);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	//sprintf(title, "D:\\JAEGON\\Projects\\VEGA\\VEGA\\geo%d.dat", id);
	sprintf(title, "geo%d.dat", id);

	fp = fopen(title, "w");
	//fprintf(fp, "variables= \"R\",\"Z\"\n");//,\"Br\",\"Bz\",\"BT\"
	fprintf(fp, "ZONE T=\"ZONE%d\", I=\t%d, J=\t%d, F=POINT\n", id, nW, nH);
	DoAllDomain(i, k) {
			//DoAllDomain(j, 2) {
				fprintf(fp, "%lf\t%lf\n", m.SM[id][i][0], m.SM[id][i][1]);
			//}
			//fprintf(fp, "\n");
	}
	fclose(fp);
}

void ONEDPRINT(SYSALL &sys, int id, int nH, int nIS, int nOS, int Nedge) {
	usealias_geo(g, sys);
	usealias_fld(f, sys);
	usealias_mes(m, sys);
	int		i, j, ii, jj, k, l;
	int		jsN = 0, nBdy;
	int		nLast = HALF + m.nD[Co] - 2;
	dmatrix	jsM, Cr, Cz;
	double	Rmin, Rmax, B0, dl, LBp, PiLoop;
	double	RLoc, ZLoc, BtLoc, BpLoc;
	char 	title[255];
	FILE	*out;

	jsM = alloc_matrix<double>(0, 10, 0, nH);
	Cr	= alloc_matrix<double>(0, nH, 0, nIS + nOS - 1);
	Cz  = alloc_matrix<double>(0, nH, 0, nIS + nOS - 1);

	nBdy = 2;

	j = 0;
	DoAllDomain(l, nOS) {
		DoAllDomain(i, nH) {
			Cr[i][j] = m.SMr[0][nLast - i][l];
			Cz[i][j] = m.SMz[0][nLast - i][l];
		}
		j++;
	}
	DoDomain(l, 1, nIS) {
		DoAllDomain(i, nH) {
			Cr[i][j] = m.SMr[1][nLast - i][nIS - 1 - l];
			Cz[i][j] = m.SMz[1][nLast - i][nIS - 1 - l];
		}
		j++;
	}

	ii = (int)((f.Op[0].r - g.r[0]) / (g.r[1] - g.r[0]));
	if (ii < 0)			ii = 0;
	if (ii > g.rN - 1)	ii = g.rN - 1;
	B0 = LINEAR1D(g.Bt0[ii], g.Bt0[ii + 1], g.r[ii], g.r[ii + 1], f.Op[0].r);

	DoAllDomain(i, nH) {
		PiLoop	= 0;
		LBp		= 0;
		DoAllDomain(j, nOS + nIS - 2) {
			RLoc = (Cr[i][j] + Cr[i][j + 1]) / 2.;
			ZLoc = (Cz[i][j] + Cz[i][j + 1]) / 2.;
			ii = (int)((RLoc - g.r[0]) / (g.r[1] - g.r[0]));
			jj = (int)((ZLoc - g.z[0]) / (g.z[1] - g.z[0]));
			if (ii < 0)			ii = 0;
			if (jj < 0)			jj = 0;
			if (ii > g.rN - 1)	ii = g.rN - 1;
			if (jj > g.zN - 1)	jj = g.zN - 1;

			BpLoc = BiliInterp(f.Bp, g.r, g.z, RLoc, ZLoc, ii, jj);
			dl = distance(Cr[i][j], Cz[i][j], Cr[i][j + 1], Cz[i][j + 1]);
			LBp += BpLoc * dl;

			BtLoc = LINEAR1D(g.Bt[ii], g.Bt[ii + 1], g.r[ii], g.r[ii + 1], RLoc);
			PiLoop += (fabs(Cz[i][j]) + fabs(Cz[i][j + 1])) / 2.
				 	 * fabs(Cr[i][j] - Cr[i][j + 1]) * BtLoc;
		}

		Rmin = Rmax = Cr[i][0];
		DoDomain(j, 1, nOS + nIS - 1) {
			Rmin = min(Rmin, Cr[i][j]);
			Rmax = max(Rmax, Cr[i][j]);
		}

		i++;
		jsM[0][i]	= sqrt(fabs(PiLoop / PI / B0));			//	rho
		jsM[1][i]	= (Rmax - Rmin) / 2.;					//	ameter
		jsM[2][i]	= 0.;									//	shift
		jsM[3][i]	= 0.;									//	elong
		jsM[4][i]	= 0.;									//	triag
		jsM[5][i]	= (Rmax + Rmin) / 2.;					//	rtor
		jsM[6][i]	= 0.;									//	g11
		jsM[7][i]	= 2 * PI * LBp;							//	g22
		jsM[8][i]	= 0.;									//	g33
		jsM[9][i]	= 4 * PI * PI * jsM[1][i] * jsM[5][i];	//	volp
		jsM[10][i]	= RLoc * BtLoc;							//	ftor		
		i--;
	}

	//	Adding values at rho == 0
	jsM[0][0]	= 0.;	jsM[1][0]	= 0.;	jsM[2][0]	= 0.;
	jsM[3][0]	= 0.;	jsM[4][0]	= 0.;	jsM[5][0]	= f.Op[0].r;	
	jsM[6][0]	= 0.;	jsM[7][0]	= 0.;	jsM[8][0]	= 0.;
	jsM[9][0]	= jsM[9][1] - (jsM[9][2] - jsM[9][1]) / (jsM[1][2] - jsM[1][1])
				* jsM[1][1] * 0.8;	//	Linear extrapolation with a under compensation(?)
	jsM[10][0]	= (jsM[10][1] * jsM[1][2] - jsM[10][2] * jsM[1][1])
				/ (jsM[1][2] - jsM[1][1]);		//	linear extrapolation

	//sprintf(title, "D:\\JAEGON\\Projects\\VEGA\\VEGA\\geo.%d.inp", id);
	sprintf(title, "geo.%d.inp", id);

	out = fopen(title, "a");	//	a : add below

	fprintf(out, "@js\n%d\n@ONED\n", jsN);

	DoAllDomain(j, nH + 1) {	//	inclunding 0 point
		DoAllDomain(i, 11)
			fprintf(out, "%e ", jsM[i][j]);
		fprintf(out, "\n");
	}
	fclose(out);

	out = fopen("J2.txt", "w");	//	a : add below
	
	DoAllDomain(j, nH + 1) {	//	inclunding 0 point
		DoAllDomain(i, 11)
			fprintf(out, "%e ", jsM[i][j]);
		fprintf(out, "\n");
	}
	fclose(out);
}

void POSTP(SYSALL &sys) {
	usealias_geo(g, sys);
	usealias_fld(f, sys);
	usealias_mes(m, sys);
	int k;
	int Nedge = m.nD[Co] / 10;
	const int KMULT = 30;
	int ni, nj;
	int top0[] = { 0 }, bot0[] = { 0 };
	int OSbot[] = { 4, 4, 1, 5, 4, 5, 1, 2, 2, 
					1, m.nD[ODt], m.nD[ODt] + (m.nD[OBdy] - 1) / 2,
					m.nD[ODt] + m.nD[OBdy] - 1, 2 * m.nD[ODt] + m.nD[OBdy] - 2, 0 },
		ISbot[] = { 4, 2, 1, 5, 2, 5, 3, 4, 2,
					1, m.nD[IDt], m.nD[IDt] + (m.nD[IBdy] - 1) / 2,
					m.nD[IDt] + m.nD[IBdy] - 1, 2 * m.nD[IDt] + m.nD[IBdy] - 2, 0 },
		UPtop[] = { 2, 1, 1, 3, 4, 1, m.nD[ODt], m.nD[ODt] + m.nD[IDt] - 1, 0},
		LPtop[] = { 2, 3, 1, 1, 4, 1, m.nD[IDt], m.nD[IDt] + m.nD[ODt] - 1, 0 },
		EDtop[] = { 4, 1, 3, 3, 2, 3, 3, 1, 2, 1, (m.nD[IBdy] + 1) / 2, m.nD[IBdy],
					m.nD[IBdy] + (m.nD[OBdy] - 1) / 2, m.nD[IBdy] + m.nD[OBdy] - 1, 0 },
		EDbot[] = { 1, 6, 1, 1, m.nD[IBdy] + m.nD[OBdy] - 1, 0 },
		COtop[] = { 1, 5, 1, 1, m.nD[IBdy] + m.nD[OBdy] - 1, 0 };

	MESHPART(sys);

	m.SM = alloc_tensor<double>(0, 9, 0, KMAX * KMULT, 0, 4);

	/*  Outer SOL	dom = 0 -> 2 + 0 + 6	*/
	k = 0;
	ASSEMBLER(sys, 1, 2, m.nD[SOL], m.nD[ODt]  - 1, k, -1, -1, 0, m.nD[ODt] - 1);
	ASSEMBLER(sys, 1, 0, m.nD[SOL], m.nD[OBdy] - 1, k, -1, -1, 0, m.nD[OBdy] - 1);
	ASSEMBLER(sys, 1, 6, m.nD[SOL], m.nD[ODt]	  , k, -1,  1, 0, 0);
	ni = m.nD[ODt] + m.nD[OBdy] + m.nD[ODt] - 2;
	nj = m.nD[SOL];
	PRINTMESH(sys, 1, 0, k, ni, nj, OSbot, top0);
	printf("(OSOL) > ");
	/*  Inner SOL	dom = 2 -> 8 + 1 + 4	*/
	k = 0;
	ASSEMBLER(sys, 3, 8, m.nD[SOL], m.nD[IDt]  - 1, k, -1, -1, 0, m.nD[IDt] - 1);
	ASSEMBLER(sys, 3, 1, m.nD[SOL], m.nD[IBdy] - 1, k, -1,  1, 0, 0);
	ASSEMBLER(sys, 3, 4, m.nD[SOL], m.nD[IDt]	  , k, -1,  1, 0, 0);
	ni = m.nD[IDt] + m.nD[IBdy] + m.nD[IDt] - 2;
	nj = m.nD[SOL];
	PRINTMESH(sys, 3, 2, k, ni, nj, ISbot, top0);
	printf("(ISOL) > ");
	/*  Upper Private	dom = 5 -> 2 + 4	*/
	k = 0;
	ASSEMBLER(sys, 4, 2, m.nD[Pr], m.nD[ODt] - 1, k, -1, -1, m.nD[Pr] - 1, m.nD[ODt] - 1);
	ASSEMBLER(sys, 4, 4, m.nD[Pr], m.nD[IDt]    , k, -1,  1, m.nD[Pr] - 1, 0);
	ni = m.nD[ODt] + m.nD[IDt] - 1;
	nj = m.nD[Pr];
	PRINTMESH(sys, 4, 5, k, ni, nj, bot0, UPtop);
	printf("(UPr) > ");
	/*  Lower Private	dom = 1 -> 8 + 6	*/
	k = 0;
	ASSEMBLER(sys, 2, 8, m.nD[Pr], m.nD[IDt] - 1, k, -1, -1, m.nD[Pr] - 1, m.nD[IDt] - 1);
	ASSEMBLER(sys, 2, 6, m.nD[Pr], m.nD[ODt]    , k, -1,  1, m.nD[Pr] - 1, 0);
	ni = m.nD[IDt] + m.nD[ODt] - 1;
	nj = m.nD[Pr];
	PRINTMESH(sys, 2, 1, k, ni, nj, bot0, LPtop);
	printf("(LPr) > ");
	/*  EDGE Region		dom = 3 -> 1 + 0	*/
	k = 0;
	ASSEMBLER(sys, 5, 0, Nedge, (m.nD[OBdy] - 1) / 2 + 1, k, -1, -1, Nedge - 1, (m.nD[OBdy] - 1) / 2);
	ASSEMBLER(sys, 5, 1, Nedge, m.nD[IBdy] - 1,			  k, -1,  1, Nedge - 1, 1);
	ASSEMBLER(sys, 5, 0, Nedge, (m.nD[OBdy] - 1) / 2,	  k, -1, -1, Nedge - 1, m.nD[OBdy] - 2);
	ni = m.nD[IBdy] + m.nD[OBdy] - 1;
	nj = Nedge;
	PRINTMESH(sys, 5, 3, k, ni, nj, EDbot, EDtop);
	printf("(Edge) > ");
	/*  CORE Region		dom = 4 -> 1 + 0	*/
	k = 0;
	ASSEMBLER(sys, 6, 0, m.nD[Co] - Nedge + 1, (m.nD[OBdy] - 1) / 2 + 1, k, -1, -1,
		m.nD[Co] - 1, (m.nD[OBdy] - 1) / 2);
	ASSEMBLER(sys, 6, 1, m.nD[Co] - Nedge + 1, m.nD[IBdy] - 1, k, -1, 1, m.nD[Co] - 1, 1);
	ASSEMBLER(sys, 6, 0, m.nD[Co] - Nedge + 1, (m.nD[OBdy] - 1) / 2, k, -1, -1,
		m.nD[Co] - 1, m.nD[OBdy] - 2);
	ni = m.nD[IBdy] + m.nD[OBdy] - 1;
	nj = m.nD[Co] - Nedge + 1;
	PRINTMESH(sys, 6, 4, k, ni, nj, bot0, COtop);
	printf("(Core) > ");
	if(m.Type.id == CDN)
		ONEDPRINT(sys, 6, nj, m.nD[IBdy], m.nD[OBdy], Nedge);
	else
		ONEDPRINT(sys, 6, nj, m.nD[IBdy] + m.nD[OBdy], 0, Nedge);
	printf("(1D Data)\tDone!\n");
}


void POSTP_SNOW(SYSALL &sys) {
	usealias_geo(g, sys);
	usealias_fld(f, sys);
	usealias_mes(m, sys);
	int i, j, k, ii, jj, kk;
	int Nedge = 4;
	const int KMULT = 30;
	int ni, nj;
	int top0[] = { 0 }, bot0[] = { 0 };
	int OSbot[] = { 4, 8, 1, 9, 4, 9, 1, 2, 2,
		1, m.nD[ODt], m.nD[ODt] + (m.nD[OBdy] - 1) / 2,
		m.nD[ODt] + m.nD[OBdy] - 1, 2 * m.nD[ODt] + m.nD[OBdy] - 2, 0 },
		ISbot[] = { 4, 4, 1, 9, 2, 9, 3, 6, 2,
		1, m.nD[IDt], m.nD[IDt] + (m.nD[IBdy] - 1) / 2,
		m.nD[IDt] + m.nD[IBdy] - 1, 2 * m.nD[IDt] + m.nD[IBdy] - 2, 0 },
		UPtop1[] = { 2, 7, 1, 5, 4, 1, m.nD[IDt], m.nD[IDt] + m.nD[IDt] - 1, 0 },
		UPbot2[] = { 2, 6, 1, 8, 2, 1, m.nD[IDt], m.nD[IDt] + m.nD[ODt] - 1, 0 },
		UPtop3[] = { 2, 1, 1, 7, 2, 1, m.nD[ODt], m.nD[ODt] + m.nD[ODt] - 1, 0 },
		LPtop1[] = { 2, 5, 1, 3, 2, 1, m.nD[IDt], m.nD[IDt] + m.nD[IDt] - 1, 0 },
		LPbot2[] = { 2, 2, 1, 4, 2, 1, m.nD[ODt], m.nD[ODt] + m.nD[IDt] - 1, 0 },
		LPtop3[] = { 2, 3, 1, 1, 4, 1, m.nD[ODt], m.nD[ODt] + m.nD[ODt] - 1, 0 },
		EDtop[] = { 4, 1, 3, 5, 2, 5, 3, 1, 2, 1, (m.nD[IBdy] + 1) / 2, m.nD[IBdy],
		m.nD[IBdy] + (m.nD[OBdy] - 1) / 2, m.nD[IBdy] + m.nD[OBdy] - 1, 0 },
		EDbot[] = { 1, 10, 1, 1, m.nD[IBdy] + m.nD[OBdy] - 1, 0 },
		COtop[] = { 1, 9, 1, 1, m.nD[IBdy] + m.nD[OBdy] - 1, 0 };

	MESHPART(sys);

	/*	0 : Outer Body, 1 : Inner Body
	2 ~ 5 : Upper Separatrix (Outer to Inner)
	6 ~ 9 : Lower Separatirx (Outer to Inner)	*/
	m.SM = alloc_tensor<double>(1, 10, 0, KMAX * KMULT, 0, 4);

	/*  Outer SOL	dom = 0 -> 2 + 0 + 6	*/
	k = 0;
	ASSEMBLER(sys, 1, 2, m.nD[SOL], m.nD[ODt]  - 1, k, -1, -1, 0, m.nD[ODt] - 1);
	ASSEMBLER(sys, 1, 0, m.nD[SOL], m.nD[OBdy] - 1, k, -1, -1, 0, m.nD[OBdy] - 1);
	ASSEMBLER(sys, 1, 6, m.nD[SOL], m.nD[ODt]     , k, -1,  1, 0, 0);
	ni = m.nD[ODt] + m.nD[OBdy] + m.nD[ODt] - 2;
	nj = m.nD[SOL];
	PRINTMESH(sys, 1, 0, k, ni, nj, OSbot, top0);
	printf("(OSOL) > ");
	/*  Inner SOL	dom = 2 -> 9 + 1 + 5	*/
	k = 0;
	ASSEMBLER(sys, 5, 8, m.nD[SOL], m.nD[IDt]  - 1, k, -1, -1, 0, m.nD[IDt] - 1);
	ASSEMBLER(sys, 5, 1, m.nD[SOL], m.nD[IBdy] - 1, k, -1,  1, 0, 0);
	ASSEMBLER(sys, 5, 4, m.nD[SOL], m.nD[IDt]     , k, -1,  1, 0, 0);
	ni = m.nD[IDt] + m.nD[IBdy] + m.nD[IDt] - 2;
	nj = m.nD[SOL];
	PRINTMESH(sys, 5, 2, k, ni, nj, ISbot, top0);
	printf("(ISOL) > ");
	/*  Upper Private	dom = 8 -> 2 + 3	*/
	k = 0;
	ASSEMBLER(sys, 8, 2, m.nD[Pr], m.nD[ODt] - 1, k, -1, -1, m.nD[Pr] - 1, m.nD[ODt] - 1);
	ASSEMBLER(sys, 8, 3, m.nD[Pr], m.nD[ODt]    , k, -1,  1, m.nD[Pr] - 1, 0);
	ni = m.nD[ODt] + m.nD[ODt] - 1;
	nj = m.nD[Pr];
	PRINTMESH(sys, 8, 5, k, ni, nj, bot0, UPtop3);
	printf("(UPr1) > ");
	/*  Upper Private	dom = 7 -> 5 + 3	*/
	k = 0;
	ASSEMBLER(sys, 7, 5, m.nD[SOL], m.nD[IDt] - 1, k, -1, -1, 0, m.nD[IDt] - 1);
	ASSEMBLER(sys, 7, 3, m.nD[SOL], m.nD[ODt]    , k, -1,  1, 0, 0);
	ni = m.nD[ODt] + m.nD[IDt] - 1;
	nj = m.nD[SOL];
	PRINTMESH(sys, 7, 5, k, ni, nj, UPbot2, top0);
	printf("(UPr2) > ");
	/*  Upper Private	dom = 6 -> 4 + 5	*/
	k = 0;
	ASSEMBLER(sys, 6, 5, m.nD[Pr], m.nD[IDt] - 1, k, -1, -1, m.nD[Pr] - 1, m.nD[IDt] - 1);
	ASSEMBLER(sys, 6, 4, m.nD[Pr], m.nD[IDt]    , k, -1,  1, m.nD[Pr] - 1, 0);
	ni = m.nD[IDt] + m.nD[IDt] - 1;
	nj = m.nD[Pr];
	PRINTMESH(sys, 6, 5, k, ni, nj, bot0, UPtop1);
	printf("(UPr3) > ");
	/*  Lower Private	dom = 4 -> 8 + 9	*/
	k = 0;
	ASSEMBLER(sys, 4, 8, m.nD[Pr], m.nD[IDt] - 1, k, -1, -1, m.nD[Pr] - 1, m.nD[IDt] - 1);
	ASSEMBLER(sys, 4, 9, m.nD[Pr], m.nD[IDt]    , k, -1,  1, m.nD[Pr] - 1, 0);
	ni = m.nD[IDt] + m.nD[IDt] - 1;
	nj = m.nD[Pr];
	PRINTMESH(sys, 4, 1, k, ni, nj, bot0, LPtop1);
	printf("(LPr1) > ");
	/*  Lower Private	dom = 3 -> 7 + 9	*/
	k = 0;
	ASSEMBLER(sys, 3, 7, m.nD[SOL], m.nD[ODt] - 1, k, -1, -1, 0, m.nD[ODt] - 1);
	ASSEMBLER(sys, 3, 9, m.nD[SOL], m.nD[IDt]    , k, -1,  1, 0, 0);
	ni = m.nD[IDt] + m.nD[ODt] - 1;
	nj = m.nD[SOL];
	PRINTMESH(sys, 3, 1, k, ni, nj, LPbot2, top0);
	printf("(LPr2) > ");
	/*  Lower Private	dom = 2 -> 7 + 6	*/
	k = 0;
	ASSEMBLER(sys, 2, 7, m.nD[Pr], m.nD[ODt] - 1, k, -1, -1, m.nD[Pr] - 1, m.nD[ODt] - 1);
	ASSEMBLER(sys, 2, 6, m.nD[Pr], m.nD[ODt]    , k, -1,  1, m.nD[Pr] - 1, 0);
	ni = m.nD[ODt] + m.nD[ODt] - 1;
	nj = m.nD[Pr];
	PRINTMESH(sys, 2, 1, k, ni, nj, bot0, LPtop3);
	printf("(LPr3) > ");
	/*  EDGE Region		dom = 9 -> 0 + 1 + 0	*/
	k = 0;
	ASSEMBLER(sys, 9, 0, Nedge, (m.nD[OBdy] + 1) / 2, k, -1, -1, Nedge - 1, (m.nD[OBdy] - 1) / 2);
	ASSEMBLER(sys, 9, 1, Nedge,  m.nD[IBdy] - 1     , k, -1,  1, Nedge - 1,  1);
	ASSEMBLER(sys, 9, 0, Nedge, (m.nD[OBdy] - 1) / 2, k, -1, -1, Nedge - 1,  m.nD[OBdy] - 2);
	ni = m.nD[IBdy] + m.nD[OBdy] - 1;
	nj = Nedge;
	PRINTMESH(sys, 9, 3, k, ni, nj, EDbot, EDtop);
	printf("(Edge) > ");
	/*  CORE Region		dom = 10 -> 0 + 1 + 0	*/
	k = 0;
	ASSEMBLER(sys, 10, 0, m.nD[Co] - Nedge + 1, (m.nD[OBdy] + 1) / 2, k, -1, -1,
		m.nD[Co] - 1, (m.nD[OBdy] - 1) / 2);
	ASSEMBLER(sys, 10, 1, m.nD[Co] - Nedge + 1,  m.nD[IBdy] - 1     , k, -1,  1, m.nD[Co] - 1, 1);
	ASSEMBLER(sys, 10, 0, m.nD[Co] - Nedge + 1, (m.nD[OBdy] - 1) / 2, k, -1, -1,
		m.nD[Co] - 1, m.nD[OBdy] - 2);
	ni = m.nD[IBdy] + m.nD[OBdy] - 1;
	nj = m.nD[Co] - Nedge + 1;
	PRINTMESH(sys, 10, 4, k, ni, nj, bot0, COtop);
	printf("(Core) > ");
	if (m.Type.id == CDN)
		ONEDPRINT(sys, 10, nj, m.nD[IBdy], m.nD[OBdy], Nedge);
	else
		ONEDPRINT(sys, 10, nj, m.nD[IBdy] + m.nD[OBdy], 0, Nedge);
	printf("(1D Data)\tDone!\n");
}


void POSTP_DDNU(SYSALL &sys) {
	usealias_geo(g, sys);
	usealias_fld(f, sys);
	usealias_mes(m, sys);
	int k, i, j;
	int Nedge = m.nD[Co] / 10;
	const int KMULT = 30;
	int ni, nj, ndIBdy, ndOBdy;
	int top0[] = { 0 }, bot0[] = { 0 };
	double Ax, Ay, Bx, By, C = 0;
	/*
	int OS1, OS2, IS1, IS2, UP1, UP2, LP1, LP2;
	int MT1, MT2, MT3, MT4, MT5, MT6, MT7, MT8, MT9, MT10, MT11, MT12;
	int MB1, MB2, MB3, MB4, MB5, MB6, MB7, MB8, MB9, MB10, MB11, MB12;
	int ED1, ED2, ED3, ED4, ED5, ED6, ED7, ED8;
	if (m.Type.z < 0) {
		OS1 = 1;
	}*/

	int OSbot[] = { 4, 5, 1, 5, 2, 5, 3, 2, 2, 1, m.nD[ODt], m.nD[ODt] + (m.nD[OBdy] - 1) / 2,
		m.nD[ODt] + m.nD[OBdy] - 1, 2 * m.nD[ODt] + m.nD[OBdy] - 2, 0 },
		ISbot[] = { 4, 2, 1, 5, 4, 5, 5, 5, 6, 1, m.nD[IDt], m.nD[IDt] + (m.nD[IBdy] - 1) / 2, 
		m.nD[IDt] + m.nD[IBdy] - 1, 2 * m.nD[IDt] + m.nD[IBdy] - 2, 0 },
		UPtop[] = { 2, 5, 1, 5, 6, 1, m.nD[ODt], m.nD[ODt] + m.nD[IDt] - 1, 0 },
		LPtop[] = { 2, 3, 1, 1, 4, 1, m.nD[IDt], m.nD[IDt] + m.nD[ODt] - 1, 0 },
		MZtop[] = { 6, 1, 1, 1, 2, 1, 3, 3, 2, 3, 3, 3, 4, 
		1, m.nD[ODt], m.nD[ODt] + (m.nD[OBdy] - 1) / 2, m.nD[ODt] + m.nD[OBdy] - 1, 
		m.nD[ODt] + m.nD[OBdy] + (m.nD[IBdy] - 1) / 2 - 1, m.nD[ODt] + m.nD[OBdy] + m.nD[IBdy] - 2,
		m.nD[ODt] + m.nD[OBdy] + m.nD[IBdy] + m.nD[IDt] - 3, 0 },
		MZbot[] = { 6, 4, 1, 6, 4, 6, 1, 6, 2, 6, 3, 4, 2, 
		1, m.nD[ODt], m.nD[ODt] + (m.nD[OBdy] - 1) / 2, m.nD[ODt] + m.nD[OBdy] - 1,
		m.nD[ODt] + m.nD[OBdy] + (m.nD[IBdy] - 1) / 2 - 1, m.nD[ODt] + m.nD[OBdy] + m.nD[IBdy] - 2,
		m.nD[ODt] + m.nD[OBdy] + m.nD[IBdy] + m.nD[IDt] - 3, 0 },
		EDtop[] = { 4, 5, 3, 5, 4, 5, 5, 5, 2, 1, (m.nD[OBdy] + 1) / 2, (m.nD[OBdy] + m.nD[IBdy]) / 2,
		m.nD[IBdy] + (m.nD[IBdy] - 1) / 2, m.nD[IBdy] + m.nD[OBdy] - 1, 0 },
		EDbot[] = { 1, 7, 1, 1, m.nD[IBdy] + m.nD[OBdy] - 1, 0 },
		COtop[] = { 1, 6, 1, 1, m.nD[IBdy] + m.nD[OBdy] - 1, 0 };
	

	MESHPART(sys);

	m.SM = alloc_tensor<double>(0, 9, 0, KMAX * KMULT, 0, 4);

	/*  Outer SOL	dom = 0 -> 2 + 0 + 6	*/
	k = 0;
	ASSEMBLER(sys, 1, 2, m.nD[SOL], m.nD[ODt] - 1, k, -1, -1, -(m.nD[AS] - 1), m.nD[ODt] - 1);
	ASSEMBLER(sys, 1, 0, m.nD[SOL], m.nD[OBdy] - 1, k, -1, 1, -(m.nD[AS] - 1), 0);
	ASSEMBLER(sys, 1, 6, m.nD[SOL], m.nD[ODt], k, -1, 1, 0, 0);
	ni = m.nD[ODt] + m.nD[OBdy] + m.nD[ODt] - 2;
	nj = m.nD[SOL];
	PRINTMESH(sys, 1, 0, k, ni, nj, OSbot, top0);
	printf("(OSOL) > ");
	/*  Inner SOL	dom = 2 -> 8 + 1 + 4	*/
	k = 0;
	ASSEMBLER(sys, 3, 8, m.nD[SOL], m.nD[IDt], k, -1, -1, 0, m.nD[IDt] - 1);
	ASSEMBLER(sys, 3, 1, m.nD[SOL], m.nD[IBdy] - 2, k, -1, -1, -(m.nD[AS] - 1), m.nD[IBdy] - 2);
	ASSEMBLER(sys, 3, 4, m.nD[SOL], m.nD[IDt], k, -1, 1, -(m.nD[AS] - 1), 0);
	ni = m.nD[IDt] + m.nD[IBdy] + m.nD[IDt] - 2;
	nj = m.nD[SOL];
	PRINTMESH(sys, 3, 2, k, ni, nj, ISbot, top0);
	printf("(ISOL) > ");
	/*  Upper Private	dom = 5 -> 2 + 4	*/
	k = 0;
	ASSEMBLER(sys, 4, 2, m.nD[APr], m.nD[ODt] - 1, k, -1, -1, m.nD[APr] - 1, m.nD[ODt] - 1);
	ASSEMBLER(sys, 4, 4, m.nD[APr], m.nD[IDt], k, -1, 1, m.nD[APr] - 1, 0);
	nj = m.nD[APr];
	ni = m.nD[ODt] + m.nD[IDt] - 1;	
	PRINTMESH(sys, 4, 5, k, ni, nj, bot0, UPtop);
	printf("(UPr) > ");
	/*  Lower Private	dom = 1 -> 8 + 6	*/
	k = 0;
	ASSEMBLER(sys, 2, 8, m.nD[Pr], m.nD[IDt] - 1, k, -1, -1, m.nD[Pr] - 1, m.nD[IDt] - 1);
	ASSEMBLER(sys, 2, 6, m.nD[Pr], m.nD[ODt], k, -1, 1, m.nD[Pr] - 1, 0);
	nj = m.nD[Pr];
	ni = m.nD[IDt] + m.nD[ODt] - 1;
	PRINTMESH(sys, 2, 1, k, ni, nj, bot0, LPtop);
	printf("(LPr) > ");
	/*  Active SOL	dom = 3 -> 1 + 0	*/
	k = 0;
	ASSEMBLER(sys, 5, 2, m.nD[AS], m.nD[ODt] - 1, k, -1, -1, 0, m.nD[ODt] - 1);
	ASSEMBLER(sys, 5, 0, m.nD[AS], m.nD[OBdy] - 1, k, -1, 1, 0, 0);
	ASSEMBLER(sys, 5, 1, m.nD[AS], m.nD[IBdy] - 1, k, -1, -1, 0, m.nD[IBdy] - 1);
	ASSEMBLER(sys, 5, 4, m.nD[AS], m.nD[IDt], k, -1, 1, 0, 0);
	ni = m.nD[IBdy] + m.nD[OBdy] + m.nD[IDt] + m.nD[ODt] - 3;
	nj = m.nD[AS];
	PRINTMESH(sys, 5, 6, k, ni, nj, MZbot, MZtop);
	printf("(ASOL) > ");
	/*  EDGE Region		dom = 3 -> 1 + 0	*/
	k = 0;
	ASSEMBLER(sys, 6, 0, Nedge, (m.nD[OBdy] - 1) / 2, k, -1, 1, Nedge - 1, (m.nD[OBdy] - 1) / 2);
	ASSEMBLER(sys, 6, 1, Nedge, m.nD[IBdy] - 1, k, -1, -1, Nedge - 1, m.nD[IBdy] - 1);
	ASSEMBLER(sys, 6, 0, Nedge, (m.nD[OBdy] + 1) / 2, k, -1, 1, Nedge - 1, 0);
	ni = m.nD[IBdy] + m.nD[OBdy] - 1;
	nj = Nedge;
	PRINTMESH(sys, 6, 3, k, ni, nj, EDbot, EDtop);
	printf("(Edge) > ");
	/*  CORE Region		dom = 4 -> 1 + 0	*/
	k = 0;
	ASSEMBLER(sys, 7, 0, m.nD[Co] - Nedge + 1, (m.nD[OBdy] - 1) / 2, 
		k, -1, 1, m.nD[Co] - 1, (m.nD[OBdy] - 1) / 2);
	ASSEMBLER(sys, 7, 1, m.nD[Co] - Nedge + 1, m.nD[IBdy] - 1,
		k, -1, -1, m.nD[Co] - 1, m.nD[IBdy] - 1);
	ASSEMBLER(sys, 7, 0, m.nD[Co] - Nedge + 1, (m.nD[OBdy] + 1) / 2, 
		k, -1, 1, m.nD[Co] - 1, 0);
	ni = m.nD[IBdy] + m.nD[OBdy] - 1;
	nj = m.nD[Co] - Nedge + 1;
	PRINTMESH(sys, 7, 4, k, ni, nj, bot0, COtop);
	printf("(Core) > ");
	ONEDPRINT(sys, 7, nj, m.nD[IBdy],m.nD[OBdy], Nedge);
	printf("(1D Data)\tDone!\n");
}


void POSTP_DDNL(SYSALL &sys) {
	usealias_geo(g, sys);
	usealias_fld(f, sys);
	usealias_mes(m, sys);
	int k, i, j;
	int Nedge = m.nD[Co] / 10;
	const int KMULT = 30;
	int ni, nj, ndIBdy, ndOBdy;
	int top0[] = { 0 }, bot0[] = { 0 };
	double Ax, Ay, Bx, By, C = 0;

	int OSbot[] = { 4, 4, 1, 5, 4, 5, 5, 5, 6, 1, m.nD[ODt], m.nD[ODt] + (m.nD[OBdy] - 1) / 2,
		m.nD[ODt] + m.nD[OBdy] - 1, 2 * m.nD[ODt] + m.nD[OBdy] - 2, 0 },
		ISbot[] = { 4, 5, 1, 5, 2, 5, 3, 5, 4, 2, m.nD[IDt], m.nD[IDt] + (m.nD[IBdy] - 1) / 2,
		m.nD[IDt] + m.nD[IBdy] - 1, 2 * m.nD[IDt] + m.nD[IBdy] - 2, 0 },
		UPtop[] = { 2, 1, 1, 3, 4, 1, m.nD[ODt], m.nD[ODt] + m.nD[IDt] - 1, 0 },
		LPtop[] = { 2, 5, 1, 5, 6, 1, m.nD[IDt], m.nD[IDt] + m.nD[ODt] - 1, 0 },
		MZtop[] = { 6, 3, 1, 3, 2, 3, 3, 1, 2, 1, 3, 1, 4,
		1, m.nD[IDt], m.nD[IDt] + (m.nD[IBdy] - 1) / 2, m.nD[IDt] + m.nD[IBdy] - 1,
		m.nD[IDt] + m.nD[IBdy] + (m.nD[OBdy] - 1) / 2 - 1, m.nD[IDt] + m.nD[IBdy] + m.nD[OBdy] - 2,
		m.nD[IDt] + m.nD[IBdy] + m.nD[OBdy] + m.nD[ODt] - 3, 0 },
		MZbot[] = { 6, 2, 1, 6, 2, 6, 3, 6, 4, 6, 1, 2, 2,
		1, m.nD[IDt], m.nD[IDt] + (m.nD[IBdy] - 1) / 2, m.nD[IDt] + m.nD[IBdy] - 1,
		m.nD[IDt] + m.nD[IBdy] + (m.nD[OBdy] - 1) / 2 - 1, m.nD[IDt] + m.nD[IBdy] + m.nD[OBdy] - 2,
		m.nD[IDt] + m.nD[IBdy] + m.nD[OBdy] + m.nD[ODt] - 3, 0 },
		EDtop[] = { 4, 5, 5, 5, 2, 5, 3, 5, 4, 1, (m.nD[OBdy] + 1) / 2, (m.nD[OBdy] + m.nD[IBdy]) / 2,
		m.nD[IBdy] + (m.nD[IBdy] - 1) / 2, m.nD[IBdy] + m.nD[OBdy] - 1, 0 },
		EDbot[] = { 1, 7, 1, 1, m.nD[IBdy] + m.nD[OBdy] - 1, 0 },
		COtop[] = { 1, 6, 1, 1, m.nD[IBdy] + m.nD[OBdy] - 1, 0 };


	MESHPART(sys);

	m.SM = alloc_tensor<double>(0, 9, 0, KMAX * KMULT, 0, 4);

	/*  Outer SOL	dom = 0 -> 2 + 0 + 6	*/
	k = 0;
	ASSEMBLER(sys, 1, 2, m.nD[SOL], m.nD[ODt], k, -1, -1, 0, m.nD[ODt] - 1);
	ASSEMBLER(sys, 1, 0, m.nD[SOL], m.nD[OBdy] - 2, k, -1, -1, -(m.nD[AS] - 1), m.nD[OBdy] - 2);
	ASSEMBLER(sys, 1, 6, m.nD[SOL], m.nD[ODt], k, -1, 1, -(m.nD[AS] - 1), 0);
	ni = m.nD[ODt] + m.nD[OBdy] + m.nD[ODt] - 2;
	nj = m.nD[SOL];
	PRINTMESH(sys, 1, 0, k, ni, nj, OSbot, top0);
	printf("(OSOL) > ");
	/*  Inner SOL	dom = 2 -> 8 + 1 + 4	*/
	k = 0;
	ASSEMBLER(sys, 3, 8, m.nD[SOL], m.nD[IDt] - 1, k, -1, -1, -(m.nD[AS] - 1), m.nD[IDt] - 1);
	ASSEMBLER(sys, 3, 1, m.nD[SOL], m.nD[IBdy] - 1, k, -1, 1, -(m.nD[AS] - 1), 0);
	ASSEMBLER(sys, 3, 4, m.nD[SOL], m.nD[IDt], k, -1, 1, 0, 0);
	ni = m.nD[IDt] + m.nD[IBdy] + m.nD[IDt] - 2;
	nj = m.nD[SOL];
	PRINTMESH(sys, 3, 2, k, ni, nj, ISbot, top0);
	printf("(ISOL) > ");
	/*  Upper Private	dom = 5 -> 2 + 4	*/
	k = 0;
	ASSEMBLER(sys, 4, 2, m.nD[Pr], m.nD[ODt] - 1, k, -1, -1, m.nD[Pr] - 1, m.nD[ODt] - 1);
	ASSEMBLER(sys, 4, 4, m.nD[Pr], m.nD[IDt], k, -1, 1, m.nD[Pr] - 1, 0);
	nj = m.nD[Pr];
	ni = m.nD[ODt] + m.nD[IDt] - 1;
	PRINTMESH(sys, 4, 5, k, ni, nj, bot0, UPtop);
	printf("(UPr) > ");
	/*  Lower Private	dom = 1 -> 8 + 6	*/
	k = 0;
	ASSEMBLER(sys, 2, 8, m.nD[APr], m.nD[IDt] - 1, k, -1, -1, m.nD[APr] - 1, m.nD[IDt] - 1);
	ASSEMBLER(sys, 2, 6, m.nD[APr], m.nD[ODt], k, -1, 1, m.nD[APr] - 1, 0);
	nj = m.nD[APr];
	ni = m.nD[IDt] + m.nD[ODt] - 1;
	PRINTMESH(sys, 2, 1, k, ni, nj, bot0, LPtop);
	printf("(LPr) > ");
	/*  Active SOL	dom = 3 -> 1 + 0	*/
	k = 0;
	ASSEMBLER(sys, 5, 8, m.nD[AS], m.nD[IDt] - 1, k, -1, -1, 0, m.nD[IDt] - 1);
	ASSEMBLER(sys, 5, 1, m.nD[AS], m.nD[IBdy] - 1, k, -1, 1, 0, 0);
	ASSEMBLER(sys, 5, 0, m.nD[AS], m.nD[OBdy] - 1, k, -1, -1, 0, m.nD[OBdy] - 1);
	ASSEMBLER(sys, 5, 6, m.nD[AS], m.nD[ODt], k, -1, 1, 0, 0);
	ni = m.nD[IBdy] + m.nD[OBdy] + m.nD[IDt] + m.nD[ODt] - 3;
	nj = m.nD[AS];
	PRINTMESH(sys, 5, 6, k, ni, nj, MZbot, MZtop);
	printf("(ASOL) > ");
	/*  EDGE Region		dom = 3 -> 1 + 0	*/
	k = 0;
	ASSEMBLER(sys, 6, 0, Nedge, (m.nD[OBdy] - 1) / 2, k, -1, -1, Nedge - 1, (m.nD[OBdy] - 1) / 2);
	ASSEMBLER(sys, 6, 1, Nedge, m.nD[IBdy] - 1, k, -1, 1, Nedge - 1, 0);
	ASSEMBLER(sys, 6, 0, Nedge, (m.nD[OBdy] + 1) / 2, k, -1, -1, Nedge - 1, m.nD[OBdy] - 1);
	ni = m.nD[IBdy] + m.nD[OBdy] - 1;
	nj = Nedge;
	PRINTMESH(sys, 6, 3, k, ni, nj, EDbot, EDtop);
	printf("(Edge) > ");
	/*  CORE Region		dom = 4 -> 1 + 0	*/
	k = 0;
	ASSEMBLER(sys, 7, 0, m.nD[Co] - Nedge + 1, (m.nD[OBdy] - 1) / 2,
			k, -1, -1, m.nD[Co] - 1, (m.nD[OBdy] - 1) / 2);
	ASSEMBLER(sys, 7, 1, m.nD[Co] - Nedge + 1, m.nD[IBdy] - 1,
			k, -1, 1, m.nD[Co] - 1, 0);
	ASSEMBLER(sys, 7, 0, m.nD[Co] - Nedge + 1, (m.nD[OBdy] + 1) / 2,
			k, -1, -1, m.nD[Co] - 1, m.nD[OBdy] - 1);
	ni = m.nD[IBdy] + m.nD[OBdy] - 1;
	nj = m.nD[Co] - Nedge + 1;
	PRINTMESH(sys, 7, 4, k, ni, nj, bot0, COtop);
	printf("(Core) > ");
	ONEDPRINT(sys, 7, nj, m.nD[IBdy], m.nD[OBdy], Nedge);
	printf("(1D Data)\tDone!\n");
}









