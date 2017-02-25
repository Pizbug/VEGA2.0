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

void PRINTMESH(SYSALL& sys, int id, int dom, int k, int nH, int nW, int *BC) {
	//usealias_mes(m, sys);
	MESHDATA&	m = sys.mes;
	int		i, j, ii;
	char 	title[256];
	FILE	*fp;

	sprintf(title, "grid.%03d.dat", id);
	fp = fopen(title, "w");
	fprintf(fp, "@DomainID\n%d\n", dom);
	fprintf(fp, "@NodeNum\n%d %d\n", nH, nW);
	fprintf(fp, "@Connect\n");
	DoAllDomain(i, 4)
		fprintf(fp, "%d ", BC[i]);
	fprintf(fp, "\n");

	fprintf(fp, "@Data\n");
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
				fprintf(fp, " %lf\t%lf\n", m.SM[id][i][0], m.SM[id][i][1]);
			//}
			//fprintf(fp, "\n");
	}
	fclose(fp);
}

void POSTP(SYSALL &sys){}
void POSTP_DDNU(SYSALL &sys){}
void POSTP_SNOW(SYSALL &sys){}
void POSTP_DDNL_lite(SYSALL &sys) {
	usealias_geo(g, sys);
	usealias_fld(f, sys);
	usealias_mes(m, sys);
	int k, i, j;
	const int KMULT = 30;
	int ni, nj, ndIBdy, ndOBdy;
	int BC[4] = { 0 };
	double Ax, Ay, Bx, By, C = 0;

	MESHPART(sys);

	m.SM = alloc_tensor<double>(0, 20, 0, KMAX * KMULT, 0, 4);

	//  DOMAIN - 1
	k = 0;
	ASSEMBLER(sys, 1, 8, m.nD[AS], m.nD[IDt], k, -1, -1, 0, m.nD[IDt] - 1);
	ni = m.nD[IDt];
	nj = m.nD[AS];
	BC[0] = -1;	BC[1] = 2;	BC[2] = 5;	BC[3] = 7;
	PRINTMESH(sys, 1, 0, k, ni, nj, BC);
	//  DOMAIN - 2
	k = 0;
	ASSEMBLER(sys, 2, 1, m.nD[AS], m.nD[IBdy], k, -1, 1, 0, 0);
	ni = m.nD[IBdy];
	nj = m.nD[AS];
	BC[0] = 1;	BC[1] = 3;	BC[2] = -3;	BC[3] = 8;
	PRINTMESH(sys, 2, 0, k, ni, nj, BC);
	//  DOMAIN - 3
	k = 0;
	ASSEMBLER(sys, 3, 0, m.nD[AS], m.nD[OBdy], k, -1, -1, 0, m.nD[OBdy] - 1);
	ni = m.nD[OBdy];
	nj = m.nD[AS];
	BC[0] = 2;	BC[1] = 4;	BC[2] = -3;	BC[3] = 11;
	PRINTMESH(sys, 3, 0, k, ni, nj, BC);
	//  DOMAIN - 4
	k = 0;
	ASSEMBLER(sys, 4, 6, m.nD[AS], m.nD[ODt], k, -1, 1, 0, 0);
	ni = m.nD[ODt];
	nj = m.nD[AS];
	BC[0] = 3;	BC[1] = -1;	BC[2] = 6;	BC[3] = 12;
	PRINTMESH(sys, 4, 0, k, ni, nj, BC);
	//  DOMAIN - 5
	k = 0;
	ASSEMBLER(sys, 5, 8, m.nD[APr], m.nD[IDt], k, -1, -1, m.nD[APr] - 1, m.nD[IDt] - 1);
	ni = m.nD[IDt];
	nj = m.nD[APr];
	BC[0] = -1;	BC[1] = 6;	BC[2] = -2;	BC[3] = 1;
	PRINTMESH(sys, 5, 1, k, ni, nj, BC);
	//  DOMAIN - 6
	k = 0;
	ASSEMBLER(sys, 6, 6, m.nD[APr], m.nD[ODt], k, -1, 1, m.nD[APr] - 1, 0);
	ni = m.nD[ODt];
	nj = m.nD[APr];
	BC[0] = 5;	BC[1] = -1;	BC[2] = -2;	BC[3] = 4;
	PRINTMESH(sys, 6, 1, k, ni, nj, BC);
	//  DOMAIN - 7
	k = 0;
	ASSEMBLER(sys, 7, 8, m.nD[SOL], m.nD[IDt], k, -1, -1, -(m.nD[AS] - 1), m.nD[IDt] - 1);
	ni = m.nD[IDt];
	nj = m.nD[SOL];
	BC[0] = -1;	BC[1] = 8;	BC[2] = 1;	BC[3] = -2;
	PRINTMESH(sys, 7, 0, k, ni, nj, BC);
	//  DOMAIN - 8
	k = 0;
	ASSEMBLER(sys, 8, 1, m.nD[SOL], m.nD[IBdy] - 1, k, -1, 1, -(m.nD[AS] - 1), 0);
	ASSEMBLER(sys, 8, 4, m.nD[SOL], 1, k, -1, 1, 0, 0);
	ni = m.nD[IBdy];
	nj = m.nD[SOL];
	BC[0] = 7;	BC[1] = 9;	BC[2] = 2;	BC[3] = -2;
	PRINTMESH(sys, 8, 0, k, ni, nj, BC);
	//  DOMAIN - 9
	k = 0;
	ASSEMBLER(sys, 9, 4, m.nD[SOL], m.nD[IDt], k, -1, 1, 0, 0);
	ni = m.nD[IDt];
	nj = m.nD[SOL];
	BC[0] = 8;	BC[1] = -1;	BC[2] = 13;	BC[3] = -2;
	PRINTMESH(sys, 9, 0, k, ni, nj, BC);
	//  DOMAIN - 10
	k = 0;
	ASSEMBLER(sys, 10, 2, m.nD[SOL], m.nD[ODt], k, -1, -1, 0, m.nD[ODt] - 1);
	ni = m.nD[ODt];
	nj = m.nD[SOL];
	BC[0] = -1;	BC[1] = 11;	BC[2] = 14;	BC[3] = -2;
	PRINTMESH(sys, 10, 0, k, ni, nj, BC);
	//  DOMAIN - 11
	k = 0;
	ASSEMBLER(sys, 11, 2, m.nD[SOL], 1, k, -1, 1, 0, 0);
	ASSEMBLER(sys, 11, 0, m.nD[SOL], m.nD[OBdy] - 1, k, -1, -1, -(m.nD[AS] - 1), m.nD[OBdy] - 2);
	ni = m.nD[OBdy];
	nj = m.nD[SOL];
	BC[0] = 10;	BC[1] = 12;	BC[2] = 3;	BC[3] = -2;
	PRINTMESH(sys, 11, 0, k, ni, nj, BC);
	//  DOMAIN - 12
	k = 0;
	ASSEMBLER(sys, 12, 6, m.nD[SOL], m.nD[ODt], k, -1, 1, -(m.nD[AS] - 1), 0);
	ni = m.nD[ODt];
	nj = m.nD[SOL];
	BC[0] = 11;	BC[1] = -1;	BC[2] = 4;	BC[3] = -2;
	PRINTMESH(sys, 12, 0, k, ni, nj, BC);
	//  DOMAIN - 13
	k = 0;
	ASSEMBLER(sys, 13, 4, m.nD[Pr], m.nD[IDt], k, -1, 1, m.nD[Pr] - 1, 0);
	ni = m.nD[IDt];
	nj = m.nD[Pr];
	BC[0] = 14;	BC[1] = -1;	BC[2] = -2;	BC[3] = 9;
	PRINTMESH(sys, 13, 1, k, ni, nj, BC);
	//  DOMAIN - 14
	k = 0;
	ASSEMBLER(sys, 14, 2, m.nD[Pr], m.nD[ODt], k, -1, -1, m.nD[Pr] - 1, m.nD[ODt] - 1);
	ni = m.nD[ODt];
	nj = m.nD[Pr];
	BC[0] = -1;	BC[1] = 13;	BC[2] = -2;	BC[3] = 10;
	PRINTMESH(sys, 14, 1, k, ni, nj, BC);
}





