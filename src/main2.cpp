/***************************************************/
/*												   */
/*		MAIN									   */
/*												   */
/*		JG LEE									   */
/*		14/JAN/2016								   */
/***************************************************/
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>

#include "parameter.h"
#include "svc.h"
#include "vega.h"
#include "memalloc.h"
#include "input.h"
#include "separatrix.h"
#include "meshgenerator.h"
#include "postprocess.h"

extern void memalloc_sys(SYSALL& sys);
void print_mesh(SYSALL& sys, char* filename);

int main(int argc, char* argv[]) {
	SYSALL sys;

	read_NUM(sys);

	memalloc_sys(sys);

	read_GEODATA(sys);

	FIELD(sys);

	BILINEAR(sys);

	GRAD0(sys);

	MULTIX(sys);

	FRAME(sys);

	MESHGEN2(sys);

	if (sys.mes.Type.id == DDN) {
		if (sys.mes.Type.z > 0)
			POSTP_DDNU(sys);
		else
			POSTP_DDNL(sys);
	}
	else if (sys.mes.Sn == 10) {
		POSTP_SNOW(sys);
	}
	else
		POSTP(sys);

	return 0;
}

void print_mesh(SYSALL& sys, char* filename) {
	usealias_geo(g, sys);
	usealias_fld(f, sys);
	usealias_mes(m, sys);
	int i, j;
	FILE *out;

	out = fopen(filename, "w");

	//fprintf(out, "Xpn = %d, Opn = %d, Sn = %d\n", f.xpn, f.opn, m.Sn);
	/*DoAllDomain(i, f.xpn) {
		fprintf(out, "X point: r = %lf, z = %lf\n", f.Xp[i].r, f.Xp[i].z);
	}
	DoAllDomain(i, f.opn) {
		fprintf(out, "O point: r = %lf, z = %lf\n", f.Op[i].r, f.Op[i].z);
	}*/
	DoAllDomain(j, NMAX){
		DoAllDomain(i, m.Sn){
			fprintf(out, "%lf\t%lf\t", m.Mr[i][j], m.Mz[i][j]);
		}
		fprintf(out, "\n");
	}
	fclose(out);
}
