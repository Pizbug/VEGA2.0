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
#include "inputeqdsk.h"
#include "separatrix.h"
#include "meshgenerator.h"
#include "postprocess.h"

extern void memalloc_sys(SYSALL& sys);
void print_mesh(SYSALL& sys, char* filename);
void vega();

int main() {
	vega();
	return 0;
}

void vega() {
	SYSALL sys;

	printf("========================================");
	printf("========================================\n");
	printf("\tVEGA2.0\n\n");

	read_par(sys);
	read_input(sys);
	
	printf("----------------------------------------");
	printf("----------------------------------------\n");
	
	FIELD(sys);

	BILINEAR(sys);

	GRAD0(sys);

	MULTIX(sys);

	FRAME(sys);

	MESHGEN2(sys);

	if (sys.mes.Type.id == DDN) {
		if (sys.mes.Type.z > 0)
			POSTP_DDNL_lite(sys);
		else
			POSTP_DDNU(sys);
	}
	else if (sys.mes.Sn == 10) {
		POSTP_SNOW(sys);
	}
	else
		POSTP(sys);

	printf("All Processes done\n");
	printf("========================================");
	printf("========================================\n");

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
