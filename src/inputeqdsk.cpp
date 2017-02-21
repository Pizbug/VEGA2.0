/***************************************************/
/*												   */
/*		INPUT FUNCTION							   */
/*												   */
/*		JG LEE									   */
/*		24/NOV/2015								   */
/***************************************************/

#define	_CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cstring>

#include "inputeqdsk.h"
#include "vega.h"
#include "svc.h"
#include "memalloc.h"

char FILE_HEADER[128];
char EQ_TITLE[128];

void get_id(FILE* in, char const* id){
	char str[256];

	rewind(in);
	do{
		fgets(str,256,in);
		if(str[0] == '@' && strncmp(str+1, id, strlen(id)) == 0) break;
	}
	while(!feof(in));

	if(feof(in)){
		printf("@%s : no such item in this file\n", id);
		exit(-1);
	}
}

void read_par(SYSALL& sys){
	usealias_mes(m, sys);
	FILE 	*in;
	char	title[128];
	int 	i;


	in = fopen("INPAR","r");
	if(in==NULL) fprintf(stderr, "INPAR ERROR");
	else 	printf("Read input parameters\n");

	printf("  Memory allocation: MES\t");
	memalloc_mes(sys.mes);
	printf("Done!\n");

	get_id(in, "FILE_HEADER");	fscanf(in, "%s", FILE_HEADER);
	printf("  FILE HEADER\t-> %s\t|", FILE_HEADER);
	get_id(in, "EQ_TITLE"); 	fscanf(in, "%s", EQ_TITLE);
	printf("  EQ title\t-> %s\n", EQ_TITLE);

	get_id(in, "GRID_Dist");
	DoAllDomain(i, MGD)	fscanf(in, "%d", &m.nD[i]);
	printf("  # distribution of grid:\t");
	printf("(Co, IB, OB, ID, OD,SOL, AS, PR,APR)\n\t\t\t\t");
	DoAllDomain(i, MGD)	printf(" %d ", m.nD[i]);
	printf("\n");

	get_id(in, "GRID_SCOPE");
	fscanf(in, "%lf\t%lf\t%lf", &m.GS[0], &m.GS[1], &m.GS[2]);
	printf("  Grid scope: %.4f, %.4f, %.4f\t(XO, SOL, Priv)\n",m.GS[0],m.GS[1],m.GS[2]);
	
	get_id(in, "EDS");
	DoAllDomain(i, MGD)
		fscanf(in, "%lf %lf %lf\n", &m.EDS[i][0], &m.EDS[i][1], &m.EDS[i][2]);
	printf("  Get E, D, S factors\n");

	fclose(in);
}

void read_input(SYSALL& sys){
	FILE	*in;
	char	str[256];
	int 	idum, i, j, rN, zN;
	double	**data, *fpol, *psi, *lim;
	double	ddum;

	//////////////////////////////////////////////////////////
	/*		READ EQDSK FILE 								*/
	//////////////////////////////////////////////////////////

	sprintf(str, "%s/%s", FILE_HEADER, EQ_TITLE);
	in = fopen(str,"r");
	if(in==NULL) fprintf(stderr, "EQDSK READ ERROR");
	else	printf("Read EQDSK file\n");

	while(fgetc(in) != '\n');
	fseek(in, -12, SEEK_CUR);

	fscanf(in, "%d %d %d\n", &idum, &sys.geo.rN, &sys.geo.zN);
	rN = sys.geo.rN;
	zN = sys.geo.zN;

	printf("  rN = %d, zN = %d, ", rN, zN);

	data = alloc_matrix<double>(0, 3, 0, 4);
	fpol = alloc_vector<double>(0, rN - 1);
	psi  = alloc_vector<double>(0, rN * zN - 1);

	DoAllDomain(i, 4)
		DoAllDomain(j, 5)
			fscanf(in, "%lf", &data[i][j]);

	DoAllDomain(i, rN)
		fscanf(in, "%lf", &fpol[i]);


	DoAllDomain(j, 3)
		DoAllDomain(i, rN)
			fscanf(in, "%lf", &ddum);

	DoAllDomain(j, zN)
		DoAllDomain(i, rN)
			fscanf(in, "%lf", &psi[i + j * rN]);

	DoAllDomain(i, rN)
		fscanf(in, "%lf", &ddum);

	fscanf(in, "%d %d", &j, &sys.geo.lnum);
	lim = alloc_vector<double>(0, sys.geo.lnum * 2 - 1);
	printf("lnum = %d\n", sys.geo.lnum);

	DoAllDomain(i, j * 2)
		fscanf(in, "%lf", &ddum);

	DoAllDomain(i, sys.geo.lnum * 2)
		fscanf(in, "%lf", &lim[i]);

	fclose(in);

	//////////////////////////////////////////////////////////
	/*		READ DIVETOR FILE 								*/
	//////////////////////////////////////////////////////////

	sprintf(str, "%s/div.dat", FILE_HEADER);
	in = fopen(str,"r");
	if(in==NULL) fprintf(stderr, "Divotor READ ERROR");
	else	printf("Read DIVERTOR file\n");

	fscanf(in, "NDIV\t%d\n", &sys.geo.ndiv);
	while(fgetc(in)!='\n') ;
	fscanf(in, "%d %d %d %d\n", &sys.geo.dnum, &idum, &idum, &idum);

	printf("  # of Divertor: %d\t", sys.geo.ndiv);
	printf("  Each divertor has (%d) points\n", sys.geo.dnum);

	//////////////////////////////////////////////////////////
	printf("Memory allocation: GEO, FLD\t");
	memalloc_geo(sys.geo);
	memalloc_fld(sys.fld, rN, zN);
	printf("Done!\n");
	//////////////////////////////////////////////////////////

	for(i=0;i<sys.geo.ndiv;i++){
		for(j=0;j<sys.geo.dnum;j++){
			fscanf(in, "%lf %lf\n", &sys.geo.Divt_r[i][j], &sys.geo.Divt_z[i][j]);
		}
	}

	fclose(in);

	//////////////////////////////////////////////////////////
    /*		CONVERTING GEO DATA 							*/
    //////////////////////////////////////////////////////////
    printf("Converting geo data\t\t");
    DoAllDomain(i, rN)
    	sys.geo.r[i] = data[0][3] + i * data[0][0] / (rN - 1.0);

    DoAllDomain(j, zN)
    	sys.geo.z[j] = data[0][4] + data[0][1] / (zN - 1.0) * (j - zN / 2);

    DoAllDomain(j, zN)
    	DoAllDomain(i, rN)
    		sys.geo.Psi[i][j] = psi[i + j * rN]; 

    DoAllDomain(i, sys.geo.lnum){
    	sys.geo.Lmt_r[i] = lim[2 * i];
    	sys.geo.Lmt_z[i] = lim[2 * i + 1];
    }
    sys.geo.Pfc_r[0] = sys.geo.Lmt_r[0];
    sys.geo.Pfc_r[1] = sys.geo.Lmt_r[1];
    sys.geo.Pfc_z[0] = sys.geo.Lmt_z[0];
    sys.geo.Pfc_z[1] = sys.geo.Lmt_z[1];

    DoAllDomain(i, rN){
    	sys.geo.Bt[i] = fpol[i] / sys.geo.r[i];
    	sys.geo.Bt0[i] = fpol[0] - i * (fpol[0] - fpol[rN - 1]) / (rN - 1);
    }

/*    free_matrix<double>(data, 0, 3, 0, 4);
    free_vector<double>(fpol, 0, rN - 1);
    free_vector<double>(psi, 0, rN * zN - 1);
    free_vector<double>(lim, 0, sys.geo.lnum * 2- 1);
*/
    printf("Done!\n");
}
