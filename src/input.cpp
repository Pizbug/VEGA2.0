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

#include "input.h"
#include "svc.h"

//#define GEO_FILE_HEADER "DDN_DIIID"
char GEO_FILE_HEADER[64];

void read_NUM(SYSALL& sys){		//	DO THIS FIRST!!
	FILE* fp;
	int		dummy;
	char	title[255];
	
	fp=fopen("GEO_FILE_HEADER","r");
	fscanf(fp,"%s",&GEO_FILE_HEADER);
	fclose(fp);

	/*	# of R and Z 	*/
	sprintf(title, "input/%s/PSI_RZ.dat", GEO_FILE_HEADER);
	if((fp = fopen(title,"r")) == NULL){
    	fprintf(stderr, "Can't read file\tread_NUM\trzpsi_DN.dat\n");
      	exit(-1);
    }
    fscanf(fp, "ZONE I=%d J=%d ", &sys.geo.zN, &sys.geo.rN);
    fclose(fp);

    /*	# of divertor and its nodes		*/
	sprintf(title, "input/%s/div.dat", GEO_FILE_HEADER);
    if((fp = fopen(title,"r")) == NULL){
    	fprintf(stderr, "Can't read file\tread_NUM\tdiv.dat\n");
      	exit(-1) ;
    }
	fscanf(fp, "NDIV\t%d\n", &sys.geo.ndiv);
	while(fgetc(fp)!='\n') ;
	fscanf(fp, "%d %d %d %d\n", &sys.geo.dnum, &dummy, &dummy, &dummy);
    fclose(fp);

    /*	# of limiter 	*/
	sprintf(title, "input/%s/limiter.in", GEO_FILE_HEADER);
    if((fp = fopen(title,"r")) == NULL){
    	fprintf(stderr, "Can't read file\tread_NUM\tlimiter.in\n");
      	exit(-1) ;
    }
	fscanf(fp, "lnum = %d\n\n", &sys.geo.lnum);
    fclose(fp);
}

/* PSI of R and Z */
void read_RZPSI(char* path, SYSALL& sys){	//	After all the memory allocated for sys
    int 	i, j;
    FILE 	*fp;

    if((fp = fopen(path,"r")) == NULL){
    	fprintf(stderr, "Can't read file\t%s\n",path);
    	exit(-1);
    }

    fscanf(fp, "ZONE I=%d J=%d", &sys.geo.zN, &sys.geo.rN);
	while(fgetc(fp)!='\n') ;

	for (i = 0; i<sys.geo.rN; i++) {
		for (j = 0; j<sys.geo.zN; j++) {
            fscanf(fp,"%lf %lf %lf",&sys.geo.r[i],&sys.geo.z[j],&(sys.geo.Psi[i][j]));
        }
    }
    fclose(fp);
}

/* Normalized PSI 
void read_NORPSI(char* path, SYSALL& sys){
	int 	i;
	FILE	*fp;

	if((fp = fopen(path,"r")) == NULL){
    	fprintf(stderr, "Can't read file\t%s\n",path);
    	exit(-1);
    }

	while(fgetc(fp)!='\n') ;

	for(i=0;i<100;i++){
		fscanf(fp,"%lf",&sys.geo.Norpsi[i]);
		while(fgetc(fp)!='\n')	;
	}
	fclose(fp);
}*/

/* Toroidal Magnetic Field */
void read_BT(char* path, SYSALL& sys){
	int 	i;
	double 	dummy[2];
	FILE*	fp;

	if((fp = fopen(path,"r")) == NULL){
    	fprintf(stderr, "Can't read file\t%s\n",path);
    	exit(-1);
    }

	sys.geo.Bt[0] = sys.geo.Bt[sys.geo.rN-1] = 0;
	sys.geo.Bt0[0] = sys.geo.Bt0[sys.geo.rN - 1] = 0;
	while(fgetc(fp)!='\n') ;

	for(i=1;i<sys.geo.rN-1;i++){
		fscanf(fp, "%lf %lf %lf", &dummy[0], &sys.geo.Bt0[i], &sys.geo.Bt[i]);
		while(fgetc(fp)!='\n')	;
	}
	fclose(fp);
}

/* Divertor Structure */
void read_DIVT(char* path, SYSALL& sys){
	int 	i,j;
	FILE*	fp;
	int		dummy[3];

	if((fp = fopen(path,"r")) == NULL){
    	fprintf(stderr, "Can't read file\t%s\n",path);
    	exit(-1);
    }

	fscanf(fp, "NDIV\t%d\n", &sys.geo.ndiv);
	while(fgetc(fp)!='\n') ;
	fscanf(fp, "%d %d %d %d\n", &sys.geo.dnum, &dummy[0], &dummy[1], &dummy[2]);

	// r -> x position of i'th divt plate
	// z -> y position of i'th divt plate
	/*	x1(1) y1(1)
		x2(1) y2(1)
		x3(1) y3(1)
		x4(1) y4(1)
		x1(2) y1(2)
		...........	*/

	for(i=0;i<sys.geo.ndiv;i++){
		for(j=0;j<sys.geo.dnum;j++){
			fscanf(fp, "%lf %lf\n", &sys.geo.Divt_r[i][j], &sys.geo.Divt_z[i][j]);
		}
	}
	fclose(fp);
}

/* Limiter Structure */
void read_LIMT(char* path, SYSALL& sys){
	int 	i;
	FILE*	fp;

	if((fp = fopen(path,"r")) == NULL){
    	fprintf(stderr, "Can't read file\t%s\n",path);
    	exit(-1);
    }

	fscanf(fp, "lnum = %d\n\n", &sys.geo.lnum);

	DoAllDomain(i,sys.geo.lnum){
		fscanf(fp, "%lf %lf\n", &sys.geo.Lmt_r[i], &sys.geo.Lmt_z[i]);
	}
	
	while(fgetc(fp) != 'C') ;
	while(fgetc(fp) != '\n') ;

	DoAllDomain(i,2){					//	YOU SHOULD CHECK THE NUMBER OF PFC!
		fscanf(fp, "%lf %lf\n", &sys.geo.Pfc_r[i],&sys.geo.Pfc_z[i]);
	}
	fclose(fp);
}

/*	READ GEODATA	*/
void read_GEODATA(SYSALL& sys){
	const int stringsize = 255;
	char 	path_rzpsi[stringsize]	,
			//path_Norpsi[stringsize]	,
			path_Bt[stringsize]		,
			path_Divt[stringsize]	,
			path_Limt[stringsize]	;

	sprintf(path_rzpsi, "input/%s/PSI_RZ.dat", GEO_FILE_HEADER)	;
	sprintf(path_Divt,	"input/%s/div.dat", GEO_FILE_HEADER)	;
	sprintf(path_Limt,	"input/%s/limiter.in", GEO_FILE_HEADER)	;
	//sprintf(path_Norpsi,"input/Profs_psi.dat")					;
	sprintf(path_Bt,	"input/%s/BT.dat", GEO_FILE_HEADER)		;

	read_RZPSI(path_rzpsi,sys)		;
	//read_NORPSI(path_Norpsi,sys)	;
	read_BT(path_Bt,sys)			;
	read_DIVT(path_Divt,sys)		;
	read_LIMT(path_Limt,sys)		;
}
