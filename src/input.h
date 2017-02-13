/***************************************************/
/*												   */
/*		INPUT HEADER							   */
/*												   */
/*		JG LEE									   */
/*		24/DEC/2015								   */
/***************************************************/

#ifndef _INPUT_H
#define _INPUT_H

#include "vega.h"

void read_NUM(SYSALL& sys);							//	rN, zN, ndiv, dnum, lnum
void read_RZPSI(char* path_rzpsi, SYSALL& sys);		//	r[rN], z[zN], Psi[rN][zN]
void read_NORPSI(char* path_Norpsi, SYSALL& sys);	//	Norpsi[100] -> Normailized Psi value
void read_BT(char* path_Bt, SYSALL& sys);			//	Bt[rN]
void read_DIVT(char* path_Divt, SYSALL& sys);		//	Divt_r[dnum][ndiv], Divt_z[][], dnum ~ ndiv ~ 4
void read_LIMT(char* path_Limt, SYSALL& sys);		//	Limt_r[lnum], Limt_z[], Pfc_r[2],Pfc_z[2], lnum ~ 55, # of Pfc ~ 2
void read_GEODATA(SYSALL& sys);						//	ALL THE PATH IN HERE!!





#endif
